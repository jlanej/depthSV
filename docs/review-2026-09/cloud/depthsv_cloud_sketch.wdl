version 1.0

# Sketch only (review scratch): the shape the cloud path needs.
#  - join as scatter over SAMPLE BLOCKS, then a per-WINDOW paste, so no task ever
#    sees the whole matrix and the multi-TB tabix file is never materialised;
#  - correct+analyze fused per window, reading the window shard (GBs, not TBs);
#  - streaming hints for both engines on the legacy whole-matrix path.

workflow depthsv_cloud {
  input {
    Array[File]   mosdepth_files          # per-sample *.regions.bed.gz (object-store URIs)
    File          windows                 # one region per line, e.g. chr1:1-10000000
    Int           block_size = 1000       # samples per join block (RAP guide: <=1000 inputs/job)
    File          pcs
    File          coverage
    File          phenotypes
    File          pheno_manifest
    Int           ndim = 16
    Int           min_obs = 100
    String        docker = "docker.io/<org>/depthsv:0.1.0"   # published, linux/amd64
    Int           block_mb = 1024         # DSV_BLOCK_BYTES for the R workers, in MB
    Int           preemptible = 3
  }

  Array[String] window_list = read_lines(windows)

  # --- join, stage 1: one task per block of samples --------------------------
  Int n_blocks = ceil(length(mosdepth_files) / (block_size * 1.0))
  scatter (b in range(n_blocks)) {
    Int lo = b * block_size
    Int hi = if (lo + block_size < length(mosdepth_files)) then lo + block_size else length(mosdepth_files)
    scatter (i in range(hi - lo)) { File block_file = mosdepth_files[lo + i] }
    call join_block {
      input: files = block_file, windows = windows, block_id = b, docker = docker, preemptible = preemptible
    }
  }

  # --- join, stage 2 + correct + analyze: one task per window ----------------
  scatter (w_idx in range(length(window_list))) {
    scatter (blk in join_block.window_shards) { File shard_for_window = blk[w_idx] }
    call window_unit {
      input:
        region = window_list[w_idx], block_shards = shard_for_window,
        pcs = pcs, coverage = coverage, phenotypes = phenotypes, pheno_manifest = pheno_manifest,
        ndim = ndim, min_obs = min_obs, block_mb = block_mb, docker = docker, preemptible = preemptible
    }
  }

  output {
    Array[File]        matrix_windows = window_unit.matrix_window
    Array[File]        corrected      = window_unit.corrected
    Array[Array[File]] sumstats       = window_unit.sumstats
    Array[File]        logs           = window_unit.log
  }
}

# Verify row counts / coordinates for a block of samples and emit one
# bgzipped column-block per window (rows of the window x samples of the block).
task join_block {
  input {
    Array[File] files
    File        windows
    Int         block_id
    String      docker
    Int         preemptible
    Int         cpu = 4
    Int         mem_gb = 8
    Int         disk_gb = 200
  }
  command <<<
    set -euo pipefail
    export TMPDIR="$PWD/tmp"; mkdir -p "$TMPDIR" out
    # (hypothetical) join.sh --mode block: extract+verify every column, then cut
    # the block by window; each window file is <window>.block~{block_id}.tsv.gz
    /opt/depthsv/scripts/join.sh --mode block --manifest ~{write_lines(files)} \
        --windows ~{windows} --block-id ~{block_id} --out out --jobs ~{cpu} --threads ~{cpu}
  >>>
  output {
    # ordered as the windows file, so the caller can index by window position
    Array[File] window_shards = read_lines("out/window_shards.list")
    File        block_meta    = "out/block.meta"
  }
  runtime {
    docker: docker
    cpu: cpu
    memory: mem_gb + " GB"
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: preemptible
    maxRetries: 1
  }
}

task window_unit {
  input {
    String      region
    Array[File] block_shards
    File        pcs
    File        coverage
    File        phenotypes
    File        pheno_manifest
    Int         ndim
    Int         min_obs
    Int         block_mb
    String      docker
    Int         preemptible
    Int         cpu = 4
    Int         mem_gb = 16
    Int         disk_gb = 200
  }
  command <<<
    set -euo pipefail
    export TMPDIR="$PWD/tmp"; mkdir -p "$TMPDIR" in out
    export DSV_BLOCK_BYTES=$(( ~{block_mb} * 1024 * 1024 ))
    # (hypothetical) paste the block shards for this window into one indexed slice
    /opt/depthsv/scripts/join.sh --mode paste --shards ~{write_lines(block_shards)} \
        --region "~{region}" --out in --threads ~{cpu}
    /opt/depthsv/scripts/correct.sh --matrix in/depth.matrix.txt.gz --pcs "~{pcs}" \
        --coverage "~{coverage}" --region "~{region}" --out out --ndim ~{ndim} \
        --jobs ~{cpu} --threads ~{cpu}
    /opt/depthsv/scripts/analyze.sh --corrected "$(ls out/corrected_ndim*.txt.gz)" \
        --pheno "~{phenotypes}" --pheno-manifest "~{pheno_manifest}" --region "~{region}" \
        --out out --jobs ~{cpu} --threads ~{cpu} -- --minObs ~{min_obs}
    cat out/*.log > "out/~{region}.log" || true
  >>>
  output {
    File        matrix_window = "in/depth.matrix.txt.gz"
    File        corrected     = glob("out/corrected_ndim*.txt.gz")[0]
    Array[File] sumstats      = glob("out/*.txt.gz")
    File        log           = "out/~{region}.log"
  }
  runtime {
    docker: docker
    cpu: cpu
    memory: mem_gb + " GB"
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: preemptible
    maxRetries: 1
  }
}

# Legacy whole-matrix path with streaming hints for both engines. The scripts
# would additionally need to accept a URL (no `[ -s ]`, header via `tabix -H`).
task correct_streamed {
  input {
    File   depth_matrix
    File   depth_matrix_index
    File   pcs
    File   coverage
    String region
    Int    ndim
    String docker
  }
  parameter_meta {
    depth_matrix:       { stream: true, localization_optional: true }
    depth_matrix_index: { stream: true, localization_optional: true }
  }
  command <<<
    set -euo pipefail
    case "~{depth_matrix}" in
      gs://*) export GCS_OAUTH_TOKEN="$(curl -s -H 'Metadata-Flavor: Google' \
                'http://metadata.google.internal/computeMetadata/v1/instance/service-accounts/default/token' \
                | sed -E 's/.*"access_token":"([^"]+)".*/\1/')" ;;
    esac
    mkdir -p out
    /opt/depthsv/scripts/correct.sh --matrix "~{depth_matrix}" --pcs "~{pcs}" \
        --coverage "~{coverage}" --region "~{region}" --out out --ndim ~{ndim}
  >>>
  output { File corrected = glob("out/corrected_ndim*.txt.gz")[0] }
  runtime {
    docker: docker
    cpu: 4
    memory: "16 GB"
    disks: "local-disk 100 HDD"
    preemptible: 3
    dx_instance_type: "mem2_ssd1_v2_x4"
  }
}
