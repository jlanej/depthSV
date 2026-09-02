version 1.0

# depthSV — the cloud shape: join by sample block, then one task per window
#
# Neither cloud platform can hand one task a 500,000-file manifest or copy a
# multi-terabyte matrix into every scatter shard. This workflow avoids both:
#
#   join_block   scripts/join.sh over a block of ~1,000 samples -> an indexed
#                matrix of every bin x that block (a few GB)
#   window_unit  scripts/join_paste.sh reads one region from every block
#                through its index, pastes the slices into the window's
#                whole-cohort matrix, then correct.sh and analyze.sh run on
#                it — the same scripts, the same unit of work as the SLURM
#                path and workflows/depthsv.wdl
#
# Finish with scripts/export.sh over the collected shards (one run-level
# step: coverage check, count suppression, empirical threshold).

workflow depthsv_blocks {
  input {
    Array[File] mosdepth_files          # per-sample *.regions.bed.gz
    File   regions_file                 # one region per line (scripts/regions.sh)
    Int    block_size = 1000            # samples per join block
    File   pcs
    File   coverage
    File   phenotypes
    File   pheno_manifest

    Int    ndim       = 16
    Int    min_obs    = 100
    Float  max_share  = 0.5
    Int    perms      = 0
    Int    perm_seed  = 1
    File?  sex
    String sex_col    = "SEX"
    File?  par
    Float  winsor_log2 = -3.0
    Int    block_mb   = 64              # DSV_BLOCK_BYTES for the R workers
    Int    chunk      = 2000            # DSV_CHUNK
    String docker     = "ghcr.io/jlanej/depthsv:0.2.0"

    Int    join_cpu     = 4
    Int    join_mem_gb  = 8
    Int    join_disk_gb = 300
    Int    unit_cpu     = 4
    Int    unit_mem_gb  = 16
    Int    unit_disk_gb = 100
    Int    preemptible  = 3
  }

  Array[String] regions = read_lines(regions_file)
  Int n_files  = length(mosdepth_files)
  Int n_blocks = ceil(n_files / (block_size * 1.0))

  scatter (b in range(n_blocks)) {
    Int lo = b * block_size
    Int hi = if (lo + block_size < n_files) then lo + block_size else n_files
    scatter (i in range(hi - lo)) { File block_file = mosdepth_files[lo + i] }
    call join_block {
      input:
        files = block_file, block_id = b,
        docker = docker, cpu = join_cpu, mem_gb = join_mem_gb, disk_gb = join_disk_gb,
        preemptible = preemptible
    }
  }

  scatter (region in regions) {
    call window_unit {
      input:
        region = region,
        blocks = join_block.matrix, block_indexes = join_block.matrix_index,
        pcs = pcs, coverage = coverage, phenotypes = phenotypes, pheno_manifest = pheno_manifest,
        ndim = ndim, min_obs = min_obs, max_share = max_share, perms = perms, perm_seed = perm_seed,
        sex = sex, sex_col = sex_col, par = par, winsor_log2 = winsor_log2,
        block_mb = block_mb, chunk = chunk,
        docker = docker, cpu = unit_cpu, mem_gb = unit_mem_gb, disk_gb = unit_disk_gb,
        preemptible = preemptible
    }
  }

  output {
    Array[File] block_matrices     = join_block.matrix
    Array[File] block_manifests    = join_block.manifest
    Array[File] corrected          = window_unit.corrected
    Array[File] corrected_index    = window_unit.corrected_index
    Array[File] correction_stats   = window_unit.stats
    Array[File] sumstats           = flatten(window_unit.sumstats)
    Array[File] permutation_maxima = flatten(window_unit.permmax)
    Array[File] logs               = flatten(window_unit.logs)
  }
}

task join_block {
  input {
    Array[File] files
    Int block_id
    String docker
    Int cpu
    Int mem_gb
    Int disk_gb
    Int preemptible
  }

  command <<<
    set -euo pipefail
    export TMPDIR="$PWD/tmp"; mkdir -p "$TMPDIR" out
    /opt/depthsv/scripts/join.sh --manifest ~{write_lines(files)} --out out \
        --jobs ~{cpu} --threads ~{cpu}
    # One name for either index flavour; the reader tells them apart by content.
    if [ -s out/depth.matrix.txt.gz.tbi ]; then cp out/depth.matrix.txt.gz.tbi out/index; else cp out/depth.matrix.txt.gz.csi out/index; fi
    echo "block ~{block_id}: $(grep -c . ~{write_lines(files)}) samples" >&2
  >>>

  output {
    File matrix       = "out/depth.matrix.txt.gz"
    File matrix_index = "out/index"
    File manifest     = "out/depth.matrix.manifest"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: mem_gb + " GB"
    disks: "local-disk " + disk_gb + " HDD"
    bootDiskSizeGb: 20
    preemptible: preemptible
    maxRetries: 1
    dx_timeout: "24H"
  }
}

task window_unit {
  input {
    String region
    Array[File] blocks
    Array[File] block_indexes
    File pcs
    File coverage
    File phenotypes
    File pheno_manifest
    Int ndim
    Int min_obs
    Float max_share
    Int perms
    Int perm_seed
    File? sex
    String sex_col
    File? par
    Float winsor_log2
    Int block_mb
    Int chunk
    String docker
    Int cpu
    Int mem_gb
    Int disk_gb
    Int preemptible
  }

  parameter_meta {
    blocks:        { stream: true }
    block_indexes: { stream: true }
  }

  command <<<
    set -euo pipefail
    export TMPDIR="$PWD/tmp"; mkdir -p "$TMPDIR" in out
    export DSV_BLOCK_BYTES=$(( ~{block_mb} * 1024 * 1024 ))
    export DSV_CHUNK=~{chunk}

    # Each block beside its index under the name tabix expects.
    blocks=(~{sep=' ' blocks}); indexes=(~{sep=' ' block_indexes})
    : > in/blocks.list
    for i in "${!blocks[@]}"; do
      ln -s "${blocks[$i]}" "in/block$i.txt.gz"
      if bgzip -dc "${indexes[$i]}" | head -c 3 | grep -q CSI; then ext=csi; else ext=tbi; fi
      ln -s "${indexes[$i]}" "in/block$i.txt.gz.$ext"
      echo "$PWD/in/block$i.txt.gz" >> in/blocks.list
    done

    /opt/depthsv/scripts/join_paste.sh --blocks in/blocks.list --region "~{region}" \
        --out in/window --threads ~{cpu}

    /opt/depthsv/scripts/correct.sh \
      --matrix in/window/depth.matrix.txt.gz \
      --pcs "~{pcs}" --coverage "~{coverage}" \
      --region "~{region}" --out out --ndim ~{ndim} \
      --winsor-log2 ~{winsor_log2} \
      ~{if defined(sex) then "--sex " + sex + " --sex-col " + sex_col else ""} \
      ~{"--par " + par} \
      --jobs ~{cpu} --threads ~{cpu}

    /opt/depthsv/scripts/analyze.sh \
      --corrected "$(ls out/corrected_ndim*.txt.gz)" \
      --pheno "~{phenotypes}" --pheno-manifest "~{pheno_manifest}" \
      --pcs "~{pcs}" --ndim ~{ndim} \
      --region "~{region}" --out out \
      --min-obs ~{min_obs} --max-share ~{max_share} \
      --perms ~{perms} --perm-seed ~{perm_seed} \
      --jobs ~{cpu} --threads ~{cpu}
  >>>

  output {
    File corrected          = glob("out/corrected_ndim*.txt.gz")[0]
    File corrected_index    = glob("out/corrected_ndim*.txt.gz.*i")[0]
    File stats              = glob("out/stats_ndim*.txt.gz")[0]
    Array[File] sumstats    = glob("out/*.*.*.txt.gz")
    Array[File] permmax     = glob("out/*.permmax.txt")
    Array[File] logs        = glob("out/*.log")
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: mem_gb + " GB"
    disks: "local-disk " + disk_gb + " HDD"
    bootDiskSizeGb: 20
    preemptible: preemptible
    maxRetries: 1
    dx_timeout: "12H"
  }
}
