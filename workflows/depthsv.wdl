version 1.0

# depthSV — correct and associate, scattered over regions
#
# This is a thin wrapper. Every task shells out to the same scripts the SLURM
# path uses, so there is one implementation of each stage and the workflow
# engine only decides what runs where. WDL is used because it runs on both
# target platforms; the SLURM path stays available and is not replaced.
#
# The join is deliberately not a task here: it reads every per-sample file at
# once and is run once per cohort, which suits a single large machine better
# than a scatter. Run scripts/join.sh, then start this workflow from the matrix.
# Every shard localises the whole matrix (the `depth_matrix` File), which is
# fine up to a few thousand samples and wasteful beyond; at biobank scale use
# workflows/depthsv_blocks.wdl, where no task ever sees the whole matrix.
#
# Finish with scripts/export.sh over the collected shards (one run-level
# step: coverage check, count suppression, empirical threshold).

workflow depthsv {
  input {
    File   depth_matrix
    File   depth_matrix_index
    File   pcs
    File   coverage
    File   phenotypes
    File   pheno_manifest
    File   regions_file                 # one region per line (scripts/regions.sh)

    Int    ndim       = 16
    Int    min_obs    = 100
    Float  max_share  = 0.5
    # Permutation maxima per shard (linear analyses) for scripts/export.sh's
    # empirical threshold; 0 disables. The export itself runs once on the
    # collected shards, outside this scatter.
    Int    perms      = 0
    Int    perm_seed  = 1
    # Ploidy model for chrX/chrY: a sex table (SAMPLE + sex_col, M/F or 1/2)
    # and the pseudo-autosomal BED for the build (conf/par.grch38.bed).
    File?  sex
    String sex_col    = "SEX"
    File?  par
    Float  winsor_log2 = -3.0
    Int    block_mb   = 64              # DSV_BLOCK_BYTES for the R workers
    Int    chunk      = 2000            # DSV_CHUNK
    # Published by CI on every release tag for linux/amd64; pin the digest
    # the release job prints for a reproducible run.
    String docker     = "ghcr.io/jlanej/depthsv:0.2.0"

    # Memory is bounded by block_mb (each worker parses one block at a time),
    # not by the cohort size; 4 workers at 64 MB blocks fit in 16 GB at any
    # width. Disk must hold the whole matrix plus one shard's outputs.
    Int    correct_cpu    = 4
    Int    correct_mem_gb = 16
    Int    analyze_cpu    = 4
    Int    analyze_mem_gb = 16
    Int    disk_gb        = 100
    # 0 keeps every task on non-preemptible instances. Region-sized units are
    # cheap to lose, so raising this is where the cost saving comes from.
    Int    preemptible    = 3
  }

  Array[String] regions = read_lines(regions_file)

  scatter (region in regions) {
    call correct {
      input:
        depth_matrix = depth_matrix,
        depth_matrix_index = depth_matrix_index,
        pcs = pcs, coverage = coverage,
        region = region, ndim = ndim,
        sex = sex, sex_col = sex_col, par = par, winsor_log2 = winsor_log2,
        block_mb = block_mb, chunk = chunk,
        docker = docker, cpu = correct_cpu, mem_gb = correct_mem_gb,
        disk_gb = disk_gb, preemptible = preemptible
    }

    call analyze {
      input:
        corrected = correct.corrected,
        corrected_index = correct.corrected_index,
        phenotypes = phenotypes,
        pheno_manifest = pheno_manifest,
        pcs = pcs, ndim = ndim,
        region = region, min_obs = min_obs, max_share = max_share,
        perms = perms, perm_seed = perm_seed,
        block_mb = block_mb, chunk = chunk,
        docker = docker, cpu = analyze_cpu, mem_gb = analyze_mem_gb,
        disk_gb = disk_gb, preemptible = preemptible
    }
  }

  output {
    Array[File] corrected        = correct.corrected
    Array[File] corrected_index  = correct.corrected_index
    Array[File] correction_stats = correct.stats
    Array[File] correction_logs  = correct.log
    # One file per (region x analysis); flattened so the caller sees a single
    # list of summary-statistic shards for scripts/export.sh.
    Array[File] sumstats           = flatten(analyze.sumstats)
    Array[File] permutation_maxima = flatten(analyze.permmax)
    Array[File] analysis_logs      = flatten(analyze.logs)
  }
}

task correct {
  input {
    File depth_matrix
    File depth_matrix_index
    File pcs
    File coverage
    String region
    Int ndim
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

  # dxCompiler streams these through dxfuse instead of copying them; Cromwell
  # (localization_optional) hands the task the gs:// URL, which htslib reads
  # through the index once it has a token — the scripts accept URLs.
  parameter_meta {
    depth_matrix:       { stream: true, localization_optional: true }
    depth_matrix_index: { stream: true, localization_optional: true }
  }

  # A localised matrix needs its index beside it under the expected name
  # (the index flavour, .tbi or .csi, is read from its content); a remote one
  # is passed through with its index expected at <url>.tbi / .csi.
  command <<<
    set -euo pipefail
    export TMPDIR="$PWD/tmp"; mkdir -p "$TMPDIR" in out
    export DSV_BLOCK_BYTES=$(( ~{block_mb} * 1024 * 1024 ))
    export DSV_CHUNK=~{chunk}
    case "~{depth_matrix}" in
      gs://*)
        export GCS_OAUTH_TOKEN="$(curl -s -H 'Metadata-Flavor: Google' \
          'http://metadata.google.internal/computeMetadata/v1/instance/service-accounts/default/token' \
          | sed -E 's/.*"access_token":"([^"]+)".*/\1/')"
        matrix="~{depth_matrix}" ;;
      *://*)
        matrix="~{depth_matrix}" ;;
      *)
        ln -s "~{depth_matrix}" in/depth.matrix.txt.gz
        if bgzip -dc "~{depth_matrix_index}" | head -c 3 | grep -q CSI; then ext=csi; else ext=tbi; fi
        ln -s "~{depth_matrix_index}" "in/depth.matrix.txt.gz.$ext"
        matrix=in/depth.matrix.txt.gz ;;
    esac

    /opt/depthsv/scripts/correct.sh \
      --matrix "$matrix" \
      --pcs "~{pcs}" \
      --coverage "~{coverage}" \
      --region "~{region}" \
      --out out \
      --ndim ~{ndim} \
      --winsor-log2 ~{winsor_log2} \
      ~{if defined(sex) then "--sex " + sex + " --sex-col " + sex_col else ""} \
      ~{"--par " + par} \
      --jobs ~{cpu} \
      --threads ~{cpu}
  >>>

  output {
    File corrected       = glob("out/corrected_ndim*.txt.gz")[0]
    File corrected_index = glob("out/corrected_ndim*.txt.gz.*i")[0]
    File stats           = glob("out/stats_ndim*.txt.gz")[0]
    File log             = glob("out/corrected_ndim*.log")[0]
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

task analyze {
  input {
    File corrected
    File corrected_index
    File phenotypes
    File pheno_manifest
    File pcs
    Int ndim
    String region
    Int min_obs
    Float max_share
    Int perms
    Int perm_seed
    Int block_mb
    Int chunk
    String docker
    Int cpu
    Int mem_gb
    Int disk_gb
    Int preemptible
  }

  # The corrected file keeps its name: analyze.sh reads ndim from it.
  command <<<
    set -euo pipefail
    export TMPDIR="$PWD/tmp"; mkdir -p "$TMPDIR" in out
    export DSV_BLOCK_BYTES=$(( ~{block_mb} * 1024 * 1024 ))
    export DSV_CHUNK=~{chunk}
    name="$(basename "~{corrected}")"
    ln -s "~{corrected}" "in/$name"
    if bgzip -dc "~{corrected_index}" | head -c 3 | grep -q CSI; then ext=csi; else ext=tbi; fi
    ln -s "~{corrected_index}" "in/$name.$ext"

    /opt/depthsv/scripts/analyze.sh \
      --corrected "in/$name" \
      --pheno "~{phenotypes}" \
      --pheno-manifest "~{pheno_manifest}" \
      --pcs "~{pcs}" --ndim ~{ndim} \
      --region "~{region}" \
      --out out \
      --min-obs ~{min_obs} --max-share ~{max_share} \
      --perms ~{perms} --perm-seed ~{perm_seed} \
      --jobs ~{cpu} \
      --threads ~{cpu}
  >>>

  output {
    Array[File] sumstats = glob("out/*.txt.gz")
    Array[File] permmax  = glob("out/*.permmax.txt")
    Array[File] logs     = glob("out/*.log")
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
