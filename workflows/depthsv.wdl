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

workflow depthsv {
  input {
    File   depth_matrix
    File   depth_matrix_index
    File   pcs
    File   coverage
    File   phenotypes
    File   pheno_manifest
    Array[String] regions

    Int    ndim       = 16
    Int    min_obs    = 100
    String docker     = "depthsv:dev"

    Int    correct_cpu    = 4
    Int    correct_mem_gb = 16
    Int    analyze_cpu    = 4
    Int    analyze_mem_gb = 16
    Int    disk_gb        = 100
    # 0 keeps every task on non-preemptible instances. Region-sized units are
    # cheap to lose, so raising this is where the cost saving comes from.
    Int    preemptible    = 3
  }

  scatter (region in regions) {
    call correct {
      input:
        depth_matrix = depth_matrix,
        depth_matrix_index = depth_matrix_index,
        pcs = pcs, coverage = coverage,
        region = region, ndim = ndim,
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
        region = region, min_obs = min_obs,
        docker = docker, cpu = analyze_cpu, mem_gb = analyze_mem_gb,
        disk_gb = disk_gb, preemptible = preemptible
    }
  }

  output {
    Array[File] corrected      = correct.corrected
    Array[File] correction_stats = correct.stats
    # One file per (region x analysis); flattened so the caller sees a single
    # list of summary-statistic shards to concatenate.
    Array[File] sumstats       = flatten(analyze.sumstats)
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
    String docker
    Int cpu
    Int mem_gb
    Int disk_gb
    Int preemptible
  }

  # The index must sit beside the matrix under its expected name; localisation
  # does not guarantee that, so link both into one directory first.
  command <<<
    set -euo pipefail
    mkdir -p in out
    ln -s "~{depth_matrix}"       in/depth.matrix.txt.gz
    ln -s "~{depth_matrix_index}" in/depth.matrix.txt.gz.tbi

    /opt/depthsv/scripts/correct.sh \
      --matrix in/depth.matrix.txt.gz \
      --pcs "~{pcs}" \
      --coverage "~{coverage}" \
      --region "~{region}" \
      --out out \
      --ndim ~{ndim} \
      --jobs ~{cpu} \
      --threads ~{cpu}
  >>>

  output {
    File corrected       = glob("out/corrected_ndim*.txt.gz")[0]
    File corrected_index = glob("out/corrected_ndim*.txt.gz.tbi")[0]
    File stats           = glob("out/stats_ndim*.txt.gz")[0]
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: mem_gb + " GB"
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: preemptible
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
    String docker
    Int cpu
    Int mem_gb
    Int disk_gb
    Int preemptible
  }

  command <<<
    set -euo pipefail
    mkdir -p in out
    ln -s "~{corrected}"       in/corrected.txt.gz
    ln -s "~{corrected_index}" in/corrected.txt.gz.tbi

    /opt/depthsv/scripts/analyze.sh \
      --corrected in/corrected.txt.gz \
      --pheno "~{phenotypes}" \
      --pheno-manifest "~{pheno_manifest}" \
      --pcs "~{pcs}" --ndim ~{ndim} \
      --region "~{region}" \
      --out out \
      --jobs ~{cpu} \
      --threads ~{cpu} \
      -- --minObs ~{min_obs}
  >>>

  output {
    Array[File] sumstats = glob("out/*.txt.gz")
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: mem_gb + " GB"
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: preemptible
  }
}
