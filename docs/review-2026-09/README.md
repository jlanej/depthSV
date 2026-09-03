# Reproduction scripts for REVIEW.md (September 2026)

Scripts the review's findings were established with, grouped by the pass that
wrote them. They were run against tree `fe768ec` on macOS (R 4.5, plink2
2.00a5, htslib 1.22, GNU parallel 2025) and write to whatever scratch
directory they are pointed at; several fetch the committed NGS-PCA 1000G
outputs from GitHub and the plink2-format chr22 genotypes from the PLINK 2.0
resources page.

- `fwl_check.R` — the independent simulation behind REVIEW §1.1 (PCs absent
  from the association model: deflation and covariate-induced bias).
- `mp_proto.R`, `lam_by_class.R` — the Marchenko-Pastur prototype and the
  per-chromosome-class lambda diagnostic.
- `stats/` — simulations on the real 3,202-sample PCs and phenotype: adjusted
  models, sex chromosomes, artefact size, Cox, logistic LRT, effective tests,
  relatedness, MP heteroskedasticity, projection shrinkage; `out_*.txt` hold
  their outputs.
- `example/` — the example/preamble pass: stale-shard, join-growth, fallback,
  glob, MP edge-case and independent-MP checks.
- `shell/` — the robustness pass: bash-3.2, spaces, resume, duplicate-ID,
  SLURM-shim and lock experiments.
- `cloud/` — module-order and chunk-size experiments, and
  `depthsv_cloud_sketch.wdl`, a miniwdl-checked sketch of the block-join
  workflow recommended in REVIEW §4.2.
- `perf/` — the throughput and fixed-cost benchmarks behind REVIEW §3.
