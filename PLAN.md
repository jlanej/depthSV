# depthSV — consolidation and migration plan

Working document.

**Status: Phases 0-2 and 4 done.** Phase 3 (the external-engine handoff) is
deferred until collaborator feedback says whether it is wanted; the seam it
would attach to — the corrected matrix — exists either way, so nothing about
the current design forecloses it.

Remaining: the split into a standalone repository.

## Scope

**In:** the analysis path only — join → correct → association summary statistics.

**Out:** downstream summarization, peak curation, result annotation. depthSV
emits summary statistics; what consumes them is a separate concern.

**Deliverable:** per-region association summary statistics (effect, standard error, test statistic,
p-value, N) for a phenotype. Everything after that is someone else's stage.

## Principles

1. **Reproduce the results, not the implementation.** The current numbers are the specification. How
   we get to them is open. Where a restructuring is provably the same estimator, it is plumbing and
   not a methods change — and it must be demonstrated to be the same estimator, not asserted.
2. **Do not over-optimize.** The realistic workflow is one or two phenotypes at a time, checked for
   sanity, then more. Optimizations that only pay off when sweeping every phenotype at once are not
   worth their complexity. Optimizations that make a single phenotype run cheaply and restartably are.
3. **HPC and cloud are both first-class.** SLURM must keep working. The way to get that is to make the
   engine take a region and emit a shard, and keep the dispatcher separate and swappable. No cloud
   concept leaks into the analysis code.
4. **Fail loudly.** Every silent-success path found so far produced a well-formed, indexable, wrong or
   empty file. A job that produces nothing must exit non-zero.

---

## What is done

Phases 0-2 are complete and are now described by the code and the test suite
rather than by this file. In summary: the repository was consolidated into a
single source of truth, five reproduced bugs were fixed, the join moved to
`bgzip`+`tabix` with contig names resolved from the index, a region became the
unit of work, and `tests/smoke_test.sh` asserts on results across 35 checks.

Since then, the join was reworked for scale — parallel extraction, batch sizes
chosen against the open-file limit, scratch bounded near two batches of
columns, and resume from the last finished batch — with output verified
byte-identical to the original implementation. Windowed region lists were made
an exact partition of the matrix (a bin ending on a window boundary previously
landed in two units), passthrough arguments now survive shell quoting on the
way to the R workers, worker diagnostics land in per-unit logs instead of
/dev/null, and BLAS threading is pinned so results do not depend on host core
count.

Both engines were verified against the originals — linear and logistic agree to
~1e-13, and the correction is exact against an independent recomputation.

Two findings from implementation worth keeping:

- **Logistic cannot be made dramatically faster while keeping identical
  results.** Building the design matrix once and swapping only the depth column
  is bit-identical but only ~1.45x: the IRLS iterations dominate, not the
  formula parsing. The large speed-up needs a score test, which is a *different
  test*. Linear keeps its ~33x because the projection is the same estimator.
- **The `#` prefix on the matrix header is load-bearing.** Without it `tabix`
  tries to parse the header as an interval and refuses to index the file.

---

### Phase 3 — external engine handoff (2–3 weeks)

Keep the simple linear and logistic sweep first-class; add a seam after the corrected matrix.

- Dosage exporter writing the corrected matrix to PGEN, with explicit rescaling and the recoverable
  scale constant recorded alongside.
- REGENIE first: it covers linear, logistic and Cox in one tool, reads PGEN natively at ~6×10⁻⁵ dosage
  resolution against 8-bit BGEN's 8×10⁻³, and its step-1 fit is reused across every depth-bin run.
  SAIGE cannot read PGEN, so a SAIGE path would need 16-bit BGEN or VCF with a dosage field.
- Wire the migrated GENESIS path in as the independent cross-check rather than retiring it.

Three traps, worst first:

- **Do not set `--minINFO`.** REGENIE derives an imputation-quality metric from dosage variance, and a
  depth predictor has far lower dosage variance than a genotype, so any threshold risks deleting the
  whole genome. It has no default, so the safe course is to leave it unset — and to check the metric's
  actual distribution on a test chromosome before considering otherwise.
- **Rescale explicitly.** PLINK2's behaviour on out-of-range dosages is undocumented. Do not discover
  it from a flat Manhattan plot.
- **Verify `--minMAC` is inert.** The default of 5 should be harmless when every "variant" has apparent
  frequency near 0.5, but confirm.

Encoding was checked locally: an affine map into [0,2] plus 16-bit quantization preserves the *t*
statistic and p-value and leaves β recoverable by the known constant. A literature search turned up
**no precedent** for writing continuous depth into a dosage container for REGENIE — the biobank CNV
studies found either discretize or use a bespoke regression — so treat the encoding as unvalidated
until a known signal at a well-characterized locus is reproduced through it.

**Gate:** the same phenotype through both paths on the same corrected matrix agrees where it should,
and a known locus is recovered at the expected strength.

### Phase 4 — portability — DONE

- `Dockerfile` pins R 4.4.1 and a dated CRAN snapshot. No `module load` baked in;
  module names now come from `DSV_MODULES` because they are site-specific.
- `workflows/depthsv.wdl` scatters correct-then-analyze over a region list;
  `workflows/slurm_array.sh` does the same as an array job. Both call the same
  stage scripts, so neither is privileged. Nextflow was not added: it is not
  enabled on RAP.
- `scripts/regions.sh` emits the region list both dispatchers consume, read
  from the tabix index so it matches the matrix rather than an assumed karyotype.

**Gate: met.** The suite passes inside the container and results are
byte-identical to a host run. The container caught one real portability
defect — `fread()` cannot read `.gz` without the `R.utils` package, so a gzipped
phenotype table failed on a clean R; both engines now decompress through the
shell.

**Not verified here:** the WDL was checked with `miniwdl check` and both task
bodies were executed in the container against the fixtures, but the workflow was
not run end to end through an engine — the local Docker swarm has an expired TLS
certificate, which is a pre-existing condition on that machine. Run
`miniwdl run workflows/depthsv.wdl` once on a host with a working swarm, or on
the target platform's own engine, before trusting the orchestration.

### Next — split into the standalone repo

- Fresh git history — the current history carries cluster paths, a username and a group name inside
  earlier revisions of the scripts, so a path-filtered split would preserve all of it.

---

## Open decisions

- **Bin size.** Currently 1 kb, giving ~3.1M regions and a multi-terabyte matrix
  at biobank scale. A published UK Biobank depth-association analysis (Garg et
  al., AJHG 2026) used 5 kb, which would cut that five-fold. This shapes
  region-list granularity and storage cost.
- **`pigz` → `bgzip` throughput.** Confirm parity on the cluster before a full run.
  `bgzip -@4` measured 181 MB in 1.23 s locally, but `pigz` was not installed
  here so the direct comparison was not made.

Resolved since: Cox stays first-class in the simple sweep — it is implemented
and covered by the suite, so removing it would now cost more than keeping it.

## Recorded

- On PC correction: `ndim=200` has been validated at scale on real data
  (thousands of bins, 100k+ samples), where copy-number variants are not a
  primary source of variation. Treat that empirical check as the justification;
  it is worth writing into the methods while it is fresh.
