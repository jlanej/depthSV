# 1000 Genomes 30x — end-to-end depthSV example

A real-data, end-to-end exercise of the whole pipeline on the 3,202
[1000 Genomes 30x high-coverage samples](https://www.internationalgenome.org/data-portal/data-collection/30x-grch38)
(NYGC, GRCh38), built on the output of the
[NGS-PCA `1000G_highcov` example](https://github.com/jlanej/NGS-PCA/tree/master/example/1000G_highcov).
It does three jobs at once:

1. **End-to-end test with built-in truth.** The upstream QC table defines
   phenotypes whose genomic architecture is known by construction, so the
   run can be *asserted on*, not just completed (see
   [What gets checked](#what-gets-checked)).
2. **Upstream fast-mode comparison, carried to the end.** The NGS-PCA
   example runs mosdepth twice per sample — standard and `--fast-mode` —
   and its own report compares the PCs, judged against a *seed-control*
   rerun of the PCA (same data, different `-randomSeed`) so the estimator's
   own noise is the yardstick. This example runs depthSV on both trees, and
   on the standard tree with the seed-control PCs, and compares the
   *association statistics* the same way — the number that actually matters
   downstream, with the same yardstick.
3. **Profiling.** Every stage invocation is timed (wall clock, peak RSS
   where recordable, plus `sacct` on SLURM), and the profile report says
   where the time went and which work units straggled.

## The phenotypes are the test

The NGS-PCA example's QC stage writes `sample_qc.tsv` per mode, and three of
its columns make real-data truth checks possible without any external
callset:

| analysis | phenotype | expectation |
|---|---|---|
| `mtdna_cn` (linear) | `MTDNA_CN` = 2 × chrM mean / HQ autosomal median | chrM bins must dominate the association; top autosomal hits are NUMT candidates, reported for inspection |
| `log2_mtdna_cn` (linear) | log₂ of the above | on a chrM bin the corrected depth *is* log₂(chrM/median), so the slope should sit near 1 — an effect-size check, not just a rank check |
| `sex_linear` (linear) | `SEX` from X/Y coverage ratios | chrX/chrY bins must dominate |
| `inferred_sex` (logistic) | same `SEX`, through the logistic engine | runs to completion; *documents* that the Wald z collapses on the sex chromosomes under complete separation (Hauck–Donner) — which is why the rank assertion lives on the linear run |
| `mtdna_cn_null` (linear) | `MTDNA_CN` permuted with a fixed seed | genomic-control λ near 1; ~5% of regions at p<0.05 |

Because phenotype, coverage medians and PCs all come from the *same
upstream run per mode*, the comparison between modes is a genuine
end-to-end one: everything downstream of the CRAMs differs only by
mosdepth's mode.

## Prerequisites

- The depthSV [requirements](../../README.md#requirements) (R with
  `optparse`/`data.table`, htslib, GNU parallel).
- A run of the NGS-PCA example with the fast-mode comparison on. Its
  README's "mosdepth fast-mode comparison" section is the reference; the
  shape, with the parts this example depends on marked:

  ```bash
  # in NGS-PCA/example/1000G_highcov, per its README
  export COMPARE_FAST_MODE=1                       # both trees (the current default)
  bash 00_setup.sh && bash 01_download_and_mosdepth.sh
  sbatch 02_run_ngspca.sh                          # standard PCs + autosomal.median.txt
  MOSDEPTH_DIR=$WORK_DIR/mosdepth_output_fast NGSPCA_OUTPUT=$WORK_DIR/ngspca_output_fast sbatch 02_run_ngspca.sh
  NGSPCA_OUTPUT=$WORK_DIR/ngspca_output_seed43 RANDOM_SEED=43 sbatch 02_run_ngspca.sh   # step 2b: seed control
  # QC AFTER each mode's PCA run, so HQ_MEDIAN_COV (hence MTDNA_CN) comes
  # from that run's autosomal.median.txt; the fast QC pass must name the
  # fast run's NGSPCA_OUTPUT or it would take the standard medians
  sbatch 03a_mosdepth_coverage_summary.sh && bash 03_collect_qc.sh
  MOSDEPTH_DIR=$WORK_DIR/mosdepth_output_fast QC_OUTPUT=$WORK_DIR/qc_output_fast \
    NGSPCA_OUTPUT=$WORK_DIR/ngspca_output_fast sbatch 03a_mosdepth_coverage_summary.sh
  MOSDEPTH_DIR=$WORK_DIR/mosdepth_output_fast QC_OUTPUT=$WORK_DIR/qc_output_fast bash 03_collect_qc.sh
  ```

  The seed control is optional: without it the `seedctl` mode is skipped
  and the fast comparison is reported without its yardstick. A single-mode
  run works too: `EX_MODES=standard`.

### Upstream checklist

What stages 0 and 1 verify — `bash run.sh --prepare-only` runs exactly
these and stops — and what to look at if they complain:

| Check | Why it matters |
|---|---|
| NGS-PCA checkout from 2026-08-28 or later (commit `7862645`) | earlier `03a` derived `HQ_MEDIAN_COV` with bedtools, or left it `NA`; now it is sourced from the PCA run's own `autosomal.median.txt` |
| each mode's `03a` ran *after* its `02`, and the fast `03a` with `NGSPCA_OUTPUT=…_fast` | otherwise `MTDNA_CN` is built on a different denominator than the one this pipeline normalises against — prepare compares the two and warns |
| `MTDNA_CN` is not `NA` across the board | prepare stops if it is (the phenotype does not exist) |
| both trees complete: `ls mosdepth_output/*.regions.bed.gz \| wc -l` equals the fast tree's count | samples processed before `COMPARE_FAST_MODE` was on are re-queued upstream; prepare reports the sample-set difference in `inputs/cross_mode_samples.txt` |
| all three PCA runs from the same container image | the report should reflect mosdepth's mode and the seed, nothing else |

No mosdepth run yet? See [Smoke mode](#smoke-mode-no-cluster-no-crams) —
the committed NGS-PCA results plus a simulated depth tree run the whole
example anywhere.

## Quick start (HPC)

```bash
cd example/1000G_highcov

# where the NGS-PCA example ran, and where this example may write
export NGSPCA_WORK_DIR=/scratch/$USER/1000G_highcov
export EX_WORK_DIR=/scratch/$USER/depthsv_1000G_highcov
# site specifics, if any: export EX_SBATCH_EXTRA="--partition=... --account=..."
# module names, if needed: export DSV_MODULES="parallel R htslib"

bash run.sh                  # 00 -> 01 -> 02: resolve, prepare, submit
```

While the upstream comparison is still finishing, `bash run.sh --prepare-only`
runs just the resolution and the checks (the [upstream checklist](#upstream-checklist))
and stops before submitting anything.

`run.sh` strings the numbered stages together; each can be run on its own.
`02_run_depthsv.sh` submits, per prepared mode: a **join** job, then a
**dispatch** job that derives the windowed region list from the joined
matrix and submits the **correct+analyze array** over it, then an
**evaluate** job. The `seedctl` mode shares the standard matrix, so it
submits no join and its dispatch waits on the standard one. Whichever
mode's evaluation finishes last also runs the cross-mode comparison and the
profile report, so one submission carries the example to its end state:

```
eval/<mode>/summary.md                per-mode truth checks (standard, fast, seedctl)
compare/summary.md                    the calibration verdict: fast vs standard, judged against seedctl
compare/standard_vs_fast/summary.md   per-analysis concordance detail (and standard_vs_seedctl/)
profile/profile_report.md             where the time went
```

Monitor with `squeue --me | grep dsvx-`. Resubmitting `02_run_depthsv.sh`
after a partial failure is safe and cheap: completed units are skipped via
the stage scripts' `.done` markers (pass `--force` to redo). If your site
does not allow `sbatch` from a compute node, run the pieces by hand in the
same order: `--stage join-exec`, then `--stage dispatch`, per mode — or use
`--runner local` on a whole node.

## Smoke mode (no cluster, no CRAMs)

The NGS-PCA repository commits its full 3,202-sample results (PCs and the
QC table — mosdepth trees are too large to commit). Smoke mode fetches
those and **simulates** a small mosdepth tree per mode that is numerically
consistent with each sample's real QC row (HQ median, X/Y ratios, chrM
ratio) over chr20/chrX/chrY/chrM slices — so every truth check still
fires, on any machine:

```bash
cd example/1000G_highcov
export EX_WORK_DIR=/tmp/depthsv-example-smoke

bash run.sh --smoke --runner local
```

Minutes end to end. The simulated depths exercise the machinery — treat the
numbers as plumbing, not results. If the fast-mode results are not yet
committed upstream, the fast tree is simulated from the standard tables
with a small deterministic perturbation, and labelled as such in
`inputs/fast/paths.env`. The seed control is never fetched — it is skipped
in smoke mode unless `EX_NGSPCA_DIR_SEEDCTL` points at a local PC table.

## Stages

`run.sh` runs stages 0, 1 and 2 in order (2 triggers 3–5 itself); the
stages exist separately so any one can be rerun alone.

| Stage | Script | What it does |
|---|---|---|
| 0 | `00_fetch_inputs.sh` | Resolve each mode's `svd.pcs.txt`, `sample_qc.tsv` and (when the run wrote one) `autosomal.median.txt`: local NGS-PCA trees first, the committed GitHub results as fallback; the seed control only locally. `EX_SMOKE=1` also simulates the mosdepth trees. Records the resolution in `inputs/<mode>/paths.env`. |
| 1 | `01_prepare_inputs.sh` | Build depthSV's input tables per mode: PC table with the `.by1000.` sample suffix stripped, `SAMPLE`/`AUTO_HQ_median` coverage from NGS-PCA's own median table (else the QC table), the phenotype table and analysis manifest above, the mosdepth manifest, and `chrom.sizes` read from the first region file. Verifies the mosdepth↔coverage ID overlap *before* hours of joining, that `MTDNA_CN` was built on the same median, and that the modes' sample sets agree. |
| 2 | `02_run_depthsv.sh` | Per prepared mode: `join` → windowed region list (`scripts/regions.sh`, filtered to `EX_CONTIG_REGEX`) → `correct` + `analyze` per region, as a SLURM chain or locally. `seedctl` reuses the standard matrix. Every stage runs under the timing recorder. |
| 3 | `03_evaluate.sh` | Truth checks per mode → `eval/<mode>/`. FAIL = machinery broken (non-zero exit); WARN = statistical expectation missed. |
| 4 | `04_compare_modes.sh` | Association concordance for standard-vs-fast and standard-vs-seedctl → `compare/<a>_vs_<b>/`, and the calibration verdict → `compare/summary.md`. |
| 5 | `05_profile.sh` | Timing + `sacct` aggregation → `profile/`. |

Stages 3–5 run automatically at the end of stage 2; they exist separately
so they can be rerun (or rerun with different thresholds) without touching
the pipeline output.

## What gets checked

`eval/<mode>/checks.tsv` — the intent behind each check:

- **`chrM_top_hit` (FAIL)** — the strongest `mtdna_cn` association must be
  a chrM bin. The phenotype is *defined* from chrM depth; if this fails,
  sample alignment, correction or the join is broken.
- **`chrM_log2_slope` (WARN)** — the `log2_mtdna_cn` estimate at the best
  chrM bin should sit near 1 (band 0.5–1.5). Catches scale errors that a
  rank check cannot.
- **`sex_top_hit` (FAIL) / `sex_top100_purity` (WARN)** — `sex_linear`
  must put chrX/chrY on top, and ≥90% of its top 100.
- **`wald_separation_note` (INFO)** — on the logistic run, small
  sex-chromosome z alongside a dominant linear t is the Hauck–Donner
  collapse under complete separation, recorded so nobody mistakes it for a
  miss. (It is also a live argument for the score/external-engine path on
  saturated binary traits.)
- **`null_lambda` / `null_frac_p05` (WARN)** — the permuted phenotype must
  stay calibrated (λ band 0.85–1.20 on real data).
- **`regions_unique` (FAIL)** — no bin tested twice: the windowed region
  list must partition the matrix.
- **`all_units_reported` (WARN)** — one output shard per work unit per
  analysis; fewer means lost array tasks.
- **`top_autosomal_hits` (INFO)** — for the mtDNA phenotypes these are
  NUMT / mito-correlated candidates. On real data, expect some.

## The mode comparison

Each pair's `compare/<a>_vs_<b>/summary.md` reports, per analysis, over
regions tested in both modes: Pearson r of estimates and of test
statistics, Spearman ρ of |stat|, sign agreement where |stat|>2, top-K
overlap, and the median relative difference of estimates — unsigned, and
signed b-vs-a (the same convention as NGS-PCA's QC concordance table, which
exposes a uniform bias even at r≈1). Everything is also broken out by
chrM / sex / autosome, because chrM — where fast mode skips the
mate-overlap correction that matters most at extreme depth — is exactly
where a divergence would hide inside a genome-wide average.

`compare/summary.md` then states the **calibration verdict**: per
analysis, the fast-vs-standard distance (1 − r(stat)) beside the
seedctl-vs-standard distance, and whether fast mode moved the statistics
more than a reseed of the PCA moved them (`EX_CALIBRATION_FACTOR` sets how
much more counts). That is the same judgement NGS-PCA's report makes about
the PCs, carried through to the association statistics. Without a seed
control the fast comparison stands alone and the verdict says so.

Upstream, the PCs were near-identical between modes; the same
normalisation argument (log₂ against each sample's own median) applies
here, so the expectation is r(stat) ≥ 0.99 per analysis and a fast
distance comparable to the seed distance. Falling short is a finding about
fast mode, not about the pipeline — that is the point of running the
comparison end to end.

## Profiling

Every stage invocation appends one record (mode, stage, unit, wall
seconds, peak RSS where GNU time exists, exit code, host, job) to
`profile/timings.d/`; `05_profile.sh` folds them into `timings.tsv`, pulls
`sacct` for every submitted job (per-array-task Elapsed/MaxRSS), and
writes `profile_report.md`: time by stage, per-unit distributions
(median/p95/max), the slowest 15 units, and the highest-memory steps.
Expected shape at 3,202 samples: the join and the logistic sweep dominate;
correct is I/O-bound; a straggling unit usually means a busy node or an
unusually dense window.

## Configuration

All knobs live in [`config.sh`](config.sh) as `EX_*` variables with
`${VAR:-default}` semantics — export before running to override. The ones
that matter most:

| Variable | Default | Meaning |
|---|---|---|
| `NGSPCA_WORK_DIR` | `/scratch/$USER/1000G_highcov` | the NGS-PCA example's `WORK_DIR` (per-mode trees derive from it) |
| `EX_WORK_DIR` | `/scratch/$USER/depthsv_1000G_highcov` | everything this example writes |
| `EX_MODES` | `standard fast seedctl` | in order; a mode without upstream inputs is skipped, so drop or keep freely |
| `EX_SEED_CONTROL_SEED` / `EX_NGSPCA_DIR_SEEDCTL` | 43 / `$NGSPCA_WORK_DIR/ngspca_output_seed43` | the seed-control PCA run (NGS-PCA step 2b) |
| `EX_MEDIAN_SOURCE` | `auto` | `auto` prefers NGS-PCA's `autosomal.median.txt`, `qc` uses `HQ_MEDIAN_COV`, `ngspca` insists on the table |
| `EX_CALIBRATION_FACTOR` | 1.5 | fast distance ≤ this × seed distance counts as "within seed noise" |
| `EX_NDIM` | 20 | PCs removed by the correction; sweep it against the per-region stats the correct stage writes |
| `EX_WINDOW` | 25000000 | work-unit size in bp (~140 units over chr1–22,X,Y,M); 0 = per contig |
| `EX_CONTIG_REGEX` | primary + chrM | which contigs get corrected/analysed (the matrix keeps everything) |
| `EX_MIN_OBS` | 100 | per-region completeness floor at the analysis stage |
| `EX_SBATCH_JOIN/UNIT/LIGHT` | see config | resources per job class |
| `EX_SBATCH_EXTRA` | empty | partition/account/QoS — anything site-specific |
| `EX_ARRAY_THROTTLE` | 100 | max concurrent array tasks |
| `EX_SAMPLE_SUFFIX` | `.by1000.` | suffix stripped from upstream PC sample IDs |

The `DSV_*` environment (see `conf/example.env`) is composed per mode by
the driver; you should not need to set it yourself here.

## Notes

- **IDs.** Matrix columns come from the mosdepth filenames
  (`HG00096.by1000.regions.bed.gz` → `HG00096`); the prepare stage strips
  the same `.by1000.` suffix from the PC table, and the QC table is
  already bare. Everything joins on bare IDs, so no `--sampleIdPattern` is
  needed. If you changed `MOSDEPTH_BIN_SIZE` upstream, set
  `EX_SAMPLE_SUFFIX` to match — and note the pipeline's own suffix
  stripping only knows `.by1000`, so other bin sizes additionally need
  `--sampleIdPattern` passed through to the R drivers.
- **Medians.** NGS-PCA's `autosomal.median.txt` (`SAMPLE`, `AUTO_HQ_median`,
  `N_BINS`, IDs with the same suffix as the PCs) is the median over exactly
  the bins PCA used — the definitional denominator, and depthSV's native
  coverage format — so prepare uses it directly whenever the run wrote one.
  The QC table's `HQ_MEDIAN_COV` should be the same number, and `MTDNA_CN`
  should equal 2 × chrM-mean / that median; prepare checks both and warns
  when they disagree, which is the signature of a QC pass run before its
  PCA run or against the other mode's. A run that reused a cached NGS-PCA
  matrix writes no median table; the QC fallback covers it.
- **Reruns resume.** Every unit writes a `.done` marker; resubmitting
  `02_run_depthsv.sh` redoes only what is missing. `--force` redoes a
  mode's units deliberately. Association output lands in
  `assoc_ndim<EX_NDIM>/`, so changing `EX_NDIM` can never silently reuse a
  previous sweep's results.
- **The PC correction absorbs globally structured phenotypes — measurably.**
  mtDNA copy number loads on the leading depth PCs (it is one of NGS-PCA's
  own QC overlays); across a 64-sample subset the first 20 PCs explain
  ~67% of log₂(MTDNA_CN) variance, and the `chrM_log2_slope` shrinks
  accordingly. That is not a defect of the example, it is the ndim
  trade-off made visible: sweep `EX_NDIM` (0 disables correction) and watch
  the slope and the null λ move. This is exactly the evidence the
  top-level README says the `--ndim` choice should rest on.
- **Regenerating upstream results** is NGS-PCA's business; this example
  only reads its outputs. Keep one upstream mosdepth version per
  comparison (their timing records note it).
- **External engines.** The per-mode corrected matrices under
  `work/<mode>/corrected/` are the seam for the REGENIE/GENESIS handoff
  described in [`PLAN.md`](../../PLAN.md); this example gives that gate
  its real-data substrate (a known chrM signal at known strength) when it
  lands.

## Troubleshooting

| Problem | Likely cause / fix |
|---|---|
| `no mosdepth sample matches the coverage table` | `EX_SAMPLE_SUFFIX` does not match the upstream naming; check one filename against the QC table's `SAMPLE_ID` |
| join dies with `open files` | raise `ulimit -n` as the error says, or pass a smaller `--batch-size` via a manual join |
| dispatch job fails with `sbatch: command not found` | your site forbids submission from compute nodes; run `--stage dispatch --mode <m>` from a login node after the join finishes |
| evaluation FAILs `chrM_top_hit` | inspect `work/<mode>/corrected/*.log` first — the `[align]` drop counts show a sample-ID mismatch immediately |
| `all_units_reported` WARN | compare `squeue`/`sacct` for the array; resubmit `02_run_depthsv.sh` — completed units are skipped |
| GitHub fetch fails | pin `EX_GITHUB_REF` to a tag/commit, or point `EX_NGSPCA_DIR_*`/`EX_QC_DIR_*` at local copies |
| prepare warns `HQ_MEDIAN_COV ... matches autosomal.median.txt for only N%` | the upstream `03a` ran before `02`, or the fast `03a` without `NGSPCA_OUTPUT=…_fast`; rerun `03a` then `03` for that mode |
| prepare stops: `MTDNA_CN is NA for every sample` | same cause, one step worse — no `HQ_MEDIAN_COV` at all; the phenotype does not exist until `03a`/`03` are rerun after `02` |
| `SKIP seedctl` | no seed-control PCA run yet; produce it with NGS-PCA's step 2b, then rerun `00` → `02` (only `seedctl` will run; the standard matrix is reused) |
| `seedctl shares the standard matrix: run --mode standard first` | `--mode seedctl` alone needs a completed standard join; use `--mode all`, or run standard first |
