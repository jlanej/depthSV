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
| `plink2` ≥ 2.00a5 reachable (`EX_PREAMBLE_MODULES`) and Dropbox reachable, for the preamble | otherwise pre-stage the genotype files under `preamble/genotypes/` |

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

sbatch preamble.sh           # once: MP-derived ndim + genotype-PC covariates (24 cores / 248 GB / <24 h)
bash run.sh                  # 00 -> 01 -> 02: resolve, prepare, submit
```

The preamble is optional but is what makes the mtDNA-CN models *reasonable*
rather than merely correct: without it the sweep runs unadjusted with
`EX_NDIM=20`; with it, the correction removes the Marchenko-Pastur-determined
number of coverage PCs and the models adjust for sex and genotype PCs. See
[Preamble](#preamble-how-many-pcs-and-which-covariates).

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

## Preamble: how many PCs, and which covariates

[`preamble.sh`](preamble.sh) is a standalone job (`sbatch preamble.sh` from
this directory; `#SBATCH` sizes it at 24 cores / 248 GB / 24 h to match the
upstream NGS-PCA stages, far more than it needs) that settles two decisions
the association models depend on, before the pipeline runs. Both are
idempotent and resumable; `bash preamble.sh --smoke` runs both on a laptop
(ndim from the committed spectrum, genotypes from chr22 only).

**A. How many coverage PCs to remove** — [`R/choose_ndim.R`](R/choose_ndim.R).
NGS-PCA reports the top 200 singular values of the samples × bins matrix.
For pure noise they follow the Marchenko–Pastur law, whose upper edge at
3,202 × 142,070 is a narrow, nearly flat plateau — so the trailing reported
values pin the noise scale down, and the count of singular values above the
edge is the number of components distinguishable from noise. The fit
iterates the noise scale and the count to a fixed point, per mode, and
`preamble/ndim.txt` is the rounded mean across modes so every tree is
corrected identically; `config.sh` picks it up as the `EX_NDIM` default. On
the committed 3,202-sample spectrum the count is **52 at the exact edge and
42 / 38 / 32 with a 1 / 2 / 5 % margin** — real coverage noise is not iid
(per-bin variance differs), so the observed plateau decays a little faster
than MP and the exact-edge count runs high. The default margin is 1 %
(`EX_MP_MARGIN`); that margin is far larger than the Tracy–Widom
fluctuation of the noise edge at this size (Johnstone's scaling puts the
edge's relative sd near 4 × 10⁻⁴, so 1 % is ~25 sd) — it absorbs the
misfit between real, heteroskedastic coverage noise and the iid MP model,
not sampling noise, and the count also depends on whether observed rank *j*
is matched to noise rank *j* or *j − k* (42 vs 48 at 1 % here). Treat the
1–2 % counts as a defensible range, not a number; `ndim/ndim_by_mode.tsv`
carries all four counts and `ndim/ndim_mp.png` the fit, so the sensitivity
is in the open rather than in a default.

**B. Genotype PCs** — the NYGC 30x GRCh38 callset for exactly these 3,202
samples is published in PLINK 2 format on the
[PLINK 2.0 resources page](https://www.cog-genomics.org/plink/2.0/resources)
(`resources/genotype_sources.tsv`: ~3.4 GB of per-chromosome `.pgen.zst`
plus `.pvar.zst` and the corrected `.psam` with pedigree, sex and
population — no VCFs, no CRAMs). The recipe is the textbook one, per
chromosome to bound disk: biallelic A/C/G/T SNPs, MAF ≥ 0.05, missingness
≤ 1 %, the Anderson et al. 2010 long-range LD regions excluded
(`resources/long_range_ld_grch38.txt`), `--indep-pairwise 1000kb 1 0.1`;
merge; KING-robust `--king-cutoff 0.0884` (2nd degree and closer removed);
`--pca 40 allele-wts` on the unrelated set; projection of every sample onto
the allele weights with the unrelated set's frequencies; then
[`R/genotype_pcs.R`](R/genotype_pcs.R) calibrates the projection against
the in-sample PCs (proportional per PC; the r² is the check), applies the
same MP-edge count to the genotype eigenvalues, and plots PC1/2, PC3/4 by
superpopulation with relatives as open points, plus the scree. Output:
`preamble/covariates.tsv` (`SAMPLE`, `GPC1..GPC40`, `GPC_PROJECTED`,
superpopulation/population), `gpc_plots.png`, `gpc_calibration.tsv`.

**How it flows into the run.** With `covariates.tsv` present, the prepare
stage merges the genotype PCs into `phenotypes.tsv` and writes the adjusted
manifest: `mtdna_cn` (unadjusted, the pure truth check) plus
`mtdna_cn_adj`, `log2_mtdna_cn_adj`, `mtdna_cn_null_adj` with
`+SEX+GPC1..GPC10` (`EX_N_GPCS`; or set `EX_COVARIATES` to any `+`-joined
list of phenotype columns, `none` for unadjusted), `sex_linear` with the
genotype PCs only, and the logistic run unchanged. The evaluation applies
each family's checks to its adjusted variant.

Requirements beyond the pipeline's: `plink2` ≥ 2.00a5 (for `--pmerge-list`
and `--pca allele-wts`; `EX_PREAMBLE_MODULES` names the module if your site
has one) and outbound HTTPS to Dropbox from the node that runs it — if that
is blocked, download `resources/genotype_sources.tsv`'s files elsewhere and
place them as `preamble/genotypes/chr<N>.pgen.zst` / `.pvar.zst` /
`samples.psam`; the script picks up whatever is already there.

## Stages

`run.sh` runs stages 0, 1 and 2 in order (2 triggers 3–5 itself); the
stages exist separately so any one can be rerun alone. The preamble sits
before all of them and runs once.

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
  stay calibrated (λ band 0.85–1.20 on real data), judged on autosomal
  bins: every chrX/Y bin carries the same sex vector, so those tests are one
  dependent draw, not thousands of independent ones. The all-bin λ is
  reported as INFO.
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
| `EX_NDIM` | `preamble/ndim.txt`, else 20 | PCs removed by the correction; the preamble's MP count when it ran, an explicit value always wins |
| `EX_MP_MARGIN` | 0.01 | relative margin above the MP edge for the PC counts (coverage and genotype) |
| `EX_N_GPCS` / `EX_COVARIATES` | 10 / `SEX+GPC1..GPC10` when covariates exist, else `none` | covariate terms of the adjusted models |
| `EX_GENO_CHROMS` | 1–22 (22 only in smoke) | chromosomes the genotype PCA uses |
| `EX_PREAMBLE_MODULES` | `plink2` | modules loaded before plink2, where `module` exists |
| `EX_WINDOW` | 10000000 | work-unit size in bp (~310 units over chr1–22,X,Y,M, ~150 s each at 3,202 samples); 0 = per contig |
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
  ~67% of log₂(MTDNA_CN) variance, and on the full 3,202 samples the first
  42 explain ~58%. What the correction removes is *power*, not the effect
  estimate: by Frisch–Waugh the chrM slope is invariant to ndim (it stays
  near 1 at every setting), while the t-statistic falls as the PCs absorb
  the depth variance (35 → 4 across ndim 0 → 40 in the smoke slice). That
  is not a defect of the example, it is the ndim trade-off made visible:
  sweep `EX_NDIM` (0 disables correction) and watch the chrM t-statistic
  and the null λ, not the slope. This is exactly the evidence the top-level
  README says the `--ndim` choice should rest on.
- **Regenerating upstream results** is NGS-PCA's business; this example
  only reads its outputs. Keep one upstream mosdepth version per
  comparison (their timing records note it).
- **External engines.** The per-mode corrected matrices under
  `work/<mode>/corrected/` are the seam for the REGENIE/GENESIS handoff
  described in [`PLAN.md`](../../PLAN.md); this example gives that gate
  its real-data substrate (a known chrM signal at known strength) when it
  lands.

## Credits

This example, its preamble and its evaluation were developed by
[Claude](https://claude.ai) (Anthropic; Claude Fable 5) with the maintainer
in Claude Code; the upstream data and PCA come from the NGS-PCA project and
the 1000 Genomes / NYGC 30x resource (Byrska-Bishop et al. 2022), and the
genotype callset in PLINK 2 format from the PLINK 2.0 resources page.

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
| preamble: `got an HTML page instead of a file` | Dropbox is blocked from that node; download the files listed in `resources/genotype_sources.tsv` elsewhere and place them under `preamble/genotypes/` |
| preamble: `projection reproduces in-sample PCs with r^2 down to …` | the projected scores should be proportional to the in-sample PCs; a low r² on a deep PC is a near-degenerate eigenvalue pair, harmless if it is beyond `EX_N_GPCS` |
| `ndim.txt not written` | the reported spectrum never reached the noise bulk; raise `NUM_PC` upstream or set `EX_NDIM` |
