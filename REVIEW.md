# depthSV — review, September 2026

Tree reviewed: `main` at `fe768ec` (the merged 1000G example and preamble), on
2026-09-02. This supersedes [`REVIEW-2026-08.md`](REVIEW-2026-08.md) (tree
`7a1e385`), whose findings are carried forward below with their current status.

**How it was done.** Five independent specialist passes — statistical validity,
biobank-scale performance, shell/R robustness, cloud/HPC deployability, and an
adversarial pass over the example and preamble — each instructed to reproduce
rather than assert, followed by my own re-verification of every finding that
carries a *blocker* or *major* label (the reproduction scripts are under
[`docs/review-2026-09/`](docs/review-2026-09/)). The core pipeline
(`scripts/`, `lib/`, `R/`, `workflows/`) has not changed since the August
review, so its findings were re-confirmed at current line numbers rather than
re-derived. Where a finding below says *verified*, I ran it; where it says
*reviewer-reproduced*, a specialist pass ran it and I checked the mechanism.

**Verdict.** The engineering is careful and the design — regions as the unit of
work, atomic outputs, one implementation behind every dispatcher — is right.
But the pipeline is not publication-grade yet, for four reasons that are each
fixable and none of which is cosmetic:

1. **The association test is mis-specified.** Depth is residualised against the
   coverage PCs, but the PCs are absent from the association model. That makes
   every test conservative by a factor of about 1 − R²(phenotype ~ PCs) — with
   the real MTDNA_CN and the preamble's 42 PCs, tests are deflated to λ ≈ 0.44 —
   and, once covariates that correlate with a PC are in the model (sex, genotype
   PCs, ancestry), it *biases* every bin that correlates with those covariates:
   ancestry-stratified deletions come out at λ ≈ 3.9 and every chrX/chrY bin at
   λ ≈ 20 in the example's adjusted models. The August review recorded the
   opposite conclusion; it was wrong (§1.1). The fix is about ten lines.
2. **A deterministic blocker in `correct.sh`** kills any work unit of ten or more
   chunks, which at 3,202 samples is every unit of the example as merged (§2.1).
   A one-line fix, plus a stopgap already applied to the example's config.
3. **The silent-success paths from August are all still open** — logistic level
   ordering, the `--minDepth` floor, parameter-blind `.done` markers, the
   swallowed worker exit status, the empty-shard commit — and the example
   inherits the marker problem in a worse form (§2, §5).
4. **It does not deploy on either cloud as shipped**: the container image is not
   published and the WDL defaults to a local-only tag; every scatter shard
   localises the whole multi-terabyte matrix; the join has no cloud shape at all
   (§4). On SLURM the documented `sbatch workflows/slurm_array.sh` cannot find
   its own scripts (§4.1).
5. **At biobank scale the design spends its time on overhead.** Measured, not
   estimated: at 500k samples the join is a 60–70 h single-node job with 7 TB of
   scratch, and correct/analyze spend > 90% of their CPU on process restarts and
   ASCII parsing (≈ 5,000 CPU-h for the example's six-phenotype sweep, 10.5 GB
   per worker at ndim 200). At 3,202 samples none of this bites (§3).

The suggested order at the end turns these into five pull requests.

---

## 1. Scientific validity

### 1.1 Coverage PCs must enter the association model — blocker (verified)

`R/correct.R:141-170` residualises each bin against `[1, PC1..PCk]`;
`R/analyze.R:175-180,238-250` builds the design `Z` from the model formula
only, so the PCs never reach the test. The August review (§7) argued this is
harmless by Frisch–Waugh. It is not: FWL requires the *phenotype* to be
residualised against the same basis, or the basis to be in the model.

Independent simulation (`docs/review-2026-09/fwl_check.R`, n = 3,000, 10
orthonormal PCs built exactly as `correct.R` builds them, phenotype 50%
explained by the PCs):

| bins | pipeline as written | same, PCs added to the model |
|---|---|---|
| technical null (x ⟂ everything) | λ = **0.513**, P(p<.05) = 0.006 | λ = 1.03 |
| null bins correlated with a covariate z that correlates with PC2 | λ = **1.81**, mean t = +0.90, 90% same sign | λ = 0.994, mean t = 0.01 |

The reviewer's run on the **real** 3,202-sample PCs and MTDNA_CN (`ndim = 42`,
covariates SEX + population; `docs/review-2026-09/stats/sim_adjusted_real.R`,
`sim_sexchrom.R`) gives the same picture at scale: technical-null λ 0.60
adjusted / 0.44 unadjusted (full model 0.99); an ancestry-stratified null
deletion λ = 3.86 with 96% of effects sharing a sign; chrY null bins under the
example's adjusted model **λ = 22.6**, chrX 20.7 (0.96 / 0.99 with the PCs in
the model). Mechanism: `x_res = x − P_PC x`; after the covariate projection,
`M_Z x_res` is no longer orthogonal to the PC space, so the omitted PC-part of
the phenotype (56% of MTDNA_CN) leaks into β̂ wherever depth correlates with a
covariate — and the sex vector, which the autosomal PCs capture at R² = 0.026,
is such a covariate on every sex-chromosome bin.

**Fix.** Give `analyze.R` `--pcs FILE --ndim N` and append `PC1..PCN` to the
covariate design whenever `ndim > 0` (equivalently, residualise the phenotype
against the same basis inside `correct.R`'s projection). Then
`M_[Z,PC] x_res = M_[Z,PC] x` and FWL holds exactly; the simulations above
reproduce λ = 0.99 and no bias. Re-run the example afterwards: its adjusted
sweep as merged would call the whole of chrX and chrY.

### 1.2 Sex chromosomes need a ploidy model — major (reviewer-reproduced)

`R/correct.R:208` normalises every contig by the autosomal median with no
ploidy branch and no PAR handling (August §2.3, still open). With real sexes
and X/Y ratios: female chrY depth is 5.7% of autosomal (log2 ≈ −4.1, a
mismapping signal, not zero), and the female half carries 82–98% of a chrY
bin's within-sex sum of squares — so with `+SEX` the chrY test is a test of
female mismapping noise. An X-linked deletion coded in log2 gives
β̂ = −3.9 per unit in females and −1.9 in males (hemizygous loss sits at the
−12.7 floor, heterozygous loss at −1). `sex_linear`/`sex_top100_purity` in the
example cannot see any of this: SEX is inferred from the same X/Y depth.
**Fix:** expected copies by sex per contig (X non-PAR 2/1, Y 0/1, PAR
autosomal), emit copies-lost or copy ratio relative to expected ploidy, test
chrY in males only, stratify or fit a sex × depth interaction on chrX.

### 1.3 The log2 + absolute-floor coding lets homozygous deletions dominate common-deletion tests — major (reviewer-reproduced; August §1.2 generalised)

`R/correct.R:39,208`. At a deletion with allele frequency 0.3 (9% CN0 at
−6…−12.7), the CN0 samples carry a median 96.5% of the depth sum of squares;
log2 coding correlates −0.78 with true copy number against −0.99 for the copy
ratio `2·2^x`. With an additive per-copy effect on the real MTDNA_CN: log2
coding power at 5e-8 = **0.17**, winsorised at −3 = 0.69, copy-ratio coding =
**0.74** (the oracle). A rare homozygous deletion carried by one to three
people is exactly August §1.2's single-outlier case (P(p<5e-8) = 1.5e-4 per
bin, 3,000× nominal, with the raw phenotype). **Fix:** test on the copy-ratio
scale (or winsorise log2 near −3); make `--minDepth` relative to the sample's
median; emit each bin's largest single-sample share of the sum of squares and
flag or drop bins above 50% (a MAC-like filter); rank-INT quantitative
phenotypes by default for the tail problem — noting INT does not touch §1.1
or §1.2.

### 1.4 The Cox path is miscalibrated by one floored sample and has never been checked numerically — major (reviewer-reproduced)

`R/analyze.R:268-276` wraps `coxph` in `try` and ignores warnings and
convergence. n = 3,199, 30% events, null: one floored sample per bin gives
P(p<.05) = **0.13**, P(p<1e-4) = **0.010** (100× nominal), min p = 0 over
2,000 bins, no warning; three floored carriers with the three earliest events
give z = −6.2, p = 0. The test suite asserts only a row count for coxph
(`tests/smoke_test.sh:150-151`). **Fix:** the §1.3 leverage handling before
the fit; sandwich or Firth Cox for sparse-carrier bins; a convergence column;
one numeric assertion in the suite.

### 1.5 No effective number of tests, no threshold — major (reviewer-reproduced)

Nothing in `scripts/`, `lib/` or `R/` reports λ, M_eff or a threshold (August
§2.5). At 3.1M independent bins, 5e-8 is a 14% family-wise error rate
(Bonferroni 1.6e-8); with neighbour correlation ρ = 0.3 / 0.6 / 0.9 / 0.97 the
effective count is 97 / 75 / 54 / 29% of M. Because the linear path is a
projection, permutation maxima are one matrix product: 400 permutations ×
30,000 bins × 2,000 samples ran in 0.3 s. **Fix:** per shard, emit max|t| over
K permuted (or structured-null, §1.6) phenotypes and the along-genome
autocorrelation of corrected depth; a run-level step turns shard maxima into an
empirical genome-wide threshold and M_eff; report Bonferroni beside it and
count the phenotype sweep in it.

### 1.6 Relatedness: λ ≈ 1 + 0.19 h² at CNV bins, invisible to genome-wide λ and to the permuted null — major (reviewer-reproduced)

613 of 3,199 example samples are related (563 trio children). A 1000G-like
pedigree simulation gives λ = 1.05 / 1.09 / 1.13 / 1.20 at CNV bins for
h² = 0.25 / 0.5 / 0.8 / 1.0 while technical bins stay at 1.00 — so a
genome-wide λ (99% technical bins) cannot show it, and the example's permuted
null has λ = 0.998 at the same bins: permutation destroys exactly the
structure that matters. The only kinship-aware path,
`R/external/analyze_nullmodel.R`, is untested and fails every region on a
homogeneous-variance null model (August §5, still open). **Fix:** run the
example's sweep on the 2,504 unrelated (or the KING-unrelated set the
preamble already computes); replace the permuted null with a structured null
`y ~ MVN(0, h²·2Φ + (1−h²)I)` from the KING kinship plus an "adversarial"
null (a coverage PC, or the release batch, as phenotype); fix and test the
GENESIS path against that null.

### 1.7 The Marchenko–Pastur ndim procedure answers the wrong question and its stated justification was wrong — major, methodological (partly verified)

`example/1000G_highcov/R/choose_ndim.R`. Three problems, in decreasing
order:

- It counts components *above noise*, which is not the number that should be
  *removed* from depth: with 30 technical factors, 12 CNV blocks and
  heteroskedastic noise the count is 58 against a true 42, and the captured
  CNVs sit on PCs 25–39 — inside the removed set (August §2.1's cliff).
  Ancestry-structured PCs remove 18.5% of an ancestry-stratified deletion's
  variance even at the "right" count.
- The 1% margin was documented as "about four Tracy–Widom sd". Verified wrong:
  Johnstone's scaling at 3,202 × 142,070 puts the edge's relative sd near
  4 × 10⁻⁴, so 1% is ≈ 25 sd. The margin is an ad-hoc allowance for
  heteroskedastic misfit (per-bin variance heterogeneity of sdlog 0.5 produces
  49 / 39 / 30 / 11 spurious components at 0 / 1 / 2 / 5% margin; column
  permutation gives 1). README corrected in this pass.
- The count depends on tuning it does not expose: gap 5–120 → k@1% = 42–47;
  observed-rank-*j* vs noise-rank-*j−k* matching → 42 vs 48; "mean across
  modes" was, in every smoke run so far, a mean over one mode.

**Fix:** standardise bins before the SVD or use parallel analysis; choose
ndim on the outcome that matters — calibration on a structured null and
recovery of the NYGC SV-callset copy numbers at known deletions as a function
of ndim (which also makes the cliff measurable); report sensitivity to the gap
and both rank models.

### 1.8 Logistic: separation reports a null; the likelihood-ratio test is free — minor (reviewer-reproduced; August §2.4)

`R/analyze.R:252-266`: a chrY-like bin vs sex gives β = 64, SE = 5,930,
z = 0.011, Wald p = 0.99, `converged = FALSE`, LRT p = 0. `dev0 − f$deviance`
costs nothing per bin. Emit the LRT p and a `converged` column. Level ordering
(August §1.1: `as.factor` makes whichever label sorts first the non-event) is
still at `R/analyze.R:218`.

### 1.9 Smaller points

- p = 0 above |t| ≈ 38 (`R/analyze.R:249`); every chrM bin in the example is
  there and `evaluate.R` silently drops p = 0 from λ. Emit `LOG10P` via
  `pt(..., log.p = TRUE)`.
- β is per log2 unit — a heterozygous deletion is −1, a duplication +0.58, a
  homozygous deletion −12.7 — neither per copy nor symmetric. State it in the
  output header and README, or report on the copy scale (§1.3).
- Genotype-PC projection (`example/.../R/genotype_pcs.R:69-86`): the
  calibration r² is 1 by construction (the in-sample projection *is* the
  in-sample PC) and cannot see out-of-sample shrinkage; at the example's
  γ ≈ 0.017 the 698 projected relatives' deep PCs (GRM eigenvalue ≤ 5) are
  shrunk 20–50% toward the origin. GPC1–4 are fine. Use a held-out calibration
  or an OADP-style correction, or PC-AiR.
- The "calibration verdict" (`R/calibration_summary.R`) is descriptive: the
  fast PCs come from their own randomised-SVD run, so its distance contains one
  reseed of noise by construction and the null expectation of the ratio is 1,
  not 0; one seed control gives no distribution. Report the excess
  d_fast − d_seed with a block bootstrap over regions and ≥ 3 seed controls.
- The example's "truth checks" are alignment/scale checks (MTDNA_CN is
  2 × chrM / median; SEX is inferred from X/Y depth) and are labelled "by
  construction" honestly; the comment calling the adjusted runs "the science
  models" (`01_prepare_inputs.sh:68`) is not. Real validation for exactly these
  samples exists and is cheap: corrected depth vs NYGC SV-callset copy numbers
  at known deletions (by ndim), cis association with tag SNPs of known common
  deletions, and superpopulation-differentiated CNVs (UGT2B17, GSTM1, RHD,
  LCE3B/C) as positive controls with genotype PCs as the confounding control.

---

## 2. Silent-success paths and correctness in the core pipeline

### 2.1 The stats merge kills every unit of ten or more chunks — blocker (verified)

`scripts/correct.sh:124-127` merges the per-chunk statistics files by a
lexicographic glob (`stats.1, stats.10, stats.11, …, stats.2`), so from ten
chunks on the merged BED is unsorted, `tabix` refuses it
(`[E::hts_idx_push] Unsorted positions`), and `set -e` kills the script
*before* `dsv_output_commit`: `.tmp.gz` left, no `.done`, deterministic on
rerun. August §3.1, still open. Verified on the smoke matrix with 12 chunks.
At 3,202 samples a row is ≈ 19.6 KB and `DSV_CHUNK` caps a chunk at 2,000
rows, so the example's 25 Mb windows (25,000 rows) were 13 chunks — every unit
of the real run would have failed. **Fix:** sort the merged stats
(`sort -k1,1 -k2,2n`) or the chunk list numerically before `bgzip`; add a
≥ 10-chunk case to the smoke test. A stopgap (10 Mb windows) is applied to
the example's `config.sh` in this pass.

### 2.2 Duplicate `SAMPLE` rows in the PC or coverage table are resolved by row order — major (verified)

`R/correct.R:96` (`merge` multiplies rows), `:116` (`match()` takes the first);
no `anyDuplicated` on either table (the phenotype side has one,
`R/analyze.R:95-101`). Fixture coverage table with `SAMPLE001 1000` appended:
SAMPLE001's corrected value 0.0164 (unchanged); the same rows prepended:
**−4.849**, identical stderr, exit 0. **Fix:** refuse duplicates in both
tables.

### 2.3 August findings still open, at current lines

All reproduced by the shell pass; none has changed.

| Aug § | finding | where |
|---|---|---|
| 1.1 | logistic level ordering silently flips the sign | `R/analyze.R:218,327-331` |
| 1.2 | absolute `--minDepth` floor manufactures significance for non-normal phenotypes | `R/correct.R:39,79,208` |
| 1.3 | `.done` keyed on a filename that omits every parameter; failed `--force` keeps the old result; response-name collision | `scripts/correct.sh:56-63`, `scripts/analyze.sh:98-105,166` |
| 1.4 | no matrix-side duplicate check after `--sampleIdPattern` | `R/correct.R:116`, `R/analyze.R:148` |
| 1.5 | join checks row counts only by default; a short column is invisible | `lib/join_extract.sh:39-46`, `scripts/join.sh:242` |
| 1.6 | `analyze.sh` worker's exit status swallowed by `| tail -n +2`; commit without a row count | `scripts/analyze.sh:129,136` |
| 1.7 | "already complete" before the manifest signature — a grown cohort is silently half-analysed | `scripts/join.sh:68-71` vs `:112-122` |
| 1.8 | `regions.sh --sizes` silently drops contigs by exact string match | `scripts/regions.sh:77-88` |
| 1.9 | env thresholds override the flag the user typed | `scripts/analyze.sh:59-60` |
| 1.10 | one NA median empties the genome, exit 0; header-only shard committed with `.done` | `R/correct.R:79`, `scripts/analyze.sh:136` |
| 3.2 | `regions.sh` needs a writable matrix directory | `scripts/regions.sh:68-69` |
| 3.3 | `.tbi` assumed, no `.csi` | `lib/common.sh:180-183` |
| 5 | GENESIS path: scalar `cholSigmaInv` unhandled, exact-zero variance, no dup check, undocumented env knob, `family` unused | `R/external/analyze_nullmodel.R:90,97,145,213,217-218` |

### 2.4 New, smaller

- `DSV_MODULES` is loaded **after** the command check in every stage
  (`scripts/correct.sh:45-46`, `analyze.sh:54-55`, `join.sh:58-59`;
  `regions.sh` never loads modules), so on a cluster where `bgzip` lives behind
  `module load` the check fails first. Verified by reading; the reviewer
  reproduced it with a shell `module` shim. Swap the two lines.
- `join.sh:85-91`'s duplicate guard dies with SIGPIPE (exit 141, no message)
  when the duplicate list is large (`… | uniq -d | head -3` under `pipefail`).
- An environment-supplied `DSV_PHENO_MANIFEST` makes `--model` on the command
  line an error (`scripts/analyze.sh:23,51-52`), contradicting "flags win".
- `tests/make_fixtures.R` crashes below 75 regions per chromosome (its
  injected-signal index is hard-coded) — verified.
- `bad_log` in `tests/smoke_test.sh:27` returns non-zero as the last command of
  a `||` list and can abort the suite under `set -e` (August §4).

---

## 3. Performance at biobank scale

Measured on an M1 Max with single-threaded BLAS as the scripts pin it
(cluster cores run 1.5–2× slower; scale up); synthetic 40-PC tables and rows
at 3,202 / 50k / 200k / 500k samples, a 3,202 × 100k-bin mosdepth tree with
the real 1000G IDs, medians, 200-PC table and phenotypes, and the real chr22
PLINK files. Scripts and logs in `docs/review-2026-09/perf/`. Row widths: raw
text 6.0 B per value (**3.0 MB per row at 500k**), corrected `%.6g` text
9.6 B (4.8 MB), so `dsv_chunk_lines` hands each worker 22 rows (correct) or
13 rows (analyze) at 500k.

### 3.1 The join is a 2.5–3 day single-node job with 7 TB of scratch, and cannot run where biobank inputs live — blocker at 100k+ samples

`scripts/join.sh:165-172,222-232,242`, `lib/common.sh:170-180`. GNU `paste`
runs at 132–142 MB/s regardless of column count, so at 500k × 3.1M bins
(9.3 TB of text) the 707 batch pastes take 18–20 h, the FIFO fan-in through
one `paste` 10–20 h, and then four **serial** verification passes over the
finished file — awk `NF` (156 MB/s, 16.6 h), `bgzip -t`, the row count and
`tabix` — another 32 h, none of them resumable (resume covers batches only).
Peak scratch ≈ 7 TB in one directory (compressed batches accumulate until the
fan-in). On RAP or AoU the 10 TB of per-sample inputs cannot be localised to
one instance at all. At 3,202 × 100k bins the whole join is 232 s.

**Fix:** a tile-parallel join with no fan-in — either a task per sample batch
streaming each file once (0.6–1.7 s per genome-sized file, coordinate
checksum in the same awk pass, which also makes `--strict-coords` free) and
writing one tile per window, or a task per window range-reading it from every
sample's csi-indexed `regions.bed.gz` (0.01 s per read; 70M reads ≈ 194 CPU-h
in total, fully parallel). Either yields the region-addressed store of §3.5
directly and removes every serial stage. Interim: fold the four verification
passes into one, use `gawk`/`mawk` explicitly (BSD awk is 4× slower), and
raise the final `bgzip` threads or use `-l 1`.

### 3.2 correct/analyze spend > 90% of CPU on process restarts and ASCII — blocker at scale

Every 64 MB chunk starts a fresh `Rscript` that re-reads the **whole** PC
table (all 200 columns, not the `ndim` used), re-factorises the design and
rebuilds `Q` (`R/correct.R:81,140-174`), or re-reads the phenotype table and
rebuilds the design (`R/analyze.R:93,180`). Fixed cost per process at 500k:
2.9 s (40-PC table, ndim 40), **22 s and 10.5 GB** (200-PC table, ndim 200);
analyze 3.0 s linear / 3.6 s logistic. Marginal cost per row at 500k: correct
0.59 s (of which `sprintf` + `paste` + `cat` 85%, the maths 5%), analyze
0.45 s (almost all the character-then-numeric parse; the fit is 4 ms/bin
blocked).

| stage at 500k × 3.1M bins | fixed | per row | total | useful |
|---|---|---|---|---|
| correct, ndim 40 | 115 CPU-h | 508 | **≈ 620 CPU-h** | 32 h (5%) |
| correct, ndim 200 (PLAN's validated setting) | 861 | 610 | **≈ 1,470 CPU-h**, 10.5 GB × jobs | ≈ 150 h |
| analyze linear, per phenotype | 197 | 388 | **≈ 585 CPU-h** | 3 h (< 1%) |
| analyze logistic, per phenotype | 236 | 1,214 | ≈ 1,450 CPU-h | see §3.3 |

The example's six-phenotype manifest at 500k is ≈ 5,000 CPU-h. At
`DSV_JOBS=12` the correct workers alone need 29 GB (ndim 40) or 126 GB
(ndim 200) — the WDL's 16 GB default and the example's 16 GB unit cannot host
them. **Fix, with measured gains:** cache `Q` once as an uncompressed `.rds`
(176 MB at 500k × 41, 0.26 s to load) and read only the `ndim` columns
(`fread(select=)`) — removes the 115–861 h of prologue; process a chunk as a
matrix (numeric `fread`, one GEMM pair, write from the matrix) — 1.8–3× on
the same chunk today; binary chunks (`readBin` float32 at 14–16 ms/row vs
190–450 ms text parse, `writeBin` 6 ms vs 481 ms) — projected 10× for correct
(≈ 50 CPU-h per pass) and ~25× for analyze-linear (≈ 17 CPU-h per phenotype).
A persistent worker over a contiguous row range removes the width dependence
entirely; a bigger block is not the lever. GNU parallel itself is not a
bottleneck (700 MB/s through `--pipe -k`), but `-L` must stay derived from row
width — a fixed `-L 2000` at 3 MB rows makes parallel regrow the block to the
whole stream.

### 3.3 Logistic: per-bin IRLS is 430–780 CPU-h per phenotype at 500k; a score test is 56–350× cheaper — major

`R/analyze.R:252-266`. `glm.fit` per bin: 2.8–5.3 ms at 3,202 (fine — the
example's logistic unit took 34 s), 77–89 ms at 50k, **503–909 ms at 500k**
(2% prevalence needs 6–7 iterations). A block score test with the null model
fitted once is 0.015 ms → 7–9 ms per bin, ≈ 7.8 CPU-h per phenotype at 500k.
For imbalanced case/control at biobank scale a saddlepoint-corrected score
(SAIGE-style) or Firth is needed for calibration, not just speed — which is
PLAN.md's Phase-3 seam: feasibility argues for finishing the REGENIE/GENESIS
handoff rather than optimising `fit_logistic`.

### 3.4 `analyze.sh` re-reads, decompresses and parses the region once per phenotype — major for any sweep

`scripts/analyze.sh:112-114,132`. At 3,202 × 25k rows the six-phenotype
manifest costs 124 s where one phenotype costs 16 s; at 500k the per-phenotype
parse alone is 388 CPU-h, so a six-phenotype sweep parses the same 15 TB six
times. **Fix:** parse once per worker and loop the manifest rows over the
parsed block (the linear block fit is 4 ms/bin per phenotype).

### 3.5 Text is the wrong storage format at this width — major

Parse 18–67 MB/s and format ~10–60 MB/s against decompression at 550–700 MB/s;
`readBin` float32 at 36–38 M values/s versus `fread` at 6.7 M. Compressed text
is no smaller than compressed float32 (a 500k × 3.1M corrected matrix is
5.6–6.4 TB either way, per ndim); the disk win needs quantisation — int16
(×1000, max error 5e-4 log2 units) is 2.4 TB, 2.3–2.5× smaller. A row-group
store (tiles of ~1,000 bins × all samples, zstd, a contig/start → offset
index) keeps `--region` as the unit of work, adds sample slicing that text
cannot do, and is cloud-native (Zarr/HDF5); for the corrected matrix
specifically, PLAN.md's Phase-3 **PGEN 16-bit dosage** export is the same
2 B/value with a variant-range index and makes the matrix directly readable
by REGENIE/PLINK. The `#` header, `tabix` and bgzf are load-bearing only
because the payload is text.

### 3.6 Smaller

- `scripts/correct.sh:78` decompresses each region twice (a `wc -l` pass for a
  log denominator and a parity check) and pays a full worker prologue for the
  header probe: 200 s of serial time per 25 Mb unit at 500k before parallel
  starts. Count inside the workers; derive the header from the sample list.
- `R/analyze.R:290-301`: `fread(colClasses = "character")` then
  `as.numeric(unlist())` is 2.4–6.2× slower than a numeric `fread` with
  `colClasses = list(character = 1:4)`, identical results — a one-line change
  worth 2–6× on the analyze stage today. The per-row loop vs one GEMM is a
  further 4–6× (max |Δβ| 5e-17), after the parse is fixed.
- Example sizing at 3,202 samples (measured on a real 25 Mb unit, 8 jobs,
  ndim 20): correct 27 s, the six-phenotype analyze 124 s, 0.65 GB peak — so
  ≈ 150 s and < 1 GB per unit against the requested 4 h / 16 GB, and the join
  ≈ 30 min / 46 GB scratch / < 1 GB RSS against 12 h / 48 GB. Over-requesting
  by 25–95× is harmless for correctness and expensive in queue time; 30–45 min
  / 2 GB per unit and 1 h / 8 GB for the join would schedule far sooner. The
  preamble's plink2 work is ≈ 3 min for the whole genome (< 1 GB); the
  serial ~4 GB Dropbox download dominates, and `--memory` at 85% of a 248 GB
  node reserves ~210 GB of workspace for nothing. 4 CPU / 8 GB / 2 h with the
  downloads prefetched in parallel would do; the 24-core node was requested
  by design, so this is information, not a defect.
- `lib/join_extract.sh:36`: one `gawk` pass can emit the column and the
  coordinate checksum together (1.7 s vs 0.6 + 2.4 s per genome-sized file),
  which is the cheap route to making `--strict-coords` the default (August
  §1.5).

---

## 4. Deployment

### 4.1 SLURM

- **`workflows/slurm_array.sh` cannot find the stage scripts under `sbatch`** —
  blocker for the documented invocation (verified by reading). `:26` resolves
  `here` from `BASH_SOURCE[0]`, but sbatch runs a spooled copy of the script, so
  `$here/../scripts/` does not exist there. The example's `lib.sh:10-13`
  documents this exact trap and works around it; the dispatcher never got the
  fix. Resolve through an exported `DSV_ROOT` (`lib/common.sh` already exports
  it once sourced) and refuse to run if `scripts/correct.sh` is not found.
- `--array=1-N` without a throttle, and N > 1,000 rejected by the default
  `MaxArraySize` (`README.md:199`; 1 Mb windows give ~3,100 units). Document
  `%100` and splitting.
- `--jobs "$cpus" --threads "$cpus"` together oversubscribe the node 2×
  (`slurm_array.sh:43,48`).
- The `DSV_MODULES` ordering bug (§2.4) means the module mechanism the README
  advertises cannot work on the cluster it was written for.
- `#SBATCH --output` in `preamble.sh:3` and `slurm_array.sh:3-4` writes into the
  submit directory — the repository checkout — and `.gitignore`'s `*.out`
  hides the litter; a read-only shared checkout kills the job at start.
- The example's chain relies on `sbatch` from a compute node, on
  `--export=ALL` being the site default, and on `--dependency=afterok` without
  `--kill-on-invalid-dep=yes` (a failed join leaves the dependents pending
  forever). `sbatch --parsable` output is used verbatim; federated clusters
  return `jobid;cluster`.

### 4.2 Both clouds: the image and the matrix

- **No pullable image — blocker on AoU and RAP.** `workflows/depthsv.wdl:26`
  defaults to `depthsv:dev`; CI builds `depthsv:ci` and never pushes; no
  registry reference exists anywhere in the repo, while `README.md` and
  `test.yml` describe the image as "the artefact the cloud platforms actually
  run". Cromwell resolves the tag to a digest before any VM starts; AoU batch
  VMs reach Docker Hub only through its Artifact Registry remote; RAP wants
  images stored in the project. The workstation that builds it is Apple
  Silicon, so an unqualified `docker build` yields `linux/arm64`. Publish
  `linux/amd64` from CI on tags, by digest, and make that the WDL default;
  document the AoU remote-repo path and the `dx upload` tarball route.
- **Every scatter shard localises the entire matrix — blocker at biobank
  scale, wasteful at any scale.** `depthsv.wdl:15-17,90-91,114`: the matrix
  and index are plain `File` inputs with a 100 GB `local-disk` and no
  `parameter_meta`. At 500k samples (≈ 9 TB uncompressed, 2–3 TB bgzipped)
  the copy fails outright; at 1000G scale 140 shards copy ~2.5 TB for ~7 GB of
  useful reads. On dxCompiler, `parameter_meta { depth_matrix: { stream: true } }`
  (dxfuse) works with the scripts unchanged; on Cromwell,
  `localization_optional: true` hands the task a `gs://` path, which needs the
  scripts to stop `[ -s ]`-testing inputs (use `tabix -l` / `tabix -H`), a
  `GCS_OAUTH_TOKEN` from the metadata server (htslib reads only that variable),
  and `curl` in the image (rocker purges it).
- **The join has no cloud shape — blocker on AoU, major on RAP.** It reads
  every per-sample file from a local manifest on one node
  (`scripts/join.sh:38-58,142-146`); on both platforms the mosdepth outputs are
  objects in a bucket, and a 500k-file `Array[File]` task is not practical.
  The redesign is already half-written: `join.sh`'s batch files
  (`:165-172`) are sample-block matrices. Make the block the unit — a
  `join_block` task per ~1,000 samples writing an indexed block, and a
  `window_unit` task that pastes a region from each block (`join.sh:222-232`
  verbatim) and runs correct + analyze on it. A miniwdl-checked sketch is in
  `docs/review-2026-09/cloud/depthsv_cloud_sketch.wdl`.
- **WDL runtime gaps.** No `maxRetries`, no `bootDiskSizeGb`, and GNU parallel's
  `--pipe` buffers land in `$TMPDIR` on the boot disk (set `TMPDIR=$PWD/tmp`);
  dxCompiler ignores `preemptible`/`maxRetries` (needs `dx_restart`,
  `dx_timeout`); outputs drop the `.tbi` indexes and the per-unit logs the
  README calls "the first thing to check"; `regions` must be hand-converted to
  JSON (`read_lines`); `DSV_BLOCK_BYTES`, `DSV_CHUNK`, `DSV_MIN_VARIANCE` are
  unreachable from the WDL; the default 16 GB / 4 CPU will OOM at ndim = 200
  and 500k samples (four concurrent R processes need ~42 GB — see §3).
- **All of Us policy.** Every output row carries `N`, `NCase`, `NControl`
  (`R/analyze.R:233,324-335`) with no minimum on the case count, and the
  correction stats carry per-bin `min`/`max` — a single participant's value.
  AoU's dissemination policy forbids publishing counts of 1–20 or anything they
  can be derived from. Add `--min-cases`, an export step that concatenates
  shards, verifies coverage against the region list and suppresses small
  counts, and keep the stats files in the workspace. PLAN.md names the policy;
  nothing implements it.
- **RAP.** UKB already provides 40 genotype PCs (Data-Field 22009), and AoU
  ships ancestry PCs with the CDR — the preamble's plink2 recipe is for cohorts
  without them; say so. Neither platform allows the example's GitHub/Dropbox
  fetches from workers.

### 4.3 "Interchangeable dispatchers" is overstated

The example never calls `workflows/slurm_array.sh`; it has its own timed unit
stage with different resources (8 CPU / 16 GB / 4 h vs 8 / 32 / 8 h), and only
one of the two survives `sbatch` (§4.1). They agree on the unit of work, which
is what matters — say that instead.

---

## 5. The example and the preamble

These were written in the last month, by me, and had not been reviewed. The
adversarial pass found one blocker and ten majors; the ones I verified are
marked.

- **B1 — every real unit dies (verified).** §2.1 at the example's default
  window; stopgap applied (`config.sh`, 10 Mb → 5 chunks). The core fix is
  still needed.
- **M1 — same analysis name, changed model → stale shards under the new
  model.** `01_prepare_inputs.sh:84-98` keeps `sex_linear` as the name for both
  the unadjusted and the covariate-adjusted sex model, so `analyze.sh` skips it
  as complete after the preamble lands (reviewer-reproduced: shard hash
  unchanged, evaluation PASS on the stale file); `EX_COVARIATES`,
  `EX_N_GPCS`, `EX_PHENO_SEED`, `EX_MIN_OBS` and `EX_MEDIAN_SOURCE` all feed
  outputs whose names do not encode them. The README's "changing `EX_NDIM` can
  never silently reuse" is true only of ndim. This is August §1.3 in a worse
  form, and the example documents `.done` reuse as its recovery procedure.
  **Fix:** freeze the resolved parameters into `inputs/<mode>/run.env` in
  stage 1, hash them into the output directory name the way `assoc_ndim`
  already is, and make every stage read that file — or port `join.sh`'s
  `resume.meta` to `correct.sh`/`analyze.sh` and fix all of §1.3 at once.
- **M2 — `EX_NDIM` and `EX_COVARIATES` are re-derived from disk inside every
  SLURM job.** `config.sh` reads `preamble/ndim.txt` and tests for
  `covariates.tsv` at source time and exports neither; the driver passes no
  `--export`; the quick start submits the preamble and the run back to back. A
  preamble finishing mid-array splits units between `assoc_ndim20` and
  `assoc_ndim42` and the evaluation reports "PASS with warnings" on a partial
  genome. Same fix as M1; refuse to submit while `dsvx-preamble` is queued.
- **M3 — a grown manifest is skipped by the join** (August §1.7 through the
  example's own "lagging fast tree" scenario; reviewer-reproduced: 65-sample
  manifest, 64-column matrix, PASS). The driver even logs the new count before
  `join.sh` skips. Compare the manifest against `depth.matrix.manifest` before
  the short-circuit.
- **M4 — the calibration headline says "within seed noise" when every row is
  undetermined** (`R/calibration_summary.R:47-64`: the headline keys on
  `n_exceed == 0`; a seed distance of exactly 0 yields `Inf` → "undetermined"
  rather than "exceeds"). Verified by inspection of my own code.
- **M5 — a stale `.finalize.lock` disables the comparison and profile on every
  later run, exit 0** (`02_run_depthsv.sh:127-139`; the eval job has a 2 h
  limit and 04 reads 3.1M-row shards). Record the owner, expire it, clean it
  in the driver.
- **M6 — a concordance FAIL is reported as "one mode has no output yet", exit
  0** (`04_compare_modes.sh:48-67`: `compare_pair` returns the Rscript status,
  so a FAIL takes the missing-directory branch; `maybe_finalize` turns a
  non-zero 04 into a WARN). Verified by inspection.
- **M7 — the real-run GitHub fallback pairs local depths with a different
  run's PCs and medians** with one log line and 0 WARN (`00_fetch_inputs.sh:176-179`;
  August §6's last bullet, still present). Make it opt-in.
- **M8 — every smoke and CI "fast" tree is synthetic and no report says so.**
  Upstream has not committed `output_fast/` or `autosomal.median.txt`, so the
  fallback (standard tables + 1% jitter) runs in every current smoke run; only
  `paths.env` and `smoke.params.txt` record it. Carry the source tag into every
  summary and name the mode `synthetic-fast`.
- **M9 — the permuted null differs between modes** whenever the QC row sets
  differ (`R/prepare_inputs.R:162-166` permutes each mode's own rows), so a
  lagging tree makes `mtdna_cn_null` "exceed seed noise" for a reason that has
  nothing to do with fast mode. Build the null once.
- **M10 — concordance PASS on two common regions** (`R/compare_modes.R:88-108`,
  no minimum `n_common`).
- **Three real samples are silently dropped (verified):** `HG02635`,
  `HG03025`, `HG03366` carry a trailing space in the upstream ID
  (`"HG02635 .by1000."` in the PC table, `"HG02635 "` in the QC table);
  `fread` trims one and not the other, and `prepare_inputs.R` notes "3199 of
  3202" without a WARN. `trimws()` every ID and warn on any mismatch.
- **Two statements I wrote were wrong and are corrected in this pass:** the
  Tracy–Widom magnitude (§1.7), and the instruction to watch the chrM *slope*
  move with ndim — by FWL the slope is invariant (0.95/0.94/0.94/0.91/0.87/0.93
  at ndim 0/2/4/10/20/40) and the t-statistic is what shrinks (35 → 4).
- Smaller: `choose_ndim.R` takes column 2 unchecked and crashes for a gap past
  the reported spectrum; `preamble.sh` republishes a stale `ndim.txt` when the
  fit is undetermined; `evaluate.R`/`compare_modes.R` stop reading
  `analyses.tsv` at the first mid-file comment; the logistic sex run is scored
  in the comparison like a real test (r(Estimate) = 0.08 under separation);
  the null bands (0.85–1.20, 0.03–0.07) are ±100 sd wide at 3M bins;
  `run.sh` prints `ndim=4 (from the preamble)` in smoke mode; hand-rerun stages
  after a smoke run resolve `assoc_ndim42` and exit 0 having evaluated nothing;
  a failed stage 1 leaves the previous run's tables in place and
  `ex_mode_ready` still says ready; `EX_JOBS`/`EX_THREADS` ignore the
  allocation; profile records accumulate skipped-unit entries on resubmission.

---

## 6. Publication-grade hygiene

- **No release, no version, no provenance.** No tags; no version string in any
  script; no output carries the tree it was made by, the parameters, or the
  R/htslib versions. Every summary-statistics header should carry
  `# depthSV <version> <commit>` plus the resolved parameters.
- **`CITATION.cff` lists "depthSV contributors"** with no ORCID, version, date
  or DOI; a paper needs real authors and a release DOI (Zenodo on tag).
- **Six quantitative claims in README have no reproducer in the repo**
  (August §6: "12–20× faster at 25,000 samples", "agrees with qr.resid to
  ~6e-14", "linear and logistic agree to ~1e-13", "exact against an
  independent recomputation", "byte-identical to the original", "byte-identical
  between container and host" — the CI jobs never compare outputs). Turn each
  into a script under `benchmarks/` or delete it.
- **The test suite cannot see the estimator** (August §4, unchanged): deleting
  the PC correction passes 35/35; five of eight estimator mutations survive;
  coxph has no numeric assertion; the λ gate cannot fail; nothing asserts a
  sign or an estimate; fixtures have no zero-depth bins, no missing values, no
  sex chromosomes. Add reference values, a case/control-string fixture, a
  ≥ 10-chunk unit, a sex-chromosome fixture, and a structured-null λ.
- **CI** depends on a third party's HEAD (`raw.githubusercontent.com/…/master`),
  asserts file existence and the evaluations' exit code only (the compare and
  calibration verdicts never gate; `grep 'verdict'` is case-sensitive and the
  calibration file writes `Verdict`), and never runs the plink2 half of the
  preamble or a macOS/bash-3.2 job although bash 3.2 is a stated requirement.
- **Docs vs code** (August §6, all still true): the 0-based `--region` form in
  README duplicates a boundary bin; the windowed partition guarantee is
  conditional on uniform bins; `--chunk` is not what README says it is; the
  scratch-footprint claim omits the accumulated batches; `--minDepth`,
  `--digits`, `--statsFile`, `--skipOutputHeader` undocumented; `cov_resids`
  appears nowhere in README; `.by1000` hard-coded in `dsv_sample_name` while
  bin size is an open decision; `USER` unbound under `env -i`.
- **Dockerfile** pins R and CRAN but not `tabix`/`parallel`; no
  `org.opencontainers` labels; no `curl` (needed for GCS token fetch and for
  the example); `PARALLEL_HOME` only silences the citation notice for root.

---

## 7. Checked and not a problem

Recorded so they are not re-litigated.

- bash 3.2: the smoke suite, the example end to end (both modes), the
  dispatcher's login-node path and every `--help` pass under `/bin/bash`
  3.2.57; no bash-4 constructs anywhere.
- Paths with spaces through the join (`parallel -q`), the example
  (`EX_WORK_DIR` with a space) and `printf %q` → `source paths.env`
  round-trips.
- `join.sh` resume: killed mid-batch, rerun resumes and the matrix is
  byte-identical to a clean join; `--force` rebuilds everything it claims to.
- EXIT traps run on SIGTERM; no `.tmp`/`.part` litter after runs.
- Shard globs: `<name>.<method>.*.txt.gz` never cross between `mtdna_cn`,
  `mtdna_cn_adj`, `mtdna_cn_null_adj`, `log2_mtdna_cn`; `.tmp.gz`/`.tbi`/`.log`
  never match.
- The plink2 recipe: pvar codes are `22` (no `chr`) and match the LD file;
  `--exclude range` behaves identically with either code; `.eigenvec.allele`
  layout is `#CHROM POS ID REF ALT A1 PC1…` so `--score … 3 6` and
  `--score-col-nums 7-46` are right; `--pmerge-list` merges; KING, PCA and
  projection run; calibration r² = 1.0000 by construction (§1.9 explains why
  that is not evidence).
- MP numbers on the committed spectrum reproduce under an independent
  implementation (52 / 42 / 38 / 32 at 0 / 1 / 2 / 5%) *given* the script's
  rank model.
- `explicit` vs `qr` projection agree, including rank truncation (August §7).
- The `--minDepth` floor does not break a *normal* phenotype under the
  *unadjusted* design (August §7) — it is §1.1 and §1.3 that do.

---

## 8. Changes made in this pass

- Credit added to `README.md` and the example README (requested).
- Example README: Tracy–Widom statement corrected; "watch the slope" corrected
  to the t-statistic; `EX_WINDOW` table entry updated.
- `example/1000G_highcov/config.sh`: `EX_WINDOW` default 25 Mb → 10 Mb with the
  reason inline (stopgap for §2.1).
- The August review preserved as `REVIEW-2026-08.md`; reproduction scripts
  under `docs/review-2026-09/`.

No pipeline code was changed: the findings above are the case for the fixes,
in the order below.

## 9. Suggested order

1. **Estimator PR** — §1.1 (PCs into the model, ~10 lines + a smoke assertion
   on λ with a PC-structured phenotype), §2.1 (sort the stats merge, one line +
   a ≥ 10-chunk test), August §1.1 (explicit 0/1 or a reference-level option),
   §2.2 and August §1.4 (duplicate-ID guards on all three tables and the
   matrix side).
2. **Silent-success PR** — parameter-aware `.done` markers via a sidecar hash
   (fixes August §1.3, M1, M2 at once); drop `| tail -n +2` and add the row
   parity check (§1.6); make "produced nothing" non-zero (§1.10 and the
   header-only shard); manifest check before the join short-circuit (§1.7,
   M3); `DSV_MODULES` ordering; `slurm_array.sh` path resolution.
3. **Coding & QC PR** — copy-ratio or winsorised scale and a relative floor
   (§1.3), per-bin max-share leverage flag, sex-chromosome ploidy (§1.2), LRT
   and `converged` for logistic (§1.8), `--min-cases`, `LOG10P`, a Cox numeric
   test (§1.4).
4. **Inference PR** — permutation maxima per shard and an export step that
   concatenates, verifies coverage, computes the empirical threshold and M_eff,
   and applies count suppression (§1.5, §4.2); structured null and the
   unrelated-set sweep in the example (§1.6); ndim chosen against the SV
   callset (§1.7).
5. **Deployment PR** — publish the image from CI by digest and default the WDL
   to it; streaming/`localization_optional` plus the `tabix -H` / token
   changes; the block-join WDL; runtime attributes; release tag, version
   stamp in headers, `CITATION.cff` authors, benchmark reproducers.

## 10. Status

Where each finding stands, by commit on the `review-2026-09` branch. "PR3"
is the coding-and-QC commit that carries this section.

| Finding | Status | Where |
|---|---|---|
| §1.1 PCs into the association model | fixed | `beaf00d` (PR1): `--pcs`/`--ndim` on the analysis stage, refusal without them, PC-null λ assertion |
| §1.2 Sex-chromosome ploidy | fixed | PR3: `--sex`/`--sex-col`/`--par` on the correction, `conf/par.grch3[78].bed`, chrY fitted on males only, example runs with it and checks that chrX carries no sex signal |
| §1.3 Floor coding / single-sample leverage | fixed | PR3: log2 winsorised at `--winsor-log2` (−3), `MAXSHARE` on every row, `--max-share` (0.5) skip |
| §1.4 Cox numerics | fixed | PR1: `coxph()` reference assertion; PR3: the floor and leverage filter above, `--robust` |
| §1.5 Effective number of tests, threshold | fixed | `d6480d8` (PR4a): `--perms` max-|t| per shard, `scripts/export.sh` with the Westfall–Young threshold, M_eff, Bonferroni and λ; the example exports every analysis |
| §1.6 Relatedness | partly | PR4b: the example runs the mtDNA phenotypes and nulls on the KING-unrelated set, draws a structured null from the kinship (`--make-king square` in the preamble) and a coverage-PC null; the GENESIS path is still untested |
| §1.7 Marchenko–Pastur ndim | partly | `dc3a42c` corrected the stated justification; PR4b adds `06_sv_recovery.sh`, recovery of NYGC-callset deletions as a function of ndim with the MP count marked and a plateau recommendation; bin standardisation before the SVD is upstream (NGS-PCA) |
| §1.8 Logistic separation | fixed | PR1: `LRT_P`, `CONVERGED`, null deviance once |
| §1.9 Smaller | partly | `--min-cases`, `LOG10P`, explicit binary coding (PR1); rank-INT, robust SEs, no mean imputation (PR3); the rest open |
| §2.1 Stats merge at ≥ 10 chunks | fixed | `dc3a42c` stopgap, `beaf00d` sorted merge + ten-chunk test |
| §2.2 Duplicate `SAMPLE` rows | fixed | `beaf00d`: refused in every table and in the matrix header |
| §2.3 August findings (markers, `tail`, parity, empty shards, manifest check, modules, `DSV_ROOT`) | fixed | `42f2823` (PR2): parameter-aware `.params` signatures with `--force`, row parity, empty-shard refusal, manifest-content check, module order, spool-safe `DSV_ROOT` |
| §2.4 New, smaller | fixed | `42f2823` |
| §3.1 Join at biobank scale | open | PR5 / later: tile-parallel and block join |
| §3.2 CPU overhead in correct/analyze | partly | PR1/PR2: basis computed once per unit, block projection, numeric `fread`, `select=`; process restarts and text parsing remain |
| §3.3 Logistic cost | open | score test not implemented |
| §3.4 Region parsed once per phenotype | open | |
| §3.5 Text storage | open | |
| §4 Deployment | partly | PR3: WDL carries the ploidy, winsor and leverage inputs; PR4a: `--min-count` suppression in the export (All of Us), permutation inputs in the WDL; image publish, streaming localisation, runtime attributes are PR5 |
| §5 Example and preamble | fixed / partly | `42f2823`: lock ownership, `.selected_modes`, run.env, provenance; PR3: ploidy wiring, INT rows, new sex-table checks; PR4b: export per mode, unrelated-set and structured-null rows, SV-callset recovery stage |
| §6 Hygiene | partly | credit and review docs landed; version stamp, `CITATION.cff` authors, release tag are PR5 |
