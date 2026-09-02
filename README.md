# depthSV

Association testing of read depth against phenotypes, at biobank scale.

depthSV takes per-sample read depth over fixed genomic bins, normalises and
corrects it for technical structure, and tests every bin against one or more
phenotypes. It emits association summary statistics and stops there.

The premise is that copy-number variation shows up as a continuous shift in
read depth, and that testing that shift directly avoids the information loss of
first discretising depth into copy-number calls.

---

## Scope

**In scope:** joining per-sample depth into a matrix, normalisation and
technical correction, per-region quality control, and association testing
(linear, logistic, Cox) producing summary statistics.

**Out of scope:** report rendering, peak curation and result annotation.
depthSV emits summary statistics; what consumes them is a separate concern.

---

## Quickstart

No real data required — the test suite generates its own.

```bash
tests/smoke_test.sh
```

That builds a synthetic cohort, runs the whole pipeline over it, and asserts on
the results: that an injected association is recovered, that a null phenotype
stays calibrated, and that the guardrails refuse malformed input.

For the same pipeline on real data — the 3,202-sample 1000 Genomes 30x
cohort, with phenotypes whose genomic architecture is known by construction,
an upstream mosdepth fast-mode comparison carried through to the association
statistics, and a runtime profile — see
[`example/1000G_highcov/`](example/1000G_highcov/). It runs as a SLURM chain
on an HPC, or anywhere in minutes via its smoke mode.

To work through it by hand:

```bash
Rscript tests/make_fixtures.R tests/fixtures 60 200

scripts/join.sh \
  --manifest tests/fixtures/mosdepth.input.txt \
  --out work/join

scripts/correct.sh \
  --matrix work/join/depth.matrix.txt.gz \
  --pcs tests/fixtures/svd.pcs.txt \
  --coverage tests/fixtures/autosomal.median.txt \
  --region chr1 --out work/corrected --ndim 4

scripts/analyze.sh \
  --corrected work/corrected/corrected_ndim4.chr1.txt.gz \
  --pheno tests/fixtures/phenotypes.tsv \
  --pheno-manifest conf/phenotypes.example.tsv \
  --pcs tests/fixtures/svd.pcs.txt --case-level case --min-cases 5 \
  --region chr1 --out work/assoc -- --minObs 30
```

---

## The three stages

### 1. Join — `scripts/join.sh`

Combines per-sample `mosdepth` region files into one matrix and indexes it.

The matrix is written with `bgzip` and a `tabix` index rather than plain gzip.
That is what lets every later stage read a region instead of scanning the whole
file, and it is what makes an arbitrary interval a valid unit of work.

Row counts and a checksum of the coordinate columns are verified inside every
extraction worker, in the same decompression pass. `paste` aligns by position
and cannot tell that one sample was processed against a different contig set
or bin order — the offending sample's depths would simply be attributed to
the wrong coordinates, and nothing downstream would notice. A mismatch
therefore stops the run at the offending sample, not at the end. A manifest
listing the same sample twice is refused for the same reason: later stages
match columns by name. A finished matrix is reused only for the same manifest
*content*; adding a sample and re-running rebuilds it.

Columns are extracted `--jobs` at a time and pasted in batches sized near √N
by default, which keeps every step clear of the open-file limit — the limit is
checked up front, with the fix in the error message. Each finished batch is
compressed and its inputs deleted, so scratch stays near two batches of
columns rather than a full uncompressed copy of the matrix, and a re-run
resumes from the last finished batch as long as the manifest and batching are
unchanged.

Writes `depth.matrix.txt.gz`, its index, and a manifest recording the sample
list, region count and coordinate source.

### 2. Correct — `scripts/correct.sh`

Converts raw depth to a log2 ratio against each sample's median autosomal
coverage, then residualises against the leading principal components of the
depth matrix to remove technical structure.

`--ndim 0` gives the log2 ratio with no correction. How many components to
remove depends on the batch structure of your cohort and is best established
empirically against your own data; the pipeline does not assume a value.

**Sex chromosomes are normalised by their expected copies** when a sex table
is given (`--sex FILE`, `--sex-col NAME`; coded M/F, male/female or 1/2 —
0/1 is refused as ambiguous), with the pseudo-autosomal regions from `--par`
(`conf/par.grch38.bed`, `conf/par.grch37.bed`): chrX outside the PARs is one
copy in males, chrY one copy in males and none in females, whose chrY depth
becomes missing rather than the log2 of nothing. Without it every
sex-chromosome bin is a sex indicator and dominates any phenotype that
correlates with sex. The log2 ratio is floored at `--winsor-log2` (default
−3) so a zero-depth sample cannot carry a region's whole leverage.

Also writes per-region pre- and post-correction summary statistics. Nothing
reads them automatically; they are the evidence for choosing the QC thresholds
applied at the next stage.

The residualisation is computed as `v − Q(Qᵀv)` from an explicitly formed
orthonormal basis. This is the same projection `qr.resid()` computes, measured
12–20× faster per region at 25,000 samples. It costs one `qr.Q()` reconstruction
up front, which at that scale paid for itself after ~9 regions at `ndim=16` and
~76 at `ndim=200` — any real work unit is thousands. Pass `--projection qr` for
the Householder path instead; the test suite asserts the two produce identical
output.

### 3. Analyze — `scripts/analyze.sh`

Tests each region against a phenotype and streams one row of summary statistics
per region.

Pass `--model` and `--method` for a single analysis, or `--pheno-manifest` for a
sweep. The manifest is tab-separated — `name`, `method`, `model`, and an
optional fourth column of flags (`rank-int`, `robust`, comma-separated) — so
adding a phenotype is a line in a file. Analyses that already completed are
skipped, so a sweep can be extended without repeating finished work.

**The PCs the correction removed are part of every model.** Residualising
depth against PC1..PCk and then testing without those PCs in the model deflates
every test by roughly 1 − R²(phenotype ~ PCs), and biases any bin that
correlates with a covariate that itself correlates with a PC — every
sex-chromosome bin once sex is a covariate. So the analysis stage appends
PC1..PCk to the covariates of every model: `k` is read from the corrected
file's name (`corrected_ndim<k>`), the table from `--pcs` or `DSV_PCS`, and a
corrected matrix with k > 0 refuses to run without them (`--ndim 0` is the
explicit opt-out). The design is still built once per worker; each region
costs one projection, computed for a whole block of regions at a time.

A binary response must be coded `0`/`1` or `TRUE`/`FALSE`, or be a two-level
text column with the case named by `--case-level`; other codings are refused
rather than guessed, because a guess once reversed effect directions. A
binary or survival phenotype with fewer cases, controls or events than
`--min-cases` (default 20) is refused up front.

A region whose depth is missing for some samples — chrY in females under the
ploidy model — is fitted on the samples that have it, with the design
refactorised for that subset (and cached, since the subset is the same for
every such region); nothing is imputed. Every row reports `MAXSHARE`, the
largest share of the residualised depth's sum of squares carried by one
sample, and a region above `--max-share` (default 0.5) is skipped: a test
carried by one participant is the read-depth analogue of a minor allele
count of one. `--rank-int` applies a rank-based inverse-normal transform to a
quantitative response and `--robust` uses HC1 standard errors (linear) or the
robust variance (Cox); both are also per-row manifest flags.

### Options

Each stage takes its inputs and output directory as flags; `--help` on any of
them prints the full usage. The rest are these:

| Flag | Stages | Meaning |
|---|---|---|
| `--jobs N` | all | parallel workers within one unit |
| `--threads N` | all | compression threads |
| `--chunk N` | correct, analyze | maximum rows handed to one worker; the actual number is whatever fits in the memory bound, so this is a cap rather than a target |
| `--force` | all | redo a unit that already completed |
| `--min-obs N` | analyze | skip a region with fewer complete observations |
| `--min-variance X` | analyze | skip a region whose depth does not vary |
| `--min-cases N` | analyze | refuse a binary or survival phenotype with fewer cases, controls or events (default 20) |
| `--max-share X` | analyze | skip a region where one sample carries more than this share of the residual depth sum of squares (default 0.5) |
| `--rank-int`, `--robust` | analyze | inverse-normal transform of the response; HC1 / robust variance. Per row in a manifest: fourth column `rank-int,robust` |
| `--pcs FILE`, `--ndim K` | analyze | the PC table and count the correction used; `k` defaults to the corrected file's name, the table to `DSV_PCS` |
| `--case-level L` | analyze | which level of a two-level text response is the case |
| `--sex FILE`, `--sex-col NAME` | correct | per-sample sex, turning on the ploidy model for chrX/chrY (`DSV_SEX`, `DSV_SEX_COL`) |
| `--par BED` | correct | pseudo-autosomal regions, diploid in both sexes (`conf/par.grch38.bed`) |
| `--winsor-log2 X` | correct | floor on the log2 ratio (default −3) |
| `--name` | analyze | label for a single-model run; defaults to the response variable |
| `--batch-size N` | join | samples per `paste` batch; the default sizes it near √N, inside the open-file limit |
| `--projection qr` | correct, analyze — after `--` | use the Householder projection instead of the default |
| `--sampleIdPattern` | correct, analyze — after `--` | PCRE to recover sample IDs from column names |

Anything after `--` is passed through to the R driver unchanged, which is how
the last two are set — the stage scripts themselves reject flags they do not
know.

### Configuration

The flags a site would fix — inputs, output directories, parallelism, QC
thresholds — have matching `DSV_*` environment variables, so they can be set
once instead of repeated. Per-run flags such as `--region` stay on the command
line:

```bash
cp conf/example.env conf/mysite.env    # conf/*.env is gitignored
$EDITOR conf/mysite.env
source conf/mysite.env
scripts/correct.sh --region chr1       # inputs and output come from the env
```

Flags always win over the environment. `conf/example.env` documents only
variables the code actually reads.

---

## Regions are the unit of work

Every stage after the join takes `--region`, which is a contig (`chr1`) or an
interval (`chr1:0-10000000`). One invocation reads that region through the index,
processes it, and writes one output shard.

`scripts/regions.sh` produces the list, read from the tabix index so it
reflects what the matrix actually contains:

```bash
scripts/regions.sh --matrix work/join/depth.matrix.txt.gz > regions.txt

# split into smaller units, so losing one to preemption costs less
scripts/regions.sh --matrix work/join/depth.matrix.txt.gz \
                   --window 10000000 --sizes hg38.chrom.sizes > regions.txt
```

A windowed list partitions the matrix: every bin lands in exactly one unit,
with windows aligned up to bin edges (the bin size is read from the first
data row, so bins are assumed uniform), so concatenating the shards yields
each region exactly once. The sizes file may use either contig naming
convention; a contig the matrix has and the sizes file lacks is reported on
stderr rather than silently left out, and no match at all is an error.
Contigs longer than 2²⁹ bp get a `.csi` index instead of `.tbi`.

That one list drives every dispatcher, which is what keeps them interchangeable:

```bash
# SLURM
sbatch --array=1-$(wc -l < regions.txt) workflows/slurm_array.sh regions.txt

# WDL, on either cloud platform
miniwdl run workflows/depthsv.wdl -i inputs.json

# anywhere, locally
parallel -j8 scripts/correct.sh --region {} ... :::: regions.txt
```

None of them is the "real" one — each calls the same stage scripts, and the
pipeline contains no scheduler-specific code. A preempted job costs one region
rather than a whole chromosome, which is what makes preemptible instances
usable.

---

## Failure behaviour

Every stage writes to a temporary file, verifies it, indexes it, and only then
moves it into place and writes a `.done` marker. A stage that dies leaves the
temporary file for inspection and no marker, so re-running redoes exactly that
unit and nothing else. A failed join additionally resumes from its last
finished batch.

A finished unit is reused only for the same inputs and parameters. Each stage
records what it computed — every input's identity (content checksum for the
sample tables, size and mtime for the matrix), every parameter, the model
formula, the driver's checksum — in a `.params` file beside the output, and a
rerun whose signature differs redoes the unit and logs which lines changed. So
a swapped coverage table, an edited model with the same name, a different
`--min-obs`, or a grown manifest never inherits an old result. `--force`
removes the markers first, so a forced run that then fails cannot leave the
previous result behind under a valid marker. Two models of the same response
under a defaulted name are refused rather than taking turns overwriting one
output; give the second one `--name`.

The correction and analysis stages keep each unit's R diagnostics in a `.log`
beside its output; the sample-alignment drop counts there are the first thing
to check when a cohort mismatch is suspected.

The pipeline refuses to proceed on: a sample whose region count disagrees with
the others, a manifest listing the same sample twice, a chromosome absent from
the input, a duplicated `SAMPLE` in the phenotype, PC or coverage table (or in
the matrix header after `--sampleIdPattern`), a missing or non-positive
coverage median, a sex table coded 0/1, a binary response whose coding is
ambiguous, too few cases,
a PC-corrected matrix analysed without its PCs, an association shard with no
rows, an output that fails its integrity check, and a correction stage that
produces fewer rows than it read.

Contig names are resolved from the tabix index rather than assumed, so a matrix
using `1` and one using `chr1` both work and a genuinely absent chromosome is an
error rather than an empty result.

---

## Input formats

| File | Required columns |
|---|---|
| per-sample depth | `mosdepth` `*.regions.bed.gz`: chrom, start, end, mean depth |
| manifest | one absolute path per line |
| principal components | `SAMPLE`, `PC1` … `PCn` |
| coverage | `SAMPLE`, `AUTO_HQ_median` |
| phenotypes | `SAMPLE` plus every term in your models; `SAMPLE` must be unique |
| phenotype manifest | `name`, `method`, `model`, optional `flags`, tab-separated, `#` for comments |
| sex table | `SAMPLE` and a sex column (M/F, male/female or 1/2); may be the phenotype table |
| PAR | BED of the pseudo-autosomal regions (`conf/par.grch38.bed`, `conf/par.grch37.bed`) |

If your matrix column names are not plain sample identifiers, pass
`--sampleIdPattern` after `--` — a PCRE whose first capture group is kept. There is no
default rewriting: a built-in pattern that silently mangles identifiers is worse
than none.

## Output

One tab-separated row per tested region, bgzip-compressed and tabix-indexed:

```
#CHROM  START  END  Region  N  NCase  NControl  <statistics>  MAXSHARE
```

`N` is the number of samples the region was tested on: a sample whose depth
is missing for the region (no expected copies — chrY in a female) is left out
of that region's fit, not imputed. `NCase`/`NControl` are meaningful for
logistic and Cox; for linear both carry the sample count. `MAXSHARE` is the
largest share of the residualised depth's sum of squares carried by one
sample (see `--max-share`).

The statistics between them depend on the method:

| Method | Statistics |
|---|---|
| linear | `BETA`, `SE`, `STAT` (t), `P`, `LOG10P` |
| logistic | `BETA`, `SE`, `STAT` (Wald z), `P`, `LOG10P`, `LRT_P`, `CONVERGED` |
| coxph | `BETA`, `HR`, `SE`, `STAT` (z), `P`, `LOG10P`, `CONVERGED` |

`BETA` is per unit of corrected depth, i.e. per log2 ratio: a heterozygous
deletion is about −1, a duplication about +0.58. `LOG10P` is −log10(P) computed
in log space, so it does not underflow where `P` prints as `0`. For logistic,
`LRT_P` is the likelihood-ratio p-value; under complete separation the Wald
statistic collapses toward zero (Hauck–Donner) while the LRT does not, and
`CONVERGED` says whether the fit converged. Regions dropped by the
quality-control filters are absent rather than reported as missing, so the row
count is normally lower than the region count; a shard with no rows at all is
an error unless `DSV_ALLOW_EMPTY=1`.

---

## Handing off to an external engine

The simple sweep above is first-class. It does not model relatedness.

For analyses that need a relatedness term, calibrated tests for imbalanced
case/control phenotypes, or whole-genome-regression calibration, the corrected
matrix is the seam. What exists today is the GENESIS score test in
[`R/external/`](R/external/), run against a null model pre-fitted with
`GENESIS::fitNullModel()`. A dosage export for REGENIE or SAIGE is designed but
not yet built. Running the same phenotype both ways on the same corrected
matrix is the strongest check either path gets.

See [`PLAN.md`](PLAN.md) for the design and the current state of that work.

---

## Container

```bash
docker build -t depthsv:dev .
docker run --rm depthsv:dev /opt/depthsv/tests/smoke_test.sh /tmp/s
```

The image pins R and a dated CRAN snapshot, so a rebuild resolves the same
package versions rather than whatever is current. It carries no site
configuration and no scheduler assumptions; the stages are separate commands
that a workflow engine calls individually.

Results are byte-identical between the container and a host run.

## Requirements

- R ≥ 4.2 with `optparse` and `data.table`; `survival` only for the Cox path
  (it ships with R)
- `htslib` — `bgzip` and `tabix`
- GNU `parallel`
- `bash` ≥ 3.2, and standard POSIX tools (`awk`, `paste`, `split`, `cut`)

`R.utils` is deliberately *not* required: gzipped input tables are decompressed
through the shell so a base R install is enough.

The scripts pin BLAS/OpenMP linear algebra to one thread per worker —
parallelism comes from the worker fan-out, and a threaded BLAS underneath it
would oversubscribe the node and make results depend on core count. Export the
usual `*_NUM_THREADS` variables yourself to override.

On a cluster, set `DSV_MODULES` to whatever your site calls these — the module
names are not guessed, and the scripts check for the commands regardless.

The external-engine path in [`R/external/`](R/external/) additionally needs
`GENESIS` (Bioconductor). It is not required for the pipeline itself.

Relatedness/GRM construction is out of scope; any standard PLINK2 recipe
produces what the external path expects.

## Testing

```bash
tests/smoke_test.sh          # full pipeline on synthetic data
```

The suite asserts on results rather than exit codes, and needs no cohort data.

## Acknowledgements

The 1000 Genomes example (`example/1000G_highcov/`, including its preamble
for choosing the number of coverage PCs and building genotype-PC covariates)
and the repository's structural and scientific reviews were developed by
[Claude](https://claude.ai) (Anthropic; Claude Fable 5) working with the
maintainer in Claude Code. Those commits carry a `Co-Authored-By` trailer.
The design decisions, and any mistakes, were reviewed and accepted by the
maintainer.

## Licence

MIT — see [LICENSE](LICENSE).
