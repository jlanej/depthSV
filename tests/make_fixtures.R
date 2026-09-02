#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# depthSV — synthetic fixture generator
#
# Produces a small, fully synthetic dataset in the exact formats the pipeline
# reads, so the test suite (and the worked example in the README) never needs
# access to real cohort data.
#
# A known association is injected at one region on each autosome so tests can
# assert that the pipeline recovers signal, not merely that it runs; a
# phenotype driven by the leading PC is included so the suite can tell whether
# the association test conditions on the PCs the correction removed. The sex
# chromosomes carry their real ploidy (chrX one copy in males outside a small
# pseudo-autosomal region, chrY one copy in males and none in females), one
# autosomal bin has zero depth in a few samples, and one has a single-sample
# outlier, so the ploidy model, the winsor and the leverage filter are all
# exercised on data whose truth is known.
#
# The random draws are ordered so that the autosomes and phenotypes come
# first and the sex chromosomes last: the autosomal fixture is the same
# realisation it was before the sex chromosomes were added.
#
#   Rscript tests/make_fixtures.R [outDir] [nSamples] [nRegionsPerChrom]
# ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(data.table))

args     <- commandArgs(trailingOnly = TRUE)
out_dir  <- if (length(args) >= 1) args[1] else "tests/fixtures"
n_sample <- if (length(args) >= 2) as.integer(args[2]) else 60L
n_region <- if (length(args) >= 3) as.integer(args[3]) else 200L
if (n_region < 125L) stop("nRegionsPerChrom must be at least 125 (the special bins sit at regions 50-120)", call. = FALSE)

set.seed(20260812)   # fixtures must be byte-reproducible

autosomes <- c("chr1", "chr2")
sex_chroms <- c("chrX", "chrY")
bin_size <- 1000L
par_end  <- 20L * bin_size       # chrX:0-20000 is pseudo-autosomal here
dir.create(file.path(out_dir, "mosdepth"), recursive = TRUE, showWarnings = FALSE)

samples <- sprintf("SAMPLE%03d", seq_len(n_sample))

# --- per-sample technical structure ----------------------------------------
# Sequencing-centre batch plus a per-sample depth scale, so the PC correction
# has something real to remove.
batch      <- factor(sample(c("siteA", "siteB", "siteC"), n_sample, replace = TRUE))
batch_load <- matrix(rnorm(n_sample * 3), n_sample, 3)
median_cov <- round(pmax(5, rnorm(n_sample, mean = 30, sd = 4)), 3)

make_regions <- function(chroms) rbindlist(lapply(chroms, function(cc) {
  data.table(CHR = cc,
             START = seq(0L, by = bin_size, length.out = n_region),
             STOP  = seq(bin_size, by = bin_size, length.out = n_region))
}))
regions_a <- make_regions(autosomes)
n_auto    <- nrow(regions_a)

# Region-level sensitivity to the batch factors.
region_load <- matrix(rnorm(n_auto * 3, sd = 0.35), n_auto, 3)

# --- the injected signal ---------------------------------------------------
# One carrier group with a real deletion at a designated region per autosome.
signal_idx <- c(which(regions_a$CHR == "chr1")[50], which(regions_a$CHR == "chr2")[75])
carrier    <- rbinom(n_sample, 1, 0.20)

# --- autosomal log2 ratios -------------------------------------------------
log2ratio_a <- tcrossprod(region_load, batch_load) + matrix(rnorm(n_auto * n_sample, sd = 0.18),
                                                            n_auto, n_sample)
for (r in signal_idx) log2ratio_a[r, ] <- log2ratio_a[r, ] - 1.0 * carrier

# --- principal components --------------------------------------------------
# Derived from the batch loadings so they genuinely capture the technical
# structure, with a little noise in the trailing components.
n_pc <- 10L
pcs  <- cbind(batch_load, matrix(rnorm(n_sample * (n_pc - 3)), n_sample, n_pc - 3))
colnames(pcs) <- paste0("PC", seq_len(n_pc))
fwrite(data.table(SAMPLE = samples, as.data.table(pcs)),
       file.path(out_dir, "svd.pcs.txt"), sep = "\t")

# --- phenotypes ------------------------------------------------------------
# quant_trait is genuinely caused by carrier status; null_trait is not, so the
# suite can check both power and calibration. pc_null is driven by PC1 — the
# structure the correction removes from depth — and by nothing else: a test
# that fails to condition on the PCs is deflated on it (lambda ~ 0.2 here).
age <- round(rnorm(n_sample, 55, 9), 1)
sex <- sample(c("M", "F"), n_sample, replace = TRUE)
case_status <- rbinom(n_sample, 1, plogis(-1.1 + 1.4 * carrier))

pheno <- data.table(
  SAMPLE      = samples,
  quant_trait = round(2.0 * carrier + 0.02 * age + rnorm(n_sample), 4),
  null_trait  = round(rnorm(n_sample), 4),
  pc_null     = round(2 * pcs[, 1] + rnorm(n_sample), 4),
  case_status = case_status,
  # The same binary phenotype as text, so the explicit case level is exercised.
  case_label  = ifelse(case_status == 1L, "case", "control"),
  # Time-to-event outcome, so the coxph path is exercised too. Carriers fail
  # faster, matching the direction of the injected deletion.
  time        = round(rexp(n_sample, rate = 0.04 * exp(0.8 * carrier)), 2),
  event       = rbinom(n_sample, 1, 0.7),
  age         = age,
  sex         = sex,
  batch       = as.character(batch)
)
for (k in 1:4) pheno[[paste0("ancestry_PC", k)]] <- round(rnorm(n_sample), 4)
fwrite(pheno, file.path(out_dir, "phenotypes.tsv"), sep = "\t")

# --- sex chromosomes -------------------------------------------------------
# Copies relative to the diploid autosomes: chrX outside the PAR is one copy
# in males; chrY is one copy in males and absent in females, whose chrY depth
# is mapping noise at ~1% of the median.
regions_s <- make_regions(sex_chroms)
n_sex     <- nrow(regions_s)
log2ratio_s <- tcrossprod(matrix(rnorm(n_sex * 3, sd = 0.35), n_sex, 3), batch_load) +
  matrix(rnorm(n_sex * n_sample, sd = 0.18), n_sex, n_sample)
x_nonpar <- regions_s$CHR == "chrX" & regions_s$START >= par_end
y_rows   <- regions_s$CHR == "chrY"
male     <- sex == "M"
log2ratio_s[x_nonpar, male]  <- log2ratio_s[x_nonpar, male] - 1
log2ratio_s[y_rows,   male]  <- log2ratio_s[y_rows,   male] - 1
log2ratio_s[y_rows,  !male]  <- log2ratio_s[y_rows,  !male] + log2(0.01)

regions   <- rbind(regions_a, regions_s)
log2ratio <- rbind(log2ratio_a, log2ratio_s)
n_total   <- nrow(regions)

# --- depth matrix ----------------------------------------------------------
# log2 ratio -> linear depth, scaled by each sample's own median coverage.
depth <- sweep(2^log2ratio, 2, median_cov, `*`)

# A zero-depth bin (three samples) and a single-sample outlier bin, both on
# chr1 away from the injected signal.
zero_idx    <- which(regions$CHR == "chr1")[100]
outlier_idx <- which(regions$CHR == "chr1")[120]
depth[zero_idx, 1:3]  <- 0
depth[outlier_idx, 5] <- depth[outlier_idx, 5] * 24

depth[depth < 0] <- 0
depth <- round(depth, 2)

# --- per-sample mosdepth regions files -------------------------------------
for (j in seq_len(n_sample)) {
  dt <- data.table(regions$CHR, regions$START, regions$STOP, depth[, j])
  f  <- file.path(out_dir, "mosdepth", sprintf("%s.by1000.regions.bed", samples[j]))
  fwrite(dt, f, sep = "\t", col.names = FALSE)
  system2("bgzip", c("-f", shQuote(f)))
}

manifest <- file.path(out_dir, "mosdepth.input.txt")
writeLines(file.path(normalizePath(file.path(out_dir, "mosdepth")),
                     sprintf("%s.by1000.regions.bed.gz", samples)), manifest)

# --- coverage stats --------------------------------------------------------
fwrite(data.table(SAMPLE = samples, AUTO_HQ_median = median_cov),
       file.path(out_dir, "autosomal.median.txt"), sep = "\t")

# --- pseudo-autosomal regions (0-based BED, like conf/par.grch38.bed) -----
writeLines(sprintf("chrX\t0\t%d\tPAR1", par_end), file.path(out_dir, "par.bed"))

# --- truth file ------------------------------------------------------------
fwrite(data.table(CHROM  = regions$CHR[signal_idx],
                  START  = regions$START[signal_idx],
                  END    = regions$STOP[signal_idx],
                  Region = sprintf("%s:%d-%d", regions$CHR[signal_idx],
                                   regions$START[signal_idx], regions$STOP[signal_idx]),
                  effect = "deletion in carriers; associated with quant_trait and case_status"),
       file.path(out_dir, "truth.tsv"), sep = "\t")

cat(sprintf(paste0("fixtures written to %s\n  %d samples (%d male) x %d regions over %s\n",
                   "  injected signal at: %s\n  zero-depth bin: chr1:%d  outlier bin: chr1:%d  PAR: chrX:0-%d\n"),
            out_dir, n_sample, sum(male), n_total, paste(c(autosomes, sex_chroms), collapse = ", "),
            paste(sprintf("%s:%d", regions$CHR[signal_idx], regions$START[signal_idx]), collapse = ", "),
            regions$START[zero_idx], regions$START[outlier_idx], par_end))
