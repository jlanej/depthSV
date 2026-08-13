#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# depthSV — synthetic fixture generator
#
# Produces a small, fully synthetic dataset in the exact formats the pipeline
# reads, so the test suite (and the worked example in the README) never needs
# access to real cohort data.
#
# A known association is injected at one region on each chromosome so tests can
# assert that the pipeline recovers signal, not merely that it runs.
#
#   Rscript tests/make_fixtures.R [outDir] [nSamples] [nRegionsPerChrom]
# ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(data.table))

args     <- commandArgs(trailingOnly = TRUE)
out_dir  <- if (length(args) >= 1) args[1] else "tests/fixtures"
n_sample <- if (length(args) >= 2) as.integer(args[2]) else 60L
n_region <- if (length(args) >= 3) as.integer(args[3]) else 200L

set.seed(20260812)   # fixtures must be byte-reproducible

chroms   <- c("chr1", "chr2")
bin_size <- 1000L
dir.create(file.path(out_dir, "mosdepth"), recursive = TRUE, showWarnings = FALSE)

samples <- sprintf("SAMPLE%03d", seq_len(n_sample))

# --- per-sample technical structure ----------------------------------------
# Sequencing-centre batch plus a per-sample depth scale, so the PC correction
# has something real to remove.
batch      <- factor(sample(c("siteA", "siteB", "siteC"), n_sample, replace = TRUE))
batch_load <- matrix(rnorm(n_sample * 3), n_sample, 3)
median_cov <- round(pmax(5, rnorm(n_sample, mean = 30, sd = 4)), 3)

regions <- rbindlist(lapply(chroms, function(cc) {
  data.table(CHR = cc,
             START = seq(0L, by = bin_size, length.out = n_region),
             STOP  = seq(bin_size, by = bin_size, length.out = n_region))
}))
n_total <- nrow(regions)

# Region-level sensitivity to the batch factors.
region_load <- matrix(rnorm(n_total * 3, sd = 0.35), n_total, 3)

# --- the injected signal ---------------------------------------------------
# One carrier group with a real deletion at a designated region per chromosome.
signal_idx <- c(which(regions$CHR == "chr1")[50], which(regions$CHR == "chr2")[75])
carrier    <- rbinom(n_sample, 1, 0.20)

# --- depth matrix ----------------------------------------------------------
# log2 ratio -> linear depth, scaled by each sample's own median coverage.
log2ratio <- tcrossprod(region_load, batch_load) + matrix(rnorm(n_total * n_sample, sd = 0.18),
                                                          n_total, n_sample)
for (r in signal_idx) log2ratio[r, ] <- log2ratio[r, ] - 1.0 * carrier

depth <- sweep(2^log2ratio, 2, median_cov, `*`)
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
# suite can check both power and calibration.
age <- round(rnorm(n_sample, 55, 9), 1)
sex <- sample(c("M", "F"), n_sample, replace = TRUE)

pheno <- data.table(
  SAMPLE      = samples,
  quant_trait = round(1.1 * carrier + 0.02 * age + rnorm(n_sample), 4),
  null_trait  = round(rnorm(n_sample), 4),
  case_status = rbinom(n_sample, 1, plogis(-1.1 + 1.4 * carrier)),
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

# --- truth file ------------------------------------------------------------
fwrite(data.table(CHROM  = regions$CHR[signal_idx],
                  START  = regions$START[signal_idx],
                  END    = regions$STOP[signal_idx],
                  Region = sprintf("%s:%d-%d", regions$CHR[signal_idx],
                                   regions$START[signal_idx], regions$STOP[signal_idx]),
                  effect = "deletion in carriers; associated with quant_trait and case_status"),
       file.path(out_dir, "truth.tsv"), sep = "\t")

cat(sprintf("fixtures written to %s\n  %d samples x %d regions over %s\n  injected signal at: %s\n",
            out_dir, n_sample, n_total, paste(chroms, collapse = ", "),
            paste(sprintf("%s:%d", regions$CHR[signal_idx], regions$START[signal_idx]),
                  collapse = ", ")))
