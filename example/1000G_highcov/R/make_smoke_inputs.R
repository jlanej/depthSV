#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# depthSV 1000G example — simulate a small mosdepth tree from a real QC table
#
# Smoke mode only. Each simulated sample's depths are anchored to its actual
# QC row from the upstream NGS-PCA run — autosomal bins at HQ_MEDIAN_COV,
# chrX/chrY scaled by the observed coverage ratios, chrM at the observed
# mitochondrial ratio — so the example's positive controls (chrM for
# MTDNA_CN, the sex chromosomes for INFERRED_SEX) hold by construction and
# the whole pipeline can be exercised without CRAMs or a cluster.
#
# The output exercises the machinery; it is not data.
#
#   Rscript make_smoke_inputs.R --qc sample_qc.tsv --out DIR \
#       --samples 64 --seed 20260818 --jitter 0
#
# --jitter > 0 adds a deterministic multiplicative perturbation (used for the
# synthetic fast tree when the real fast outputs are unavailable).
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

option_list <- list(
  make_option("--qc",      type = "character", help = "sample_qc.tsv from the NGS-PCA example"),
  make_option("--out",     type = "character", help = "output directory for *.by1000.regions.bed.gz"),
  make_option("--samples", type = "integer", default = 64L),
  make_option("--seed",    type = "integer", default = 20260818L),
  make_option("--jitter",  type = "double",  default = 0)
)
opt <- parse_args(OptionParser(option_list = option_list))
for (req in c("qc", "out")) if (is.null(opt[[req]])) stop("--", req, " is required", call. = FALSE)

needed <- c("SAMPLE_ID", "HQ_MEDIAN_COV", "X_COV_RATIO", "Y_COV_RATIO",
            "MITO_COV_RATIO", "MEAN_AUTOSOMAL_COV", "INFERRED_SEX")
qc <- fread(opt$qc, select = needed)
qc <- qc[complete.cases(qc) & HQ_MEDIAN_COV > 0 & INFERRED_SEX %in% c("M", "F")]
if (!nrow(qc)) stop("no usable rows in ", opt$qc, call. = FALSE)
if (nrow(qc) < opt$samples) {
  message(sprintf("[warn] only %d usable samples (requested %d)", nrow(qc), opt$samples))
}
qc <- qc[seq_len(min(opt$samples, nrow(qc)))]
n  <- nrow(qc)

# --- bins ------------------------------------------------------------------
# A deliberately small genome: one autosomal slice, non-PAR slices of X and
# Y, and all of chrM. 1 kb bins to match the upstream convention (and the
# .by1000 suffix the sample-name stripping expects).

bin <- 1000L
mk <- function(chrom, from, to) {
  s <- seq(from, to - 1L, by = bin)
  data.table(CHR = chrom, START = s, STOP = pmin(s + bin, to))
}
regions <- rbindlist(list(
  mk("chr20", 0L,        2000000L),   # autosomal null background
  mk("chrX",  3000000L,  5000000L),   # non-PAR: males ~0.5x
  mk("chrY",  3000000L,  4000000L),   # non-PAR: females ~0
  mk("chrM",  0L,        16569L)      # the MTDNA_CN numerator
))
n_bins <- nrow(regions)

# --- depth matrix ----------------------------------------------------------

set.seed(opt$seed)

target <- matrix(0, n_bins, n)
auto_rows <- regions$CHR == "chr20"
x_rows    <- regions$CHR == "chrX"
y_rows    <- regions$CHR == "chrY"
m_rows    <- regions$CHR == "chrM"
# The sex chromosomes follow the inferred sex exactly — males at half the
# autosomal median on X and Y, females at the median on X and ~1% on Y —
# rather than the observed X/Y ratios, so the ploidy model is exact on this
# tree and any sex signal left on chrX after correction is a wiring fault.
male <- as.character(qc$INFERRED_SEX) == "M"
for (j in seq_len(n)) {
  target[auto_rows, j] <- qc$HQ_MEDIAN_COV[j]
  target[x_rows,    j] <- qc$HQ_MEDIAN_COV[j] * if (male[j]) 0.5 else 1.0
  target[y_rows,    j] <- qc$HQ_MEDIAN_COV[j] * if (male[j]) 0.5 else 0.01
  target[m_rows,    j] <- qc$MEAN_AUTOSOMAL_COV[j] * qc$MITO_COV_RATIO[j]
}

# A shared per-bin profile (GC-like) plus independent sampling noise. The
# noise is kept below the cohort's log2(MTDNA_CN) spread (~0.2) so the chrM
# signal survives the PC correction's absorption of that phenotype.
bin_effect <- exp(rnorm(n_bins, sd = 0.08))
depth <- target * bin_effect * matrix(exp(rnorm(n_bins * n, sd = 0.035)), n_bins, n)

# Known deletions on the autosomal slice, for the SV-recovery stage: twelve
# deletions of 5-30 kb with carrier frequencies from 5% to 40%, one copy
# lost in carriers (a fifth of them lose both). Drawn after everything
# above, so the rest of the tree is unchanged by their presence; the same
# seed gives the fast tree the same carriers.
set.seed(opt$seed + 2000L)
auto_start <- min(regions$START[auto_rows]); auto_end <- max(regions$STOP[auto_rows])
del_len  <- c(5, 8, 10, 15, 20, 30, 5, 10, 15, 20, 8, 12) * 1000L
del_freq <- c(0.05, 0.10, 0.20, 0.30, 0.40, 0.10, 0.20, 0.05, 0.30, 0.15, 0.40, 0.25)
slot <- (auto_end - auto_start) %/% length(del_len)
calls <- list()
for (i in seq_along(del_len)) {
  s <- as.integer((auto_start + (i - 1L) * slot + 10000L) %/% bin * bin)   # inside its slot, bin-aligned
  e <- s + del_len[i]
  rows <- which(auto_rows & regions$START >= s & regions$STOP <= e)
  n_car <- max(2L, round(del_freq[i] * n))
  car <- sort(sample.int(n, n_car))
  hom <- car[seq_len(round(0.2 * n_car))]
  depth[rows, car] <- depth[rows, car] * 0.5
  if (length(hom)) depth[rows, hom] <- depth[rows, hom] * 0.04
  calls[[i]] <- data.table(CHROM = "chr20", START = s, END = e, ID = sprintf("smokeDEL%02d", i),
                           N_CARRIERS = n_car, CARRIERS = paste(qc$SAMPLE_ID[car], collapse = ","))
}
calls <- rbindlist(calls)

if (opt$jitter > 0) {
  set.seed(opt$seed + 1000L)
  depth <- depth * exp(matrix(rnorm(n_bins * n, sd = opt$jitter), n_bins, n)) * 1.002
}
depth <- round(depth, 2)

# --- write -----------------------------------------------------------------

dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
for (j in seq_len(n)) {
  f <- file.path(opt$out, sprintf("%s.by1000.regions.bed", qc$SAMPLE_ID[j]))
  fwrite(data.table(regions$CHR, regions$START, regions$STOP, depth[, j]),
         f, sep = "\t", col.names = FALSE)
  if (system2("bgzip", c("-f", shQuote(f))) != 0L) stop("bgzip failed on ", f, call. = FALSE)
}

fwrite(calls, file.path(opt$out, "sv_calls.tsv"), sep = "\t")

writeLines(c(sprintf("generated\t%s", format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")),
             sprintf("qc\t%s", opt$qc),
             sprintf("samples\t%d", n),
             sprintf("bins\t%d", n_bins),
             sprintf("seed\t%d", opt$seed),
             sprintf("jitter\t%g", opt$jitter),
             "note\tsimulated smoke data - exercises the machinery, not results"),
           file.path(opt$out, "smoke.params.txt"))

cat(sprintf("simulated %d samples x %d bins into %s\n", n, n_bins, opt$out))
