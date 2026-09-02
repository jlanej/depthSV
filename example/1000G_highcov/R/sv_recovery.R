#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# depthSV 1000G example — recovery of known deletions as a function of ndim
#
#   --stage select    pick the deletions: autosomal, on the matrix's contigs,
#                     within the length and carrier-frequency bands, spread
#                     over chromosomes; writes selected.tsv and units.txt
#   --stage evaluate  per deletion and ndim, the corrected depth (mean over
#                     the bins inside the deletion) of carriers against
#                     non-carriers: AUC, log2 shift, Welch t; a per-ndim
#                     summary, a plot, and the recommended ndim
#
# The recommendation is the smallest ndim whose median AUC is within 0.01
# of the best: the start of the plateau, before the correction begins to
# absorb the deletions themselves.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

option_list <- list(
  make_option("--stage",     type = "character", help = "select | evaluate"),
  make_option("--calls",     type = "character", default = NULL, help = "CHROM START END ID N_CARRIERS CARRIERS"),
  make_option("--contigs",   type = "character", default = "", help = "comma-separated contigs of the matrix"),
  make_option("--samples",   type = "character", default = NULL, help = "samples.txt: the matrix columns"),
  make_option("--max-dels",  type = "integer",   dest = "max_dels", default = 200L),
  make_option("--min-len",   type = "integer",   dest = "min_len",  default = 5000L),
  make_option("--min-af",    type = "double",    dest = "min_af",   default = 0.02),
  make_option("--max-af",    type = "double",    dest = "max_af",   default = 0.5),
  make_option("--selected",  type = "character", default = NULL),
  make_option("--corrected", type = "character", default = NULL, help = "directory of corrected_ndim<k>.<slug>.txt.gz"),
  make_option("--ndims",     type = "character", default = ""),
  make_option("--mp",        type = "character", default = "NA", help = "the preamble's Marchenko-Pastur count"),
  make_option("--ndim",      type = "character", default = "NA", help = "the run's EX_NDIM"),
  make_option("--mode",      type = "character", default = "standard"),
  make_option("--smoke",     type = "character", default = "0"),
  make_option("--out",       type = "character", help = "output directory")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$stage) || !opt$stage %in% c("select", "evaluate")) stop("--stage must be select or evaluate", call. = FALSE)
if (is.null(opt$out)) stop("--out is required", call. = FALSE)
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)

slug_of <- function(region) gsub(":", "_", region, fixed = TRUE)

# --- select ---------------------------------------------------------------------

if (opt$stage == "select") {
  for (req in c("calls", "samples")) if (is.null(opt[[req]])) stop("--", req, " is required", call. = FALSE)
  calls <- fread(opt$calls, colClasses = list(character = c("CHROM", "ID", "CARRIERS")))
  for (col in c("CHROM", "START", "END", "ID", "N_CARRIERS", "CARRIERS")) {
    if (!col %in% names(calls)) stop(opt$calls, " lacks a ", col, " column", call. = FALSE)
  }
  samples <- trimws(readLines(opt$samples)); samples <- samples[nzchar(samples)]
  contigs <- strsplit(opt$contigs, ",", fixed = TRUE)[[1]]
  calls[, LEN := END - START]
  calls[, AUTOSOME := grepl("^(chr)?[0-9]+$", CHROM)]
  # Carriers among the matrix's samples: the callset covers more samples than
  # a smoke tree simulates, and the frequency band is judged in-matrix.
  calls[, N_IN := vapply(strsplit(CARRIERS, ",", fixed = TRUE), function(x) sum(x %in% samples), integer(1))]
  n_s <- length(samples)
  sel <- calls[AUTOSOME & CHROM %in% contigs & LEN >= opt$min_len & LEN <= 5e6 &
               N_IN >= 2 & N_IN / n_s >= opt$min_af & N_IN / n_s <= opt$max_af]
  message(sprintf("[select] %d of %d deletions are autosomal, on the matrix, %d bp-5 Mb long, with %.0f-%.0f%% carriers among %d samples",
                  nrow(sel), nrow(calls), opt$min_len, 100 * opt$min_af, 100 * opt$max_af, n_s))
  if (nrow(sel)) {
    # Spread over chromosomes: the longest of each chromosome first, round robin.
    setorder(sel, CHROM, -LEN)
    sel[, rank_in_chrom := seq_len(.N), by = CHROM]
    setorder(sel, rank_in_chrom, -LEN)
    sel <- sel[seq_len(min(opt$max_dels, nrow(sel)))]
    sel[, rank_in_chrom := NULL]
    setorder(sel, CHROM, START)
  }
  fwrite(sel[, .(CHROM, START, END, ID, LEN, N_CARRIERS, N_IN, CARRIERS)], file.path(opt$out, "selected.tsv"), sep = "\t")
  writeLines(sprintf("%s:%d-%d", sel$CHROM, sel$START + 1L, sel$END), file.path(opt$out, "units.txt"))
  message(sprintf("[select] %d deletions -> selected.tsv", nrow(sel)))
  quit(status = 0)
}

# --- evaluate ---------------------------------------------------------------------

for (req in c("selected", "corrected")) if (is.null(opt[[req]])) stop("--", req, " is required", call. = FALSE)
sel <- fread(opt$selected, colClasses = list(character = c("CHROM", "ID", "CARRIERS")))
ndims <- as.integer(strsplit(trimws(opt$ndims), "[ ,]+")[[1]])
ndims <- sort(unique(ndims[is.finite(ndims)]))
if (!length(ndims)) stop("--ndims is empty", call. = FALSE)
mp       <- suppressWarnings(as.integer(opt$mp))
run_ndim <- suppressWarnings(as.integer(opt$ndim))

auc_lower <- function(x_car, x_non) {
  # P(carrier value < non-carrier value): 1 when every carrier sits below.
  x_car <- x_car[is.finite(x_car)]; x_non <- x_non[is.finite(x_non)]
  if (!length(x_car) || !length(x_non)) return(NA_real_)
  r <- rank(c(x_car, x_non))
  u <- sum(r[seq_along(x_car)]) - length(x_car) * (length(x_car) + 1) / 2
  1 - u / (length(x_car) * length(x_non))
}

rows <- list()
for (i in seq_len(nrow(sel))) {
  region <- sprintf("%s:%d-%d", sel$CHROM[i], sel$START[i] + 1L, sel$END[i])
  carriers <- strsplit(sel$CARRIERS[i], ",", fixed = TRUE)[[1]]
  for (k in ndims) {
    f <- file.path(opt$corrected, sprintf("corrected_ndim%d.%s.txt.gz", k, slug_of(region)))
    if (!file.exists(f)) { message(sprintf("[warn] missing %s", basename(f))); next }
    d <- fread(cmd = paste("gzip -cd", shQuote(f)), na.strings = c("NA", ""))
    if (!nrow(d)) next
    setnames(d, 1:4, c("CHROM", "START", "END", "Region"))
    inside <- d[START >= sel$START[i] & END <= sel$END[i]]
    if (!nrow(inside)) next
    M <- as.matrix(inside[, -(1:4)]); storage.mode(M) <- "double"
    x <- colMeans(M, na.rm = TRUE)
    is_car <- names(x) %in% carriers
    x_car <- x[is_car]; x_non <- x[!is_car]
    if (length(x_car) < 2 || length(x_non) < 2) next
    tt <- tryCatch(t.test(x_car, x_non)$statistic, error = function(e) NA_real_)
    rows[[length(rows) + 1L]] <- data.table(
      ndim = k, ID = sel$ID[i], CHROM = sel$CHROM[i], START = sel$START[i], END = sel$END[i],
      LEN = sel$END[i] - sel$START[i], n_bins = nrow(inside), n_carriers = length(x_car),
      n_noncarriers = length(x_non), auc = auc_lower(x_car, x_non),
      shift = mean(x_car) - mean(x_non), t = as.numeric(tt))
  }
}
res <- rbindlist(rows)
if (!nrow(res)) stop("no deletion could be evaluated (no corrected units found under ", opt$corrected, ")", call. = FALSE)
fwrite(res, file.path(opt$out, "sv_recovery.tsv"), sep = "\t")

summ <- res[, .(n_deletions = .N,
                median_auc = median(auc, na.rm = TRUE),
                frac_auc_ge_0.95 = mean(auc >= 0.95, na.rm = TRUE),
                median_shift = median(shift, na.rm = TRUE),
                median_abs_t = median(abs(t), na.rm = TRUE)), by = ndim][order(ndim)]
summ[, is_mp := ndim %in% mp]
summ[, is_run := ndim %in% run_ndim]
fwrite(summ, file.path(opt$out, "sv_recovery_summary.tsv"), sep = "\t")

best <- max(summ$median_auc, na.rm = TRUE)
rec  <- summ[median_auc >= best - 0.01][order(ndim)]$ndim[1]
writeLines(as.character(rec), file.path(opt$out, "recommended_ndim.txt"))
# Row-wise tags (by = ndim makes each group one row, so the ifs see scalars).
summ[, tag := paste(c(if (is_mp) "MP", if (is_run) "run", if (ndim == rec) "recommended"), collapse = ", "), by = ndim]

# --- plot -------------------------------------------------------------------------
ok_png <- tryCatch({
  png(file.path(opt$out, "sv_recovery.png"), width = 1400, height = 600, res = 130); TRUE
}, error = function(e) FALSE)
if (!ok_png) pdf(file.path(opt$out, "sv_recovery.pdf"), width = 11, height = 4.5)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
xs <- seq_along(summ$ndim)
plot(xs, summ$median_auc, type = "b", pch = 19, xaxt = "n", xlab = "coverage PCs removed (ndim)",
     ylab = "median AUC (carriers below non-carriers)", ylim = c(min(0.5, min(summ$median_auc, na.rm = TRUE)), 1),
     main = sprintf("deletion recovery, %s mode (%d deletions)", opt$mode, max(summ$n_deletions)))
axis(1, at = xs, labels = summ$ndim)
if (is.finite(mp) && mp %in% summ$ndim) abline(v = match(mp, summ$ndim), lty = 2, col = "grey40")
if (is.finite(run_ndim) && run_ndim %in% summ$ndim) abline(v = match(run_ndim, summ$ndim), lty = 3, col = "firebrick")
abline(v = match(rec, summ$ndim), lty = 1, col = "steelblue")
legend("bottomleft", bty = "n", lty = c(2, 3, 1), col = c("grey40", "firebrick", "steelblue"),
       legend = c(sprintf("MP count (%s)", opt$mp), sprintf("run ndim (%s)", opt$ndim), sprintf("recommended (%d)", rec)))
plot(xs, summ$median_shift, type = "b", pch = 19, xaxt = "n", xlab = "coverage PCs removed (ndim)",
     ylab = "median carrier shift (log2)", main = "effect size retained (one lost copy is about -1)")
axis(1, at = xs, labels = summ$ndim)
abline(h = -1, lty = 2, col = "grey60")
dev.off()

# --- summary ------------------------------------------------------------------------
md <- c(sprintf("# SV-callset recovery — %s mode", opt$mode), "",
        if (opt$smoke == "1") "**SMOKE**: the deletions are the ones the simulated tree carries; the curve tests the machinery, not the callset." else NULL,
        "",
        sprintf("- deletions evaluated: %d (selected.tsv); ndims: %s", max(summ$n_deletions), paste(summ$ndim, collapse = ", ")),
        sprintf("- Marchenko-Pastur count: %s; this run's ndim: %s", opt$mp, opt$ndim),
        sprintf("- **recommended ndim: %d** — the smallest ndim within 0.01 of the best median AUC (%.3f). Informational; set EX_NDIM to adopt it.", rec, best),
        "",
        "| ndim | deletions | median AUC | AUC >= 0.95 | median shift (log2) | median abs t | |",
        "|---|---|---|---|---|---|---|",
        summ[, sprintf("| %d | %d | %.3f | %.0f%% | %.3f | %.1f | %s |", ndim, n_deletions, median_auc,
                       100 * frac_auc_ge_0.95, median_shift, median_abs_t, tag)],
        "",
        "AUC is the probability that a carrier's corrected depth over the deletion sits below a",
        "non-carrier's; the shift is the carrier mean minus the non-carrier mean in log2 units (one",
        "lost copy is about -1). A curve that rises then falls shows the correction absorbing the",
        "deletions themselves; the recommendation is the start of the plateau.")
writeLines(md, file.path(opt$out, "summary.md"))
message(sprintf("[sv] %d deletions x %d ndims; best median AUC %.3f; recommended ndim %d -> %s",
                max(summ$n_deletions), length(ndims), best, rec, file.path(opt$out, "summary.md")))
