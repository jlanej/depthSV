#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# depthSV 1000G example — how many coverage PCs to remove (Marchenko-Pastur)
#
# NGS-PCA reports the top singular values of the median-centred log2-ratio
# matrix (samples x bins). For a pure-noise matrix with n rows, p columns
# and entry variance sigma^2, those values follow the Marchenko-Pastur law:
# the bulk lies in sigma * [sqrt(p) - sqrt(n), sqrt(p) + sqrt(n)], and a
# component is only distinguishable from noise above the upper edge
# sigma * (sqrt(p) + sqrt(n)). At 3,202 x 142,070 that band is narrow
# (edge / centre ~ 1.15), so the trailing reported singular values sit on a
# nearly flat plateau that pins sigma down well.
#
# sigma is fitted by matching the trailing reported ranks to the MP
# quantiles of the corresponding ranks (only ranks past the current signal
# estimate plus a gap), the count above the edge is recomputed, and the two
# are iterated to a fixed point. Real coverage noise is not iid — per-bin
# variance differs — so the observed plateau decays a little faster than MP
# and the exact-edge count runs high; a margin above the edge (default 1%,
# about four Tracy-Widom standard deviations at this size) is the sensible
# operating point, and the count at 0/1/2/5% is reported so the sensitivity
# is visible rather than hidden in a default.
#
# One number is chosen for every mode — the rounded mean of the determined
# per-mode counts — so the fast-vs-standard comparison corrects both trees
# identically.
#
#   Rscript choose_ndim.R --runs standard=/path/ngspca_output,fast=/path/ngspca_output_fast \
#       --margin 0.01 --gap 20 --out DIR
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

option_list <- list(
  make_option("--runs",   type = "character", help = "comma-separated label=dir, each dir an NGS-PCA output"),
  make_option("--margin", type = "double",  default = 0.01, help = "relative margin above the MP edge [default %default]"),
  make_option("--gap",    type = "integer", default = 20L, help = "ranks skipped between the signal estimate and the sigma fit [default %default]"),
  make_option("--out",    type = "character", help = "output directory")
)
opt <- parse_args(OptionParser(option_list = option_list))
for (req in c("runs", "out")) if (is.null(opt[[req]])) stop("--", req, " is required", call. = FALSE)
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)

count_lines <- function(f) {
  n <- 0L; con <- file(f, "r"); on.exit(close(con))
  while (length(chunk <- readLines(con, n = 100000L))) n <- n + length(chunk)
  n
}

# MP quantiles for the eigenvalues of (1/p) X X', unit variance, ratio gamma = n/p.
mp_quantile_fn <- function(gamma) {
  a <- (1 - sqrt(gamma))^2; b <- (1 + sqrt(gamma))^2
  x <- seq(a, b, length.out = 20001)
  dens <- sqrt(pmax(0, (b - x) * (x - a))) / (2 * pi * gamma * x)
  cdf <- cumsum(c(0, (dens[-1] + dens[-length(dens)]) / 2 * diff(x)))
  cdf <- cdf / max(cdf)
  function(q) approx(cdf, x, xout = q, ties = "ordered")$y
}

fit_one <- function(sv, n, p, margin, gap) {
  k_rep <- length(sv)
  gamma <- n / p
  qf <- mp_quantile_fn(gamma)
  # rank j of n bulk eigenvalues sits at CDF 1 - (j - 0.5)/n
  s_unit <- sqrt(p * qf(1 - (seq_len(k_rep) - 0.5) / n))   # noise singular value per rank, sigma = 1
  edge_unit <- sqrt(p) + sqrt(n)
  step <- function(k, m) {
    # Fewer than ten ranks left past the gap: the reported spectrum does not
    # reach the bulk (or the gap is too wide), so sigma cannot be fitted.
    if (k + gap + 1 > k_rep - 9) return(NULL)
    ranks <- (k + gap + 1):k_rep
    sigma <- exp(mean(log(sv[ranks] / s_unit[ranks])))
    edge <- sigma * edge_unit * (1 + m)
    list(sigma = sigma, edge = edge, k = sum(sv > edge))
  }
  solve_k <- function(m) {
    k <- 0L
    for (i in 1:100) {
      f <- step(k, m)
      if (is.null(f)) return(list(k = NA_integer_, sigma = NA_real_, edge = NA_real_, status = "undetermined"))
      if (f$k == k) return(c(f, status = "ok"))
      k <- f$k
    }
    list(k = NA_integer_, sigma = NA_real_, edge = NA_real_, status = "no-convergence")
  }
  chosen <- solve_k(margin)
  sens <- sapply(c(0, 0.01, 0.02, 0.05), function(m) solve_k(m)$k)
  list(n = n, p = p, k_reported = k_rep, gamma = gamma, s_unit = s_unit,
       sigma = chosen$sigma, edge = chosen$edge, k = chosen$k, status = chosen$status,
       k_m0 = sens[1], k_m1 = sens[2], k_m2 = sens[3], k_m5 = sens[4],
       plateau_obs = if (k_rep >= 200) sv[150] / sv[200] else NA_real_,
       plateau_mp  = if (k_rep >= 200) s_unit[150] / s_unit[200] else NA_real_)
}

runs <- strsplit(opt$runs, ",", fixed = TRUE)[[1]]
rows <- list(); fits <- list()
for (r in runs) {
  kv <- strsplit(r, "=", fixed = TRUE)[[1]]
  label <- kv[1]; dir <- kv[2]
  svf <- file.path(dir, "svd.singularvalues.txt")
  if (!file.exists(svf)) { message("[ndim] ", label, ": no ", svf, " - skipped"); next }
  svt <- fread(svf)
  if (!"SINGULAR_VALUES" %in% names(svt)) {
    stop(svf, " has no SINGULAR_VALUES column (columns: ", paste(names(svt), collapse = ", "), ")", call. = FALSE)
  }
  sv <- as.numeric(svt[["SINGULAR_VALUES"]])
  if (is.unsorted(rev(sv))) stop(svf, ": singular values are not in decreasing order", call. = FALSE)
  n <- count_lines(file.path(dir, "svd.samples.txt"))
  p <- count_lines(file.path(dir, "svd.bins.txt"))
  f <- fit_one(sv, n, p, opt$margin, opt$gap)
  fits[[label]] <- c(f, list(sv = sv))
  rows[[label]] <- data.table(mode = label, n_samples = n, n_bins = p, k_reported = f$k_reported,
                              sigma = f$sigma, edge = f$edge,
                              k_margin_0 = f$k_m0, k_margin_1pct = f$k_m1,
                              k_margin_2pct = f$k_m2, k_margin_5pct = f$k_m5,
                              k_chosen = f$k, margin = opt$margin, status = f$status,
                              plateau_ratio_obs = f$plateau_obs, plateau_ratio_mp = f$plateau_mp)
  message(sprintf("[ndim] %s: n=%d p=%d sigma=%.4f edge=%.2f k=%s (0/1/2/5%%: %s/%s/%s/%s) %s",
                  label, n, p, f$sigma, f$edge, f$k, f$k_m0, f$k_m1, f$k_m2, f$k_m5, f$status))
}
if (!length(rows)) stop("no NGS-PCA run had a svd.singularvalues.txt", call. = FALSE)
tab <- rbindlist(rows)
fwrite(tab, file.path(opt$out, "ndim_by_mode.tsv"), sep = "\t")

# --- plots -----------------------------------------------------------------

png(file.path(opt$out, "ndim_mp.png"), width = 700 * length(fits), height = 1100, res = 130)
par(mfrow = c(2, length(fits)), mar = c(4.5, 4.5, 3, 1))
for (label in names(fits)) {
  f <- fits[[label]]; sv <- f$sv; k_rep <- length(sv)
  plot(seq_len(k_rep), sv, log = "y", pch = 16, cex = 0.6, xlab = "rank", ylab = "singular value (log)",
       main = sprintf("%s: spectrum vs MP bulk", label))
  if (!is.na(f$sigma)) {
    lines(seq_len(k_rep), f$sigma * f$s_unit, col = "steelblue", lwd = 2)
    abline(h = f$edge, col = "firebrick", lty = 2); abline(v = f$k + 0.5, col = "firebrick", lty = 3)
  }
  legend("topright", c("observed", "MP bulk (fitted sigma)", sprintf("edge (+%.0f%%)", 100 * opt$margin)),
         col = c("black", "steelblue", "firebrick"), pch = c(16, NA, NA), lty = c(NA, 1, 2), bty = "n")
}
for (label in names(fits)) {
  f <- fits[[label]]; sv <- f$sv; k_rep <- length(sv)
  lo <- max(1, (if (is.na(f$k)) 1 else f$k) - 15)
  r <- lo:k_rep
  plot(r, sv[r], pch = 16, cex = 0.6, xlab = "rank", ylab = "singular value",
       main = sprintf("%s: tail, k = %s", label, f$k))
  if (!is.na(f$sigma)) { lines(r, f$sigma * f$s_unit[r], col = "steelblue", lwd = 2); abline(h = f$edge, col = "firebrick", lty = 2) }
}
dev.off()

# --- the one number --------------------------------------------------------

ok <- tab[status == "ok" & !is.na(k_chosen)]
md <- c("# Coverage PCs to remove (Marchenko-Pastur)", "",
        sprintf("- margin above the MP edge: %.0f%%; sigma fitted on ranks past k + %d", 100 * opt$margin, opt$gap),
        "", "| mode | n | bins | reported | sigma | edge | k @0% | k @1% | k @2% | k @5% | chosen | status |",
        "|---|---|---|---|---|---|---|---|---|---|---|---|",
        tab[, sprintf("| %s | %d | %d | %d | %.4f | %.2f | %s | %s | %s | %s | %s | %s |",
                      mode, n_samples, n_bins, k_reported, sigma, edge,
                      k_margin_0, k_margin_1pct, k_margin_2pct, k_margin_5pct, k_chosen, status)])
if (nrow(ok)) {
  ndim <- as.integer(round(mean(ok$k_chosen)))
  writeLines(as.character(ndim), file.path(opt$out, "ndim.txt"))
  md <- c(md, "", sprintf("**ndim = %d** — the rounded mean over %s, applied to every mode.", ndim,
                          paste(ok$mode, collapse = ", ")))
  message(sprintf("[ndim] chosen ndim = %d (mean over %s)", ndim, paste(ok$mode, collapse = ", ")))
} else {
  md <- c(md, "", "No mode could be determined: the reported spectrum does not reach the noise bulk.",
          "Raise NUM_PC upstream (NGS-PCA computes 200 by default) or set EX_NDIM by hand.")
  message("[ndim] undetermined for every mode; ndim.txt not written")
}
md <- c(md, "",
        "`k @0%` is the exact-edge count; the plateau in real coverage data decays a little",
        "faster than iid Marchenko-Pastur (per-bin variance differs), which makes it run high —",
        sprintf("observed s[150]/s[200] = %s vs MP %s.",
                paste(sprintf("%.4f", tab$plateau_ratio_obs), collapse = "/"),
                paste(sprintf("%.4f", tab$plateau_ratio_mp), collapse = "/")),
        "The margin absorbs that; treat the 1-2% counts as the defensible range and the",
        "diagnostic plot (`ndim_mp.png`) as the evidence. Override with EX_NDIM at any time.")
writeLines(md, file.path(opt$out, "summary.md"))
