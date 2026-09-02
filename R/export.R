#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# depthSV — summary of one exported analysis
#
# Called by scripts/export.sh on the concatenated, count-suppressed table.
# Reports the genomic-control lambda, the Bonferroni threshold and, when the
# analysis stage ran with --perms, the empirical family-wise threshold from
# the folded max-|t| distribution (Westfall-Young):
#
#   threshold  the k-th largest permutation maximum, k = floor(alpha (B+1)),
#              so that P(max |t| >= threshold) <= alpha under the null;
#              B must be at least ceil(1/alpha) - 1 (19 at alpha = 0.05)
#   M_eff      the number of independent tests a Sidak correction would need
#              to reach that threshold's p: log(1 - alpha) / log(1 - p_thr)
#   P_ADJ      (1 + #{maxima >= |t|}) / (B + 1) for every hit
#
# The permutation statistic is the classical t of the linear model, so the
# threshold applies to STAT as written by the linear method; --robust
# changes SE and STAT but not the permutation distribution.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

option_list <- list(
  make_option("--input",   type = "character", help = "exported table (bgzip)"),
  make_option("--permmax", type = "character", default = NULL, help = "folded permutation maxima"),
  make_option("--counts",  type = "character", help = "rows_in / rows_suppressed from the concatenation"),
  make_option("--alpha",   type = "double",    default = 0.05),
  make_option("--minCount", type = "integer",  default = 20L),
  make_option("--name",    type = "character", default = ""),
  make_option("--method",  type = "character", default = ""),
  make_option("--regionsListed", type = "integer", default = NA_integer_),
  make_option("--shards",  type = "integer",   default = NA_integer_),
  make_option("--out",     type = "character", help = "summary tsv"),
  make_option("--hits",    type = "character", help = "hits tsv")
)
opt <- parse_args(OptionParser(option_list = option_list))
for (req in c("input", "counts", "out", "hits")) if (is.null(opt[[req]])) stop("--", req, " is required", call. = FALSE)
if (!(opt$alpha > 0 && opt$alpha < 1)) stop("--alpha must be in (0, 1)", call. = FALSE)

cnt <- fread(opt$counts, header = FALSE, col.names = c("key", "value"))
rows_in <- cnt[key == "rows_in"]$value
rows_sup <- cnt[key == "rows_suppressed"]$value

d <- fread(cmd = paste("gzip -cd", shQuote(opt$input)), na.strings = c("NA", ""))
M <- nrow(d)
if (M) setnames(d, 1:7, c("CHROM", "START", "END", "Region", "N", "NCase", "NControl"))

lambda_gc <- function(p) {
  p <- p[is.finite(p) & p > 0 & p <= 1]
  if (length(p) < 50) return(NA_real_)
  median(qchisq(p, 1, lower.tail = FALSE)) / qchisq(0.5, 1)
}

p <- if (M) as.numeric(d$P) else numeric(0)
stat <- if (M) as.numeric(d$STAT) else numeric(0)
bonf <- if (M) opt$alpha / M else NA_real_
hits_bonf <- sum(p <= bonf, na.rm = TRUE)

# --- permutation threshold ---------------------------------------------------
B <- 0L; df <- NA_integer_; thr <- NA_real_; thr_p <- NA_real_; m_eff <- NA_real_
hits_emp <- NA_integer_; perm_note <- "no permutation maxima (analyze.sh --perms 0)"
mx <- numeric(0)
if (!is.null(opt$permmax) && file.exists(opt$permmax)) {
  lines <- readLines(opt$permmax)
  hdr <- lines[grepl("^#", lines)]
  getv <- function(k) { v <- sub(paste0("^#", k, "="), "", hdr[grepl(paste0("^#", k, "="), hdr)]); if (length(v)) v[1] else NA }
  df <- as.integer(getv("df"))
  body <- lines[!grepl("^#", lines) & !grepl("^perm\t", lines) & nzchar(lines)]
  if (length(body)) {
    pm <- fread(text = body, header = FALSE, col.names = c("perm", "max_abs_stat"))
    mx <- as.numeric(pm$max_abs_stat)
    B <- length(mx)
  }
  k <- floor(opt$alpha * (B + 1))
  if (B > 0 && k >= 1 && is.finite(df)) {
    thr <- sort(mx, decreasing = TRUE)[k]
    thr_p <- 2 * pt(-thr, df)
    m_eff <- log(1 - opt$alpha) / log(1 - thr_p)
    hits_emp <- sum(abs(stat) >= thr, na.rm = TRUE)
    perm_note <- sprintf("Westfall-Young max-|t| over %d permutations; k=%d", B, k)
    if (is.finite(m_eff) && m_eff > M) {
      perm_note <- paste0(perm_note, "; M_eff exceeds the region count: the permutation null is heavier-tailed",
                          " than the t reference (small n, non-normal response), so the empirical threshold",
                          " is stricter than Bonferroni")
    }
  } else if (B > 0) {
    perm_note <- sprintf("%d permutations are too few for alpha=%g: need at least %d", B, opt$alpha,
                         ceiling(1 / opt$alpha) - 1L)
  }
}
p_adj <- function(t_abs) if (B > 0) (1 + sum(mx >= t_abs)) / (B + 1) else NA_real_

# --- hits ----------------------------------------------------------------------
if (M) {
  sel <- (p <= bonf) | (is.finite(thr) & abs(stat) >= thr)
  sel[is.na(sel)] <- FALSE
  h <- d[sel]
  h[, P_ADJ := if (nrow(h)) vapply(abs(as.numeric(STAT)), p_adj, numeric(1)) else numeric(0)]
  h[, PASSES := ifelse(as.numeric(P) <= bonf, "bonferroni", "") ]
  if (is.finite(thr)) h[abs(as.numeric(STAT)) >= thr, PASSES := ifelse(nzchar(PASSES), paste0(PASSES, ",empirical"), "empirical")]
  setorder(h, P)
  fwrite(h, opt$hits, sep = "\t")
} else {
  writeLines("#CHROM\tSTART\tEND\tRegion\tN\tNCase\tNControl\tP_ADJ\tPASSES", opt$hits)
}

# --- summary ---------------------------------------------------------------------
fmt <- function(x) if (is.na(x)) "NA" else formatC(x, digits = 6, format = "g")
# (`key` is an argument of data.table() itself, hence the neutral names.)
summary <- data.table(k = c(
  "analysis", "method", "regions_listed", "shards", "rows_in", "rows_suppressed", "min_count",
  "regions_exported", "samples_max", "lambda_gc", "alpha", "bonferroni_p", "hits_bonferroni",
  "perms", "perm_df", "perm_threshold_stat", "perm_threshold_p", "m_eff", "hits_empirical", "perm_note"),
  v = c(
  opt$name, opt$method, opt$regionsListed, opt$shards, rows_in, rows_sup, opt$minCount,
  M, if (M) max(as.numeric(d$N), na.rm = TRUE) else 0, fmt(lambda_gc(p)), opt$alpha, fmt(bonf), hits_bonf,
  B, df, fmt(thr), fmt(thr_p), fmt(m_eff), hits_emp, perm_note))
fwrite(summary, opt$out, sep = "\t", col.names = FALSE)

message(sprintf("[export] %s.%s: %d regions (%d suppressed at min_count=%d), lambda=%s, Bonferroni p<=%s: %d hit(s)%s",
                opt$name, opt$method, M, rows_sup, opt$minCount, fmt(lambda_gc(p)), fmt(bonf), hits_bonf,
                if (is.finite(thr)) sprintf("; empirical |t|>=%s (p<=%s, M_eff=%s): %d hit(s)", fmt(thr), fmt(thr_p), fmt(m_eff), hits_emp) else ""))
