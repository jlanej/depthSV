#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# depthSV 1000G example — compare association results between mosdepth modes
#
# The question this answers end to end: does running mosdepth in --fast-mode
# upstream change what the depth-association pipeline concludes? Upstream,
# NGS-PCA's own evaluation found the two modes' PCs near-identical because
# fast mode's effects are largely per-sample multiplicative shifts that the
# log2-vs-own-median normalisation absorbs; the same argument applies to the
# corrected depths here, so near-perfect concordance of the association
# statistics is the expectation — and anything less is a finding.
#
# Per analysis, on the regions tested in both modes: Pearson r of estimates
# and test statistics, rank correlation, sign agreement where |stat| > 2,
# top-K overlap, and the median relative difference of the estimates (the
# NGS-PCA QC comparison's convention, which exposes a uniform bias even
# where correlation is perfect). Broken out by chrM / sex / autosome so a
# class-specific divergence (chrM is where fast mode's skipped mate-overlap
# correction bites hardest) cannot hide in a genome-wide average.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

option_list <- list(
  make_option("--a-dir",    type = "character", dest = "a_dir", help = "assoc dir, reference mode"),
  make_option("--b-dir",    type = "character", dest = "b_dir", help = "assoc dir, comparison mode"),
  make_option("--a-name",   type = "character", dest = "a_name", default = "standard"),
  make_option("--b-name",   type = "character", dest = "b_name", default = "fast"),
  make_option("--analyses", type = "character", help = "analyses.tsv (either mode's)"),
  make_option("--top-k",    type = "integer",   dest = "top_k", default = 100L),
  make_option("--profile",  type = "character", default = "real", help = "real | smoke"),
  make_option("--timings",  type = "character", default = NULL, help = "optional timings.tsv for the runtime table"),
  make_option("--a-source", type = "character", dest = "a_source", default = "", help = "provenance of mode a's inputs"),
  make_option("--b-source", type = "character", dest = "b_source", default = "", help = "provenance of mode b's inputs"),
  make_option("--out",      type = "character", help = "output directory")
)
opt <- parse_args(OptionParser(option_list = option_list))
for (req in c("a_dir", "b_dir", "analyses", "out")) {
  if (is.null(opt[[req]])) stop("--", gsub("_", "-", req), " is required", call. = FALSE)
}
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)

# Same reader as R/evaluate.R; kept local so each script stands alone.
read_shards <- function(dir, name, method) {
  files <- Sys.glob(file.path(dir, sprintf("%s.%s.*.txt.gz", name, method)))
  files <- files[!grepl("\\.tmp\\.gz$", files)]
  if (!length(files)) return(NULL)
  dt <- rbindlist(lapply(files, function(f) fread(cmd = paste("gzip -cd", shQuote(f)))), fill = TRUE)
  setnames(dt, 1:7, c("CHROM", "START", "END", "Region", "N", "NCase", "NControl"))
  for (col in c("BETA", "STAT", "P")) {
    if (!col %in% names(dt)) stop("shards for ", name, " lack a ", col, " column", call. = FALSE)
  }
  dt[, .(CHROM, Region, N,
         Estimate = as.numeric(BETA), STAT = as.numeric(STAT), P = as.numeric(P))]
}

chrom_class <- function(x) {
  b <- toupper(sub("^chr", "", x))
  fifelse(b %in% c("M", "MT"), "chrM", fifelse(b %in% c("X", "Y"), "sex", "autosome"))
}

lambda_gc <- function(p) {
  p <- p[is.finite(p) & p > 0 & p <= 1]
  if (length(p) < 50) return(NA_real_)
  median(qchisq(p, 1, lower.tail = FALSE)) / qchisq(0.5, 1)
}

# Comment lines may sit anywhere in the manifest (the sweep skips them
# wherever they are), so filter before parsing rather than letting fread stop
# at the first one.
mf_lines <- readLines(opt$analyses)
mf_lines <- mf_lines[!grepl("^\\s*#", mf_lines) & nzchar(trimws(mf_lines))]
manifest <- fread(text = mf_lines, header = FALSE, sep = "\t", col.names = c("name", "method", "model"))

r_pass <- if (opt$profile == "smoke") 0.95 else 0.99

rows <- list(); verdicts <- character(0)
for (i in seq_len(nrow(manifest))) {
  nm <- manifest$name[i]; method <- manifest$method[i]
  a <- read_shards(opt$a_dir, nm, method)
  b <- read_shards(opt$b_dir, nm, method)
  if (is.null(a) || is.null(b)) {
    rows[[nm]] <- data.table(analysis = nm, status = "FAIL", n_common = 0L,
                             detail = sprintf("missing results (%s: %s, %s: %s)",
                                              opt$a_name, !is.null(a), opt$b_name, !is.null(b)))
    verdicts <- c(verdicts, "FAIL")
    next
  }
  j <- merge(a, b, by = "Region", suffixes = c("_a", "_b"))
  j <- j[is.finite(STAT_a) & is.finite(STAT_b)]
  if (!nrow(j)) {
    rows[[nm]] <- data.table(analysis = nm, status = "FAIL", n_common = 0L,
                             detail = "no region tested in both modes")
    verdicts <- c(verdicts, "FAIL")
    next
  }
  j[, class := chrom_class(CHROM_a)]

  strong <- j[abs(STAT_a) > 2 | abs(STAT_b) > 2]
  k <- min(opt$top_k, nrow(j))
  top_a <- j[order(-abs(STAT_a))]$Region[seq_len(k)]
  top_b <- j[order(-abs(STAT_b))]$Region[seq_len(k)]

  by_class <- j[, .(n = .N, r_stat = if (.N >= 10) suppressWarnings(cor(STAT_a, STAT_b)) else NA_real_),
                by = class]
  cls <- function(cl) by_class[class == cl]$r_stat[1]

  r_stat <- suppressWarnings(cor(j$STAT_a, j$STAT_b))
  # A correlation over a handful of regions says nothing: require a real
  # overlap before calling the pair concordant.
  min_common <- max(100L, as.integer(0.5 * min(nrow(a), nrow(b))))
  status <- if (nrow(j) < min_common) "WARN" else if (!is.finite(r_stat)) "WARN" else if (r_stat >= r_pass) "PASS" else "WARN"
  verdicts <- c(verdicts, status)

  rows[[nm]] <- data.table(
    analysis        = nm,
    status          = status,
    n_a             = nrow(a),
    n_b             = nrow(b),
    n_common        = nrow(j),
    r_estimate      = suppressWarnings(cor(j$Estimate_a, j$Estimate_b)),
    r_stat          = r_stat,
    rho_abs_stat    = suppressWarnings(cor(abs(j$STAT_a), abs(j$STAT_b), method = "spearman")),
    sign_agree_gt2  = if (nrow(strong)) mean(sign(strong$STAT_a) == sign(strong$STAT_b)) else NA_real_,
    topk_jaccard    = length(intersect(top_a, top_b)) / length(union(top_a, top_b)),
    med_rel_dEst    = if (nrow(strong)) {
                        median(abs(strong$Estimate_a - strong$Estimate_b) /
                               pmax(abs(strong$Estimate_a), abs(strong$Estimate_b)))
                      } else NA_real_,
    # Signed, b relative to a: a uniform bias (fast mode's skipped mate-overlap
    # correction shows as one) has a sign, and NGS-PCA's QC concordance
    # reports it the same way.
    med_signed_rel_dEst = if (nrow(strong)) {
                        median((strong$Estimate_b - strong$Estimate_a) / abs(strong$Estimate_a))
                      } else NA_real_,
    max_abs_dStat   = max(abs(j$STAT_a - j$STAT_b)),
    r_stat_chrM     = cls("chrM"),
    r_stat_sex      = cls("sex"),
    r_stat_autosome = cls("autosome"),
    lambda_a        = lambda_gc(a$P),
    lambda_b        = lambda_gc(b$P),
    detail          = if (nrow(j) < min_common) sprintf("only %d regions in common (need %d)", nrow(j), min_common)
                      else sprintf("top-%d overlap on |stat|", k)
  )
}

conc <- rbindlist(rows, fill = TRUE)
fwrite(conc, file.path(opt$out, "concordance.tsv"), sep = "\t")

# --- runtime side-by-side --------------------------------------------------

timing_md <- character(0)
if (!is.null(opt$timings) && file.exists(opt$timings)) {
  tm <- fread(opt$timings)
  if (nrow(tm)) {
    agg <- tm[, .(units = .N, total_s = sum(elapsed_s), max_s = max(elapsed_s)),
              by = .(mode, stage)][order(stage, mode)]
    timing_md <- c("", "## depthSV runtime by mode", "",
                   "Depth of the depthSV stages themselves — expected to be mode-independent.",
                   "The upstream mosdepth speedup is measured by NGS-PCA's own",
                   "`04_fast_mode_eval.sh`, not here.", "",
                   "| stage | mode | units | total (s) | max (s) |", "|---|---|---|---|---|",
                   agg[, sprintf("| %s | %s | %d | %d | %d |", stage, mode, units, total_s, max_s)])
  }
}

# --- report ----------------------------------------------------------------

overall <- if ("FAIL" %in% verdicts) "FAIL" else if ("WARN" %in% verdicts) "PASS with warnings" else "PASS"
fmt <- function(x, d = 4) ifelse(is.na(x), "NA", formatC(x, digits = d, format = "f"))

synthetic <- grepl("synthetic", paste(opt$a_source, opt$b_source), ignore.case = TRUE)
md <- c(sprintf("# depthSV 1000G example — %s vs %s", opt$a_name, opt$b_name), "",
        sprintf("- verdict: **%s** (PASS needs r(stat) >= %.2f per analysis; `%s` profile)",
                overall, r_pass, opt$profile),
        if (nzchar(opt$a_source)) sprintf("- %s inputs: `%s`", opt$a_name, opt$a_source) else character(0),
        if (nzchar(opt$b_source)) sprintf("- %s inputs: `%s`", opt$b_name, opt$b_source) else character(0),
        if (synthetic) paste("- **SYNTHETIC COMPARISON.** One side's depths were simulated from the other's",
                             "tables with a small perturbation (smoke mode, no upstream fast-mode",
                             "outputs). This exercises the machinery and says nothing about mosdepth.") else character(0),
        "",
        "| analysis | status | common | r(Estimate) | r(stat) | rho(|stat|) | sign agree (|stat|>2) | top-K Jaccard | med rel dEst | r(stat) chrM | r(stat) sex | r(stat) auto |",
        "|---|---|---|---|---|---|---|---|---|---|---|---|",
        conc[, sprintf("| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |",
                       analysis, status, n_common,
                       fmt(r_estimate), fmt(r_stat), fmt(rho_abs_stat),
                       fmt(sign_agree_gt2), fmt(topk_jaccard), fmt(med_rel_dEst),
                       fmt(r_stat_chrM), fmt(r_stat_sex), fmt(r_stat_autosome))],
        "",
        sprintf("Null calibration: lambda(%s) / lambda(%s) per analysis are in `concordance.tsv`.",
                opt$a_name, opt$b_name),
        timing_md)
writeLines(md, file.path(opt$out, "summary.md"))

cat(sprintf("[compare] %s vs %s: %s -> %s\n", opt$a_name, opt$b_name, overall,
            file.path(opt$out, "summary.md")))
if (overall == "FAIL") quit(status = 1)
