#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# depthSV 1000G example — the calibration verdict
#
# Puts the fast-mode comparison beside the seed-control comparison and says,
# per analysis, whether mosdepth --fast-mode moved the association
# statistics more than the randomized PCA estimator moves them by itself.
# The distance used is 1 - r(stat) over the regions tested in both modes;
# "within seed noise" means the fast distance is at most --factor times the
# seed-control distance. Top-K overlap and the largest single-region shift
# are shown alongside so a verdict never rests on one number.
#
# This mirrors how NGS-PCA's own report judges the PCs (fast-vs-normal
# against seed-vs-seed) — carried through to the end of the pipeline.
#
#   Rscript calibration_summary.R --primary standard_vs_fast/concordance.tsv \
#       [--control standard_vs_seedctl/concordance.tsv] --factor 1.5 --out DIR
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

option_list <- list(
  make_option("--primary", type = "character", help = "concordance.tsv of standard vs fast"),
  make_option("--control", type = "character", default = NULL, help = "concordance.tsv of standard vs seedctl"),
  make_option("--factor",  type = "double",    default = 1.5),
  make_option("--out",     type = "character", help = "output directory (compare/)")
)
opt <- parse_args(OptionParser(option_list = option_list))
for (req in c("primary", "out")) if (is.null(opt[[req]])) stop("--", req, " is required", call. = FALSE)

p <- fread(opt$primary)
ctl <- if (!is.null(opt$control) && file.exists(opt$control)) fread(opt$control) else NULL

fmt <- function(x, d = 4) ifelse(is.na(x), "NA", formatC(x, digits = d, format = "f"))
dist <- function(r) 1 - r

tab <- p[, .(analysis, status_fast = status, n_common_fast = n_common,
             r_stat_fast = r_stat, topk_fast = topk_jaccard, max_dstat_fast = max_abs_dStat)]

if (!is.null(ctl)) {
  c2 <- ctl[, .(analysis, status_seed = status, n_common_seed = n_common,
                r_stat_seed = r_stat, topk_seed = topk_jaccard, max_dstat_seed = max_abs_dStat)]
  tab <- merge(tab, c2, by = "analysis", all.x = TRUE, sort = FALSE)
  tab[, ratio := dist(r_stat_fast) / dist(r_stat_seed)]
  tab[, verdict := fifelse(!is.finite(ratio), "undetermined",
                   fifelse(ratio <= opt$factor, "within seed noise", "exceeds seed noise"))]
} else {
  tab[, `:=`(status_seed = NA_character_, n_common_seed = NA_integer_, r_stat_seed = NA_real_,
             topk_seed = NA_real_, max_dstat_seed = NA_real_, ratio = NA_real_,
             verdict = "no seed control")]
}
fwrite(tab, file.path(opt$out, "calibration.tsv"), sep = "\t")

n_exceed <- sum(tab$verdict == "exceeds seed noise")
overall <- if (is.null(ctl)) {
  "fast-mode comparison without a calibration yardstick (no seed-control run)"
} else if (n_exceed == 0) {
  "fast mode moves the association statistics no more than a reseed of the PCA does"
} else {
  sprintf("fast mode moves the statistics beyond seed noise for %d of %d analyses", n_exceed, nrow(tab))
}

md <- c("# depthSV 1000G example — mode comparison", "",
        sprintf("**Verdict:** %s.", overall), "",
        "Distances are 1 - r(stat) over regions tested in both modes; `within seed noise`",
        sprintf("means the fast distance is at most %.2fx the seed-control distance.", opt$factor), "",
        "| analysis | fast r(stat) | seed r(stat) | fast/seed distance | fast top-K | seed top-K | fast max dstat | seed max dstat | verdict |",
        "|---|---|---|---|---|---|---|---|---|",
        tab[, sprintf("| %s | %s | %s | %s | %s | %s | %s | %s | %s |",
                      analysis, fmt(r_stat_fast), fmt(r_stat_seed), fmt(ratio, 2),
                      fmt(topk_fast, 3), fmt(topk_seed, 3),
                      fmt(max_dstat_fast, 2), fmt(max_dstat_seed, 2), verdict)],
        "",
        "Per-pair detail, including the chrM / sex / autosome breakdown and the runtime table:",
        "",
        sprintf("- `%s`", file.path(dirname(opt$primary), "summary.md")),
        if (!is.null(ctl)) sprintf("- `%s`", file.path(dirname(opt$control), "summary.md")) else
          "- no seed-control pair: produce one with NGS-PCA's step 2b (a reseeded 02_run_ngspca.sh) and rerun 00 -> 02 for the `seedctl` mode",
        "",
        "Reading it: on chrM — where fast mode skips the mate-overlap correction that",
        "matters most at extreme depth — a fast distance well above the seed distance is",
        "a finding about fast mode; on the null phenotype both distances should be tiny.")
writeLines(md, file.path(opt$out, "summary.md"))

cat(sprintf("[calibration] %s -> %s\n", overall, file.path(opt$out, "summary.md")))
