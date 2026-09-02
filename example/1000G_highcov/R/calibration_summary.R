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
# The headline is only as strong as its rows: it is stated only when at
# least one analysis could be judged, a reseed that moved nothing (seed
# distance 0) makes any fast movement "exceeds", and an analysis that could
# not be compared stays "undetermined". Analyses named by --exclude (the
# logistic sex run, whose Wald statistics are unstable under separation) are
# reported but never judged.
#
# This mirrors how NGS-PCA's own report judges the PCs (fast-vs-normal
# against seed-vs-seed) — carried through to the end of the pipeline. It is
# descriptive: the fast PCs come from their own randomized run, so the fast
# distance carries one reseed of noise by construction and the expected
# ratio under "no effect" is about 1, not 0.
#
#   Rscript calibration_summary.R --primary standard_vs_fast/concordance.tsv \
#       [--control standard_vs_seedctl/concordance.tsv] --factor 1.5 \
#       [--exclude inferred_sex] [--note "..."] --out DIR
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

option_list <- list(
  make_option("--primary", type = "character", help = "concordance.tsv of standard vs fast"),
  make_option("--control", type = "character", default = NULL, help = "concordance.tsv of standard vs seedctl"),
  make_option("--factor",  type = "double",    default = 1.5),
  make_option("--exclude", type = "character", default = "", help = "comma-separated analysis names reported but not judged"),
  make_option("--note",    type = "character", default = "", help = "provenance line printed under the verdict"),
  make_option("--out",     type = "character", help = "output directory (compare/)")
)
opt <- parse_args(OptionParser(option_list = option_list))
for (req in c("primary", "out")) if (is.null(opt[[req]])) stop("--", req, " is required", call. = FALSE)
excluded <- setdiff(trimws(strsplit(opt$exclude, ",", fixed = TRUE)[[1]]), "")

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
  tab[, d_fast := dist(r_stat_fast)]
  tab[, d_seed := dist(r_stat_seed)]
  tab[, ratio := d_fast / d_seed]
  tab[, verdict := fifelse(analysis %in% excluded, "not judged",
                   fifelse(!is.finite(d_fast) | !is.finite(d_seed), "undetermined",
                   fifelse(d_fast <= 0 & d_seed <= 0, "within seed noise",
                   fifelse(d_seed <= 0, "exceeds seed noise",
                   fifelse(ratio <= opt$factor, "within seed noise", "exceeds seed noise")))))]
} else {
  tab[, `:=`(status_seed = NA_character_, n_common_seed = NA_integer_, r_stat_seed = NA_real_,
             topk_seed = NA_real_, max_dstat_seed = NA_real_, d_fast = dist(r_stat_fast),
             d_seed = NA_real_, ratio = NA_real_,
             verdict = fifelse(analysis %in% excluded, "not judged", "no seed control"))]
}
fwrite(tab, file.path(opt$out, "calibration.tsv"), sep = "\t")

judged   <- tab[verdict %in% c("within seed noise", "exceeds seed noise")]
n_exceed <- sum(judged$verdict == "exceeds seed noise")
overall <- if (is.null(ctl)) {
  "fast-mode comparison without a calibration yardstick (no seed-control run)"
} else if (!nrow(judged)) {
  "undetermined: no analysis could be compared against the seed control"
} else if (n_exceed == 0) {
  sprintf("fast mode moves the association statistics no more than a reseed of the PCA does (%d of %d analyses judged)",
          nrow(judged), nrow(tab))
} else {
  sprintf("fast mode moves the statistics beyond seed noise for %d of %d judged analyses", n_exceed, nrow(judged))
}

md <- c("# depthSV 1000G example — mode comparison", "",
        sprintf("**Verdict:** %s.", overall), "",
        if (nzchar(opt$note)) c(opt$note, "") else character(0),
        "Distances are 1 - r(stat) over regions tested in both modes; `within seed noise`",
        sprintf("means the fast distance is at most %.2fx the seed-control distance. A reseed", opt$factor),
        "that moved nothing makes any fast movement `exceeds`; an analysis without a usable",
        "pair stays `undetermined`; excluded analyses are `not judged`.", "",
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
        "a finding about fast mode; on the null phenotype both distances should be tiny.",
        "The verdict is descriptive, not a test: one seed control gives no distribution",
        "for the ratio, and the fast PCs carry a reseed of noise of their own.")
writeLines(md, file.path(opt$out, "summary.md"))

cat(sprintf("[calibration] %s -> %s\n", overall, file.path(opt$out, "summary.md")))
