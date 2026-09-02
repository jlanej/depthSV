#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# depthSV 1000G example — turn the timing records into a profile report
#
# Input is the recorder's table (one row per timed stage invocation) and,
# when available, sacct's view of the same jobs. Output is what a person
# deciding "what do we speed up next?" needs: where the wall time went, how
# it distributes over work units, the stragglers, and peak memory.
#
#   Rscript profile_report.R --timings timings.tsv [--sacct sacct.tsv] --out DIR
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

option_list <- list(
  make_option("--timings", type = "character"),
  make_option("--sacct",   type = "character", default = NULL),
  make_option("--out",     type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))
for (req in c("timings", "out")) if (is.null(opt[[req]])) stop("--", req, " is required", call. = FALSE)
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)

tm <- fread(opt$timings)
if (!nrow(tm)) stop("no timing records in ", opt$timings, call. = FALSE)
tm[, elapsed_s := as.numeric(elapsed_s)]
tm[, rss_mb := suppressWarnings(as.numeric(max_rss_kb)) / 1024]

# Failed invocations are part of the profile: a stage that died after two
# hours is a finding, not a row to drop. They are flagged, not filtered.
tm[, ok := exit == 0]

# --- by stage --------------------------------------------------------------

q95 <- function(x) as.numeric(quantile(x, 0.95, names = FALSE, type = 7))
stage_summary <- tm[, .(
  invocations = .N,
  failed      = sum(!ok),
  total_s     = round(sum(elapsed_s)),
  mean_s      = round(mean(elapsed_s), 1),
  median_s    = round(median(elapsed_s), 1),
  p95_s       = round(q95(elapsed_s), 1),
  max_s       = max(elapsed_s),
  max_rss_mb  = if (all(is.na(rss_mb))) NA_real_ else round(max(rss_mb, na.rm = TRUE))
), by = .(mode, stage)][order(-total_s)]
fwrite(stage_summary, file.path(opt$out, "stage_summary.tsv"), sep = "\t")

grand_total <- sum(tm$elapsed_s)
stage_share <- tm[, .(total_s = sum(elapsed_s)), by = stage][order(-total_s)]
stage_share[, share := sprintf("%.1f%%", 100 * total_s / grand_total)]

slowest <- tm[order(-elapsed_s)][seq_len(min(15L, nrow(tm))),
              .(mode, stage, unit, elapsed_s, rss_mb = round(rss_mb), exit, host, jobid)]
fwrite(slowest, file.path(opt$out, "slowest_units.tsv"), sep = "\t")

# --- sacct (optional) ------------------------------------------------------

sacct_md <- character(0)
if (!is.null(opt$sacct) && file.exists(opt$sacct)) {
  sa <- tryCatch(fread(opt$sacct, sep = "|"), error = function(e) NULL)
  if (!is.null(sa) && nrow(sa) && "MaxRSS" %in% names(sa)) {
    # MaxRSS lives on the .batch step rows; --noconvert leaves it in bytes.
    steps <- sa[grepl("\\.batch$", JobID)]
    if (nrow(steps)) {
      steps[, rss_gb := suppressWarnings(as.numeric(sub("K$", "e3", sub("M$", "e6", sub("G$", "e9", MaxRSS)))) / 1e9)]
      top_mem <- steps[order(-rss_gb)][seq_len(min(10L, nrow(steps))),
                                       .(JobID, Elapsed, MaxRSS_GB = round(rss_gb, 2))]
      sacct_md <- c("", "## Scheduler accounting (sacct)", "",
                    sprintf("Full table: `%s`. Highest-memory steps:", basename(opt$sacct)), "",
                    "| JobID | Elapsed | MaxRSS (GB) |", "|---|---|---|",
                    top_mem[, sprintf("| %s | %s | %s |", JobID, Elapsed, MaxRSS_GB)])
    }
  }
}

# --- report ----------------------------------------------------------------

hms <- function(s) sprintf("%dh%02dm%02ds", s %/% 3600, (s %% 3600) %/% 60, round(s %% 60))
na_dash <- function(x) ifelse(is.na(x), "-", as.character(x))

md <- c("# depthSV 1000G example — runtime profile", "",
        sprintf("- timed invocations: %d (%d failed)", nrow(tm), sum(!tm$ok)),
        sprintf("- recorded wall time: %s (sum over invocations; concurrent units overlap in real time)",
                hms(grand_total)),
        "",
        "## Where the time went", "",
        "| stage | total | share |", "|---|---|---|",
        stage_share[, sprintf("| %s | %s | %s |", stage, hms(total_s), share)],
        "",
        "## By mode and stage", "",
        "| mode | stage | n | failed | total (s) | median (s) | p95 (s) | max (s) | max RSS (MB) |",
        "|---|---|---|---|---|---|---|---|---|",
        stage_summary[, sprintf("| %s | %s | %d | %d | %d | %s | %s | %s | %s |",
                                mode, stage, invocations, failed, total_s,
                                median_s, p95_s, max_s, na_dash(max_rss_mb))],
        "",
        "## Slowest work units", "",
        "| mode | stage | unit | s | RSS (MB) | exit | host |", "|---|---|---|---|---|---|---|",
        slowest[, sprintf("| %s | %s | %s | %d | %s | %d | %s |",
                          mode, stage, unit, elapsed_s, na_dash(rss_mb), exit, host)],
        sacct_md, "",
        "Peak RSS comes from GNU time where it exists; `-` means it was not",
        "recordable there (sacct fills that gap on SLURM). A failed invocation",
        "keeps its row - a stage that dies late is exactly what this table is for.")
writeLines(md, file.path(opt$out, "profile_report.md"))

cat(sprintf("[profile] %d records, %s recorded wall time -> %s\n",
            nrow(tm), hms(grand_total), file.path(opt$out, "profile_report.md")))
