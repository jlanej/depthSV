#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# depthSV 1000G example — evaluate one mode's association results
#
# The example's phenotypes carry their own truth, so this asserts on results
# rather than exit codes (the same stance as tests/smoke_test.sh):
#
#   mtdna_cn / log2_mtdna_cn   MTDNA_CN is defined upstream as
#                              2 x chrM-mean / HQ-median, so chrM bins must
#                              dominate the association. The log2 phenotype
#                              additionally pins the effect size: on a chrM
#                              bin the corrected depth IS log2(chrM/median),
#                              so the slope should sit near 1.
#   inferred_sex               inferred from X/Y coverage ratios upstream, so
#                              the sex chromosomes must dominate.
#   mtdna_cn_null              a seeded permutation; must stay calibrated.
#
# Autosomal hits for the mtDNA phenotypes are reported, not judged — real
# NUMT regions are expected to appear there on real data.
#
# Severity: FAIL = the machinery is wrong (exit 1). WARN = a statistical
# expectation missed, worth a look. INFO = context.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

option_list <- list(
  make_option("--assoc",    type = "character", help = "association output directory of one mode"),
  make_option("--analyses", type = "character", help = "analyses.tsv used for the sweep"),
  make_option("--regions",  type = "character", default = NULL, help = "the region list the sweep ran over"),
  make_option("--mode",     type = "character", default = "?"),
  make_option("--profile",  type = "character", default = "real", help = "threshold profile: real | smoke"),
  make_option("--out",      type = "character", help = "output directory")
)
opt <- parse_args(OptionParser(option_list = option_list))
for (req in c("assoc", "analyses", "out")) if (is.null(opt[[req]])) stop("--", req, " is required", call. = FALSE)
if (!opt$profile %in% c("real", "smoke")) stop("--profile must be real or smoke", call. = FALSE)
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
smoke <- opt$profile == "smoke"

checks <- data.table(analysis = character(), check = character(), status = character(),
                     value = character(), detail = character())
add <- function(analysis, check, status, value, detail = "") {
  checks <<- rbind(checks, data.table(analysis = analysis, check = check, status = status,
                                      value = as.character(value), detail = detail))
  invisible(NULL)
}

# --- readers ---------------------------------------------------------------

read_shards <- function(dir, name, method) {
  files <- Sys.glob(file.path(dir, sprintf("%s.%s.*.txt.gz", name, method)))
  # .tmp.gz files are failed partial units; never read them as results.
  files <- files[!grepl("\\.tmp\\.gz$", files)]
  if (!length(files)) return(list(files = character(0), dt = NULL))
  dt <- rbindlist(lapply(files, function(f) fread(cmd = paste("gzip -cd", shQuote(f)))))
  setnames(dt, 1:7, c("CHROM", "START", "END", "Region", "N", "NCase", "NControl"))
  ncol_stat <- ncol(dt) - 7L
  stat_idx <- if (method == "coxph") 11L else 10L   # z for coxph, t/z otherwise
  setnames(dt, ncol(dt), "P")
  setnames(dt, 8L, "Estimate")
  setnames(dt, stat_idx, "STAT")
  dt[, `:=`(P = as.numeric(P), Estimate = as.numeric(Estimate), STAT = as.numeric(STAT))]
  list(files = files, dt = dt)
}

chrom_class <- function(x) {
  b <- toupper(sub("^chr", "", x))
  fifelse(b %in% c("M", "MT"), "chrM",
  fifelse(b %in% c("X", "Y"), "sex", "autosome"))
}

lambda_gc <- function(p) {
  p <- p[is.finite(p) & p > 0 & p <= 1]
  if (length(p) < 50) return(NA_real_)
  median(qchisq(p, 1, lower.tail = FALSE)) / qchisq(0.5, 1)
}

manifest <- fread(opt$analyses, header = FALSE, sep = "\t",
                  col.names = c("name", "method", "model"))
manifest <- manifest[!grepl("^#", name) & nzchar(name)]

n_units <- NA_integer_
if (!is.null(opt$regions) && file.exists(opt$regions)) {
  n_units <- length(readLines(opt$regions))
}

# --- per-analysis checks ---------------------------------------------------

for (i in seq_len(nrow(manifest))) {
  nm <- manifest$name[i]; method <- manifest$method[i]
  res <- read_shards(opt$assoc, nm, method)

  if (is.null(res$dt) || !nrow(res$dt)) {
    add(nm, "results_present", "FAIL", "0 rows",
        sprintf("no output shards matched %s.%s.*.txt.gz in %s", nm, method, opt$assoc))
    next
  }
  d <- res$dt
  add(nm, "results_present", "PASS", sprintf("%d regions in %d shards", nrow(d), length(res$files)))

  if (!is.na(n_units) && length(res$files) != n_units) {
    add(nm, "all_units_reported", "WARN",
        sprintf("%d shards vs %d work units", length(res$files), n_units),
        "some array tasks may not have completed (or were QC-skipped whole)")
  } else if (!is.na(n_units)) {
    add(nm, "all_units_reported", "PASS", sprintf("%d shards", n_units))
  }

  ndup <- sum(duplicated(d$Region))
  add(nm, "regions_unique", if (ndup == 0) "PASS" else "FAIL",
      sprintf("%d duplicated", ndup),
      if (ndup) "a bin was tested in two work units; the windowed region list must partition the matrix" else "")

  d[, class := chrom_class(CHROM)]
  d[, abs_stat := abs(STAT)]
  top <- d[order(-abs_stat)]

  if (nm %in% c("mtdna_cn", "log2_mtdna_cn")) {
    n_m <- sum(d$class == "chrM")
    add(nm, "chrM_bins_tested", if (n_m >= 10) "PASS" else "WARN", n_m,
        "17 x 1 kb chrM bins exist; fewer tested means QC dropped some")
    if (n_m == 0) {
      add(nm, "chrM_top_hit", "FAIL", "no chrM bins",
          "the phenotype is built from chrM depth; its bins must be present")
    } else {
      top_is_m <- top$class[1] == "chrM"
      add(nm, "chrM_top_hit", if (top_is_m) "PASS" else "FAIL",
          sprintf("top |stat| at %s (stat=%.1f)", top$Region[1], top$STAT[1]),
          "MTDNA_CN is 2 x chrM-mean / HQ-median by construction")
      frac_sig <- mean(d[class == "chrM"]$P < 1e-6)
      add(nm, "chrM_bins_significant", if (frac_sig >= 0.8) "PASS" else "WARN",
          sprintf("%.0f%% at p<1e-6", 100 * frac_sig))
      if (nm == "log2_mtdna_cn") {
        best <- d[class == "chrM"][which.max(abs_stat)]
        slope_ok <- best$Estimate > 0.5 && best$Estimate < 1.5
        add(nm, "chrM_log2_slope", if (slope_ok) "PASS" else "WARN",
            sprintf("%.3f at %s", best$Estimate, best$Region),
            "log2 phenotype vs log2-ratio depth: slope ~1 expected at the strongest chrM bin")
      }
      auto_top <- d[class == "autosome"][order(-abs_stat)][seq_len(min(10, sum(d$class == "autosome")))]
      if (nrow(auto_top)) {
        add(nm, "top_autosomal_hits", "INFO",
            sprintf("best %s (stat=%.1f, p=%.2g)", auto_top$Region[1], auto_top$STAT[1], auto_top$P[1]),
            "on real data these are NUMT / mito-copy-number-correlated candidates, not errors")
      }
    }
    fwrite(rbind(top[seq_len(min(25, nrow(top)))], d[class == "chrM"])[!duplicated(Region)],
           file.path(opt$out, sprintf("top_hits.%s.tsv", nm)), sep = "\t")
  }

  if (nm == "sex_linear") {
    n_sex <- sum(d$class == "sex")
    if (n_sex == 0) {
      add(nm, "sex_top_hit", "FAIL", "no chrX/chrY bins tested", "")
    } else {
      add(nm, "sex_top_hit", if (top$class[1] == "sex") "PASS" else "FAIL",
          sprintf("top |stat| at %s (stat=%.1f)", top$Region[1], top$STAT[1]),
          "INFERRED_SEX comes from X/Y coverage ratios; the sex chromosomes must dominate")
      k <- min(100L, nrow(top))
      frac_top_sex <- mean(top$class[seq_len(k)] == "sex")
      add(nm, "sex_top100_purity", if (frac_top_sex >= 0.9) "PASS" else "WARN",
          sprintf("%.0f%% of top %d on chrX/Y", 100 * frac_top_sex, k))
    }
    fwrite(top[seq_len(min(25, nrow(top)))],
           file.path(opt$out, "top_hits.sex_linear.tsv"), sep = "\t")
  }

  if (nm == "inferred_sex") {
    # Engine coverage, not a rank assertion: sex separates the corrected
    # depth on chrX/Y completely, and under complete separation the
    # logistic Wald z COLLAPSES (Hauck-Donner) — the strongest real effects
    # report the weakest z. The linear run above carries the truth check;
    # here the collapse itself is the observation worth recording.
    k <- min(100L, nrow(top))
    frac_top_sex <- mean(top$class[seq_len(k)] == "sex")
    med_sex_z  <- suppressWarnings(median(d[class == "sex"]$abs_stat))
    med_auto_z <- suppressWarnings(median(d[class == "autosome"]$abs_stat))
    add(nm, "logistic_ran", "PASS", sprintf("%d regions", nrow(d)),
        "the logistic engine completed over the sweep")
    add(nm, "wald_separation_note", "INFO",
        sprintf("median |z|: sex=%.2f autosome=%.2f; %.0f%% of top %d on chrX/Y",
                med_sex_z, med_auto_z, 100 * frac_top_sex, k),
        "small sex-chromosome z under complete separation is Hauck-Donner, not a miss")
    add(nm, "case_control_counts", "INFO",
        sprintf("NCase=%d NControl=%d", top$NCase[1], top$NControl[1]),
        "NCase should equal the male count from the QC table")
    fwrite(top[seq_len(min(25, nrow(top)))],
           file.path(opt$out, "top_hits.inferred_sex.tsv"), sep = "\t")
  }

  if (nm == "mtdna_cn_null") {
    lam <- lambda_gc(d$P)
    band <- if (smoke) c(0.6, 1.6) else c(0.85, 1.20)
    add(nm, "null_lambda", if (!is.na(lam) && lam >= band[1] && lam <= band[2]) "PASS" else "WARN",
        sprintf("%.3f", lam), sprintf("acceptance band %.2f-%.2f (%s profile)", band[1], band[2], opt$profile))
    frac05 <- mean(d$P < 0.05)
    # The smoke band reaches lower: removing ndim real-PC dimensions from 64
    # samples deflates the residual variance the test never learns about.
    fb <- if (smoke) c(0.015, 0.10) else c(0.03, 0.07)
    add(nm, "null_frac_p05", if (frac05 >= fb[1] && frac05 <= fb[2]) "PASS" else "WARN",
        sprintf("%.3f", frac05), "fraction of regions at p<0.05; ~0.05 when calibrated")
    add(nm, "null_min_p", "INFO", sprintf("%.2g over %d regions", min(d$P), nrow(d)))
  }
}

# --- write -----------------------------------------------------------------

fwrite(checks, file.path(opt$out, "checks.tsv"), sep = "\t")

status_order <- c(FAIL = 1L, WARN = 2L, PASS = 3L, INFO = 4L)
n_fail <- sum(checks$status == "FAIL"); n_warn <- sum(checks$status == "WARN")

md <- c(sprintf("# depthSV 1000G example — evaluation: %s mode", opt$mode),
        "",
        sprintf("- profile: `%s`", opt$profile),
        sprintf("- verdict: **%s** (%d FAIL, %d WARN, %d PASS)",
                if (n_fail) "FAIL" else if (n_warn) "PASS with warnings" else "PASS",
                n_fail, n_warn, sum(checks$status == "PASS")),
        "",
        "| analysis | check | status | value | detail |",
        "|---|---|---|---|---|")
ord <- checks[order(status_order[status], analysis)]
md <- c(md, ord[, sprintf("| %s | %s | %s | %s | %s |", analysis, check, status, value, detail)],
        "", "Top-hit tables sit beside this file; the per-region correction",
        "statistics are under the mode's `corrected/` directory.")
writeLines(md, file.path(opt$out, "summary.md"))

cat(sprintf("[evaluate] %s: %d FAIL, %d WARN -> %s\n", opt$mode, n_fail, n_warn,
            file.path(opt$out, "summary.md")))
if (n_fail > 0) quit(status = 1)
