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
#   sex_linear                 SEX was inferred from X/Y coverage upstream.
#                              Under the ploidy model chrX is normalised by
#                              its expected copies, so the share of SEX's
#                              variance that corrected chrX depth explains
#                              must be small: a misaligned sex table gives
#                              R^2 ~ 0.5 there. chrY is fitted on males
#                              only, so N on chrY must equal the male count.
#   inferred_sex               the same response through the logistic engine.
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
  make_option("--source",   type = "character", default = "", help = "provenance of this mode's inputs (from paths.env)"),
  make_option("--pheno",    type = "character", default = NULL, help = "phenotypes.tsv, for the male count chrY tests must equal"),
  make_option("--samples",  type = "character", default = NULL, help = "samples.txt: the matrix columns, so that count covers tested samples only"),
  make_option("--ploidy",   type = "character", default = "1", help = "1 if the correction ran with the ploidy model (EX_PLOIDY)"),
  make_option("--out",      type = "character", help = "output directory")
)
opt <- parse_args(OptionParser(option_list = option_list))
for (req in c("assoc", "analyses", "out")) if (is.null(opt[[req]])) stop("--", req, " is required", call. = FALSE)
if (!opt$profile %in% c("real", "smoke")) stop("--profile must be real or smoke", call. = FALSE)
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
smoke  <- opt$profile == "smoke"
ploidy <- opt$ploidy == "1"

# Males among the samples the sweep could test: in the matrix (samples.txt)
# and with the phenotype. The phenotype table carries every upstream sample,
# which is more than a smoke tree simulates.
n_male <- NA_integer_
if (!is.null(opt$pheno) && file.exists(opt$pheno)) {
  ph <- fread(opt$pheno)
  if ("SEX_MF" %in% names(ph)) {
    tested <- rep(TRUE, nrow(ph))
    if (!is.null(opt$samples) && file.exists(opt$samples)) {
      tested <- tested & ph$SAMPLE %in% trimws(readLines(opt$samples))
    }
    if ("MTDNA_CN" %in% names(ph)) tested <- tested & is.finite(ph$MTDNA_CN)
    n_male <- sum(tested & ph$SEX_MF == "M", na.rm = TRUE)
  }
}

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
  dt <- rbindlist(lapply(files, function(f) fread(cmd = paste("gzip -cd", shQuote(f)))), fill = TRUE)
  setnames(dt, 1:7, c("CHROM", "START", "END", "Region", "N", "NCase", "NControl"))
  for (col in c("BETA", "STAT", "P")) {
    if (!col %in% names(dt)) stop("shards for ", name, " lack a ", col, " column", call. = FALSE)
  }
  dt[, `:=`(P = as.numeric(P), Estimate = as.numeric(BETA), STAT = as.numeric(STAT))]
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

# Comment lines may sit anywhere in the manifest; filter before parsing.
mf_lines <- readLines(opt$analyses)
mf_lines <- mf_lines[!grepl("^\\s*#", mf_lines) & nzchar(trimws(mf_lines))]
manifest <- fread(text = mf_lines, header = FALSE, sep = "\t", fill = TRUE)
setnames(manifest, 1:3, c("name", "method", "model"))

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

  # Checks key on the analysis family: the covariate-adjusted variants
  # (preamble.sh ran) and the inverse-normal-transformed ones carry the same
  # truth as their plain namesakes.
  fam <- sub("_int$", "", sub("_adj$", "", nm))

  if (fam %in% c("mtdna_cn", "log2_mtdna_cn")) {
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
      if (fam == "log2_mtdna_cn") {
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

    # Under the ploidy model chrY has no expected copies in females, so every
    # chrY region is fitted on the males alone: N there must equal the male
    # count. (A linear test of SEX itself has no variance within the males,
    # so this lives on the mtDNA family.)
    if (ploidy && fam == "mtdna_cn") {
      dy <- d[toupper(sub("^chr", "", CHROM)) == "Y"]
      if (!nrow(dy)) {
        add(nm, "chrY_tested_in_males", "WARN", "no chrY bins tested", "")
      } else if (is.na(n_male)) {
        add(nm, "chrY_tested_in_males", "INFO", sprintf("N on chrY: %s", paste(unique(dy$N), collapse = ",")),
            "no SEX_MF column to compare against")
      } else {
        ok_n <- all(dy$N == n_male)
        add(nm, "chrY_tested_in_males", if (ok_n) "PASS" else "FAIL",
            sprintf("N on chrY: %s; males in the sex table: %d", paste(unique(dy$N), collapse = ","), n_male),
            "chrY regions are fitted on the samples with an expected copy — the males")
      }
    }
  }

  if (fam == "sex_linear") {
    # Under the ploidy model chrX is normalised by its expected copies, so
    # the corrected chrX depth should explain little of SEX. The linear
    # t-statistic on a binary response converts to that R^2 exactly
    # (R^2 = t^2 / (t^2 + df)); a misaligned sex table leaves each bin at
    # R^2 ~ 0.5, while real residual sex effects sit far below.
    dx <- d[toupper(sub("^chr", "", CHROM)) == "X"]
    if (!ploidy) {
      add(nm, "sex_top_hit", if (top$class[1] == "sex") "PASS" else "FAIL",
          sprintf("top |stat| at %s (stat=%.1f)", top$Region[1], top$STAT[1]),
          "without the ploidy model (EX_PLOIDY=0) the sex chromosomes must dominate")
    } else if (!nrow(dx)) {
      add(nm, "sex_signal_removed", "WARN", "no chrX bins tested", "")
    } else {
      dx[, df := N - 2L]
      dx[, r2 := STAT^2 / (STAT^2 + df)]
      med_r2 <- median(dx$r2, na.rm = TRUE); max_r2 <- max(dx$r2, na.rm = TRUE)
      status <- if (med_r2 > 0.35) "FAIL" else if (med_r2 > 0.15) "WARN" else "PASS"
      add(nm, "sex_signal_removed", status,
          sprintf("median R^2(SEX ~ chrX depth) = %.3f, max %.3f over %d bins", med_r2, max_r2, nrow(dx)),
          "ploidy model applied: a misaligned sex table gives ~0.5; residual per-bin sex effects stay well below")
      k <- min(100L, nrow(top))
      add(nm, "sex_top100_share", "INFO",
          sprintf("%.0f%% of top %d on chrX/Y; top |stat| at %s (stat=%.1f)",
                  100 * mean(top$class[seq_len(k)] == "sex"), k, top$Region[1], top$STAT[1]))
    }
    fwrite(top[seq_len(min(25, nrow(top)))],
           file.path(opt$out, "top_hits.sex_linear.tsv"), sep = "\t")
  }

  if (fam == "inferred_sex") {
    # Engine coverage on a binary response. Without the ploidy model the
    # sexes are completely separated on chrX/Y and the Wald z collapses
    # (Hauck-Donner); with it there is nothing to separate.
    k <- min(100L, nrow(top))
    frac_top_sex <- mean(top$class[seq_len(k)] == "sex")
    med_sex_z  <- suppressWarnings(median(d[class == "sex"]$abs_stat))
    med_auto_z <- suppressWarnings(median(d[class == "autosome"]$abs_stat))
    add(nm, "logistic_ran", "PASS", sprintf("%d regions", nrow(d)),
        "the logistic engine completed over the sweep")
    add(nm, "logistic_sex_note", "INFO",
        sprintf("median |z|: sex=%.2f autosome=%.2f; %.0f%% of top %d on chrX/Y",
                med_sex_z, med_auto_z, 100 * frac_top_sex, k),
        if (ploidy) "ploidy model applied; chrX/Y are expected to look like autosomes here"
        else "small sex-chromosome z under complete separation is Hauck-Donner, not a miss")
    add(nm, "case_control_counts", "INFO",
        sprintf("NCase=%d NControl=%d", top$NCase[1], top$NControl[1]),
        "NCase should equal the male count from the QC table")
    fwrite(top[seq_len(min(25, nrow(top)))],
           file.path(opt$out, "top_hits.inferred_sex.tsv"), sep = "\t")
  }

  if (fam == "mtdna_cn_null") {
    # Calibration is judged on the autosomes. Sex-chromosome bins all carry
    # the same sex vector, so their tests are one dependent draw rather
    # than thousands of independent ones — a majority of the bins in the
    # smoke slice, a few percent at genome scale, and misleading either way.
    dn <- d[class == "autosome"]
    scope <- "autosomal"
    if (nrow(dn) < 50) { dn <- d; scope <- "all" }
    lam <- lambda_gc(dn$P)
    band <- if (smoke) c(0.6, 1.6) else c(0.85, 1.20)
    add(nm, "null_lambda", if (!is.na(lam) && lam >= band[1] && lam <= band[2]) "PASS" else "WARN",
        sprintf("%.3f (%s bins)", lam, scope),
        sprintf("acceptance band %.2f-%.2f (%s profile)", band[1], band[2], opt$profile))
    frac05 <- mean(dn$P < 0.05)
    # The smoke band reaches lower: removing ndim real-PC dimensions from 64
    # samples deflates the residual variance the test never learns about.
    fb <- if (smoke) c(0.015, 0.10) else c(0.03, 0.07)
    add(nm, "null_frac_p05", if (frac05 >= fb[1] && frac05 <= fb[2]) "PASS" else "WARN",
        sprintf("%.3f (%s bins)", frac05, scope), "fraction of regions at p<0.05; ~0.05 when calibrated")
    add(nm, "null_lambda_all_bins", "INFO", sprintf("%.3f", lambda_gc(d$P)),
        "includes chrX/Y/M; moves with a single sex draw in small cohorts")
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
        if (nzchar(opt$source)) sprintf("- inputs: `%s`%s", opt$source,
                                        if (grepl("synthetic", opt$source, ignore.case = TRUE))
                                          " — **SYNTHETIC**: simulated depths, a machinery check only" else "")
        else character(0),
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
