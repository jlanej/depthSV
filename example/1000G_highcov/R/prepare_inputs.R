#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# depthSV 1000G example — build depthSV input tables from NGS-PCA outputs
#
# From one mode's upstream results this writes the three tables the pipeline
# reads (see the top-level README's input-format table):
#
#   svd.pcs.txt          SAMPLE, PC1..PCn — the upstream PC table with the
#                        NGS-PCA sample suffix (.by1000.) stripped, so IDs
#                        match the depth-matrix columns
#   autosomal.median.txt SAMPLE, AUTO_HQ_median — from NGS-PCA's own
#                        autosomal.median.txt when the run wrote one (the
#                        median over exactly the bins PCA used), else from
#                        the QC table's HQ_MEDIAN_COV
#   phenotypes.tsv       SAMPLE, MTDNA_CN, LOG2_MTDNA_CN, MTDNA_CN_NULL, SEX,
#                        plus covariate candidates carried through verbatim
#
# MTDNA_CN_NULL is MTDNA_CN permuted across samples with a fixed seed: same
# distribution, no genomic structure, so the association sweep gets a
# calibration control alongside the real phenotype.
#
# Two consistency checks guard the one assumption the whole test rests on —
# that the phenotype's denominator is the same median this pipeline
# normalises against: the QC table's HQ_MEDIAN_COV against the median table,
# and MTDNA_CN against its own recomputation. Either disagreeing means the
# upstream QC pass ran before its NGS-PCA run, or against the other mode's.
#
#   Rscript prepare_inputs.R --qc sample_qc.tsv --pcs svd.pcs.txt \
#       [--median autosomal.median.txt] --suffix .by1000. --seed 20260818 --out DIR
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

option_list <- list(
  make_option("--qc",     type = "character", help = "sample_qc.tsv from the NGS-PCA example"),
  make_option("--pcs",    type = "character", help = "svd.pcs.txt from the NGS-PCA example"),
  make_option("--median", type = "character", default = NULL,
              help = "NGS-PCA's autosomal.median.txt (SAMPLE, AUTO_HQ_median, N_BINS); optional"),
  make_option("--covariates", type = "character", default = NULL,
              help = "preamble covariates.tsv (SAMPLE, GPC1..); merged into the phenotype table"),
  make_option("--null-from", type = "character", dest = "null_from", default = NULL,
              help = "phenotypes.tsv of another mode whose MTDNA_CN_NULL is copied by sample ID"),
  make_option("--suffix", type = "character", default = ".by1000.",
              help = "trailing suffix to strip from upstream sample IDs [default %default]"),
  make_option("--seed",   type = "integer",   default = 20260818L),
  make_option("--out",    type = "character", help = "output directory")
)
opt <- parse_args(OptionParser(option_list = option_list))
for (req in c("qc", "pcs", "out")) if (is.null(opt[[req]])) stop("--", req, " is required", call. = FALSE)
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)

note <- function(...) message(sprintf(...))
warned <- 0L
warn <- function(...) { warned <<- warned + 1L; message("[warn] ", sprintf(...)) }

# Replace an output only when its content changed: the pipeline keys finished
# work on its inputs, and a byte-identical table rewritten with a new mtime
# should not look like a new input.
write_if_changed <- function(dt, path) {
  tmp <- paste0(path, ".tmp")
  fwrite(dt, tmp, sep = "\t")
  if (file.exists(path) && isTRUE(unname(tools::md5sum(tmp)) == unname(tools::md5sum(path)))) {
    unlink(tmp)
    return(invisible(FALSE))
  }
  file.rename(tmp, path)
  invisible(TRUE)
}

strip_suffix <- function(ids, suffix) {
  # Upstream tables have carried IDs with a stray space BEFORE the suffix
  # ("HG02635 .by1000."); fread trims one table's IDs and not the other's, so
  # trim both before and after the suffix comes off.
  ids <- trimws(as.character(ids))
  hit <- endsWith(ids, suffix)
  ids[hit] <- substr(ids[hit], 1L, nchar(ids[hit]) - nchar(suffix))
  ids <- trimws(ids)
  attr(ids, "stripped") <- sum(hit)
  ids
}

# Fraction of paired values agreeing to within a relative tolerance.
agree_frac <- function(a, b, tol = 0.005) {
  ok <- is.finite(a) & is.finite(b) & b != 0
  if (!any(ok)) return(NA_real_)
  mean(abs(a[ok] - b[ok]) / abs(b[ok]) <= tol)
}

# --- QC table --------------------------------------------------------------

qc <- fread(opt$qc)
for (col in c("SAMPLE_ID", "HQ_MEDIAN_COV", "MTDNA_CN", "INFERRED_SEX")) {
  if (!col %in% names(qc)) stop(opt$qc, " lacks the ", col, " column", call. = FALSE)
}
qc[, SAMPLE_ID := trimws(as.character(SAMPLE_ID))]
if (anyDuplicated(qc$SAMPLE_ID)) stop("duplicated SAMPLE_ID in ", opt$qc, call. = FALSE)
qc[, HQ_MEDIAN_COV := suppressWarnings(as.numeric(HQ_MEDIAN_COV))]
qc[, MTDNA_CN := suppressWarnings(as.numeric(MTDNA_CN))]
note("[qc] %d samples in %s", nrow(qc), opt$qc)

# --- coverage medians ------------------------------------------------------

if (!is.null(opt$median)) {
  med <- fread(opt$median)
  for (col in c("SAMPLE", "AUTO_HQ_median")) {
    if (!col %in% names(med)) stop(opt$median, " lacks the ", col, " column", call. = FALSE)
  }
  ids <- strip_suffix(med$SAMPLE, opt$suffix)
  if (anyDuplicated(ids)) stop("stripping '", opt$suffix, "' left duplicated IDs in ", opt$median, call. = FALSE)
  med[, SAMPLE := ids]
  med[, AUTO_HQ_median := suppressWarnings(as.numeric(AUTO_HQ_median))]
  # NGS-PCA does not floor the median: a failed or empty sample reads as 0.
  bad <- med[!is.finite(AUTO_HQ_median) | AUTO_HQ_median <= 0]
  if (nrow(bad)) warn("%d samples with a zero/missing median in %s are dropped (e.g. %s)",
                      nrow(bad), basename(opt$median), paste(utils::head(bad$SAMPLE, 3), collapse = ", "))
  cov <- med[is.finite(AUTO_HQ_median) & AUTO_HQ_median > 0, .(SAMPLE, AUTO_HQ_median)]
  note("[coverage] %d samples from NGS-PCA's %s (suffix stripped from %d IDs) -> autosomal.median.txt",
       nrow(cov), basename(opt$median), attr(ids, "stripped"))

  # The QC pass is supposed to have taken HQ_MEDIAN_COV from this very
  # table. If it did not, MTDNA_CN was built on a different denominator
  # than the one this pipeline normalises against.
  chk <- merge(cov, qc[, .(SAMPLE = SAMPLE_ID, HQ_MEDIAN_COV)], by = "SAMPLE")
  f <- agree_frac(chk$HQ_MEDIAN_COV, chk$AUTO_HQ_median)
  if (is.na(f)) {
    warn("QC table has no usable HQ_MEDIAN_COV to check against the median table")
  } else if (f < 0.9) {
    warn(paste("HQ_MEDIAN_COV in the QC table matches %s for only %.0f%% of samples:",
               "the upstream 03a pass ran before its NGS-PCA run, or with NGSPCA_OUTPUT",
               "pointing at another mode's run. Rerun 03a then 03 for this mode."),
         basename(opt$median), 100 * f)
  } else {
    note("[check] QC HQ_MEDIAN_COV agrees with %s for %.1f%% of %d samples",
         basename(opt$median), 100 * f, nrow(chk))
  }
} else {
  cov <- qc[is.finite(HQ_MEDIAN_COV) & HQ_MEDIAN_COV > 0,
            .(SAMPLE = SAMPLE_ID, AUTO_HQ_median = HQ_MEDIAN_COV)]
  if (nrow(cov) < nrow(qc)) {
    note("[coverage] %d samples dropped for missing/zero HQ_MEDIAN_COV", nrow(qc) - nrow(cov))
  }
  note("[coverage] %d samples from the QC table's HQ_MEDIAN_COV -> autosomal.median.txt", nrow(cov))
}
if (!nrow(cov)) stop("no usable coverage medians; the correction stage needs them", call. = FALSE)
write_if_changed(cov, file.path(opt$out, "autosomal.median.txt"))

# --- phenotypes ------------------------------------------------------------

pheno <- data.table(SAMPLE = qc$SAMPLE_ID, MTDNA_CN = qc$MTDNA_CN)
n_cn <- sum(is.finite(pheno$MTDNA_CN))
if (n_cn == 0) {
  stop(paste("MTDNA_CN is NA for every sample in", opt$qc, "- upstream needs HQ_MEDIAN_COV,",
             "which 03a takes from NGS-PCA's autosomal.median.txt (run 03a after 02) or from bedtools."),
       call. = FALSE)
}
if (n_cn < nrow(pheno)) note("[phenotypes] MTDNA_CN missing for %d samples", nrow(pheno) - n_cn)

# MTDNA_CN should equal 2 x chrM-mean / median with the median this
# pipeline uses; MITO_COV_RATIO x MEAN_AUTOSOMAL_COV is chrM-mean.
if (all(c("MITO_COV_RATIO", "MEAN_AUTOSOMAL_COV") %in% names(qc))) {
  re <- merge(qc[, .(SAMPLE = SAMPLE_ID, MTDNA_CN,
                     chrM = suppressWarnings(as.numeric(MITO_COV_RATIO)) *
                            suppressWarnings(as.numeric(MEAN_AUTOSOMAL_COV)))],
              cov, by = "SAMPLE")
  re[, recomputed := 2 * chrM / AUTO_HQ_median]
  f <- agree_frac(re$MTDNA_CN, re$recomputed)
  if (!is.na(f) && f < 0.9) {
    warn(paste("MTDNA_CN agrees with 2 x chrM-mean / median for only %.0f%% of samples",
               "(median relative difference %.2f%%): the phenotype and this pipeline's",
               "normalisation use different denominators. See the HQ_MEDIAN_COV check above."),
         100 * f, 100 * median(abs(re$MTDNA_CN - re$recomputed) / re$recomputed, na.rm = TRUE))
  } else if (!is.na(f)) {
    note("[check] MTDNA_CN agrees with 2 x chrM-mean / AUTO_HQ_median for %.1f%% of %d samples",
         100 * f, nrow(re))
  }
}

pheno[, LOG2_MTDNA_CN := ifelse(is.finite(MTDNA_CN) & MTDNA_CN > 0, log2(MTDNA_CN), NA_real_)]

# Permuted null: shuffle the observed values among the samples that have one
# — or, for a second mode, copy the first mode's permutation by sample ID so
# every mode tests the same null values (a tree with a different sample set
# would otherwise get an unrelated permutation, and the cross-mode
# comparison of the null would measure that, not the mode).
pheno[, MTDNA_CN_NULL := NA_real_]
if (!is.null(opt$null_from) && file.exists(opt$null_from)) {
  src <- fread(opt$null_from, select = c("SAMPLE", "MTDNA_CN_NULL"))
  src[, SAMPLE := trimws(as.character(SAMPLE))]
  pheno[, MTDNA_CN_NULL := src$MTDNA_CN_NULL[match(SAMPLE, src$SAMPLE)]]
  note("[phenotypes] MTDNA_CN_NULL copied from %s for %d of %d samples", basename(opt$null_from),
       sum(is.finite(pheno$MTDNA_CN_NULL)), nrow(pheno))
} else {
  set.seed(opt$seed)
  have <- which(is.finite(pheno$MTDNA_CN))
  pheno$MTDNA_CN_NULL[have] <- sample(pheno$MTDNA_CN[have])
}

pheno[, SEX := fifelse(qc$INFERRED_SEX == "M", 1L,
                fifelse(qc$INFERRED_SEX == "F", 0L, NA_integer_))]
# The same sex as M/F, for the correction stage's ploidy model (which refuses
# a 0/1 coding as ambiguous).
pheno[, SEX_MF := fifelse(as.character(qc$INFERRED_SEX) %in% c("M", "F"),
                          as.character(qc$INFERRED_SEX), NA_character_)]

# Covariate candidates for models a user might add; unused by the shipped
# analyses.tsv but free to carry.
for (extra in c("MEAN_AUTOSOMAL_COV", "SUPERPOPULATION", "POPULATION")) {
  if (extra %in% names(qc)) pheno[[extra]] <- qc[[extra]]
}

# Genotype PCs from the preamble, when it ran. A left join: a sample without
# covariates keeps its phenotype and drops out of the adjusted models only.
if (!is.null(opt$covariates)) {
  cv <- fread(opt$covariates)
  if (!"SAMPLE" %in% names(cv)) stop(opt$covariates, " lacks a SAMPLE column", call. = FALSE)
  cv[, SAMPLE := as.character(SAMPLE)]
  keep <- c("SAMPLE", grep("^GPC[0-9]+$", names(cv), value = TRUE), intersect("GPC_PROJECTED", names(cv)))
  cv <- cv[, keep, with = FALSE]
  n_hit <- sum(pheno$SAMPLE %in% cv$SAMPLE)
  pheno <- merge(pheno, cv, by = "SAMPLE", all.x = TRUE, sort = FALSE)
  note("[covariates] %d genotype PCs for %d of %d samples merged from %s",
       length(keep) - 1L - ("GPC_PROJECTED" %in% keep), n_hit, nrow(pheno), basename(opt$covariates))
  if (n_hit < 0.9 * nrow(pheno)) warn("fewer than 90%% of samples have genotype PCs; adjusted models shrink accordingly")
}

write_if_changed(pheno, file.path(opt$out, "phenotypes.tsv"))
note("[phenotypes] %d samples (MTDNA_CN present for %d; SEX: %d M / %d F) -> phenotypes.tsv",
     nrow(pheno), n_cn, sum(pheno$SEX == 1L, na.rm = TRUE), sum(pheno$SEX == 0L, na.rm = TRUE))

# --- PC table --------------------------------------------------------------

pcs <- fread(opt$pcs)
if (!"SAMPLE" %in% names(pcs)) stop(opt$pcs, " lacks a SAMPLE column", call. = FALSE)
ids <- strip_suffix(pcs$SAMPLE, opt$suffix)
if (anyDuplicated(ids)) stop("stripping '", opt$suffix, "' left duplicated sample IDs", call. = FALSE)
pcs[, SAMPLE := as.character(ids)]
write_if_changed(pcs, file.path(opt$out, "svd.pcs.txt"))
note("[pcs] %d samples, %d PCs; suffix '%s' stripped from %d IDs -> svd.pcs.txt",
     nrow(pcs), sum(grepl("^PC[0-9]+$", names(pcs))), opt$suffix, attr(ids, "stripped"))

# --- consistency -----------------------------------------------------------
# The correction stage inner-joins PCs and coverage; a thin overlap there
# silently shrinks the cohort, so surface it here where the cause (an ID
# suffix mismatch) is still obvious.

overlap <- length(intersect(pcs$SAMPLE, cov$SAMPLE))
note("[align] PC/coverage overlap: %d of %d PC samples", overlap, nrow(pcs))
if (overlap == 0) stop("no overlap between PC and coverage sample IDs; check --suffix", call. = FALSE)
if (overlap < nrow(pcs)) {
  only_pc <- setdiff(pcs$SAMPLE, cov$SAMPLE)
  warn("%d PC sample(s) have no coverage median and will be dropped by the correction stage: %s",
       length(only_pc), paste(utils::head(only_pc, 5), collapse = ", "))
}

if (warned > 0L) note("[prepare] finished with %d warning(s) - read them before running", warned)
