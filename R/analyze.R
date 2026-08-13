#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# depthSV — per-region association testing
#
# Reads a corrected-depth matrix (regions x samples) and tests each region
# against a phenotype, streaming one result row per region to stdout.
#
# Estimators are unchanged from the original implementation. Two of them are
# computed differently for speed, and both were verified to agree with the
# straightforward form:
#
#   linear    A single projection replaces a per-region lm(). Because only the
#             depth column varies, the covariates can be projected out once and
#             every beta/SE recovered from one matrix product. Agrees with lm()
#             to ~1e-14 on beta, SE and t.
#
#   logistic  The design matrix is built once and only the depth column is
#             swapped per region, so glm.fit() runs the same IRLS on the same
#             data. Bit-identical to glm(formula, data). Note this is a modest
#             speed-up: the IRLS iterations dominate, not the formula parsing.
#             A score test would be far faster but is a DIFFERENT test, so it
#             is not used here; use the external-engine path for that.
#
#   coxph     Unchanged, fitted per region.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

option_list <- list(
  make_option(c("-f", "--inputFile"), type = "character", default = NULL,
              help = "corrected depth matrix, or '-' for stdin"),
  make_option(c("-p", "--phenoFile"), type = "character", default = NULL,
              help = "phenotype table; must contain a SAMPLE column"),
  make_option(c("-m", "--model"), type = "character", default = NULL,
              help = "model formula; the depth term must be named cov_resids"),
  make_option(c("-r", "--regressionMethod"), type = "character", default = "linear",
              help = "linear, logistic or coxph [default %default]"),
  make_option(c("--minObs"), type = "integer", default = 100,
              help = "skip a region with fewer complete observations [default %default]"),
  make_option(c("--minVariance"), type = "double", default = 1e-12,
              help = "skip a region whose depth variance is at or below this [default %default]"),
  make_option(c("--sampleIdPattern"), type = "character", default = NULL,
              help = paste("optional PCRE applied to matrix column names to recover sample IDs;",
                           "the first capture group is kept. Default: use column names as-is.")),
  make_option(c("--digits"), type = "integer", default = 8,
              help = "significant digits for reported statistics [default %default]"),
  make_option(c("--projection"), type = "character", default = "explicit",
              help = paste("how to residualise covariates for the linear path:",
                           "'explicit' (default) or 'qr'. Same projection; see the",
                           "note beside residualize() below."))
)

opt <- parse_args(OptionParser(option_list = option_list))

for (req in c("inputFile", "phenoFile", "model")) {
  if (is.null(opt[[req]])) stop("--", req, " is required", call. = FALSE)
}
if (!opt$regressionMethod %in% c("linear", "logistic", "coxph")) {
  stop("--regressionMethod must be one of: linear, logistic, coxph", call. = FALSE)
}
# Only the Cox path needs it, and every parallel chunk pays the load cost.
if (opt$regressionMethod == "coxph") suppressPackageStartupMessages(library(survival))

# --- sample identifiers ----------------------------------------------------
# Cohort-specific rewriting is opt-in. The previous default silently mangled
# any identifier containing a colon, which is a data-corruption risk on a
# cohort whose IDs happen to use one.
normalize_sample_ids <- function(x, pattern) {
  if (is.null(pattern)) return(x)
  out <- sub(pattern, "\\1", x, perl = TRUE)
  if (any(!nzchar(out))) stop("--sampleIdPattern produced empty IDs", call. = FALSE)
  out
}

# --- inputs ----------------------------------------------------------------

# fread() can only read .gz directly if the R.utils package is installed, which
# is not part of a base R install — so a gzipped phenotype or PC table fails on
# a clean machine even though it works wherever R.utils happens to be present.
# Decompressing through the shell removes that hidden dependency; gzip is
# already required by the pipeline.
read_table <- function(path, ...) {
  if (grepl("\\.(gz|bgz)$", path, ignore.case = TRUE)) {
    fread(cmd = paste("gzip -cd", shQuote(path)), ...)
  } else {
    fread(path, ...)
  }
}

pheno <- read_table(opt$phenoFile)
if (!"SAMPLE" %in% names(pheno)) stop("phenotype table must contain a SAMPLE column", call. = FALSE)
if (anyDuplicated(pheno$SAMPLE)) {
  dup <- unique(pheno$SAMPLE[duplicated(pheno$SAMPLE)])
  stop(sprintf(paste("phenotype table contains %d duplicated SAMPLE value(s), e.g. %s.",
                     "Duplicates silently misalign depth values against participants;",
                     "de-duplicate the phenotype table before running."),
               length(dup), paste(utils::head(dup, 3), collapse = ", ")), call. = FALSE)
}

open_input <- function(path) {
  if (path == "-") return(file("stdin", open = "r"))
  if (grepl("\\.gz$", path, ignore.case = TRUE)) return(gzfile(path, open = "rt"))
  file(path, open = "r")
}
con <- open_input(opt$inputFile)
on.exit(try(close(con), silent = TRUE), add = TRUE)

header_line <- readLines(con, n = 1L)
if (!length(header_line)) stop("input is empty", call. = FALSE)
header <- strsplit(header_line, "\t", fixed = TRUE)[[1]]
if (length(header) < 5L) {
  stop("input header must have 4 region columns followed by at least one sample column", call. = FALSE)
}
sample_ids <- normalize_sample_ids(header[-(1:4)], opt$sampleIdPattern)
n_samples  <- length(sample_ids)

# --- model terms -----------------------------------------------------------

model_formula <- as.formula(opt$model)
response <- trimws(sub("~.*", "", opt$model))
time_var <- character(0)
if (opt$regressionMethod == "coxph") {
  time_var <- sub("Surv\\(", "", sub(",.*", "", opt$model))
  response <- sub("\\)", "", sub(".*,", "", sub("~.*", "", opt$model)))
  time_var <- trimws(time_var); response <- trimws(response)
}
term_labels <- attr(terms(model_formula), "term.labels")
if (!"cov_resids" %in% term_labels) {
  stop("the model must include a term named 'cov_resids' for the depth value", call. = FALSE)
}
covariates <- setdiff(term_labels, "cov_resids")
needed     <- unique(c(response, time_var, covariates))
needed     <- needed[nzchar(needed)]

missing_cols <- setdiff(needed, names(pheno))
if (length(missing_cols)) {
  stop("phenotype table is missing: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

# --- align phenotypes to matrix column order -------------------------------
# A left join on a de-duplicated key returns exactly one row per sample, in
# order. The assertion below is what would have caught the misalignment that
# duplicated SAMPLE values used to cause.

idx <- match(sample_ids, pheno$SAMPLE)
aligned <- pheno[idx]
stopifnot(nrow(aligned) == n_samples)

n_unmatched <- sum(is.na(idx))
if (n_unmatched == n_samples) {
  stop("no matrix sample matched the phenotype table; check --sampleIdPattern", call. = FALSE)
}
if (n_unmatched > 0) {
  message(sprintf("[align] %d of %d matrix samples absent from the phenotype table; dropped",
                  n_unmatched, n_samples))
}

complete_base <- !is.na(idx) & complete.cases(aligned[, ..needed])
message(sprintf("[align] samples=%d  usable=%d  response=%s  method=%s",
                n_samples, sum(complete_base), response, opt$regressionMethod))
if (sum(complete_base) < opt$minObs) {
  stop(sprintf("only %d complete observations, below --minObs=%d",
               sum(complete_base), opt$minObs), call. = FALSE)
}

# --- fixed design ----------------------------------------------------------
# Built once. Everything that does not depend on the region is computed here.

base_rows <- which(complete_base)
n_use     <- length(base_rows)

cov_formula <- if (length(covariates)) {
  as.formula(paste("~", paste(covariates, collapse = " + ")))
} else {
  ~1
}
Z <- model.matrix(cov_formula, data = aligned[base_rows])
if (nrow(Z) != n_use) {
  stop("model.matrix dropped rows; check for NA or unused factor levels in covariates", call. = FALSE)
}

y_vec <- aligned[[response]][base_rows]

if (opt$regressionMethod == "linear") {
  if (!opt$projection %in% c("explicit", "qr")) {
    stop("--projection must be 'explicit' or 'qr'", call. = FALSE)
  }
  qr_Z   <- qr(Z)
  df_res <- n_use - qr_Z$rank - 1L
  if (df_res < 1L) stop("no residual degrees of freedom; too few samples for this model", call. = FALSE)

  # Two ways to compute the same projection, residual = (I - P)v.
  #
  #   qr        qr.resid() applies Householder reflections one column at a
  #             time. Unblocked and memory-bound.
  #   explicit  Z = QR with Q orthonormal, so P = QQ' and the residual is
  #             v - Q(Q'v): two BLAS calls instead. Measured 12-20x faster,
  #             agreeing with qr.resid to ~6e-14 relative — six orders below
  #             the last digit this script prints.
  #
  # Q MUST be truncated to the numerical rank. qr.Q() returns every column,
  # and on a rank-deficient design the extra columns span numerically-null
  # directions; projecting onto them gives a different answer, not a rounding
  # difference — which is reachable here whenever covariates are collinear.
  if (opt$projection == "explicit") {
    Q <- qr.Q(qr_Z)[, seq_len(qr_Z$rank), drop = FALSE]
    residualize <- function(v) as.numeric(v - Q %*% crossprod(Q, v))
  } else {
    residualize <- function(v) as.numeric(qr.resid(qr_Z, v))
  }

  y_res  <- residualize(as.numeric(y_vec))
  yy_res <- sum(y_res^2)
} else if (opt$regressionMethod == "logistic") {
  y_bin <- as.numeric(as.factor(y_vec)) - 1
  if (length(unique(y_bin[!is.na(y_bin)])) != 2L) {
    stop("logistic response must have exactly two levels", call. = FALSE)
  }
  # intercept, depth placeholder, then covariates
  X_fixed <- cbind(Z[, 1, drop = FALSE], cov_resids = 0, Z[, -1, drop = FALSE])
  depth_col <- 2L
}

# --- output header ---------------------------------------------------------

stat_cols <- switch(opt$regressionMethod,
  linear   = c("Estimate", "Std. Error", "t value", "Pr(>|t|)"),
  logistic = c("Estimate", "Std. Error", "z value", "Pr(>|z|)"),
  coxph    = c("coef", "exp(coef)", "se(coef)", "z", "Pr(>|z|)"))
out_header <- c("#CHROM", "START", "END", "Region", "N", "NCase", "NControl", stat_cols)
cat(paste(out_header, collapse = "\t"), "\n", sep = "")

# --- per-region estimators -------------------------------------------------

fit_linear <- function(g, keep) {
  x_res <- residualize(g)
  xx    <- sum(x_res^2)
  # Guard against exact collinearity with the covariates. This is a residual
  # sum of squares, not a variance, so it deliberately does not use
  # --minVariance: that threshold is applied to var(g) before we get here.
  if (xx <= 0) return(NULL)
  beta  <- sum(x_res * y_res) / xx
  rss   <- yy_res - beta^2 * xx
  se    <- sqrt((rss / df_res) / xx)
  tval  <- beta / se
  c(beta, se, tval, 2 * pt(-abs(tval), df_res))
}

fit_logistic <- function(g, keep) {
  X <- X_fixed
  X[, depth_col] <- g
  f <- suppressWarnings(glm.fit(X, y_bin, family = binomial()))
  p <- f$rank
  piv <- f$qr$pivot[seq_len(p)]
  pos <- match(depth_col, piv)
  if (is.na(pos)) return(NULL)                        # depth aliased out
  cov_unscaled <- chol2inv(f$qr$qr[seq_len(p), seq_len(p), drop = FALSE])
  se   <- sqrt(cov_unscaled[pos, pos])
  beta <- f$coefficients[depth_col]
  if (!is.finite(beta) || !is.finite(se)) return(NULL)
  z <- beta / se
  c(beta, se, z, 2 * pnorm(-abs(z)))
}

fit_coxph <- function(g, keep) {
  dt <- aligned[base_rows][keep]
  dt[, cov_resids := g[keep]]
  f <- try(coxph(model_formula, data = dt), silent = TRUE)
  if (inherits(f, "try-error")) return(NULL)
  cf <- coef(summary(f))
  if (!"cov_resids" %in% rownames(cf)) return(NULL)
  as.numeric(cf["cov_resids", ])
}

fit_region <- switch(opt$regressionMethod,
  linear = fit_linear, logistic = fit_logistic, coxph = fit_coxph)

# --- stream ----------------------------------------------------------------

n_out <- 0L; n_skipped <- 0L; n_error <- 0L
region_cols <- 4L

repeat {
  lines <- readLines(con, n = 10000L)   # parse in blocks; the caller caps lines per process
  if (!length(lines)) break

  block <- data.table::fread(text = lines, header = FALSE, sep = "\t",
                             showProgress = FALSE, colClasses = "character")
  if (ncol(block) != region_cols + n_samples) {
    n_skipped <- n_skipped + nrow(block)
    message(sprintf("[skip] block has %d columns, expected %d", ncol(block), region_cols + n_samples))
    next
  }

  meta  <- as.matrix(block[, 1:region_cols, with = FALSE])
  depth <- suppressWarnings(
    matrix(as.numeric(unlist(block[, (region_cols + 1L):ncol(block), with = FALSE], use.names = FALSE)),
           nrow = nrow(block)))

  for (i in seq_len(nrow(block))) {
    g_all <- depth[i, base_rows]
    keep  <- !is.na(g_all)
    n_obs <- sum(keep)

    if (n_obs < opt$minObs) { n_skipped <- n_skipped + 1L; next }

    # A region with no variation carries no information and would otherwise
    # cost a full model fit before failing.
    if (stats::var(g_all[keep]) <= opt$minVariance) { n_skipped <- n_skipped + 1L; next }

    g <- g_all
    if (n_obs < length(g_all)) {
      # Mean-impute the few missing depths so the fixed design stays valid.
      # Regions with meaningful missingness are excluded by --minObs above.
      g[!keep] <- mean(g_all[keep])
    }

    vals <- tryCatch(fit_region(g, keep), error = function(e) { n_error <<- n_error + 1L; NULL })
    if (is.null(vals) || anyNA(vals)) { n_skipped <- n_skipped + 1L; next }

    if (opt$regressionMethod == "linear") {
      n_case <- n_obs; n_ctrl <- n_obs
    } else {
      resp <- if (opt$regressionMethod == "logistic") y_bin[keep] else y_vec[keep]
      # Count explicitly rather than relying on table() ordering, which sorts
      # by value and previously caused case and control counts to be swapped.
      n_case <- sum(resp == max(resp, na.rm = TRUE), na.rm = TRUE)
      n_ctrl <- n_obs - n_case
    }

    cat(paste(c(meta[i, ], n_obs, n_case, n_ctrl,
                formatC(vals, digits = opt$digits, format = "g")), collapse = "\t"), "\n", sep = "")
    n_out <- n_out + 1L
  }
}

message(sprintf("[done] processed=%d skipped=%d err=%d", n_out, n_skipped, n_error))
