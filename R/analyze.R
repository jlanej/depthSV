#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# depthSV — per-region association testing
#
# Reads a corrected-depth matrix (regions x samples) and tests each region
# against a phenotype, streaming one result row per region to stdout.
#
# The design is built ONCE per process and only the depth column changes per
# region:
#
#   linear    One projection replaces a per-region lm(). The covariates —
#             which include the coverage PCs the correction stage removed —
#             are projected out of the phenotype once, and a whole block of
#             regions is projected and fitted with two matrix products.
#             Agrees with lm() to ~1e-14 on beta, SE and t.
#
#   logistic  The design is built once and only the depth column is swapped,
#             so glm.fit() runs the same IRLS on the same data as glm(). The
#             Wald test collapses under complete separation (Hauck-Donner),
#             so the likelihood-ratio p-value and the convergence flag are
#             reported beside it; the null deviance is computed once.
#
#   coxph     Fitted per region; a convergence flag is reported.
#
# Why the PCs must be in the model. correct.R residualises depth against
# PC1..PCk. If those PCs are absent from the association model, the test is
# deflated by roughly 1 - R^2(phenotype ~ PCs), and any bin correlated with a
# covariate that itself correlates with a PC is BIASED (every sex-chromosome
# bin, once sex is a covariate). With the PCs in the model the projection
# is the same on both sides and Frisch-Waugh holds exactly. Hence --pcs and
# --ndim: pass the table and count the correction used.
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
  make_option(c("-i", "--pcs"), type = "character", default = NULL,
              help = "PC table (SAMPLE, PC1..PCn) the correction stage used"),
  make_option(c("-d", "--ndim"), type = "integer", default = 0L,
              help = "PCs removed by the correction; PC1..PCndim join the covariates [default %default]"),
  make_option("--caseLevel", type = "character", default = NULL,
              help = "for a character/factor binary response, the level that is the case"),
  make_option("--minObs", type = "integer", default = 100,
              help = "skip a region with fewer complete observations [default %default]"),
  make_option("--minCases", type = "integer", default = 20,
              help = "refuse a binary or survival phenotype with fewer cases or events, or controls [default %default]"),
  make_option("--minVariance", type = "double", default = 1e-12,
              help = "skip a region whose depth variance is at or below this [default %default]"),
  make_option("--sampleIdPattern", type = "character", default = NULL,
              help = paste("optional PCRE applied to matrix column names to recover sample IDs;",
                           "the first capture group is kept. Default: use column names as-is.")),
  make_option("--skipOutputHeader", action = "store_true", default = FALSE,
              help = "suppress the output header (for parallel chunks)"),
  make_option("--digits", type = "integer", default = 8,
              help = "significant digits for reported statistics [default %default]"),
  make_option("--projection", type = "character", default = "explicit",
              help = "how to residualise covariates for the linear path: 'explicit' (default) or 'qr'")
)

opt <- parse_args(OptionParser(option_list = option_list))

for (req in c("inputFile", "phenoFile", "model")) {
  if (is.null(opt[[req]])) stop("--", req, " is required", call. = FALSE)
}
if (!opt$regressionMethod %in% c("linear", "logistic", "coxph")) {
  stop("--regressionMethod must be one of: linear, logistic, coxph", call. = FALSE)
}
if (opt$ndim < 0) stop("--ndim must be >= 0", call. = FALSE)
if (opt$ndim > 0 && is.null(opt$pcs)) {
  stop("--ndim ", opt$ndim, " needs --pcs: the correction removed those PCs from the depth, ",
       "so the test must condition on them", call. = FALSE)
}
if (!opt$projection %in% c("explicit", "qr")) stop("--projection must be 'explicit' or 'qr'", call. = FALSE)
# Only the Cox path needs it, and every parallel chunk pays the load cost.
if (opt$regressionMethod == "coxph") suppressPackageStartupMessages(library(survival))

normalize_sample_ids <- function(x, pattern) {
  if (is.null(pattern)) return(x)
  out <- sub(pattern, "\\1", x, perl = TRUE)
  if (any(!nzchar(out))) stop("--sampleIdPattern produced empty IDs", call. = FALSE)
  out
}

# fread() can only read .gz directly if the R.utils package is installed, which
# is not part of a base R install. Decompress through the shell instead.
read_table <- function(path, ...) {
  if (grepl("\\.(gz|bgz)$", path, ignore.case = TRUE)) {
    fread(cmd = paste("gzip -cd", shQuote(path)), ...)
  } else {
    fread(path, ...)
  }
}

refuse_duplicates <- function(ids, what) {
  if (anyDuplicated(ids)) {
    dup <- unique(ids[duplicated(ids)])
    stop(sprintf(paste("%s contains %d duplicated sample ID(s), e.g. %s.",
                       "Duplicates silently misalign depth values against participants."),
                 what, length(dup), paste(utils::head(dup, 3), collapse = ", ")), call. = FALSE)
  }
}

# --- inputs ----------------------------------------------------------------

pheno <- read_table(opt$phenoFile)
if (!"SAMPLE" %in% names(pheno)) stop("phenotype table must contain a SAMPLE column", call. = FALSE)
pheno[, SAMPLE := as.character(SAMPLE)]
refuse_duplicates(pheno$SAMPLE, "the phenotype table")

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
# match() is many-to-one: two matrix columns mapping to one participant would
# both be kept and that participant counted twice.
refuse_duplicates(sample_ids, "the matrix header (after --sampleIdPattern)")

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

# The coverage PCs join the covariates. They are matched to the matrix
# columns like the phenotypes, and a sample without them drops out.
pc_cols <- character(0)
if (opt$ndim > 0) {
  pc_cols <- paste0("PC", seq_len(opt$ndim))
  pcs_dt <- read_table(opt$pcs, select = c("SAMPLE", pc_cols))
  if (!"SAMPLE" %in% names(pcs_dt)) stop("PC table must contain a SAMPLE column", call. = FALSE)
  missing_pc <- setdiff(pc_cols, names(pcs_dt))
  if (length(missing_pc)) {
    stop(sprintf("PC table has %d of the %d requested components", opt$ndim - length(missing_pc), opt$ndim),
         call. = FALSE)
  }
  pcs_dt[, SAMPLE := as.character(SAMPLE)]
  refuse_duplicates(pcs_dt$SAMPLE, "the PC table")
  clash <- intersect(pc_cols, names(pheno))
  if (length(clash)) {
    stop("phenotype table already has column(s) ", paste(clash, collapse = ", "),
         "; the coverage PCs use those names", call. = FALSE)
  }
  covariates <- c(covariates, pc_cols)
  if (opt$regressionMethod == "coxph") {
    model_formula <- as.formula(paste(opt$model, "+", paste(pc_cols, collapse = " + ")))
  }
}

needed <- unique(c(response, time_var, covariates))
needed <- needed[nzchar(needed)]
missing_cols <- setdiff(setdiff(needed, pc_cols), names(pheno))
if (length(missing_cols)) {
  stop("phenotype table is missing: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

# --- align phenotypes (and PCs) to matrix column order ---------------------

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

if (length(pc_cols)) {
  pidx <- match(sample_ids, pcs_dt$SAMPLE)
  n_nopc <- sum(is.na(pidx) & !is.na(idx))
  if (n_nopc > 0) message(sprintf("[align] %d matrix samples lack PCs; dropped", n_nopc))
  for (pc in pc_cols) aligned[[pc]] <- pcs_dt[[pc]][pidx]
}

complete_base <- !is.na(idx) & complete.cases(aligned[, ..needed])
message(sprintf("[align] samples=%d  usable=%d  response=%s  method=%s  pcs=%d",
                n_samples, sum(complete_base), response, opt$regressionMethod, length(pc_cols)))
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

# A binary response is accepted only when its coding is unambiguous: numeric
# 0/1, logical, or a character/factor with --caseLevel naming the case.
# Anything else used to be recoded by sort order, which made whichever label
# sorted first the non-event and silently reversed effect directions.
encode_binary <- function(y, what) {
  if (is.logical(y)) return(as.integer(y))
  if (is.numeric(y)) {
    vals <- sort(unique(as.numeric(y[!is.na(y)])))
    if (identical(vals, c(0, 1))) return(as.integer(y))
    if (length(vals) == 1L) stop(what, " takes a single value (", vals, ")", call. = FALSE)
    stop(what, " is numeric but not 0/1 (values: ", paste(utils::head(vals, 5), collapse = ", "),
         "); recode it as 0 = control, 1 = case", call. = FALSE)
  }
  y <- as.character(y)
  lv <- sort(unique(y[!is.na(y)]))
  if (length(lv) != 2L) stop(what, " must have exactly two levels; found: ", paste(lv, collapse = ", "), call. = FALSE)
  if (is.null(opt$caseLevel)) {
    stop(what, " is coded ", paste(lv, collapse = "/"), "; pass --caseLevel <label> to say which is the case",
         call. = FALSE)
  }
  if (!opt$caseLevel %in% lv) stop("--caseLevel '", opt$caseLevel, "' is not a level of ", what, call. = FALSE)
  as.integer(y == opt$caseLevel)
}

if (opt$regressionMethod == "linear") {
  qr_Z   <- qr(Z)
  df_res <- n_use - qr_Z$rank - 1L
  if (df_res < 1L) stop("no residual degrees of freedom; too few samples for this model", call. = FALSE)

  # Two ways to compute the same projection, residual = (I - P)v: 'qr' uses
  # qr.resid() (Householder, unblocked); 'explicit' forms the orthonormal Q
  # and computes v - Q(Q'v), 12-20x faster and agreeing to ~6e-14. Q MUST be
  # truncated to the numerical rank: qr.Q() returns every column, and on a
  # rank-deficient design the extra columns span numerically-null directions.
  Q <- qr.Q(qr_Z)[, seq_len(qr_Z$rank), drop = FALSE]
  residualize <- if (opt$projection == "explicit") {
    function(v) as.numeric(v - Q %*% crossprod(Q, v))
  } else {
    function(v) as.numeric(qr.resid(qr_Z, v))
  }
  # The same projection applied to a block of regions at once (rows = regions).
  residualize_block <- if (opt$projection == "explicit") {
    function(M) M - (M %*% Q) %*% t(Q)
  } else {
    function(M) t(apply(M, 1L, function(v) qr.resid(qr_Z, v)))
  }

  y_res  <- residualize(as.numeric(y_vec))
  yy_res <- sum(y_res^2)
} else if (opt$regressionMethod == "logistic") {
  y_bin <- encode_binary(y_vec, response)
  n_case <- sum(y_bin == 1L); n_ctrl <- sum(y_bin == 0L)
  if (min(n_case, n_ctrl) < opt$minCases) {
    stop(sprintf("%s has %d cases and %d controls; --minCases=%d", response, n_case, n_ctrl, opt$minCases),
         call. = FALSE)
  }
  # intercept, depth placeholder, then covariates
  X_fixed <- cbind(Z[, 1, drop = FALSE], cov_resids = 0, Z[, -1, drop = FALSE])
  depth_col <- 2L
  # The null deviance (no depth term) is the same for every region.
  null_fit <- suppressWarnings(glm.fit(X_fixed[, -depth_col, drop = FALSE], y_bin, family = binomial()))
  dev_null <- null_fit$deviance
} else {
  ev <- aligned[[response]][base_rows]
  ev_bin <- encode_binary(ev, response)
  if (min(sum(ev_bin == 1L), sum(ev_bin == 0L)) < opt$minCases) {
    stop(sprintf("%s has %d events and %d censored; --minCases=%d", response,
                 sum(ev_bin == 1L), sum(ev_bin == 0L), opt$minCases), call. = FALSE)
  }
  aligned[[response]] <- NA_integer_
  aligned[[response]][base_rows] <- ev_bin
}

# --- output header ---------------------------------------------------------

stat_cols <- switch(opt$regressionMethod,
  linear   = c("BETA", "SE", "STAT", "P", "LOG10P"),
  logistic = c("BETA", "SE", "STAT", "P", "LOG10P", "LRT_P", "CONVERGED"),
  coxph    = c("BETA", "HR", "SE", "STAT", "P", "LOG10P", "CONVERGED"))
out_header <- c("#CHROM", "START", "END", "Region", "N", "NCase", "NControl", stat_cols)
if (!opt$skipOutputHeader) cat(paste(out_header, collapse = "\t"), "\n", sep = "")

fmt <- function(x) formatC(x, digits = opt$digits, format = "g")
log10p_t <- function(tval, df) -(pt(-abs(tval), df, log.p = TRUE) + log(2)) / log(10)
log10p_z <- function(z) -(pnorm(-abs(z), log.p = TRUE) + log(2)) / log(10)

# --- per-region estimators (logistic and Cox) ------------------------------

fit_logistic <- function(g) {
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
  lrt <- max(0, dev_null - f$deviance)
  c(beta, se, z, 2 * pnorm(-abs(z)), log10p_z(z), pchisq(lrt, 1, lower.tail = FALSE), as.integer(f$converged))
}

fit_coxph <- function(g, keep) {
  dt <- aligned[base_rows][keep]
  dt[, cov_resids := g[keep]]
  converged <- 1L
  f <- withCallingHandlers(
    try(coxph(model_formula, data = dt), silent = TRUE),
    warning = function(w) { converged <<- 0L; invokeRestart("muffleWarning") })
  if (inherits(f, "try-error")) return(NULL)
  cf <- coef(summary(f))
  if (!"cov_resids" %in% rownames(cf)) return(NULL)
  v <- as.numeric(cf["cov_resids", ])               # coef exp(coef) se(coef) z p
  c(v[1], v[2], v[3], v[4], v[5], log10p_z(v[4]), converged)
}

# --- stream ----------------------------------------------------------------

n_out <- 0L; n_skipped <- 0L; n_error <- 0L
region_cols <- 4L

repeat {
  lines <- readLines(con, n = 10000L)   # parse in blocks; the caller caps lines per process
  if (!length(lines)) break

  block <- data.table::fread(text = lines, header = FALSE, sep = "\t", showProgress = FALSE,
                             colClasses = list(character = seq_len(region_cols)))
  if (ncol(block) != region_cols + n_samples) {
    n_skipped <- n_skipped + nrow(block)
    message(sprintf("[skip] block has %d columns, expected %d", ncol(block), region_cols + n_samples))
    next
  }

  meta  <- as.matrix(block[, seq_len(region_cols), with = FALSE])
  depth <- as.matrix(block[, (region_cols + 1L):ncol(block), with = FALSE])
  storage.mode(depth) <- "double"
  depth <- depth[, base_rows, drop = FALSE]        # usable samples only, in design order

  # Per-region QC: enough complete observations, and some variation. Missing
  # depths are mean-imputed so the fixed design stays valid; regions with
  # meaningful missingness are excluded by --minObs.
  n_obs <- rowSums(!is.na(depth))
  row_mean <- rowMeans(depth, na.rm = TRUE)
  na_pos <- which(is.na(depth), arr.ind = TRUE)
  if (nrow(na_pos)) depth[na_pos] <- row_mean[na_pos[, 1]]
  row_var <- apply(depth, 1L, stats::var)
  ok <- n_obs >= opt$minObs & is.finite(row_var) & row_var > opt$minVariance
  n_skipped <- n_skipped + sum(!ok)

  if (opt$regressionMethod == "linear") {
    if (any(ok)) {
      Dr   <- residualize_block(depth[ok, , drop = FALSE])
      xx   <- rowSums(Dr^2)
      beta <- as.numeric(Dr %*% y_res) / xx
      rss  <- yy_res - beta^2 * xx
      se   <- sqrt((rss / df_res) / xx)
      tval <- beta / se
      pval <- 2 * pt(-abs(tval), df_res)
      lp   <- log10p_t(tval, df_res)
      good <- xx > 0 & is.finite(tval)             # exact collinearity with the covariates
      n_skipped <- n_skipped + sum(!good)
      rows <- which(ok)[good]
      out <- cbind(meta[rows, , drop = FALSE], n_obs[rows], n_obs[rows], n_obs[rows],
                   fmt(beta[good]), fmt(se[good]), fmt(tval[good]), fmt(pval[good]), fmt(lp[good]))
      if (nrow(out)) writeLines(apply(out, 1L, paste, collapse = "\t"))
      n_out <- n_out + nrow(out)
    }
  } else {
    for (i in which(ok)) {
      g <- depth[i, ]
      keep <- rep(TRUE, length(g))
      vals <- tryCatch(
        if (opt$regressionMethod == "logistic") fit_logistic(g) else fit_coxph(g, keep),
        error = function(e) { n_error <<- n_error + 1L; NULL })
      if (is.null(vals) || anyNA(vals)) { n_skipped <- n_skipped + 1L; next }
      resp <- if (opt$regressionMethod == "logistic") y_bin else aligned[[response]][base_rows]
      nc <- sum(resp == 1L); nk <- n_obs[i] - nc
      cat(paste(c(meta[i, ], n_obs[i], nc, nk, fmt(vals[-length(vals)]), as.integer(vals[length(vals)])),
                collapse = "\t"), "\n", sep = "")
      n_out <- n_out + 1L
    }
  }
}

message(sprintf("[done] processed=%d skipped=%d err=%d", n_out, n_skipped, n_error))
