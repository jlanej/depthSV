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
#
# A region whose depth is missing for some samples (no expected copies:
# chrY in females) is fitted on the samples that have it, with the design
# refactorised for that subset; the subset is the same for every such
# region, so the refactorisation is cached. Nothing is imputed.
#
# Per-region QC, beyond --minObs and --minVariance: MAXSHARE is the largest
# share of the residualised depth's sum of squares carried by one sample. A
# region above --maxShare is a test of one participant — the read-depth
# analogue of a minor allele count of one — and is skipped.
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
  make_option("--rankInt", action = "store_true", default = FALSE,
              help = "rank-based inverse-normal transform of a quantitative response (linear only)"),
  make_option("--robust", action = "store_true", default = FALSE,
              help = "heteroskedasticity-robust (HC1) SE for linear; robust variance for coxph"),
  make_option("--perms", type = "integer", default = 0L,
              help = "permutations of the response for a max-|t| null over this input (linear only) [default %default]"),
  make_option("--permSeed", type = "integer", default = 1L,
              help = "seed of the permutations; every shard of a run must use the same [default %default]"),
  make_option("--permOut", type = "character", default = NULL,
              help = "file for the per-permutation max |t| (with --perms)"),
  make_option("--minObs", type = "integer", default = 100,
              help = "skip a region with fewer complete observations [default %default]"),
  make_option("--minCases", type = "integer", default = 20,
              help = "refuse a binary or survival phenotype with fewer cases or events, or controls [default %default]"),
  make_option("--minVariance", type = "double", default = 1e-12,
              help = "skip a region whose depth variance is at or below this [default %default]"),
  make_option("--maxShare", type = "double", default = 0.5,
              help = "skip a region where one sample carries more than this share of the residual depth SS [default %default]"),
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
if (opt$rankInt && opt$regressionMethod != "linear") stop("--rankInt applies to linear models only", call. = FALSE)
if (opt$perms < 0) stop("--perms must be >= 0", call. = FALSE)
if (opt$perms > 0 && opt$regressionMethod != "linear") stop("--perms is available for linear models only", call. = FALSE)
if (opt$perms > 0 && is.null(opt$permOut)) stop("--perms needs --permOut", call. = FALSE)
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
  for (pc in pc_cols) set(aligned, j = pc, value = pcs_dt[[pc]][pidx])
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

rank_int <- function(y) qnorm((rank(y, ties.method = "average") - 0.5) / length(y))

# The projection of a vector or a block of rows on the orthogonal complement
# of the covariate design. Q MUST be truncated to the numerical rank.
make_projector <- function(Zm) {
  q <- qr(Zm)
  Qm <- qr.Q(q)[, seq_len(q$rank), drop = FALSE]
  list(rank = q$rank, qr = q,
       vec   = if (opt$projection == "explicit") function(v) as.numeric(v - Qm %*% crossprod(Qm, v))
               else function(v) as.numeric(qr.resid(q, v)),
       block = if (opt$projection == "explicit") function(M) M - (M %*% Qm) %*% t(Qm)
               else function(M) t(apply(M, 1L, function(v) qr.resid(q, v))),
       cols  = if (opt$projection == "explicit") function(M) M - Qm %*% crossprod(Qm, M)
               else function(M) qr.resid(q, M))
}
proj <- make_projector(Z)

if (opt$regressionMethod == "linear") {
  if (opt$rankInt) y_vec <- rank_int(as.numeric(y_vec))
  df_res <- n_use - proj$rank - 1L
  if (df_res < 1L) stop("no residual degrees of freedom; too few samples for this model", call. = FALSE)
  y_res  <- proj$vec(as.numeric(y_vec))
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

# --- permutations (linear only) -------------------------------------------
# Freedman-Lane: the covariate-adjusted response is permuted, projected
# again, and the largest |t| over this input's regions is kept per
# permutation. The permutations come from --permSeed over the usable
# samples in matrix order, so every shard of a run permutes identically and
# the per-shard maxima combine into one genome-wide max-T distribution
# (scripts/export.sh). Computed on the classical t, also under --robust.
n_perm <- opt$perms
perm_idx <- NULL
if (n_perm > 0) {
  set.seed(opt$permSeed)
  perm_idx <- vapply(seq_len(n_perm), function(b) sample.int(n_use), integer(n_use))
}
perm_max <- rep(0, n_perm)
perm_columns <- function(y_res_s, proj_s, idx) {  # permuted, re-projected responses (n x B)
  Yp <- proj_s$cols(matrix(y_res_s[idx], nrow(idx), n_perm))
  list(Yp = Yp, yyp = colSums(Yp^2))
}

# Subset designs, for regions with missing samples. One-entry cache: the
# missing set is the same for every region of the same kind (chrY).
sub_cache <- list(mask = NULL)
subset_design <- function(mask) {
  if (!is.null(sub_cache$mask) && identical(sub_cache$mask, mask)) return(sub_cache)
  s <- list(mask = mask, n = sum(mask), proj = make_projector(Z[mask, , drop = FALSE]))
  if (opt$regressionMethod == "linear") {
    s$df <- s$n - s$proj$rank - 1L
    s$y_res <- s$proj$vec(as.numeric(y_vec[mask]))
    s$yy <- sum(s$y_res^2)
    if (n_perm > 0 && s$df >= 1L) {
      # The global permutation restricted to the subset, as ranks, is a
      # permutation of the subset — the same one in every shard.
      idx_s <- apply(perm_idx[mask, , drop = FALSE], 2L, rank, ties.method = "first")
      s[c("Yp", "yyp")] <- perm_columns(s$y_res, s$proj, idx_s)
    }
  } else if (opt$regressionMethod == "logistic") {
    s$X_fixed <- X_fixed[mask, , drop = FALSE]
    s$y <- y_bin[mask]
    s$dev_null <- suppressWarnings(glm.fit(s$X_fixed[, -depth_col, drop = FALSE], s$y, family = binomial()))$deviance
  }
  sub_cache <<- s
  s
}
full_design <- list(mask = rep(TRUE, n_use), n = n_use, proj = proj)
if (opt$regressionMethod == "linear") {
  full_design$df <- df_res; full_design$y_res <- y_res; full_design$yy <- yy_res
  if (n_perm > 0) full_design[c("Yp", "yyp")] <- perm_columns(y_res, proj, perm_idx)
} else if (opt$regressionMethod == "logistic") {
  full_design$X_fixed <- X_fixed; full_design$y <- y_bin
  full_design$dev_null <- suppressWarnings(glm.fit(X_fixed[, -depth_col, drop = FALSE], y_bin, family = binomial()))$deviance
}

# --- output header ---------------------------------------------------------

stat_cols <- switch(opt$regressionMethod,
  linear   = c("BETA", "SE", "STAT", "P", "LOG10P"),
  logistic = c("BETA", "SE", "STAT", "P", "LOG10P", "LRT_P", "CONVERGED"),
  coxph    = c("BETA", "HR", "SE", "STAT", "P", "LOG10P", "CONVERGED"))
out_header <- c("#CHROM", "START", "END", "Region", "N", "NCase", "NControl", stat_cols, "MAXSHARE")
if (!opt$skipOutputHeader) cat(paste(out_header, collapse = "\t"), "\n", sep = "")

fmt <- function(x) formatC(x, digits = opt$digits, format = "g")
log10p_t <- function(tval, df) -(pt(-abs(tval), df, log.p = TRUE) + log(2)) / log(10)
log10p_z <- function(z) -(pnorm(-abs(z), log.p = TRUE) + log(2)) / log(10)

# --- per-region estimators -------------------------------------------------

# Linear fit of one region on a design s; g has no NA within s$mask.
fit_linear_one <- function(g, s) {
  x_res <- s$proj$vec(g)
  xx <- sum(x_res^2)
  if (xx <= 0) return(NULL)
  share <- max(x_res^2) / xx
  beta <- sum(x_res * s$y_res) / xx
  e <- s$y_res - beta * x_res
  se <- if (opt$robust) sqrt(sum(x_res^2 * e^2) / xx^2 * s$n / s$df) else sqrt((sum(e^2) / s$df) / xx)
  tval <- beta / se
  pm <- NULL
  if (n_perm > 0 && !is.null(s$Yp)) {
    bp <- as.numeric(crossprod(s$Yp, x_res)) / xx
    pm <- abs(bp) / sqrt(((s$yyp - bp^2 * xx) / s$df) / xx)
  }
  list(vals = c(beta, se, tval, 2 * pt(-abs(tval), s$df), log10p_t(tval, s$df)), share = share, perm = pm)
}

fit_logistic_one <- function(g, s) {
  x_res <- s$proj$vec(g); xx <- sum(x_res^2)
  if (xx <= 0) return(NULL)
  share <- max(x_res^2) / xx
  X <- s$X_fixed
  X[, depth_col] <- g
  f <- suppressWarnings(glm.fit(X, s$y, family = binomial()))
  p <- f$rank
  piv <- f$qr$pivot[seq_len(p)]
  pos <- match(depth_col, piv)
  if (is.na(pos)) return(NULL)                        # depth aliased out
  cov_unscaled <- chol2inv(f$qr$qr[seq_len(p), seq_len(p), drop = FALSE])
  se   <- sqrt(cov_unscaled[pos, pos])
  beta <- f$coefficients[depth_col]
  if (!is.finite(beta) || !is.finite(se)) return(NULL)
  z <- beta / se
  lrt <- max(0, s$dev_null - f$deviance)
  list(vals = c(beta, se, z, 2 * pnorm(-abs(z)), log10p_z(z), pchisq(lrt, 1, lower.tail = FALSE),
                as.integer(f$converged)), share = share)
}

fit_coxph_one <- function(g, s) {
  x_res <- s$proj$vec(g); xx <- sum(x_res^2)
  if (xx <= 0) return(NULL)
  share <- max(x_res^2) / xx
  dt <- aligned[base_rows][s$mask]
  dt[, cov_resids := g]
  converged <- 1L
  f <- withCallingHandlers(
    try(coxph(model_formula, data = dt, robust = opt$robust), silent = TRUE),
    warning = function(w) { converged <<- 0L; invokeRestart("muffleWarning") })
  if (inherits(f, "try-error")) return(NULL)
  cf <- coef(summary(f))
  if (!"cov_resids" %in% rownames(cf)) return(NULL)
  v <- as.numeric(cf["cov_resids", ])
  # coef exp(coef) se(coef) [robust se] z p — the z and p are the last two
  z <- v[length(v) - 1L]; pv <- v[length(v)]; se <- v[3L + as.integer(opt$robust)]
  list(vals = c(v[1], v[2], se, z, pv, log10p_z(z), converged), share = share)
}

fit_one <- switch(opt$regressionMethod, linear = fit_linear_one, logistic = fit_logistic_one, coxph = fit_coxph_one)

emit <- function(meta_row, n_obs, s, r) {
  if (opt$regressionMethod == "linear") { nc <- n_obs; nk <- n_obs }
  else {
    resp <- if (opt$regressionMethod == "logistic") s$y else aligned[[response]][base_rows][s$mask]
    nc <- sum(resp == 1L); nk <- n_obs - nc
  }
  vals <- r$vals
  conv <- if (opt$regressionMethod == "linear") character(0) else as.character(as.integer(vals[length(vals)]))
  if (opt$regressionMethod != "linear") vals <- vals[-length(vals)]
  cat(paste(c(meta_row, n_obs, nc, nk, fmt(vals), conv, fmt(r$share)), collapse = "\t"), "\n", sep = "")
}

# --- stream ----------------------------------------------------------------

n_out <- 0L; n_skipped <- 0L; n_share <- 0L; n_error <- 0L
region_cols <- 4L

repeat {
  lines <- readLines(con, n = 10000L)   # parse in blocks; the caller caps lines per process
  if (!length(lines)) break

  block <- data.table::fread(text = lines, header = FALSE, sep = "\t", showProgress = FALSE,
                             colClasses = list(character = seq_len(region_cols)), na.strings = c("NA", ""))
  if (ncol(block) != region_cols + n_samples) {
    n_skipped <- n_skipped + nrow(block)
    message(sprintf("[skip] block has %d columns, expected %d", ncol(block), region_cols + n_samples))
    next
  }

  meta  <- as.matrix(block[, seq_len(region_cols), with = FALSE])
  depth <- as.matrix(block[, (region_cols + 1L):ncol(block), with = FALSE])
  storage.mode(depth) <- "double"
  depth <- depth[, base_rows, drop = FALSE]        # usable samples only, in design order

  n_obs <- rowSums(!is.na(depth))
  row_var <- apply(depth, 1L, function(v) stats::var(v[!is.na(v)]))
  ok <- n_obs >= opt$minObs & is.finite(row_var) & row_var > opt$minVariance
  n_skipped <- n_skipped + sum(!ok)
  complete <- ok & n_obs == n_use

  # Complete regions of a linear model: one block projection for all of them.
  if (opt$regressionMethod == "linear" && any(complete)) {
    rows <- which(complete)
    Dr   <- proj$block(depth[rows, , drop = FALSE])
    xx   <- rowSums(Dr^2)
    share <- apply(Dr^2, 1L, max) / xx
    beta <- as.numeric(Dr %*% y_res) / xx
    E    <- matrix(y_res, nrow = length(rows), ncol = n_use, byrow = TRUE) - beta * Dr
    se   <- if (opt$robust) sqrt(rowSums(Dr^2 * E^2) / xx^2 * n_use / df_res)
            else sqrt((rowSums(E^2) / df_res) / xx)
    tval <- beta / se
    pval <- 2 * pt(-abs(tval), df_res)
    lp   <- log10p_t(tval, df_res)
    good <- xx > 0 & is.finite(tval)
    over <- good & share > opt$maxShare
    n_share <- n_share + sum(over)
    n_skipped <- n_skipped + sum(!good)
    keepi <- which(good & !over)
    if (length(keepi)) {
      out <- cbind(meta[rows[keepi], , drop = FALSE], n_obs[rows[keepi]], n_obs[rows[keepi]], n_obs[rows[keepi]],
                   fmt(beta[keepi]), fmt(se[keepi]), fmt(tval[keepi]), fmt(pval[keepi]), fmt(lp[keepi]),
                   fmt(share[keepi]))
      writeLines(apply(out, 1L, paste, collapse = "\t"))
      n_out <- n_out + length(keepi)
      if (n_perm > 0) {
        # One product for every region x permutation; the tested regions only.
        bp   <- (Dr[keepi, , drop = FALSE] %*% full_design$Yp) / xx[keepi]
        rssp <- matrix(full_design$yyp, nrow = length(keepi), ncol = n_perm, byrow = TRUE) - bp^2 * xx[keepi]
        tp   <- abs(bp) / sqrt((rssp / df_res) / xx[keepi])
        perm_max <- pmax(perm_max, apply(tp, 2L, max))
      }
    }
    todo <- which(ok & !complete)
  } else {
    todo <- which(ok)
  }

  # Everything else — logistic and Cox regions, and regions with missing
  # samples — one at a time on the matching design.
  for (i in todo) {
    g_all <- depth[i, ]
    mask <- !is.na(g_all)
    s <- if (all(mask)) full_design else subset_design(mask)
    if (!all(mask) && (opt$regressionMethod != "linear" && s$n < opt$minObs)) { n_skipped <- n_skipped + 1L; next }
    if (opt$regressionMethod == "linear" && !all(mask) && s$df < 1L) { n_skipped <- n_skipped + 1L; next }
    r <- tryCatch(fit_one(g_all[mask], s), error = function(e) { n_error <<- n_error + 1L; NULL })
    if (is.null(r) || anyNA(r$vals)) { n_skipped <- n_skipped + 1L; next }
    if (r$share > opt$maxShare) { n_share <- n_share + 1L; next }
    emit(meta[i, ], n_obs[i], s, r)
    n_out <- n_out + 1L
    if (n_perm > 0 && !is.null(r$perm)) perm_max <- pmax(perm_max, r$perm)
  }
}

if (n_perm > 0) {
  writeLines(c(sprintf("#perms=%d", n_perm), sprintf("#seed=%d", opt$permSeed), sprintf("#df=%d", df_res),
               sprintf("#regions=%d", n_out), "#stat=classical |t| under Freedman-Lane permutation of the response",
               "perm\tmax_abs_stat",
               sprintf("%d\t%s", seq_len(n_perm), formatC(perm_max, digits = 8, format = "g"))),
             opt$permOut)
}

message(sprintf("[done] processed=%d skipped=%d single_sample=%d err=%d", n_out, n_skipped, n_share, n_error))
