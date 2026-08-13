#!/usr/bin/env Rscript
# =============================================================================
# analyze.nullmodel.R
# -----------------------------------------------------------------------------
# Exact GENESIS single-variant score test, applied to depth-based "genotype"
# vectors (one corrected-coverage value per sample, per region), using a null
# model object previously saved by GENESIS::fitNullModel(), passed via
# --nullModelFile. By convention these live at
#   <nullmodel_root>/<TRAIT>/analysis_null_model_invnorm.RData   (Gaussian)
#   <nullmodel_root>/<TRAIT>/analysis_null_model.RData           (binary)
#
# Why this script rather than the plain analyze.R fixed-effects fit?
# -----------------------------------------------------------------------------
# A plain OLS / GLM against a phenotype table can mimic the *fixed effect
# covariates* only. It can NOT reproduce:
#   (a) the 4th-degree sparse empirical kinship matrix (PC-Relate), and
#   (b) the heterogeneous-by-(study x ancestry) residual variance,
#   (c) the rank-based inverse-normal transform of stage-1 marginal residuals
#       rescaled by the pre-transform variance.
# All three are baked into nullmod$cholSigmaInv / nullmod$fit$resid.PY /
# nullmod${CX,CXCXI,RSS0}. Re-using that object gives an association test that
# is mathematically identical to the stage-2 score test the null model was
# fitted for -- the only thing that changes is the column of "genotype" plugged
# in (here: a batch-corrected coverage residual instead of an SV dosage).
#
# Math (lifted verbatim from GENESIS R/testGeno.R::.testGenoSingleVarScore and
# R/nullModelTestPrep.R::calcGtilde):
#   CG      = t(cholSigmaInv) %*% G                   # n x 1, sigma^{-1/2}G
#   Gtilde  = CG - CXCXI %*% (t(CG) %*% CX)           # residualise G in Sigma^{-1} metric
#   GPG     = sum(Gtilde^2)                           # G' P G
#   score   = sum(G * resid.PY)                       # G' P Y, resid.PY precomputed
#   Stat    = score / sqrt(GPG)
#   pval    = pchisq(Stat^2, df = 1, lower.tail=FALSE)
#   Est     = score / GPG                             # beta hat
#   Est.SE  = 1 / sqrt(GPG)
#   PVE     = Stat^2 / RSS0
#
# These formulas are valid for BOTH the LMM (Gaussian, two-stage INT) and the
# GMMAT logistic case (binary BASO outcome). In both cases GENESIS' fitNullModel
# returns the same set of test-prep fields, and the score statistic
# distribution under the null is chi^2_1.
#
# Missing-G handling: per the manuscript ("missing SV genotype calls were
# imputed to the mean before performing the association tests"), we mean-impute
# missing coverage values across the nullmod sample set. Regions with too few
# observed samples are skipped.
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(Matrix)
})

option_list <- list(
  make_option(c("-f", "--inputFile"), type = "character",
              help = "corrected-depth file (- for stdin)", default = "-"),
  make_option(c("-n", "--nullModelFile"), type = "character",
              help = "RData containing GENESIS nullmod (must contain object 'nullmod')",
              default = NULL),
  make_option(c("-o", "--outcome"), type = "character",
              help = "optional outcome label to emit in the OUTCOME column",
              default = NA_character_),
  make_option(c("--minObs"), type = "integer",
              help = "skip region if fewer than this many non-missing nullmod samples have coverage",
              default = 100)
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$nullModelFile)) stop("--nullModelFile is required")
if (!file.exists(opt$nullModelFile)) stop("nullModelFile not found: ", opt$nullModelFile)

# -------------------------------------------------------------------------
# 1. Load null model and unpack test-prep matrices.
# -------------------------------------------------------------------------
env <- new.env()
load(opt$nullModelFile, envir = env)
if (!"nullmod" %in% ls(env)) {
  stop("Expected an object named 'nullmod' inside ", opt$nullModelFile,
       "; found: ", paste(ls(env), collapse = ","))
}
nullmod <- env$nullmod

required <- c("cholSigmaInv", "CX", "CXCXI", "RSS0", "fit")
miss <- setdiff(required, names(nullmod))
if (length(miss)) stop("nullmod is missing required GENESIS fields: ", paste(miss, collapse = ","))
if (!all(c("resid.PY", "sample.id") %in% colnames(nullmod$fit)))
  stop("nullmod$fit must contain columns 'resid.PY' and 'sample.id' (GENESIS >= 2.18)")

C        <- nullmod$cholSigmaInv         # n x n sparse triangular (dtCMatrix)
# Pre-compute t(C) ONCE here.  Inside the per-region loop we used to call
# crossprod(C, G), which is equivalent to `t(C) %*% G`. For a dtCMatrix that
# allocates a brand-new transposed sparse matrix (~hundreds of MB) on every
# region. Under GNU parallel with 24 workers that allocation rate outruns R's
# GC and the cgroup OOM-kills workers (and that is what was causing the
# OUT_OF_MEMORY SLURM failures while partial output was still produced).
tC       <- t(C)
CX       <- nullmod$CX                   # n x k
CXCXI    <- nullmod$CXCXI                # n x k
RSS0     <- as.numeric(nullmod$RSS0)
resid.PY <- as.numeric(nullmod$fit$resid.PY)
nm_ids   <- as.character(nullmod$fit$sample.id)
N        <- length(nm_ids)
family   <- if (!is.null(nullmod$model$family$family)) nullmod$model$family$family else "gaussian"
hetResid <- isTRUE(nullmod$model$hetResid)

outcome_label <- if (is.na(opt$outcome) && !is.null(nullmod$model$outcome))
                   nullmod$model$outcome else opt$outcome
if (is.na(outcome_label)) outcome_label <- "NA"

write(sprintf("[nullmod] file=%s  N=%d  k=%d  family=%s  hetResid=%s  outcome=%s",
              opt$nullModelFile, N, ncol(CX), family, hetResid, outcome_label),
      stderr())

# Index lookup: nullmod-sample-id -> position in nullmod vectors.
nm_pos <- setNames(seq_len(N), nm_ids)

# -------------------------------------------------------------------------
# 2. Open the coverage stream.  Header == 4 region columns then sample IDs.
# -------------------------------------------------------------------------
open_input <- function(path) {
  # NOTE: it is CRUCIAL to pass open="r" (or "rt") here. With the default
  # open="", R opens the connection inside each readLines() call and then
  # closes it again. For unseekable streams like stdin that close() drops
  # R's read-ahead buffer, so the FIRST body line after the header was
  # being silently skipped (one row lost per Rscript invocation, i.e. one
  # per parallel chunk in the production pipeline). Keeping the connection
  # explicitly open avoids that data loss.
  if (path == "-") return(file("stdin", open = "r"))
  if (grepl("\\.gz$", path, ignore.case = TRUE)) return(gzfile(path, open = "rt"))
  file(path, open = "r")
}
con <- open_input(opt$inputFile)
on.exit(try(close(con), silent = TRUE), add = TRUE)

header_line <- readLines(con, n = 1)
if (!length(header_line)) stop("empty input")
header <- strsplit(header_line, "\t", fixed = TRUE)[[1]]
if (length(header) < 5) stop("input header must have >=4 region cols + sample cols")

# Sample-ID rewriting is opt-in via --sampleIdPattern (a PCRE whose first
# capture group is kept). The previous hardcoded rules were cohort-specific and
# one of them silently mangled any identifier containing a colon.
cov_sample_ids_raw <- header[-(1:4)]
.id_pattern <- Sys.getenv("DSV_SAMPLE_ID_PATTERN", unset = "")
if (nzchar(.id_pattern)) {
  cov_sample_ids <- sub(.id_pattern, "\\1", cov_sample_ids_raw, perl = TRUE)
} else {
  cov_sample_ids <- cov_sample_ids_raw
}
M <- length(cov_sample_ids)

# Build mapping: for each coverage column j, which nullmod row (or NA) does it
# correspond to? This is the only sample-alignment work; done once.
nm_row_for_cov <- nm_pos[cov_sample_ids]   # length M, NA where the cov sample
                                           # is not in the nullmod sample set
covered <- !is.na(nm_row_for_cov)
n_overlap <- sum(covered)
write(sprintf("[align] coverage samples=%d  nullmod samples=%d  overlap=%d",
              M, N, n_overlap), stderr())
if (n_overlap < opt$minObs)
  stop("Insufficient sample overlap (", n_overlap, ") between coverage file and nullmod")

# Pre-compute the destination row indices used at every region.
nm_target <- as.integer(nm_row_for_cov[covered])
cov_src   <- which(covered)

# Pre-compute G' * resid.PY contribution for missing samples (mean-imputed G
# contributes mean(G_obs) * sum(resid.PY[missing])). Same trick for G' P G.
all_idx <- seq_len(N)
miss_idx_default <- setdiff(all_idx, nm_target)   # never observed in cov file

# -------------------------------------------------------------------------
# 3. Emit output header.
# -------------------------------------------------------------------------
out_header <- c("#CHROM", "START", "END", "Region", "OUTCOME",
                "n.obs", "n.imp", "Est", "Est.SE",
                "Score", "Score.SE", "Score.Stat", "Score.pval", "PVE")
writeLines(paste(out_header, collapse = "\t"), stdout())

# -------------------------------------------------------------------------
# 4. Stream regions and compute score statistic.
# -------------------------------------------------------------------------
chunk_size <- 10000L
nProc <- 0L; nSkipped <- 0L; nErr <- 0L
expected_cols <- M + 4L

repeat {
  lines <- readLines(con, n = chunk_size)
  if (!length(lines)) break

  for (line in lines) {
    parts <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (length(parts) < expected_cols) { nSkipped <- nSkipped + 1L; next }

    chr <- parts[1]; start <- parts[2]; end <- parts[3]; region <- parts[4]
    cov_vals <- suppressWarnings(as.numeric(parts[5:expected_cols]))

    # Build G aligned to nullmod rows (length N).
    G <- rep(NA_real_, N)
    G[nm_target] <- cov_vals[cov_src]

    # Mean-impute missing (manuscript: "missing SV genotype calls were imputed
    # to the mean before performing the association tests").
    obs <- !is.na(G)
    n.obs <- sum(obs)
    if (n.obs < opt$minObs) { nSkipped <- nSkipped + 1L; next }
    Gmean <- mean(G[obs])
    n.imp <- N - n.obs
    if (n.imp) G[!obs] <- Gmean

    # Drop constant G (e.g. a region with the same value everywhere -> GPG=0).
    if (n.obs > 1L && var(G[obs]) == 0) { nSkipped <- nSkipped + 1L; next }

    fit_res <- tryCatch({
      # Use pre-transposed tC; G is a plain numeric vector — no Matrix() wrap.
      CG     <- tC %*% G                                  # n x 1
      Gtilde <- CG - tcrossprod(CXCXI, crossprod(CG, CX)) # n x 1
      GPG    <- sum(as.numeric(Gtilde)^2)
      if (!is.finite(GPG) || GPG <= 0)
        stop("non-positive G'PG (", GPG, ")")
      score   <- sum(G * resid.PY)
      score.SE <- sqrt(GPG)
      Stat    <- score / score.SE
      pval    <- pchisq(Stat^2, df = 1, lower.tail = FALSE)
      list(Est = score / GPG, Est.SE = 1 / score.SE,
           Score = score, Score.SE = score.SE,
           Score.Stat = Stat, Score.pval = pval,
           PVE = (Stat^2) / RSS0)
    }, error = function(e) {
      nErr <<- nErr + 1L
      if (nErr <= 10L) write(sprintf("[err] region=%s : %s", region, e$message), stderr())
      NULL
    })
    if (is.null(fit_res)) { nSkipped <- nSkipped + 1L; next }

    out <- c(chr, start, end, region, outcome_label,
             as.character(n.obs), as.character(n.imp),
             formatC(fit_res$Est,        digits = 10, format = "g"),
             formatC(fit_res$Est.SE,     digits = 10, format = "g"),
             formatC(fit_res$Score,      digits = 10, format = "g"),
             formatC(fit_res$Score.SE,   digits = 10, format = "g"),
             formatC(fit_res$Score.Stat, digits = 10, format = "g"),
             formatC(fit_res$Score.pval, digits = 10, format = "g"),
             formatC(fit_res$PVE,        digits = 10, format = "g"))
    writeLines(paste(out, collapse = "\t"), stdout())
    nProc <- nProc + 1L
    if (nProc %% 1000L == 0L)
      write(sprintf("[progress] processed=%d skipped=%d err=%d", nProc, nSkipped, nErr), stderr())
  }
}

write(sprintf("[done] processed=%d skipped=%d err=%d", nProc, nSkipped, nErr), stderr())


