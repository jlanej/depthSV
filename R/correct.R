#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# depthSV — depth normalisation and technical correction
#
# Converts raw per-region read depth into a log2 ratio against each sample's
# median autosomal coverage, then optionally residualises against the leading
# principal components of the depth matrix to remove technical structure.
#
# Input : CHR START STOP <sample columns...>   (the joined depth matrix)
# Output: #CHROM START END Region <sample columns...>
#
# The sample alignment and the projection basis Q are built once per region
# by the driver (--saveBasis) and loaded by every parallel worker
# (--loadBasis): at biobank width, re-reading the PC table and refactorising
# the design in every worker cost more than the correction itself.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

option_list <- list(
  make_option(c("-i", "--inputPCs"), type = "character", default = NULL,
              help = "PC table; must contain SAMPLE and PC1..PCn"),
  make_option(c("-f", "--inputFile"), type = "character", default = NULL,
              help = "joined depth matrix, or '-' for stdin"),
  make_option(c("-d", "--ndim"), type = "integer", default = 16,
              help = "number of PCs to residualise against; 0 = log2 ratio only [default %default]"),
  make_option(c("-c", "--coverageStats"), type = "character", default = NULL,
              help = "per-sample coverage; must contain SAMPLE and AUTO_HQ_median"),
  make_option("--saveBasis", type = "character", default = NULL,
              help = "read only the input header, write the alignment and projection basis here, and exit"),
  make_option("--loadBasis", type = "character", default = NULL,
              help = "reuse a basis written by --saveBasis instead of reading the PC and coverage tables"),
  make_option("--skipOutputHeader", action = "store_true", default = FALSE,
              help = "suppress the output header (for parallel chunks)"),
  make_option(c("-s", "--statsFile"), type = "character", default = NULL,
              help = "optional per-region pre/post-correction summary statistics"),
  make_option("--sampleIdPattern", type = "character", default = NULL,
              help = "optional PCRE applied to column names; first capture group is kept"),
  make_option("--minDepth", type = "double", default = 0.005,
              help = "floor applied to depth before log2 [default %default]"),
  make_option("--digits", type = "integer", default = 6,
              help = "significant digits for corrected values [default %default]"),
  make_option("--projection", type = "character", default = "explicit",
              help = "how to residualise: 'explicit' (default) or 'qr'")
)

opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$inputFile)) stop("--inputFile is required", call. = FALSE)
if (is.null(opt$loadBasis)) {
  for (req in c("inputPCs", "coverageStats")) {
    if (is.null(opt[[req]])) stop("--", req, " is required (or --loadBasis)", call. = FALSE)
  }
}
if (opt$ndim < 0) stop("--ndim must be >= 0", call. = FALSE)
if (!opt$projection %in% c("explicit", "qr")) stop("--projection must be 'explicit' or 'qr'", call. = FALSE)

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
                       "Duplicates silently assign one sample's values to another."),
                 what, length(dup), paste(utils::head(dup, 3), collapse = ", ")), call. = FALSE)
  }
}

# --- input stream and header ----------------------------------------------

open_input <- function(path) {
  if (path == "-") return(file("stdin", open = "r"))
  if (grepl("\\.gz$", path, ignore.case = TRUE)) return(gzfile(path, open = "rt"))
  file(path, open = "r")
}
con <- open_input(opt$inputFile)
on.exit(try(close(con), silent = TRUE), add = TRUE)

header <- strsplit(readLines(con, n = 1L), "\t", fixed = TRUE)[[1]]
if (length(header) < 4L) {
  stop("input header must have 3 region columns followed by at least one sample column", call. = FALSE)
}
header_samples <- normalize_sample_ids(header[-(1:3)], opt$sampleIdPattern)
n_in <- length(header_samples)
# match() is many-to-one: two columns mapping to one participant would both
# be kept under one name and every later stage would count them twice.
refuse_duplicates(header_samples, "the matrix header (after --sampleIdPattern)")

# --- alignment and projection basis ----------------------------------------

build_basis <- function() {
  cov_dt <- read_table(opt$coverageStats, select = c("SAMPLE", "AUTO_HQ_median"))
  for (col in c("SAMPLE", "AUTO_HQ_median")) {
    if (!col %in% names(cov_dt)) stop("coverage stats must contain a ", col, " column", call. = FALSE)
  }
  cov_dt <- cov_dt[, .(SAMPLE = as.character(SAMPLE), AUTO_HQ_median = as.numeric(AUTO_HQ_median))]
  refuse_duplicates(cov_dt$SAMPLE, "the coverage table")
  # A missing or zero median would put NA or -Inf into every value of that
  # sample and, through the projection, into every sample; refuse it here.
  bad <- cov_dt[!is.finite(AUTO_HQ_median) | AUTO_HQ_median <= 0]
  if (nrow(bad)) {
    stop(sprintf("coverage table has a missing or non-positive AUTO_HQ_median for %d sample(s), e.g. %s",
                 nrow(bad), paste(utils::head(bad$SAMPLE, 3), collapse = ", ")), call. = FALSE)
  }

  pc_cols <- if (opt$ndim > 0) paste0("PC", seq_len(opt$ndim)) else character(0)
  # Only the requested components are parsed: the upstream table may carry
  # hundreds of columns, and every worker used to pay for all of them.
  pcs_dt <- read_table(opt$inputPCs, select = c("SAMPLE", pc_cols))
  if (!"SAMPLE" %in% names(pcs_dt)) stop("PC table must contain a SAMPLE column", call. = FALSE)
  missing_pc <- setdiff(pc_cols, names(pcs_dt))
  if (length(missing_pc)) {
    stop(sprintf("PC table has %d of the %d requested components; missing e.g. %s",
                 opt$ndim - length(missing_pc), opt$ndim,
                 paste(utils::head(missing_pc, 3), collapse = ", ")), call. = FALSE)
  }
  pcs_dt[, SAMPLE := as.character(SAMPLE)]
  refuse_duplicates(pcs_dt$SAMPLE, "the PC table")

  ref <- merge(pcs_dt[, c("SAMPLE", pc_cols), with = FALSE], cov_dt, by = "SAMPLE")
  if (!nrow(ref)) stop("no samples common to the PC table and the coverage stats", call. = FALSE)

  idx <- match(header_samples, ref$SAMPLE)
  keep <- !is.na(idx)
  if (!any(keep)) {
    stop("no matrix sample matched the PC/coverage tables; check --sampleIdPattern", call. = FALSE)
  }
  if (any(!keep)) {
    message(sprintf("[align] %d of %d samples lack PCs or coverage and are dropped", sum(!keep), n_in))
  }
  ref_aligned <- ref[idx[keep]]

  Q <- NULL; rank <- 0L
  if (opt$ndim > 0) {
    X <- cbind(Intercept = 1, as.matrix(ref_aligned[, pc_cols, with = FALSE]))
    qr_X <- qr(X)
    rank <- qr_X$rank
    if (rank < ncol(X)) {
      message(sprintf("[warn] PC design is rank deficient (%d of %d); correction proceeds on the reduced basis",
                      rank, ncol(X)))
    }
    if (sum(keep) <= rank) {
      stop(sprintf("cannot residualise %d samples against a rank-%d design", sum(keep), rank), call. = FALSE)
    }
    # Q truncated to the numerical rank: qr.Q() returns every column, and on a
    # rank-deficient design the extra columns span numerically-null directions.
    Q <- qr.Q(qr_X)[, seq_len(rank), drop = FALSE]
    if (opt$projection == "qr") attr(Q, "qr") <- qr_X
  }
  list(header = header_samples, keep = keep, kept_samples = header_samples[keep],
       median_vec = ref_aligned$AUTO_HQ_median, Q = Q, rank = rank, ndim = opt$ndim,
       projection = opt$projection)
}

if (!is.null(opt$loadBasis)) {
  basis <- readRDS(opt$loadBasis)
  if (!identical(basis$header, header_samples)) {
    stop("the loaded basis was built for a different matrix header", call. = FALSE)
  }
  if (basis$ndim != opt$ndim) stop("the loaded basis was built for ndim ", basis$ndim, call. = FALSE)
} else {
  basis <- build_basis()
  if (!is.null(opt$saveBasis)) {
    saveRDS(basis, opt$saveBasis, compress = FALSE)
    message(sprintf("[basis] samples=%d ndim=%d rank=%d -> %s",
                    length(basis$kept_samples), basis$ndim, basis$rank, opt$saveBasis))
    quit(status = 0)
  }
}

keep         <- basis$keep
kept_samples <- basis$kept_samples
median_vec   <- basis$median_vec
n_keep       <- length(kept_samples)
use_resid    <- basis$ndim > 0
Q            <- basis$Q

# Two ways to compute the same projection, residual = (I - P)v: 'explicit'
# forms P = QQ' and computes v - Q(Q'v) (two BLAS calls, 12-20x faster);
# 'qr' applies qr.resid(). They agree to ~6e-14 relative.
residualize <- NULL
if (use_resid) {
  residualize <- if (basis$projection == "explicit") {
    function(v) as.numeric(v - Q %*% crossprod(Q, v))
  } else {
    function(v) as.numeric(qr.resid(attr(Q, "qr"), v))
  }
}

out_header <- c("#CHROM", "START", "END", "Region", kept_samples)
if (!opt$skipOutputHeader) cat(paste(out_header, collapse = "\t"), "\n", sep = "")

stats_con <- NULL
if (!is.null(opt$statsFile)) {
  stats_con <- file(opt$statsFile, "w")
  on.exit(try(close(stats_con), silent = TRUE), add = TRUE)
  writeLines(paste(c("#CHROM", "START", "END", "Region",
                     paste0("pre_",  c("mean","median","sd","mad","min","max")),
                     paste0("post_", c("mean","median","sd","mad","min","max"))),
                   collapse = "\t"), stats_con)
}

# --- stream ----------------------------------------------------------------

summarise <- function(v) {
  m <- median(v)
  c(mean(v), m, sd(v), mad(v, center = m), min(v), max(v))
}

fmt <- paste0("%.", opt$digits, "g")   # loop-invariant
n_proc <- 0L; n_bad <- 0L
while (length(line <- readLines(con, n = 1L)) > 0L) {
  fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
  if (length(fields) != n_in + 3L) { n_bad <- n_bad + 1L; next }

  chrom <- fields[1]; start <- fields[2]; end <- fields[3]
  depth <- suppressWarnings(as.numeric(fields[-(1:3)]))[keep]
  if (anyNA(depth)) { n_bad <- n_bad + 1L; next }

  log2_ratio <- log2(pmax(opt$minDepth, depth) / median_vec)
  vals <- if (use_resid) residualize(log2_ratio) else log2_ratio
  region_id <- sprintf("%s:%s-%s", chrom, start, end)

  cat(paste(c(chrom, start, end, region_id, sprintf(fmt, vals)), collapse = "\t"), "\n", sep = "")

  if (!is.null(stats_con)) {
    pre  <- summarise(log2_ratio)
    post <- if (use_resid) summarise(vals) else pre
    writeLines(paste(c(chrom, start, end, region_id, sprintf("%.6f", c(pre, post))), collapse = "\t"),
               stats_con)
  }
  n_proc <- n_proc + 1L
}

if (n_bad > 0L) message(sprintf("[warn] %d malformed rows skipped", n_bad))
message(sprintf("[done] processed=%d samples=%d ndim=%d", n_proc, n_keep, basis$ndim))
