#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# depthSV — depth normalisation and technical correction
#
# Converts raw per-region read depth into a log2 ratio against each sample's
# EXPECTED depth — its median autosomal coverage scaled by the copies the
# region should carry — then optionally residualises against the leading
# principal components of the depth matrix to remove technical structure.
#
# Input : CHR START STOP <sample columns...>   (the joined depth matrix)
# Output: #CHROM START END Region <sample columns...>
#
# Expected copies are 2 everywhere unless --sex is given: then X outside the
# pseudo-autosomal regions is 2 in females and 1 in males, and Y outside them
# 1 in males and NONE in females (their value is NA, so a chrY region is
# tested in males only). Without a sex table every sex-chromosome bin is a
# near-perfect sex indicator; the pipeline says so once per region.
#
# The log2 ratio is winsorised at --winsorLog2 (default -3, an eighth of the
# expected depth). A zero-depth bin otherwise sits at log2(floor/median), a
# value that encodes the sample's sequencing depth rather than its copy
# number and carries almost all of a region's leverage; the winsor keeps a
# homozygous deletion the most extreme value while bounding its weight.
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
  make_option("--sex", type = "character", default = NULL,
              help = "table with SAMPLE and a sex column (M/F, male/female, or 1=male/2=female) for the ploidy model"),
  make_option("--sexCol", type = "character", default = "SEX",
              help = "name of the sex column in --sex [default %default]"),
  make_option("--par", type = "character", default = NULL,
              help = "BED of pseudo-autosomal regions (treated as diploid in both sexes)"),
  make_option("--winsorLog2", type = "double", default = -3,
              help = "floor on the log2 ratio; -Inf disables [default %default]"),
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
              help = "floor applied to raw depth before the log (guards log(0) only) [default %default]"),
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

# Sex codes accepted: M/F, male/female (any case), PLINK 1/2. 0/1 is refused
# because either sex could be the 1.
decode_sex <- function(x) {
  s <- toupper(trimws(as.character(x)))
  out <- rep(NA_character_, length(s))
  out[s %in% c("M", "MALE", "1")]   <- "M"
  out[s %in% c("F", "FEMALE", "2")] <- "F"
  bad <- unique(s[is.na(out) & !s %in% c("", "NA")])
  if (length(bad)) {
    stop("unrecognised sex code(s): ", paste(utils::head(bad, 5), collapse = ", "),
         " (use M/F, male/female, or 1=male/2=female)", call. = FALSE)
  }
  out
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
  kept <- header_samples[keep]

  sex <- NULL
  # Exact indexing: with no --sex given, `opt$sex` partial-matches --sexCol.
  if (!is.null(opt[["sex"]])) {
    sx <- read_table(opt[["sex"]], select = c("SAMPLE", opt$sexCol))
    if (!opt$sexCol %in% names(sx)) stop("--sex table has no column ", opt$sexCol, call. = FALSE)
    sx[, SAMPLE := as.character(SAMPLE)]
    refuse_duplicates(sx$SAMPLE, "the sex table")
    sex <- decode_sex(sx[[opt$sexCol]])[match(kept, sx$SAMPLE)]
    message(sprintf("[sex] %d male, %d female, %d unknown of %d kept samples",
                    sum(sex == "M", na.rm = TRUE), sum(sex == "F", na.rm = TRUE), sum(is.na(sex)), length(kept)))
  }

  X <- NULL; Q <- NULL; rank <- 0L
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
  }
  list(header = header_samples, keep = keep, kept_samples = kept,
       median_vec = ref_aligned$AUTO_HQ_median, sex = sex, X = X, Q = Q, rank = rank,
       ndim = opt$ndim, projection = opt$projection)
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
    message(sprintf("[basis] samples=%d ndim=%d rank=%d sex=%s -> %s",
                    length(basis$kept_samples), basis$ndim, basis$rank,
                    if (is.null(basis$sex)) "none" else "yes", opt$saveBasis))
    quit(status = 0)
  }
}

keep         <- basis$keep
kept_samples <- basis$kept_samples
median_vec   <- basis$median_vec
sex          <- basis$sex
n_keep       <- length(kept_samples)
use_resid    <- basis$ndim > 0
Q            <- basis$Q
X            <- basis$X

# --- pseudo-autosomal regions and contig classes ---------------------------

par <- NULL
if (!is.null(opt$par)) {
  par <- fread(opt$par, header = FALSE, select = 1:3, col.names = c("chrom", "start", "end"))
  par[, chrom := toupper(sub("^chr", "", as.character(chrom)))]
}
contig_class <- function(chrom) {
  c0 <- toupper(sub("^chr", "", chrom))
  if (c0 == "X") "X" else if (c0 == "Y") "Y" else "A"
}
# Base subsetting on purpose: inside a data.table `[` the arguments would
# resolve to the PAR table's own start/end columns.
in_par <- function(chrom, s, e) {
  if (is.null(par)) return(FALSE)
  c0 <- toupper(sub("^chr", "", chrom))
  any(par$chrom == c0 & par$end > s & par$start < e)
}
# Expected copies per kept sample for a region; NA where there are none.
expected_copies <- function(cls, is_par) {
  if (cls == "A" || is_par) return(rep(2, n_keep))
  if (is.null(sex)) return(rep(2, n_keep))
  e <- rep(NA_real_, n_keep)
  if (cls == "X") { e[sex == "F"] <- 2; e[sex == "M"] <- 1 }
  else            { e[sex == "M"] <- 1 }          # Y: none in females (NA)
  e
}
warned_sex <- FALSE

# --- projection ------------------------------------------------------------
# Two ways to compute the same projection, residual = (I - P)v: 'explicit'
# forms P = QQ' and computes v - Q(Q'v) (two BLAS calls, 12-20x faster);
# 'qr' applies qr.resid(). They agree to ~6e-14 relative. A vector with
# missing values (no expected copies for some samples) is residualised on
# the complete subset with a basis refactorised for that subset; the mask
# is the same for every chrY region, so the refactorisation is cached.

qr_full <- if (use_resid && basis$projection == "qr") qr(X) else NULL
sub_cache <- list(mask = NULL, Q = NULL)
residualize <- function(v) {
  if (!use_resid) return(v)
  mask <- !is.na(v)
  if (all(mask)) {
    if (basis$projection == "explicit") return(as.numeric(v - Q %*% crossprod(Q, v)))
    return(as.numeric(qr.resid(qr_full, v)))
  }
  if (is.null(sub_cache$mask) || !identical(sub_cache$mask, mask)) {
    qr_s <- qr(X[mask, , drop = FALSE])
    if (sum(mask) <= qr_s$rank) return(rep(NA_real_, length(v)))
    sub_cache <<- list(mask = mask, Q = qr.Q(qr_s)[, seq_len(qr_s$rank), drop = FALSE])
  }
  out <- rep(NA_real_, length(v))
  vs <- v[mask]
  out[mask] <- as.numeric(vs - sub_cache$Q %*% crossprod(sub_cache$Q, vs))
  out
}

out_header <- c("#CHROM", "START", "END", "Region", kept_samples)
if (!opt$skipOutputHeader) cat(paste(out_header, collapse = "\t"), "\n", sep = "")

stats_con <- NULL
if (!is.null(opt$statsFile)) {
  stats_con <- file(opt$statsFile, "w")
  on.exit(try(close(stats_con), silent = TRUE), add = TRUE)
  writeLines(paste(c("#CHROM", "START", "END", "Region", "n",
                     paste0("pre_",  c("mean","median","sd","mad","min","max")),
                     paste0("post_", c("mean","median","sd","mad","min","max"))),
                   collapse = "\t"), stats_con)
}

# --- stream ----------------------------------------------------------------

summarise <- function(v) {
  v <- v[!is.na(v)]
  if (!length(v)) return(rep(NA_real_, 6))
  m <- median(v)
  c(mean(v), m, sd(v), mad(v, center = m), min(v), max(v))
}

fmt <- paste0("%.", opt$digits, "g")   # loop-invariant
fmt_vals <- function(v) { s <- sprintf(fmt, v); s[is.na(v)] <- "NA"; s }
n_proc <- 0L; n_bad <- 0L
while (length(line <- readLines(con, n = 1L)) > 0L) {
  fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
  if (length(fields) != n_in + 3L) { n_bad <- n_bad + 1L; next }

  chrom <- fields[1]; start <- fields[2]; end <- fields[3]
  depth <- suppressWarnings(as.numeric(fields[-(1:3)]))[keep]
  if (anyNA(depth)) { n_bad <- n_bad + 1L; next }

  cls <- contig_class(chrom)
  if (cls != "A" && is.null(sex) && !warned_sex) {
    message("[warn] sex-chromosome regions present but no --sex table: every chrX/chrY value is a sex indicator")
    warned_sex <- TRUE
  }
  expected <- expected_copies(cls, cls != "A" && in_par(chrom, as.numeric(start), as.numeric(end)))
  log2_ratio <- log2(pmax(opt$minDepth, depth) / (median_vec * expected / 2))
  log2_ratio <- pmax(log2_ratio, opt$winsorLog2)
  vals <- residualize(log2_ratio)
  region_id <- sprintf("%s:%s-%s", chrom, start, end)

  cat(paste(c(chrom, start, end, region_id, fmt_vals(vals)), collapse = "\t"), "\n", sep = "")

  if (!is.null(stats_con)) {
    pre  <- summarise(log2_ratio)
    post <- if (use_resid) summarise(vals) else pre
    writeLines(paste(c(chrom, start, end, region_id, sum(!is.na(vals)),
                       sprintf("%.6f", c(pre, post))), collapse = "\t"), stats_con)
  }
  n_proc <- n_proc + 1L
}

if (n_bad > 0L) message(sprintf("[warn] %d malformed rows skipped", n_bad))
message(sprintf("[done] processed=%d samples=%d ndim=%d", n_proc, n_keep, basis$ndim))
