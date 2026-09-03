#!/usr/bin/env Rscript
# Replicates R/correct.R's per-process prologue step by step and times each
# piece: package load, PC table read, coverage read, merge, header match, qr,
# qr.Q. Args: <pcs> <coverage> <ndim> <headerFile>
t_all <- proc.time()
tm <- function(label, expr) {
  t0 <- proc.time(); v <- expr; t1 <- proc.time()
  cat(sprintf("%-28s %7.3f s wall %7.3f s cpu\n", label, (t1 - t0)[["elapsed"]], (t1 - t0)[["user.self"]] + (t1 - t0)[["sys.self"]]))
  invisible(v)
}
tm("library(optparse,data.table)", suppressPackageStartupMessages({library(optparse); library(data.table)}))
args <- commandArgs(trailingOnly = TRUE)
pcs_f <- args[1]; cov_f <- args[2]; ndim <- as.integer(args[3]); hdr_f <- args[4]
cov_dt <- tm("fread coverage", fread(cov_f, select = c("SAMPLE", "AUTO_HQ_median")))
cov_dt[, SAMPLE := as.character(SAMPLE)]
pcs_dt <- tm("fread PCs (all columns)", fread(pcs_f))
pcs_dt[, SAMPLE := as.character(SAMPLE)]
pc_cols <- paste0("PC", seq_len(ndim))
ref <- tm("merge", merge(pcs_dt[, c("SAMPLE", pc_cols), with = FALSE], cov_dt, by = "SAMPLE"))
header <- tm("read+split header", strsplit(readLines(hdr_f, n = 1L), "\t", fixed = TRUE)[[1]])
hs <- header[-(1:3)]
idx <- tm("match samples", match(hs, ref$SAMPLE))
keep <- !is.na(idx)
ref_aligned <- ref[idx[keep]]
X <- tm("build X", cbind(Intercept = 1, as.matrix(ref_aligned[, pc_cols, with = FALSE])))
qr_X <- tm("qr(X)", qr(X))
Q <- tm("qr.Q + truncate", qr.Q(qr_X)[, seq_len(qr_X$rank), drop = FALSE])
cat(sprintf("TOTAL prologue %.3f s wall, n=%d ndim=%d, maxRSS %.0f MB\n",
            (proc.time() - t_all)[["elapsed"]], nrow(X), ndim, as.numeric(system2("ps", c("-o", "rss=", "-p", Sys.getpid()), stdout = TRUE)) / 1024))
