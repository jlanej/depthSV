#!/usr/bin/env Rscript
# Synthetic reference tables and a few wide matrix rows at biobank widths.
# For each N: svd.pcs_N.txt (N x 40), autosomal.median_N.txt, phenotypes_N.tsv,
# raw_N.txt (matrix header + 40 raw-depth rows, mosdepth-style 2 decimals),
# corr_N.txt (corrected header + 40 rows at %.6g).
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
out <- args[1]
Ns  <- as.integer(strsplit(args[2], ",")[[1]])
nrow_out <- if (length(args) >= 3) as.integer(args[3]) else 40L
set.seed(1)
for (N in Ns) {
  t0 <- proc.time()[["elapsed"]]
  ids <- sprintf("S%07d", seq_len(N))
  pcs <- matrix(rnorm(N * 40, sd = 0.02), N, 40)
  colnames(pcs) <- paste0("PC", 1:40)
  fwrite(data.table(SAMPLE = ids, as.data.table(pcs)), sprintf("%s/svd.pcs_%d.txt", out, N), sep = "\t")
  med <- round(pmax(5, rnorm(N, 30, 4)), 4)
  fwrite(data.table(SAMPLE = ids, AUTO_HQ_median = med), sprintf("%s/autosomal.median_%d.txt", out, N), sep = "\t")
  ph <- data.table(SAMPLE = ids, y = round(rnorm(N), 4), ybin = rbinom(N, 1, 0.5),
                   ybin_rare = rbinom(N, 1, 0.02),
                   age = round(rnorm(N, 55, 9), 1), sex = sample(0:1, N, TRUE))
  for (k in 1:10) ph[[paste0("gpc", k)]] <- round(rnorm(N, sd = 0.01), 6)
  fwrite(ph, sprintf("%s/phenotypes_%d.tsv", out, N), sep = "\t")
  # wide rows: log2 ratio with PC structure + noise, then raw depth from the medians
  load <- matrix(rnorm(nrow_out * 20, sd = 3), nrow_out, 20)
  l2   <- load %*% t(pcs[, 1:20]) + matrix(rnorm(nrow_out * N, sd = 0.18), nrow_out, N)
  raw  <- sweep(2^l2, 2, med, `*`)
  raw_txt <- matrix(sprintf("%.2f", raw), nrow_out, N)
  starts <- seq(0L, by = 1000L, length.out = nrow_out)
  con <- file(sprintf("%s/raw_%d.txt", out, N), "w")
  writeLines(paste(c("#CHR", "START", "STOP", ids), collapse = "\t"), con)
  for (i in seq_len(nrow_out))
    writeLines(paste(c("chr1", starts[i], starts[i] + 1000L, raw_txt[i, ]), collapse = "\t"), con)
  close(con)
  corr_txt <- matrix(sprintf("%.6g", l2 - rowMeans(l2)), nrow_out, N)
  con <- file(sprintf("%s/corr_%d.txt", out, N), "w")
  writeLines(paste(c("#CHROM", "START", "END", "Region", ids), collapse = "\t"), con)
  for (i in seq_len(nrow_out))
    writeLines(paste(c("chr1", starts[i], starts[i] + 1000L,
                       sprintf("chr1:%d-%d", starts[i], starts[i] + 1000L), corr_txt[i, ]), collapse = "\t"), con)
  close(con)
  cat(sprintf("N=%d done in %.1fs\n", N, proc.time()[["elapsed"]] - t0))
}
