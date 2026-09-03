#!/usr/bin/env Rscript
# A 3,202-sample mosdepth tree over nbins 1 kb bins on one contig, using the
# REAL 1000G sample IDs, medians and PCs, so the real PC table and phenotype
# table apply unchanged. Depth = median * 2^(PC structure + noise), 2 decimals
# like mosdepth. Written as <ID>.by1000.regions.bed.gz (bgzip).
#   Rscript gen_mosdepth.R <medianTable> <pcTable> <outDir> <nbins> [<contig>]
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
med_f <- args[1]; pc_f <- args[2]; out <- args[3]; nbins <- as.integer(args[4])
contig <- if (length(args) >= 5) args[5] else "chr1"
set.seed(7)
med <- fread(med_f); pcs <- fread(pc_f)
setkey(med, SAMPLE); setkey(pcs, SAMPLE)
ids <- intersect(med$SAMPLE, pcs$SAMPLE)
med <- med[ids]; pcs <- pcs[ids]
n <- length(ids)
P <- as.matrix(pcs[, paste0("PC", 1:20), with = FALSE])
P <- sweep(P, 2, apply(P, 2, sd), `/`)          # unit-variance factors
dir.create(out, showWarnings = FALSE, recursive = TRUE)
starts <- seq(0L, by = 1000L, length.out = nbins)
ends   <- starts + 1000L
load   <- matrix(rnorm(nbins * 20, sd = 0.06), nbins, 20)   # per-bin loadings
# one deletion carried by 20% at bin 500
carrier <- rbinom(n, 1, 0.2)
step <- 200L
for (j0 in seq(1L, n, by = step)) {
  jj <- j0:min(n, j0 + step - 1L)
  l2 <- load %*% t(P[jj, , drop = FALSE]) + matrix(rnorm(nbins * length(jj), sd = 0.15), nbins, length(jj))
  l2[500, ] <- l2[500, ] - carrier[jj]
  d  <- sweep(2^l2, 2, med$AUTO_HQ_median[jj], `*`)
  for (k in seq_along(jj)) {
    f <- file.path(out, sprintf("%s.by1000.regions.bed", ids[jj[k]]))
    fwrite(data.table(contig, starts, ends, sprintf("%.2f", d[, k])), f, sep = "\t", col.names = FALSE, quote = FALSE)
    system2("bgzip", c("-f", shQuote(f)))
  }
  cat(sprintf("%d/%d\n", max(jj), n))
}
writeLines(file.path(normalizePath(out), sprintf("%s.by1000.regions.bed.gz", ids)), file.path(out, "..", "manifest.txt"))
cat(sprintf("wrote %d samples x %d bins to %s\n", n, nbins, out))
