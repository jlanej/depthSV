#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
S <- "/private/tmp/claude-501/-Users-Kitty-git-depthSV/1db72424-c4a1-45f4-b81e-5b466f3ec8a4/scratchpad/review_example/spectrum"
q <- fread(file.path(S, "std.sample_qc.tsv"))
p <- fread(file.path(S, "svd.pcs.txt"))
p[, SAMPLE := sub("[.]by1000[.]$", "", SAMPLE)]
cat("PC-table IDs absent from the QC table:", paste(setdiff(p$SAMPLE, q$SAMPLE_ID), collapse = ", "), "\n")
cat("QC IDs absent from the PC table:", paste(setdiff(q$SAMPLE_ID, p$SAMPLE), collapse = ", "), "\n")
m <- merge(q[, .(SAMPLE = SAMPLE_ID, MTDNA_CN)], p, by = "SAMPLE", sort = FALSE)
m <- m[is.finite(MTDNA_CN) & MTDNA_CN > 0]
y <- log2(m$MTDNA_CN)
r2 <- function(k, rows) { X <- as.matrix(m[rows, paste0("PC", seq_len(k)), with = FALSE]); summary(lm(y[rows] ~ X))$r.squared }
sub64 <- match(q$SAMPLE_ID[1:64], m$SAMPLE); sub64 <- sub64[!is.na(sub64)]
cat(sprintf("smoke subset (first %d QC rows): R2 of log2(MTDNA_CN) on PC1..4 = %.2f, PC1..20 = %.2f, PC1..42 = %.2f\n",
            length(sub64), r2(4, sub64), r2(20, sub64), r2(42, sub64)))
all <- seq_len(nrow(m))
cat(sprintf("all %d samples:                 R2 on PC1..4 = %.2f, PC1..20 = %.2f, PC1..42 = %.2f, PC1..52 = %.2f, PC1..200 = %.2f\n",
            length(all), r2(4, all), r2(20, all), r2(42, all), r2(52, all), r2(200, all)))
cat("(README Notes: 'across a 64-sample subset the first 20 PCs explain ~67%'; config.sh: 'R^2 ~0.67 at ndim=20, ~0.29 at ndim=4' and 'at full scale ... dilutes over 3,202 samples')\n")
