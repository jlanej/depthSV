#!/usr/bin/env Rscript
# Prototype of the remedy: one chunk processed as a block. (a) correct: numeric
# fread of the raw chunk, one GEMM against a cached Q, signif(6)+fwrite.
# (b) analyze-linear: numeric fread of a corrected chunk, one GEMM, stats.
# Compares against the measured per-process cost of the shipped scripts.
# Args: <rawChunk> <corrChunk> <pcs> <coverage> <ndim> <qCache.rds>
suppressPackageStartupMessages(library(data.table))
setDTthreads(1)
args <- commandArgs(trailingOnly = TRUE)
raw_f <- args[1]; corr_f <- args[2]; pcs_f <- args[3]; cov_f <- args[4]; ndim <- as.integer(args[5]); qf <- args[6]
tm <- function(expr) { t0 <- proc.time()[["elapsed"]]; force(expr); proc.time()[["elapsed"]] - t0 }
# --- prologue: build Q once and cache it as an uncompressed .rds ---------------
if (!file.exists(qf)) {
  t_build <- tm({
    pcs <- fread(pcs_f, select = c("SAMPLE", paste0("PC", seq_len(ndim)))); cov <- fread(cov_f)
    ref <- merge(pcs, cov, by = "SAMPLE")
    X <- cbind(1, as.matrix(ref[, paste0("PC", seq_len(ndim)), with = FALSE]))
    q <- qr(X); Q <- qr.Q(q)[, seq_len(q$rank), drop = FALSE]
    saveRDS(list(samples = ref$SAMPLE, Q = Q, med = pmax(0.005, ref$AUTO_HQ_median)), qf, compress = FALSE)
  })
  cat(sprintf("build+cache Q (%d x %d): %.2f s, %.0f MB on disk\n", nrow(Q), ncol(Q), t_build, file.size(qf) / 1e6))
}
t_load <- tm(cache <- readRDS(qf))
Q <- cache$Q; med <- cache$med
cat(sprintf("load cached Q: %.2f s\n", t_load))
# --- (a) correct one chunk as a block ------------------------------------------
t_read <- tm({
  blk <- fread(raw_f, header = FALSE, skip = 1L, sep = "\t", colClasses = list(character = 1:3), showProgress = FALSE)
  D <- as.matrix(blk[, -(1:3), with = FALSE]); storage.mode(D) <- "double"
})
t_math <- tm({
  L2 <- log2(pmax(0.005, D) / rep(med, each = nrow(D)))
  V  <- L2 - (L2 %*% Q) %*% t(Q)
})
t_write <- tm({
  out <- data.table(blk[[1]], blk[[2]], blk[[3]], paste0(blk[[1]], ":", blk[[2]], "-", blk[[3]]))
  out <- cbind(out, as.data.table(signif(V, 6)))
  fwrite(out, "proto_corr_out.txt", sep = "\t", col.names = FALSE)
})
sz <- file.size("proto_corr_out.txt")
cat(sprintf("CORRECT block %d rows x %d: fread %.2f s  math %.2f s  signif+fwrite %.2f s  TOTAL %.2f s (%.0f ms/row); output %.1f B/value\n",
            nrow(D), ncol(D), t_read, t_math, t_write, t_read + t_math + t_write, 1e3 * (t_read + t_math + t_write) / nrow(D), sz / (nrow(D) * ncol(D))))
# --- (b) analyze-linear one chunk as a block ------------------------------------
set.seed(4); n <- ncol(D)
Z <- cbind(1, matrix(rnorm(n * 12), n, 12)); y <- rnorm(n)
qz <- qr(Z); Qz <- qr.Q(qz)[, seq_len(qz$rank), drop = FALSE]; df_res <- n - qz$rank - 1L
y_res <- as.numeric(y - Qz %*% crossprod(Qz, y)); yy <- sum(y_res^2)
t_read2 <- tm({
  cb <- fread(corr_f, header = FALSE, skip = 1L, sep = "\t", colClasses = list(character = 1:4), showProgress = FALSE)
  G <- as.matrix(cb[, -(1:4), with = FALSE]); storage.mode(G) <- "double"
})
t_fit <- tm({
  Gr <- G - (G %*% Qz) %*% t(Qz); xx <- rowSums(Gr * Gr)
  beta <- as.numeric(Gr %*% y_res) / xx; rss <- yy - beta^2 * xx; se <- sqrt((rss / df_res) / xx); tv <- beta / se
  res <- data.table(cb[[1]], cb[[2]], cb[[3]], cb[[4]], n, beta, se, tv, 2 * pt(-abs(tv), df_res))
})
cat(sprintf("ANALYZE-linear block %d rows x %d: fread %.2f s  fit %.3f s  TOTAL %.2f s (%.0f ms/row)\n",
            nrow(G), ncol(G), t_read2, t_fit, t_read2 + t_fit, 1e3 * (t_read2 + t_fit) / nrow(G)))
# --- (c) the same from a float32 binary chunk -------------------------------------
con <- file("proto_chunk.f32", "wb"); writeBin(as.numeric(t(G)), con, size = 4); close(con)
t_bin <- tm({ con <- file("proto_chunk.f32", "rb"); Gb <- matrix(readBin(con, "numeric", n = length(G), size = 4), nrow(G), byrow = TRUE); close(con) })
cat(sprintf("float32 chunk read (%d x %d, %.0f MB): %.2f s (%.0f ms/row)\n", nrow(G), ncol(G), 4 * length(G) / 1e6, t_bin, 1e3 * t_bin / nrow(G)))
invisible(file.remove("proto_corr_out.txt", "proto_chunk.f32"))
