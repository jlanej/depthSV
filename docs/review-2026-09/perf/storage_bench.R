#!/usr/bin/env Rscript
# Text vs binary: footprint and parse/format throughput for one 500k-wide row
# set (40 rows). Writes float32 / int16 files for the shell to compress.
suppressPackageStartupMessages(library(data.table))
out <- "tables"
con <- file(file.path(out, "corr_500000.txt"), "r"); header <- readLines(con, 1); lines <- readLines(con); close(con)
N <- 500000L; R <- length(lines)
tm <- function(expr) { t0 <- proc.time()[["elapsed"]]; force(expr); proc.time()[["elapsed"]] - t0 }
# --- parse: per-row strsplit + as.numeric (correct.R) --------------------------
t_split <- tm(for (l in lines) suppressWarnings(as.numeric(strsplit(l, "\t", fixed = TRUE)[[1]][-(1:4)])))
# --- parse: fread numeric of the whole block (data.table, threads as default) ---
t_fread <- tm(blk <- fread(text = lines, header = FALSE, sep = "\t", colClasses = list(character = 1:4), showProgress = FALSE))
t_fread1 <- { setDTthreads(1); tm(blk <- fread(text = lines, header = FALSE, sep = "\t", colClasses = list(character = 1:4), showProgress = FALSE)) }
setDTthreads(0)
M <- as.matrix(blk[, -(1:4), with = FALSE]); storage.mode(M) <- "double"
bytes_txt <- sum(nchar(lines, type = "bytes"))
cat(sprintf("text block: %d rows x %d = %.0f MB; strsplit+as.numeric %.2f s (%.0f MB/s); fread numeric %.2f s (%.0f MB/s, %d threads) / %.2f s (1 thread, %.0f MB/s)\n",
            R, N, bytes_txt / 1e6, t_split, bytes_txt / 1e6 / t_split, t_fread, bytes_txt / 1e6 / t_fread, getDTthreads(), t_fread1, bytes_txt / 1e6 / t_fread1))
# --- format: sprintf %.6g (correct.R) vs fwrite vs writeBin -------------------
sink_con <- file("/dev/null", "wb")
t_sprintf <- tm(for (i in seq_len(R)) cat(paste(c("chr1", "0", "1000", "id", sprintf("%.6g", M[i, ])), collapse = "\t"), "\n", sep = "", file = sink_con))
t_fwrite  <- tm(fwrite(as.data.table(M), "/dev/null", sep = "\t", col.names = FALSE))
t_fwrite1 <- { setDTthreads(1); tm(fwrite(as.data.table(M), "/dev/null", sep = "\t", col.names = FALSE)) }
setDTthreads(0)
f32 <- file.path(out, "corr_500000.f32"); f16 <- file.path(out, "corr_500000.i16"); f64 <- file.path(out, "corr_500000.f64")
t_f32 <- tm({ con <- file(f32, "wb"); writeBin(as.numeric(t(M)), con, size = 4); close(con) })
t_i16 <- tm({ con <- file(f16, "wb"); writeBin(as.integer(round(pmax(-32, pmin(32, t(M))) * 1000)), con, size = 2); close(con) })
con <- file(f64, "wb"); writeBin(as.numeric(t(M)), con, size = 8); close(con)
cat(sprintf("format %d rows: sprintf+paste+cat %.2f s (%.0f ms/row)  fwrite %.2f s (%d thr) / %.2f s (1 thr)  writeBin f32 %.2f s  writeBin i16 %.2f s\n",
            R, t_sprintf, 1e3 * t_sprintf / R, t_fwrite, getDTthreads(), t_fwrite1, t_f32, t_i16))
# --- read back binary ------------------------------------------------------------
t_r32 <- tm({ con <- file(f32, "rb"); v <- readBin(con, "numeric", n = N * R, size = 4); close(con) })
t_r16 <- tm({ con <- file(f16, "rb"); v <- readBin(con, "integer", n = N * R, size = 2) / 1000; close(con) })
cat(sprintf("readBin %d x %d: float32 %.3f s (%.0f MB/s)   int16 %.3f s\n", R, N, t_r32, 4 * N * R / 1e6 / t_r32, t_r16))
# --- quantisation error of int16 (x1000) and of %.6g -----------------------------
q16 <- round(M * 1000) / 1000
cat(sprintf("int16 x1000 max abs err %.2e (values range %.2f..%.2f); %%.6g max rel err %.1e\n",
            max(abs(q16 - M)), min(M), max(M), max(abs(as.numeric(sprintf("%.6g", M)) - M) / pmax(abs(M), 1e-300))))
# --- raw depth rows too ---------------------------------------------------------------
con <- file(file.path(out, "raw_500000.txt"), "r"); hdr <- readLines(con, 1); rl <- readLines(con); close(con)
rb <- fread(text = rl, header = FALSE, sep = "\t", colClasses = list(character = 1:3), showProgress = FALSE)
RM <- as.matrix(rb[, -(1:3), with = FALSE]); storage.mode(RM) <- "double"
con <- file(file.path(out, "raw_500000.f32"), "wb"); writeBin(as.numeric(t(RM)), con, size = 4); close(con)
con <- file(file.path(out, "raw_500000.u16"), "wb"); writeBin(as.integer(pmin(65535, round(t(RM) * 100))), con, size = 2); close(con)
cat(sprintf("raw depth: %.0f MB text for %d x %d; u16 x100 max abs err %.3f\n", sum(nchar(rl, type = "bytes")) / 1e6, nrow(RM), N, max(abs(pmin(655.35, round(RM * 100) / 100) - RM))))
