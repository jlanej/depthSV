#!/usr/bin/env Rscript
# Per-row cost of R/correct.R's streaming loop at width N, broken into parse,
# log2, residualise (explicit Q), format, write. Also a block (matrix) variant.
# Args: <raw_N.txt> <ndim>
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
f <- args[1]; ndim <- as.integer(args[2])
con <- file(f, "r"); header <- strsplit(readLines(con, 1), "\t", fixed = TRUE)[[1]]
lines <- readLines(con); close(con)
N <- length(header) - 3L; R <- length(lines)
set.seed(1)
X <- cbind(1, matrix(rnorm(N * ndim, sd = 0.02), N, ndim))
Q <- qr.Q(qr(X))
med <- pmax(5, rnorm(N, 30, 4))
fmt <- "%.6g"
sink_con <- file("/dev/null", "w")
tm <- function(expr) { t0 <- proc.time()[["elapsed"]]; v <- expr; c(v = list(v), t = proc.time()[["elapsed"]] - t0) }
tp <- ta <- tl <- tr <- tf <- tw <- 0
for (line in lines) {
  r <- tm(strsplit(line, "\t", fixed = TRUE)[[1]]); tp <- tp + r$t; fields <- r$v
  r <- tm(suppressWarnings(as.numeric(fields[-(1:3)]))); ta <- ta + r$t; depth <- r$v
  r <- tm(log2(pmax(0.005, depth) / med)); tl <- tl + r$t; l2 <- r$v
  r <- tm(as.numeric(l2 - Q %*% crossprod(Q, l2))); tr <- tr + r$t; vals <- r$v
  r <- tm(sprintf(fmt, vals)); tf <- tf + r$t; s <- r$v
  r <- tm(cat(paste(c(fields[1:3], "id", s), collapse = "\t"), "\n", sep = "", file = sink_con)); tw <- tw + r$t
}
cat(sprintf("N=%d rows=%d ndim=%d  per-row ms: split %.1f  as.numeric %.1f  log2 %.1f  resid %.1f  sprintf %.1f  paste+cat %.1f  TOTAL %.1f\n",
            N, R, ndim, 1e3 * tp / R, 1e3 * ta / R, 1e3 * tl / R, 1e3 * tr / R, 1e3 * tf / R, 1e3 * tw / R,
            1e3 * (tp + ta + tl + tr + tf + tw) / R))
# block variant: parse all rows with fread, one GEMM, one sprintf, one writeLines
t0 <- proc.time()[["elapsed"]]
blk <- fread(text = lines, header = FALSE, sep = "\t", showProgress = FALSE)
t1 <- proc.time()[["elapsed"]]
D <- as.matrix(blk[, -(1:3), with = FALSE])
L2 <- log2(sweep(pmax(0.005, D), 2, med, `/`))
Vals <- L2 - (L2 %*% Q) %*% t(Q)
t2 <- proc.time()[["elapsed"]]
S <- matrix(sprintf(fmt, Vals), nrow(Vals))
out <- paste(blk[[1]], blk[[2]], blk[[3]], "id", apply(S, 1, paste, collapse = "\t"), sep = "\t")
writeLines(out, sink_con)
t3 <- proc.time()[["elapsed"]]
cat(sprintf("block of %d rows: fread %.1f ms/row  math %.1f ms/row  format+write %.1f ms/row  TOTAL %.1f ms/row\n",
            R, 1e3 * (t1 - t0) / R, 1e3 * (t2 - t1) / R, 1e3 * (t3 - t2) / R, 1e3 * (t3 - t0) / R))
