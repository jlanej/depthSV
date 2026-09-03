#!/usr/bin/env Rscript
# analyze.R linear path: per-row loop (as written) vs one block matrix op.
# Also the parse path: fread(text=, colClasses="character") + as.numeric(unlist())
# versus a numeric fread. Args: <n> [<m bins>=2000] [<ncov>=12]
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
n <- as.integer(args[1]); m <- if (length(args) >= 2) as.integer(args[2]) else 2000L
p <- if (length(args) >= 3) as.integer(args[3]) else 12L
set.seed(2)
Z <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))
y <- rnorm(n)
qr_Z <- qr(Z); Q <- qr.Q(qr_Z)[, seq_len(qr_Z$rank), drop = FALSE]
df_res <- n - qr_Z$rank - 1L
residualize <- function(v) as.numeric(v - Q %*% crossprod(Q, v))
y_res <- residualize(y); yy_res <- sum(y_res^2)
G <- matrix(rnorm(m * n, sd = 0.2), m, n)              # bins x samples, as depth[i, base_rows]
base_rows <- seq_len(n)
minObs <- 100L; minVar <- 1e-12

fit_linear <- function(g) {
  x_res <- residualize(g); xx <- sum(x_res^2)
  if (xx <= 0) return(NULL)
  beta <- sum(x_res * y_res) / xx
  rss  <- yy_res - beta^2 * xx
  se   <- sqrt((rss / df_res) / xx)
  tval <- beta / se
  c(beta, se, tval, 2 * pt(-abs(tval), df_res))
}
# --- loop, as in analyze.R lines 303-337 (minus the cat) ----------------------
t0 <- proc.time()[["elapsed"]]
res_loop <- matrix(NA_real_, m, 4)
for (i in seq_len(m)) {
  g_all <- G[i, base_rows]; keep <- !is.na(g_all); n_obs <- sum(keep)
  if (n_obs < minObs) next
  if (stats::var(g_all[keep]) <= minVar) next
  g <- g_all
  if (n_obs < length(g_all)) g[!keep] <- mean(g_all[keep])
  vals <- tryCatch(fit_linear(g), error = function(e) NULL)
  if (is.null(vals) || anyNA(vals)) next
  res_loop[i, ] <- vals
}
t_loop <- proc.time()[["elapsed"]] - t0
# --- block ---------------------------------------------------------------------
t0 <- proc.time()[["elapsed"]]
Gr   <- G - (G %*% Q) %*% t(Q)                          # one GEMM pair
xx   <- rowSums(Gr * Gr)
beta <- as.numeric(Gr %*% y_res) / xx
rss  <- yy_res - beta^2 * xx
se   <- sqrt((rss / df_res) / xx)
tval <- beta / se
pv   <- 2 * pt(-abs(tval), df_res)
vr   <- rowSums((G - rowMeans(G))^2) / (n - 1)          # the var() QC, vectorised
res_blk <- cbind(beta, se, tval, pv)
t_blk <- proc.time()[["elapsed"]] - t0
cat(sprintf("n=%d m=%d p=%d  loop %.2f s (%.2f ms/bin)  block %.2f s (%.3f ms/bin)  speedup %.1fx  max|diff| beta %.1e t %.1e\n",
            n, m, p, t_loop, 1e3 * t_loop / m, t_blk, 1e3 * t_blk / m, t_loop / t_blk,
            max(abs(res_loop[, 1] - res_blk[, 1])), max(abs(res_loop[, 3] - res_blk[, 3]))))
# --- parse path ----------------------------------------------------------------
lines <- paste("chr1", seq_len(m) * 1000L, seq_len(m) * 1000L + 1000L, "chr1:x-y",
               apply(matrix(sprintf("%.6g", G), m, n), 1, paste, collapse = "\t"), sep = "\t")
t0 <- proc.time()[["elapsed"]]
block <- fread(text = lines, header = FALSE, sep = "\t", showProgress = FALSE, colClasses = "character")
depth <- suppressWarnings(matrix(as.numeric(unlist(block[, 5:ncol(block), with = FALSE], use.names = FALSE)), nrow = nrow(block)))
t_parse_chr <- proc.time()[["elapsed"]] - t0
t0 <- proc.time()[["elapsed"]]
block2 <- fread(text = lines, header = FALSE, sep = "\t", showProgress = FALSE,
                colClasses = list(character = 1:4))
depth2 <- as.matrix(block2[, 5:ncol(block2), with = FALSE])
t_parse_num <- proc.time()[["elapsed"]] - t0
cat(sprintf("parse %d x %d (%.0f MB text): character+as.numeric %.2f s   numeric fread %.2f s   (%.1fx)\n",
            m, n, sum(nchar(lines)) / 1e6, t_parse_chr, t_parse_num, t_parse_chr / t_parse_num))
