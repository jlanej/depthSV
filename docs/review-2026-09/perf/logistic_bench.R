#!/usr/bin/env Rscript
# analyze.R logistic path: glm.fit per bin as in fit_logistic, at width n, with
# balanced and 2% prevalence responses; a block score test for comparison; and
# logistf (Firth) per bin where n is small enough. Args: <n> <bins> [<ncov>=12]
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
n <- as.integer(args[1]); m <- as.integer(args[2]); p <- if (length(args) >= 3) as.integer(args[3]) else 12L
set.seed(3)
Z <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))
X_fixed <- cbind(Z[, 1, drop = FALSE], cov_resids = 0, Z[, -1, drop = FALSE]); depth_col <- 2L
G <- matrix(rnorm(m * n, sd = 0.2), m, n)
fit_logistic <- function(g, y_bin) {
  X <- X_fixed; X[, depth_col] <- g
  f <- suppressWarnings(glm.fit(X, y_bin, family = binomial()))
  pr <- f$rank; piv <- f$qr$pivot[seq_len(pr)]; pos <- match(depth_col, piv)
  if (is.na(pos)) return(NULL)
  cov_unscaled <- chol2inv(f$qr$qr[seq_len(pr), seq_len(pr), drop = FALSE])
  se <- sqrt(cov_unscaled[pos, pos]); beta <- f$coefficients[depth_col]
  z <- beta / se; c(beta, se, z, 2 * pnorm(-abs(z)), f$iter)
}
for (prev in c(0.5, 0.02)) {
  y <- rbinom(n, 1, prev)
  t0 <- proc.time()[["elapsed"]]
  iters <- numeric(m)
  for (i in seq_len(m)) { v <- fit_logistic(G[i, ], y); iters[i] <- v[5] }
  t_glm <- proc.time()[["elapsed"]] - t0
  # block score test: null model once, then per block one GEMM
  t0 <- proc.time()[["elapsed"]]
  f0 <- glm.fit(Z, y, family = binomial()); mu <- f0$fitted.values; w <- mu * (1 - mu)
  r  <- y - mu
  ZtWZ_inv <- chol2inv(chol(crossprod(Z * w, Z)))
  t_null <- proc.time()[["elapsed"]] - t0
  t0 <- proc.time()[["elapsed"]]
  U   <- as.numeric(G %*% r)                         # score
  GW  <- G * rep(w, each = m)                        # G diag(w)
  GWZ <- GW %*% Z                                    # m x p
  V   <- rowSums(GW * G) - rowSums((GWZ %*% ZtWZ_inv) * GWZ)
  Tsc <- U^2 / V
  t_score <- proc.time()[["elapsed"]] - t0
  cat(sprintf("n=%d p=%d prev=%.2f bins=%d  glm.fit %.1f ms/bin (mean iter %.1f)  | score: null %.2f s once, %.3f ms/bin  -> %.0fx\n",
              n, p, prev, m, 1e3 * t_glm / m, mean(iters), t_null, 1e3 * t_score / m, (t_glm / m) / (t_score / m)))
}
if (n <= 5000 && requireNamespace("logistf", quietly = TRUE)) {
  y <- rbinom(n, 1, 0.5)
  dat <- data.frame(y = y, Z[, -1], cov_resids = 0)
  k <- min(m, 20L)
  t0 <- proc.time()[["elapsed"]]
  for (i in seq_len(k)) { dat$cov_resids <- G[i, ]; invisible(logistf::logistf(y ~ ., data = dat, pl = FALSE)) }
  cat(sprintf("n=%d  logistf (Firth, no profile CI) %.0f ms/bin\n", n, 1e3 * (proc.time()[["elapsed"]] - t0) / k))
}
