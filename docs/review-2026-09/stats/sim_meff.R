# Effective number of tests for correlated adjacent bins, and the cost of a
# permutation-based genome-wide threshold in the linear path.
set.seed(7)
n <- 2000; M <- 30000; R <- 400
lambda <- function(p) median(qchisq(p, 1, lower.tail = FALSE)) / qchisq(0.5, 1)
Y <- matrix(rnorm(n * R), n, R); Y <- sweep(Y, 2, colMeans(Y))
Y <- sweep(Y, 2, sqrt(colSums(Y^2)), "/")
for (rho in c(0, 0.3, 0.6, 0.9, 0.97)) {
  E <- matrix(rnorm(n * M), n, M)
  X <- E
  if (rho > 0) for (j in 2:M) X[, j] <- rho * X[, j - 1] + sqrt(1 - rho^2) * E[, j]
  X <- sweep(X, 2, colMeans(X)); X <- sweep(X, 2, sqrt(colSums(X^2)), "/")
  t0 <- proc.time()[3]
  Rm <- crossprod(X, Y)                               # M x R correlations = all permutation statistics at once
  el <- proc.time()[3] - t0
  tt <- Rm * sqrt((n - 2) / (1 - Rm^2)); P <- 2 * pt(-abs(tt), n - 2)
  minp <- apply(P, 2, min)
  thr <- quantile(minp, 0.05)
  meff <- 0.05 / thr
  fwer_5e8_scaled <- mean(minp < 5e-8 * (M / 3.1e6))   # not meaningful at M=30000; instead report FWER at Bonferroni and at M_eff
  cat(sprintf("rho=%.2f  neighbour r(t)=%.3f  M=%d  5%%-quantile(min p)=%.2e  M_eff=%.0f (%.0f%% of M)  FWER at 0.05/M = %.3f  lambda(all)=%.3f  [%d perms in %.1fs]\n",
              rho, cor(tt[-1, 1], tt[-M, 1]), M, thr, meff, 100 * meff / M, mean(minp < 0.05 / M), lambda(P[, 1]), R, el))
}
cat("\nScaling: at 3.1M bins x 3202 samples the same product for 1000 permutations is ~2e13 flop, i.e. minutes on one node,\n")
cat("and it is a by-product of the projection analyze.R already computes (x_res' y for one y becomes x_res' Y for 1000 y).\n")
cat(sprintf("\nFor reference, with M=3.1e6 INDEPENDENT tests: P(any p < 5e-8) = %.3f ; Bonferroni 0.05/M = %.2e\n", 1 - (1 - 5e-8)^3.1e6, 0.05 / 3.1e6))
