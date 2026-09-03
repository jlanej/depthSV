# Relatedness: 1000G-like structure (2504 unrelated founders, 563 trio children whose
# parents are among the founders). Linear test of a heritable phenotype against
# (a) a CNV-like bin that segregates in families, (b) a technical bin.
set.seed(3)
nf <- 2504; nc <- 563
par_f <- sample.int(nf, nc); par_m <- sample(setdiff(seq_len(nf), par_f), nc)
n <- nf + nc
# kinship-derived prediction: sum_{i!=j} K_ij^2 / n with K=0.5 for parent-offspring
sumK2 <- 2 * (2 * nc) * 0.25
cat(sprintf("n=%d (%d founders + %d children); sum_{i!=j} K_ij^2 / n = %.3f -> predicted lambda = 1 + %.3f h2 at genetic bins\n\n", n, nf, nc, sumK2 / n, sumK2 / n))

lambda <- function(p) median(qchisq(p, 1, lower.tail = FALSE)) / qchisq(0.5, 1)
lin <- function(X, y) {
  xr <- sweep(X, 2, colMeans(X)); yr <- y - mean(y)
  xx <- colSums(xr^2); beta <- colSums(xr * yr) / xx
  rss <- sum(yr^2) - beta^2 * xx; se <- sqrt(rss / (n - 2) / xx); t <- beta / se
  2 * pt(-abs(t), n - 2)
}
B <- 20000
mk_cnv <- function(f) {
  gf <- matrix(rbinom(nf * B, 2, f), nf, B)                     # founders
  # child: one allele from each parent, transmitted with prob = parent count / 2
  a1 <- matrix(rbinom(nc * B, 1, gf[par_f, ] / 2), nc, B)
  a2 <- matrix(rbinom(nc * B, 1, gf[par_m, ] / 2), nc, B)
  rbind(gf, a1 + a2)
}
Gc <- mk_cnv(0.3)
Xg <- -Gc + rnorm(n * B, sd = 0.2)                             # copy-number-like depth at a CNV bin
Xt <- matrix(rnorm(n * B), n, B)                                # technical bin
cat(sprintf("%-6s %-22s %-22s %-14s\n", "h2", "genetic bin lambda", "technical bin lambda", "P(p<1e-4) gen"))
for (h2 in c(0, 0.25, 0.5, 0.8, 1.0)) {
  set.seed(100 + round(h2 * 100))
  g <- rnorm(nf, sd = sqrt(h2))
  gc <- (g[par_f] + g[par_m]) / 2 + rnorm(nc, sd = sqrt(h2 / 2))
  y <- c(g, gc) + rnorm(n, sd = sqrt(1 - h2))
  pg <- lin(Xg, y); pt_ <- lin(Xt, y)
  cat(sprintf("%-6.2f %-22s %-22s %-14s\n", h2, sprintf("%.3f (pred %.3f)", lambda(pg), 1 + h2 * sumK2 / n), sprintf("%.3f", lambda(pt_)), sprintf("%.5f (nominal 1e-4)", mean(pg < 1e-4))))
}
cat("\nsame at h2 = 0.5 but restricted to the 2504 founders (the KING-unrelated design):\n")
set.seed(150); g <- rnorm(nf, sd = sqrt(0.5)); y <- g + rnorm(nf, sd = sqrt(0.5))
xr <- sweep(Xg[1:nf, ], 2, colMeans(Xg[1:nf, ])); yr <- y - mean(y); xx <- colSums(xr^2); beta <- colSums(xr * yr) / xx
rss <- sum(yr^2) - beta^2 * xx; t <- beta / sqrt(rss / (nf - 2) / xx)
cat(sprintf("founders only: lambda = %.3f\n", lambda(2 * pt(-abs(t), nf - 2))))
cat("\nthe example's permuted null (MTDNA_CN_NULL) is a permutation across samples: it destroys the family\n")
set.seed(151); yp <- sample(c(g, gc <- (g[par_f] + g[par_m]) / 2 + rnorm(nc, sd = 0.5)) + rnorm(n, sd = sqrt(0.5)))
cat(sprintf("structure by construction; genetic-bin lambda under the permuted phenotype = %.3f (h2 = 0.5 before permutation)\n", lambda(lin(Xg, yp))))
