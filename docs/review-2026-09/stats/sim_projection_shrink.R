# Genotype-PC projection: in-sample PCs vs projected held-out samples at gamma = n/m of the
# 1000G preamble (2500 unrelated / ~150k pruned SNPs = 0.0167). The preamble's calibration
# regresses the in-sample projection on the in-sample PC, which is exactly proportional by
# construction and cannot see out-of-sample shrinkage.
set.seed(9)
ntr <- 250; nte <- 250; m <- 15000; gamma <- ntr / m
ell <- c(100, 30, 10, 5, 3, 2, 1.5)                     # target GRM eigenvalues (noise units, edge = (1+sqrt(gamma))^2)
r <- length(ell)
V <- qr.Q(qr(matrix(rnorm(m * r), m, r)))
a <- sqrt((ell - 1) * m / ntr)                           # per-SNP spike amplitude giving GRM eigenvalue ~ ell
mk <- function(nn) { U <- matrix(rnorm(nn * r), nn, r); list(U = U, X = U %*% (diag(a) %*% t(V)) + matrix(rnorm(nn * m), nn, m)) }
tr <- mk(ntr); te <- mk(nte)
Xc <- scale(tr$X, scale = FALSE); mu <- attr(Xc, "scaled:center")
s <- svd(Xc, nu = r, nv = r)
lam <- s$d[1:r]^2 / m                                    # GRM eigenvalues in noise units (as plink2 reports)
score_tr <- Xc %*% s$v[, 1:r]                            # in-sample scores == u * d
score_te <- sweep(te$X, 2, mu) %*% s$v[, 1:r]            # projection with the training loadings
cal <- sapply(1:r, function(k) { x <- score_tr[, k]; y <- s$u[, k] * s$d[k]; c(scale = sum(x * y) / sum(x * x), r2 = cor(x, y)^2) })
# the truth for a held-out sample: its population score a_k * U_k; how well does the projection recover it?
cat(sprintf("gamma = %.4f, MP edge = %.3f (noise units)\n", gamma, (1 + sqrt(gamma))^2))
cat(sprintf("%-4s %-10s %-12s %-14s %-20s %-22s\n", "PC", "target eig", "sample eig", "calib r^2", "sd(test)/sd(train)", "cor(test proj, truth)"))
for (k in 1:r) {
  cat(sprintf("%-4d %-10.1f %-12.2f %-14.4f %-20.3f %-22.3f\n", k, ell[k], lam[k], cal["r2", k], sd(score_te[, k]) / sd(score_tr[, k]), abs(cor(score_te[, k], te$U[, k]))))
}
cat("\nReading: calibration r^2 is 1.000 for every PC (the in-sample projection is the in-sample PC), while the\n")
cat("held-out projections are shrunk by the factor in the last two columns; the shrinkage is what the 698 projected\n")
cat("relatives carry into the covariate table. Shrinkage is negligible for eigenvalues >> the MP edge and large near it.\n")
