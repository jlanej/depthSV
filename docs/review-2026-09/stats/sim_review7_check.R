# Independent check of REVIEW.md section 7's "lambda_pipeline 0.983 when the phenotype is 49% explained by PC1".
set.seed(31)
n <- 2000; B <- 20000; K <- 16
P <- matrix(rnorm(n * K), n, K)
Q <- qr.Q(qr(cbind(1, P)))
M <- function(v) v - Q %*% crossprod(Q, v)
lambda <- function(p) median(qchisq(p, 1, lower.tail = FALSE)) / qchisq(0.5, 1)
lin <- function(xm, y, Z) { qz <- qr(Z); Qz <- qr.Q(qz); yr <- as.numeric(y - Qz %*% crossprod(Qz, y)); xr <- xm - Qz %*% crossprod(Qz, xm)
  xx <- colSums(xr^2); b <- as.numeric(crossprod(xr, yr)) / xx; df <- n - qz$rank - 1; se <- sqrt((sum(yr^2) - b^2 * xx) / df / xx); 2 * pt(-abs(b / se), df) }
for (R2 in c(0, 0.25, 0.49, 0.75)) {
  y <- sqrt(R2) * P[, 1] + sqrt(1 - R2) * rnorm(n)
  X <- matrix(rnorm(n * B), n, B)
  cat(sprintf("phenotype R^2 on PC1 = %.2f : lambda pipeline (y ~ x_res) = %.3f   full model (y ~ x + PCs) = %.3f   analytic 1-R^2 = %.2f\n",
              R2, lambda(lin(M(X), y, matrix(1, n, 1))), lambda(lin(X, y, cbind(1, P))), 1 - R2))
}
