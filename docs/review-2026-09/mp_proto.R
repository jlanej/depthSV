suppressPackageStartupMessages(library(data.table))
sv <- fread("svd.singularvalues.txt")$SINGULAR_VALUES
n <- 3202; p <- 142070; k_reported <- length(sv)
gamma <- n / p
# MP law for eigenvalues x of (1/p) X X^T, unit variance: support [a, b]
a <- (1 - sqrt(gamma))^2; b <- (1 + sqrt(gamma))^2
grid <- seq(a, b, length.out = 20001)
dens <- sqrt(pmax(0, (b - grid) * (grid - a))) / (2 * pi * gamma * grid)
cdf <- cumsum(c(0, (dens[-1] + dens[-length(dens)]) / 2 * diff(grid))); cdf <- cdf / max(cdf)
mp_quantile <- function(q) approx(cdf, grid, xout = q, ties = "ordered")$y   # eigenvalue at CDF q, unit sigma
rank_q <- 1 - (seq_len(k_reported) - 0.5) / n
s_unit <- sqrt(p * mp_quantile(rank_q))      # noise singular value expected at each rank, sigma=1
edge_unit <- sqrt(p) * (1 + sqrt(gamma))     # = sqrt(p)+sqrt(n)
fit <- function(k, gap = 20, margin = 0) {
  ranks <- (k + gap + 1):k_reported
  sigma <- exp(mean(log(sv[ranks] / s_unit[ranks])))
  edge <- sigma * edge_unit * (1 + margin)
  list(sigma = sigma, edge = edge, k = sum(sv > edge))
}
for (margin in c(0, 0.01, 0.02, 0.05)) {
  k <- 0
  for (it in 1:50) { f <- fit(k, margin = margin); if (f$k == k) break; k <- f$k }
  cat(sprintf("margin=%.2f  sigma=%.4f  edge=%.2f  k=%d  (s_k=%.2f  s_k+1=%.2f)\n",
              margin, f$sigma, f$edge, f$k, sv[f$k], sv[f$k + 1]))
}
f <- fit(0); k0 <- 0; for (it in 1:50) { f <- fit(k0); if (f$k == k0) break; k0 <- f$k }
png("mp_proto.png", width = 1400, height = 600, res = 130)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
plot(seq_len(k_reported), sv, log = "y", pch = 16, cex = 0.6, xlab = "rank", ylab = "singular value (log)",
     main = "NGS-PCA spectrum vs Marchenko-Pastur bulk")
lines(seq_len(k_reported), f$sigma * s_unit, col = "steelblue", lwd = 2)
abline(h = f$edge, col = "firebrick", lty = 2); abline(v = f$k + 0.5, col = "firebrick", lty = 3)
legend("topright", c("observed", "MP bulk (fitted sigma)", "MP edge"), col = c("black", "steelblue", "firebrick"),
       pch = c(16, NA, NA), lty = c(NA, 1, 2), bty = "n")
r <- (f$k - 10):k_reported; r <- r[r >= 1]
plot(r, sv[r], pch = 16, cex = 0.6, xlab = "rank", ylab = "singular value", main = sprintf("tail: k_MP = %d", f$k))
lines(r, f$sigma * s_unit[r], col = "steelblue", lwd = 2); abline(h = f$edge, col = "firebrick", lty = 2)
dev.off()
cat(sprintf("\nplateau check: s[150]/s[200] = %.4f ; predicted unit ratio = %.4f\n", sv[150]/sv[200], s_unit[150]/s_unit[200]))
