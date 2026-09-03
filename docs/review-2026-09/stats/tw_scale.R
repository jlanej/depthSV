n <- 3202; p <- 142070
mu <- (sqrt(n-1) + sqrt(p))^2
sig <- (sqrt(n-1) + sqrt(p)) * (1/sqrt(n-1) + 1/sqrt(p))^(1/3)
cat(sprintf("eigenvalue scale: mu=%.1f sigma_TW=%.3f  sigma/mu=%.3e\n", mu, sig, sig/mu))
cat(sprintf("singular-value relative TW scale = sigma/(2 mu) = %.3e\n", sig/(2*mu)))
cat(sprintf("a 1%% margin on the singular value = %.1f TW scale units\n", 0.01/(sig/(2*mu))))
cat(sprintf("TW1 sd ~1.27 => 1%% margin = %.1f TW1 standard deviations\n", 0.01/(sig/(2*mu))/1.268))
# empirical check at the same gamma with n=800
set.seed(1); n2 <- 800; p2 <- round(p*n2/n)
sv <- numeric(20)
for (r in 1:20) { X <- matrix(rnorm(n2*p2), n2, p2); sv[r] <- svd(X, nu=0, nv=0)$d[1] }
edge <- sqrt(n2)+sqrt(p2)
tw_pred <- ((sqrt(n2-1)+sqrt(p2))*(1/sqrt(n2-1)+1/sqrt(p2))^(1/3))/(2*(sqrt(n2-1)+sqrt(p2))^2)*1.268
cat(sprintf("empirical (n=%d,p=%d, 20 reps): mean top sv/edge = %.5f, sd/edge = %.3e ; TW prediction sd/edge = %.3e\n",
  n2, p2, mean(sv)/edge, sd(sv)/edge, tw_pred))
cat(sprintf("at n=800 a 1%% margin = %.1f empirical sd of the top noise singular value\n", 0.01/(sd(sv)/edge)))
