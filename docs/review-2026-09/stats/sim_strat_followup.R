# Follow-up: (A) mechanism of the inflation at ancestry-stratified null bins in the
# adjusted pipeline model: systematic bias -x_res' P_Z v / ||x_tilde||^2 (v = PC-part of y)
# vs heteroskedasticity (HC3 sandwich) vs non-normality (rank-INT).
# (B) plateau ratio s150/s200 of heteroskedastic pure noise vs the real spectrum (1.0099).
# (C) discreteness of the real PC scores (a CNV-driven PC has genotype-like clusters).
suppressPackageStartupMessages({ library(data.table); library(sandwich) })
S <- "/private/tmp/claude-501/-Users-Kitty-git-depthSV/1db72424-c4a1-45f4-b81e-5b466f3ec8a4/scratchpad/review_stats"
d <- fread(file.path(S, "real_pheno.tsv"))
pcs <- fread("/Users/Kitty/git/NGS-PCA/example/1000G_highcov/output/ngspca_output/svd.pcs.txt")
pcs[, SAMPLE := sub("\\.by1000\\.$", "", SAMPLE)]
d <- merge(d, pcs, by = "SAMPLE")
n <- nrow(d); K <- 42
P <- as.matrix(d[, paste0("PC", 1:K), with = FALSE]); Qpc <- qr.Q(qr(cbind(1, P)))
Mpc <- function(v) v - Qpc %*% crossprod(Qpc, v)
y <- d$MTDNA_CN; y_int <- qnorm((rank(y) - 0.5) / n)
Zadj <- model.matrix(~ SEX + POPULATION, d)
qz <- qr(Zadj); Qz <- qr.Q(qz)[, seq_len(qz$rank)]
Mz <- function(v) v - Qz %*% crossprod(Qz, v); Pz <- function(v) Qz %*% crossprod(Qz, v)
v <- Qpc %*% crossprod(Qpc, y)                          # PC-part of y (the omitted variable)
cat(sprintf("A. PC-part of MTDNA_CN: R^2 = %.3f ; of that, the share explained by SEX+pop (||P_Z v||^2/||v||^2) = %.3f\n",
            1 - sum((y - v)^2) / sum((y - mean(y))^2), sum(Pz(v - mean(v))^2) / sum((v - mean(v))^2)))
set.seed(21); B <- 2000
f <- ifelse(d$SUPERPOPULATION == "AFR", 0.55, 0.15)
G <- matrix(rbinom(n * B, 2, f), n, B)
X <- matrix(0, n, B); X[G == 1] <- -1; X[G == 2] <- runif(sum(G == 2), log2(0.005 / 33), -6); X <- X + rnorm(n * B, sd = 0.1)
Xres <- Mpc(X); Xt <- Mz(Xres)
yr <- as.numeric(Mz(y)); df <- n - qz$rank - 1L
xx <- colSums(Xt^2); beta <- as.numeric(crossprod(Xt, yr)) / xx
rss <- sum(yr^2) - beta^2 * xx; se <- sqrt(rss / df / xx); t <- beta / se
bias_pred <- -as.numeric(crossprod(Xres, Pz(v))) / xx     # the analytic omitted-PC bias
cat(sprintf("   stratified null deletion, adjusted pipeline model: mean beta = %.3f (sd %.3f); fraction beta>0 = %.3f\n", mean(beta), sd(beta), mean(beta > 0)))
cat(sprintf("   analytic bias -x_res'P_Z v/||x~||^2: mean %.3f ; cor(bias_pred, beta) = %.3f ; mean(beta - bias_pred) = %.3f\n", mean(bias_pred), cor(bias_pred, beta), mean(beta - bias_pred)))
cat(sprintf("   mean t = %.2f (bias in SE units); lambda = %.2f\n", mean(t), median(t^2) / qchisq(0.5, 1)))
# HC3 sandwich on 300 bins, and rank-INT phenotype
hc <- sapply(1:300, function(b) { fit <- lm(y ~ Xres[, b] + Zadj - 1); c(ols = coef(summary(fit))[1, 3], hc3 = coef(fit)[1] / sqrt(vcovHC(fit, type = "HC3")[1, 1])) })
cat(sprintf("   300 bins: lambda OLS = %.2f, HC3 sandwich = %.2f  (a sandwich cannot remove a bias)\n", median(hc[1, ]^2) / qchisq(0.5, 1), median(hc[2, ]^2) / qchisq(0.5, 1)))
yr2 <- as.numeric(Mz(y_int)); beta2 <- as.numeric(crossprod(Xt, yr2)) / xx; rss2 <- sum(yr2^2) - beta2^2 * xx; t2 <- beta2 / sqrt(rss2 / df / xx)
cat(sprintf("   rank-INT phenotype: lambda = %.2f (the PC-part of INT(y) is R^2 = %.3f)\n", median(t2^2) / qchisq(0.5, 1), 1 - sum((y_int - Qpc %*% crossprod(Qpc, y_int))^2) / sum(y_int^2)))
# and with the coverage PCs added to the association model
Zf <- cbind(Zadj, P); qf <- qr(Zf); Qf <- qr.Q(qf)[, seq_len(qf$rank)]
Xf <- X - Qf %*% crossprod(Qf, X); yf <- as.numeric(y - Qf %*% crossprod(Qf, y)); dff <- n - qf$rank - 1L
xxf <- colSums(Xf^2); bf <- as.numeric(crossprod(Xf, yf)) / xxf; tf <- bf / sqrt((sum(yf^2) - bf^2 * xxf) / dff / xxf)
cat(sprintf("   PCs added to the association model: mean beta = %.3f, lambda = %.2f\n", mean(bf), median(tf^2) / qchisq(0.5, 1)))
# a non-stratified deletion (same AF everywhere): does the bias vanish?
G0 <- matrix(rbinom(n * B, 2, 0.3), n, B); X0 <- matrix(0, n, B); X0[G0 == 1] <- -1; X0[G0 == 2] <- runif(sum(G0 == 2), log2(0.005 / 33), -6); X0 <- X0 + rnorm(n * B, sd = 0.1)
Xt0 <- Mz(Mpc(X0)); xx0 <- colSums(Xt0^2); b0 <- as.numeric(crossprod(Xt0, yr)) / xx0; t0 <- b0 / sqrt((sum(yr^2) - b0^2 * xx0) / df / xx0)
cat(sprintf("   NON-stratified deletion (AF 0.3 everywhere), adjusted pipeline model: lambda = %.2f, P(p<.05) = %.3f\n\n", median(t0^2) / qchisq(0.5, 1), mean(2 * pt(-abs(t0), df) < .05)))

# B. plateau ratio of heteroskedastic noise -------------------------------------------
cat("B. plateau ratio s[150]/s[200] (real spectrum: 1.0099; iid MP at 3202x142070: 1.0056):\n")
n2 <- 800; p2 <- 35495; set.seed(22)
cov <- d$HQ_MEDIAN_COV; sd_sample <- sqrt(median(cov) / sample(cov, n2))
E <- matrix(rnorm(n2 * p2), n2, p2)
mp_ratio <- { g <- n2 / p2; a <- (1 - sqrt(g))^2; b <- (1 + sqrt(g))^2; x <- seq(a, b, length.out = 20001)
  dens <- sqrt(pmax(0, (b - x) * (x - a))) / (2 * pi * g * x); cdf <- cumsum(c(0, (dens[-1] + dens[-length(dens)]) / 2 * diff(x))); cdf <- cdf / max(cdf)
  q <- approx(cdf, x, xout = 1 - (c(150, 200) - 0.5) / n2, ties = "ordered")$y; sqrt(q[1] / q[2]) }
cat(sprintf("   iid MP prediction at n=%d: %.4f\n", n2, mp_ratio))
for (s in c(0, 0.25, 0.4, 0.5, 0.6)) {
  sdb <- exp(rnorm(p2, sd = s)); X <- sweep(E * sd_sample, 2, sdb, "*"); sv <- svd(X, nu = 0, nv = 0)$d
  cat(sprintf("   per-bin sdlog = %.2f (+ real per-sample coverage): s150/s200 = %.4f ; ratio to MP = %.4f\n", s, sv[150] / sv[200], (sv[150] / sv[200]) / mp_ratio))
}
cat(sprintf("   real: ratio to MP = %.4f\n\n", 1.0099 / 1.0056))

# C. discreteness of real PC scores ---------------------------------------------------------
cat("C. real coverage PCs: bimodality coefficient (>0.555 suggests clusters), max |score|/sd (single-sample driven), R^2 on batch\n")
bc <- function(x) { m <- length(x); s <- (x - mean(x)) / sd(x); g <- mean(s^3); k <- mean(s^4) - 3; (g^2 + 1) / (k + 3 * (m - 1)^2 / ((m - 2) * (m - 3))) }
batch <- factor(d$RELATEDNESS)
out <- sapply(1:60, function(i) { x <- P[, min(i, K)]; if (i > K) x <- d[[paste0("PC", i)]]; c(bc = bc(x), maxz = max(abs(x - mean(x))) / sd(x), r2batch = summary(lm(x ~ batch))$r.squared) })
cat(sprintf("   PC: %s\n", paste(sprintf("%4d", 1:20), collapse = "")))
cat(sprintf("   BC: %s\n", paste(sprintf("%4.2f", out["bc", 1:20]), collapse = "")))
cat(sprintf("   mz: %s\n", paste(sprintf("%4.0f", out["maxz", 1:20]), collapse = "")))
cat(sprintf("   PCs 21-60 with BC > 0.555: %s ; PCs 1-60 with max|z| > 8 (one-sample PCs): %s\n",
            paste(which(out["bc", ] > 0.555 & seq_len(60) > 20), collapse = ","), paste(which(out["maxz", ] > 8), collapse = ",")))
