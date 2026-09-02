# Adjusted-model calibration and estimand, using the REAL 1000G coverage PCs,
# the REAL MTDNA_CN phenotype and the real sex / population labels.
suppressPackageStartupMessages(library(data.table))
S <- "/private/tmp/claude-501/-Users-Kitty-git-depthSV/1db72424-c4a1-45f4-b81e-5b466f3ec8a4/scratchpad/review_stats"
d <- fread(file.path(S, "real_pheno.tsv"))
pcs <- fread("/Users/Kitty/git/NGS-PCA/example/1000G_highcov/output/ngspca_output/svd.pcs.txt")
pcs[, SAMPLE := sub("\\.by1000\\.$", "", SAMPLE)]
d <- merge(d, pcs, by = "SAMPLE")
n <- nrow(d); K <- 42
P <- as.matrix(d[, paste0("PC", 1:K), with = FALSE])
Qpc <- qr.Q(qr(cbind(1, P)))
Mpc <- function(v) v - Qpc %*% crossprod(Qpc, v)          # correct.R residualize()
y_raw <- d$MTDNA_CN
y_log <- log2(y_raw)
y_int <- qnorm((rank(y_raw) - 0.5) / n)
Z1 <- matrix(1, n, 1)
Zadj <- model.matrix(~ SEX + POPULATION, d)              # SEX + population dummies (GPC stand-in)
Zfull <- cbind(Zadj, P)                                   # the proper adjusted model: coverage PCs in the model

lin_test <- function(xmat, y, Z) {                       # analyze.R fit_linear, vectorised over bins
  qz <- qr(Z); Qz <- qr.Q(qz)[, seq_len(qz$rank), drop = FALSE]
  yr <- as.numeric(y - Qz %*% crossprod(Qz, y))
  xr <- xmat - Qz %*% crossprod(Qz, xmat)
  xx <- colSums(xr^2); beta <- as.numeric(crossprod(xr, yr)) / xx
  df <- n - qz$rank - 1L
  rss <- sum(yr^2) - beta^2 * xx
  se <- sqrt(rss / df / xx); t <- beta / se
  list(beta = beta, se = se, t = t, p = 2 * pt(-abs(t), df))
}
lambda <- function(p) median(qchisq(p, 1, lower.tail = FALSE)) / qchisq(0.5, 1)
summ <- function(label, r) cat(sprintf("%-46s lambda=%.3f  P(p<.05)=%.4f  P(p<1e-3)=%.5f  P(p<1e-4)=%.5f  P(p<5e-8)=%.2e  medSE=%.4g\n",
  label, lambda(r$p), mean(r$p < .05), mean(r$p < 1e-3), mean(r$p < 1e-4), mean(r$p < 5e-8), median(r$se)))

set.seed(11)
B <- 20000
cat(sprintf("n=%d samples, K=%d coverage PCs removed, Z = SEX + %d population dummies, %d null bins\n\n", n, K, ncol(Zadj) - 2, B))

# --- 1. null technical bins ------------------------------------------------
X <- matrix(rnorm(n * B, sd = 0.1), n, B)
Xres <- Mpc(X)
cat("1. NULL technical bins (iid noise), phenotype = real MTDNA_CN\n")
summ("pipeline adjusted:   y ~ x_res + SEX + pop", lin_test(Xres, y_raw, Zadj))
summ("pipeline unadjusted: y ~ x_res", lin_test(Xres, y_raw, Z1))
summ("full model:          y ~ x + PC1..42 + SEX + pop", lin_test(X, y_raw, Zfull))
summ("pipeline adjusted, y = rank-INT(MTDNA_CN)", lin_test(Xres, y_int, Zadj))
summ("full model, y = rank-INT(MTDNA_CN)", lin_test(X, y_int, Zfull))
r1 <- lin_test(Xres, y_raw, Zadj); r2 <- lin_test(X, y_raw, Zfull)
cat(sprintf("   median SE ratio pipeline-adjusted / full = %.3f ; median |t| ratio = %.3f\n\n", median(r1$se / r2$se), median(abs(r1$t)) / median(abs(r2$t))))
rm(X, Xres, r1, r2); invisible(gc())

# --- 2. one zero-depth sample per bin (REVIEW 1.2) with the real phenotype ---
cat("2. one zero-depth sample per bin (log2(0.005/33) = -12.7), rest N(0, 0.1); real phenotype; adjusted model\n")
floor_val <- log2(0.005 / 33)
X2 <- matrix(rnorm(n * B, sd = 0.1), n, B)
idx <- sample.int(n, B, replace = TRUE)
X2[cbind(idx, seq_len(B))] <- floor_val
X2res <- Mpc(X2)
h <- 1 - colSums(Qpc[idx, , drop = FALSE]^2)              # leverage of the floored sample after PC projection
cat(sprintf("   floored sample's share of the depth SS after correction: median %.3f (min %.3f)\n", median(X2res[cbind(idx, seq_len(B))]^2 / colSums(X2res^2)), min(X2res[cbind(idx, seq_len(B))]^2 / colSums(X2res^2))))
summ("y = MTDNA_CN (raw)      adjusted", lin_test(X2res, y_raw, Zadj))
summ("y = log2(MTDNA_CN)      adjusted", lin_test(X2res, y_log, Zadj))
summ("y = rank-INT(MTDNA_CN)  adjusted", lin_test(X2res, y_int, Zadj))
summ("y = MTDNA_CN (raw)      unadjusted", lin_test(X2res, y_raw, Z1))
rm(X2, X2res); invisible(gc())
cat("   3 zero-depth samples per bin (a rare homozygous deletion carried by 3 people):\n")
X3 <- matrix(rnorm(n * B, sd = 0.1), n, B)
for (b in seq_len(B)) X3[sample.int(n, 3), b] <- floor_val + runif(3, 0, 4)
X3res <- Mpc(X3)
summ("y = MTDNA_CN (raw)      adjusted", lin_test(X3res, y_raw, Zadj))
summ("y = log2(MTDNA_CN)      adjusted", lin_test(X3res, y_log, Zadj))
summ("y = rank-INT(MTDNA_CN)  adjusted", lin_test(X3res, y_int, Zadj))
rm(X3, X3res); invisible(gc())
cat("\n")

# --- 3. population-stratified common deletion with NO effect on y ----------
cat("3. population-stratified deletion (AF: AFR 0.55, others 0.15), NO effect on MTDNA_CN; 2000 replicate bins\n")
B3 <- 2000
f <- ifelse(d$SUPERPOPULATION == "AFR", 0.55, 0.15)
G <- matrix(rbinom(n * B3, 2, f), n, B3)                   # deletion allele count
log2_code <- function(G) { x <- matrix(0, nrow(G), ncol(G)); x[G == 1] <- -1; x[G == 2] <- runif(sum(G == 2), floor_val, -6); x + rnorm(length(x), sd = 0.1) }
ratio_code <- function(x) 2 * 2^x                          # estimated copy number
X3l <- log2_code(G)
summ("unadjusted 'truth-check' model, log2 coding", lin_test(Mpc(X3l), y_raw, Z1))
summ("adjusted (SEX+pop), log2 coding", lin_test(Mpc(X3l), y_raw, Zadj))
summ("full model (PCs+SEX+pop), log2 coding", lin_test(X3l, y_raw, Zfull))
summ("unadjusted, ratio (copy-number) coding", lin_test(Mpc(ratio_code(X3l)), y_raw, Z1))
cat(sprintf("   share of the log2-coded depth SS carried by homozygous-deleted samples (%.1f%% of samples): median %.3f\n\n",
            100 * mean(G == 2), median(colSums((X3l * (G == 2))^2) / colSums(X3l^2))))

# --- 4. estimand: additive per-copy effect on y ----------------------------
cat("4. deletion (AF 0.3, HWE, not stratified) with an additive effect of delta per deleted copy on y; 2000 bins\n")
delta <- 0.15 * sd(y_raw)
G4 <- matrix(rbinom(n * B3, 2, 0.3), n, B3)
X4l <- log2_code(G4)
X4r <- ratio_code(X4l)
est <- function(label, r, scale = 1) cat(sprintf("%-52s median beta/delta = %6.3f   median |t| = %6.2f   P(p<5e-8) = %.3f\n", label, median(r$beta) / delta * scale, median(abs(r$t)), mean(r$p < 5e-8)))
for (b in seq_len(B3)) { }
Y4 <- matrix(y_raw, n, B3) + delta * G4                    # deleted copies = G (0/1/2)
tt <- function(xm, ym, Z) {                                # per-bin y differs: loop
  qz <- qr(Z); Qz <- qr.Q(qz)[, seq_len(qz$rank), drop = FALSE]; df <- n - qz$rank - 1L
  xr <- xm - Qz %*% crossprod(Qz, xm); yr <- ym - Qz %*% crossprod(Qz, ym)
  xx <- colSums(xr^2); beta <- colSums(xr * yr) / xx
  rss <- colSums(yr^2) - beta^2 * xx; se <- sqrt(rss / df / xx); t <- beta / se
  list(beta = beta, se = se, t = t, p = 2 * pt(-abs(t), df))
}
est("pipeline log2 coding (beta per log2 unit; x -1 = het)", tt(Mpc(X4l), Y4, Zadj), scale = -1)
est("pipeline log2 coding winsorised at -3", tt(Mpc(pmax(X4l, -3)), Y4, Zadj), scale = -1)
est("pipeline ratio coding (beta per copy)", tt(Mpc(X4r), Y4, Zadj), scale = -1)
est("oracle: true copy count as x", tt(Mpc(-G4 + rnorm(n * B3, sd = 0.1)), Y4, Zadj), scale = -1)
cat(sprintf("   correlation between log2 coding and true copy count: median %.3f ; ratio coding: %.3f ; winsorised(-3): %.3f\n\n",
            median(sapply(1:200, function(b) cor(X4l[, b], G4[, b]))), median(sapply(1:200, function(b) cor(X4r[, b], G4[, b]))),
            median(sapply(1:200, function(b) cor(pmax(X4l[, b], -3), G4[, b])))))

# --- 5. estimand under Z correlated with removed PCs (FWL breaks) ----------
cat("5. estimand when the covariate correlates with removed coverage PCs: stratified deletion, additive effect, RATIO coding\n")
G5 <- matrix(rbinom(n * B3, 2, f), n, B3)                  # AFR-enriched deletion
X5 <- -G5 + rnorm(n * B3, sd = 0.15)                      # clean copy-number-like predictor (no floor issue)
Y5 <- matrix(y_raw, n, B3) + delta * G5
est("pipeline adjusted (SEX+pop), x_res", tt(Mpc(X5), Y5, Zadj), scale = -1)
est("pipeline unadjusted, x_res", tt(Mpc(X5), Y5, Z1), scale = -1)
est("full model (PCs + SEX + pop), raw x", tt(X5, Y5, Zfull), scale = -1)
est("oracle: x with no PC removal, SEX+pop", tt(X5, Y5, Zadj), scale = -1)
cat(sprintf("   variance of the stratified deletion vector removed by the PC projection: median %.3f\n",
            median(1 - colSums(Mpc(X5)^2) / colSums(scale(X5, scale = FALSE)^2))))
cat("   same, with y additionally carrying a superpopulation effect not in y already (adds 0.5 sd for AFR):\n")
Y5b <- Y5 + 0.5 * sd(y_raw) * (d$SUPERPOPULATION == "AFR")
est("pipeline adjusted (SEX+pop), x_res", tt(Mpc(X5), Y5b, Zadj), scale = -1)
est("full model (PCs + SEX + pop), raw x", tt(X5, Y5b, Zfull), scale = -1)
