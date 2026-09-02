# Sex chromosomes under the pipeline's normalisation (every contig / autosomal median,
# no ploidy model), with the example's adjusted design (SEX + population covariates,
# real 1000G sexes, real coverage PCs removed, real MTDNA_CN as the phenotype).
suppressPackageStartupMessages(library(data.table))
S <- "/private/tmp/claude-501/-Users-Kitty-git-depthSV/1db72424-c4a1-45f4-b81e-5b466f3ec8a4/scratchpad/review_stats"
d <- fread(file.path(S, "real_pheno.tsv"))
pcs <- fread("/Users/Kitty/git/NGS-PCA/example/1000G_highcov/output/ngspca_output/svd.pcs.txt")
pcs[, SAMPLE := sub("\\.by1000\\.$", "", SAMPLE)]
d <- merge(d, pcs, by = "SAMPLE")
n <- nrow(d); K <- 42
P <- as.matrix(d[, paste0("PC", 1:K), with = FALSE]); Qpc <- qr.Q(qr(cbind(1, P)))
Mpc <- function(v) v - Qpc %*% crossprod(Qpc, v)
male <- d$SEX == 1
y_raw <- d$MTDNA_CN; y_int <- qnorm((rank(y_raw) - 0.5) / n)
Zadj <- model.matrix(~ SEX + POPULATION, d)
Zpop <- model.matrix(~ POPULATION, d)
med <- d$HQ_MEDIAN_COV
floor_of <- log2(0.005 / med)
lambda <- function(p) median(qchisq(p, 1, lower.tail = FALSE)) / qchisq(0.5, 1)
tt <- function(xm, ym, Z, rows = seq_len(nrow(xm))) {
  xm <- xm[rows, , drop = FALSE]; ym <- ym[rows, , drop = FALSE]; Z <- Z[rows, , drop = FALSE]
  qz <- qr(Z); Qz <- qr.Q(qz)[, seq_len(qz$rank), drop = FALSE]; df <- length(rows) - qz$rank - 1L
  xr <- xm - Qz %*% crossprod(Qz, xm); yr <- ym - Qz %*% crossprod(Qz, ym)
  xx <- colSums(xr^2); beta <- colSums(xr * yr) / xx
  rss <- colSums(yr^2) - beta^2 * xx; se <- sqrt(rss / df / xx); t <- beta / se
  list(beta = beta, se = se, t = t, p = 2 * pt(-abs(t), df), xr = xr)
}
summ <- function(label, r) cat(sprintf("%-60s lambda=%.3f  P(p<.05)=%.4f  P(p<1e-4)=%.5f  P(p<5e-8)=%.2e\n", label, lambda(r$p), mean(r$p < .05), mean(r$p < 1e-4), mean(r$p < 5e-8)))
set.seed(12)
B <- 8000

# --- chrY bins, null phenotype (real MTDNA_CN), SEX in the model ------------
# males: chrY depth ~ 0.41 x median with Poisson-ish noise (log2 sd 0.1)
# females: mismapped reads; per-sample rate from the real Y_COV_RATIO (median 0.057, up to 0.39),
#          per-bin Poisson at ~13 reads/kb -> log2 sd ~0.4, and 10% of bins with essentially no reads (floor)
cat("chrY bins (null phenotype = real MTDNA_CN):\n")
mk_y <- function(zero_frac) {
  X <- matrix(0, n, B)
  X[male, ] <- log2(0.41) + rnorm(sum(male) * B, sd = 0.1)
  reads <- matrix(rpois(sum(!male) * B, lambda = d$Y_COV_RATIO[!male] * med[!male] * 1000 / 150), sum(!male), B)
  zero <- matrix(runif(sum(!male) * B) < zero_frac, sum(!male), B)
  reads[zero] <- 0
  X[!male, ] <- log2(pmax(reads * 150 / 1000, 0.005) / med[!male])
  X
}
Zfull <- cbind(Zadj, P)                                   # coverage PCs in the association model
Xy <- mk_y(0.0)
# vectorise y: same y for all bins
tt1 <- function(xm, y, Z, rows = seq_len(n)) { qz <- qr(Z[rows, , drop = FALSE]); Qz <- qr.Q(qz)[, seq_len(qz$rank), drop = FALSE]; df <- length(rows) - qz$rank - 1L
  xr <- xm[rows, , drop = FALSE] - Qz %*% crossprod(Qz, xm[rows, , drop = FALSE]); yr <- as.numeric(y[rows] - Qz %*% crossprod(Qz, y[rows]))
  xx <- colSums(xr^2); beta <- as.numeric(crossprod(xr, yr)) / xx; rss <- sum(yr^2) - beta^2 * xx; se <- sqrt(rss / df / xx); t <- beta / se
  list(beta = beta, se = se, t = t, p = 2 * pt(-abs(t), df), xr = xr) }
r <- tt1(Mpc(Xy), y_raw, Zadj)
fem_share <- colSums(r$xr[!male, ]^2) / colSums(r$xr^2)
cat(sprintf("   no zero-read bins: share of the within-sex depth SS carried by FEMALES: median %.3f\n", median(fem_share)))
summ("   y raw, pooled with SEX+pop covariates", r)
summ("   y rank-INT, pooled with SEX+pop", tt1(Mpc(Xy), y_int, Zadj))
summ("   y raw, MALES ONLY (pop covariates)", tt1(Mpc(Xy), y_raw, Zpop, rows = which(male)))
summ("   y raw, pooled, coverage PCs IN the model (raw x)", tt1(Xy, y_raw, Zfull))
summ("   y raw, unadjusted pipeline model (no covariates)", tt1(Mpc(Xy), y_raw, matrix(1, n, 1)))
art <- -2.84 * (Qpc %*% crossprod(Qpc, as.numeric(male)))  # the fixed artefact -c * P_PC(sex) carried by every chrY bin
cat(sprintf("   fixed artefact -c*P_PC(sex) in x_res: ||.||^2 = %.1f vs within-sex noise SS per bin ~%.1f ; cor(artefact, PC-part of y) = %.3f\n",
            sum(art^2), median(colSums(r$xr^2)) , cor(art, Qpc %*% crossprod(Qpc, y_raw))))
Xy2 <- mk_y(0.10)
r <- tt1(Mpc(Xy2), y_raw, Zadj)
cat(sprintf("   10%% of female bin-values at zero reads (floor): female share of SS median %.3f; max single-sample share median %.3f\n",
            median(colSums(r$xr[!male, ]^2) / colSums(r$xr^2)), median(apply(r$xr^2, 2, max) / colSums(r$xr^2))))
summ("   y raw, pooled with SEX+pop covariates", r)
summ("   y rank-INT, pooled with SEX+pop", tt1(Mpc(Xy2), y_int, Zadj))
summ("   y raw, MALES ONLY (pop covariates)", tt1(Mpc(Xy2), y_raw, Zpop, rows = which(male)))
summ("   y raw, pooled, coverage PCs IN the model (raw x)", tt1(Xy2, y_raw, Zfull))

# --- chrX non-PAR NULL bins (no CNV): males -1, females 0, noise 0.1 ----------
cat("\nchrX non-PAR null bins (males log2 -1, females 0, noise sd 0.1), phenotype = real MTDNA_CN:\n")
Xx <- matrix(ifelse(male, -1, 0), n, B) + rnorm(n * B, sd = 0.1)
summ("   y raw, pooled with SEX+pop covariates", tt1(Mpc(Xx), y_raw, Zadj))
summ("   y rank-INT, pooled with SEX+pop", tt1(Mpc(Xx), y_int, Zadj))
summ("   y raw, unadjusted pipeline model", tt1(Mpc(Xx), y_raw, matrix(1, n, 1)))
summ("   y raw, pooled, coverage PCs IN the model (raw x)", tt1(Xx, y_raw, Zfull))
summ("   y raw, FEMALES ONLY (pop covariates)", tt1(Mpc(Xx), y_raw, Zpop, rows = which(!male)))
rm(Xx, Xy, Xy2); invisible(gc())

# --- chrX non-PAR deletion: same per-copy effect, sex-dependent log2 coding ---
cat("\nchrX non-PAR bin with a deletion at allele frequency 0.10; per-deleted-copy effect delta on y:\n")
B2 <- 2000; delta <- 0.15 * sd(y_raw); floor_v <- floor_of
mkx <- function(af) {
  G <- matrix(0L, n, B2)
  G[!male, ] <- rbinom(sum(!male) * B2, 2, af)            # females: 2 copies
  G[male, ]  <- rbinom(sum(male) * B2, 1, af)             # males: 1 copy
  G
}
G <- mkx(0.10)
cn_expected <- ifelse(male, 1, 2)
cn <- cn_expected - G                                     # copies remaining
X <- log2(pmax(cn * med / 2 * exp(rnorm(n * B2, sd = 0.07)), 0.005) / med)   # depth = cn/2 x median (+noise), floored
dim(X) <- c(n, B2)
Y <- matrix(y_raw, n, B2) + delta * G
cat(sprintf("   log2 value by class: female CN2 %.2f, CN1 %.2f, CN0 %.1f ; male CN1 %.2f, CN0 %.1f\n",
            0, -1, median(floor_v), -1, median(floor_v)))
for (lab in c("pooled, SEX+pop covariates", "females only", "males only")) {
  rows <- switch(lab, "pooled, SEX+pop covariates" = seq_len(n), "females only" = which(!male), "males only" = which(male))
  Z <- if (lab == "pooled, SEX+pop covariates") Zadj else Zpop
  r <- tt(Mpc(X), Y, Z, rows)
  cat(sprintf("   %-32s median beta = %8.3f (delta per copy = %.2f)  median |t| = %6.2f  P(p<5e-8) = %.3f\n", lab, median(r$beta), delta, median(abs(r$t)), mean(r$p < 5e-8)))
}
cat("   the same with a ploidy-aware coding (copies relative to the sex's expected copies: 2*2^x/expected):\n")
Xcn <- 2 * 2^X / cn_expected                              # fraction of expected copies (1 = normal)
for (lab in c("pooled, SEX+pop covariates", "females only", "males only")) {
  rows <- switch(lab, "pooled, SEX+pop covariates" = seq_len(n), "females only" = which(!male), "males only" = which(male))
  Z <- if (lab == "pooled, SEX+pop covariates") Zadj else Zpop
  r <- tt(Mpc(Xcn), Y, Z, rows)
  cat(sprintf("   %-32s median beta = %8.3f  median |t| = %6.2f  P(p<5e-8) = %.3f\n", lab, median(r$beta), median(abs(r$t)), mean(r$p < 5e-8)))
}
cat("\nchrX non-PAR, NULL (no effect), pooled with SEX+pop, real MTDNA_CN: is the sex-dependent coding a calibration problem?\n")
Y0 <- matrix(y_raw, n, B2)
r <- tt(Mpc(X), Y0, Zadj); summ("   log2 coding, y raw", r)
r <- tt(Mpc(X), matrix(y_int, n, B2), Zadj); summ("   log2 coding, y rank-INT", r)
cat("\nsex_linear check (SEX ~ cov_resids): the chrY/chrX |t| it asserts on is the sex dichotomy itself:\n")
r <- tt1(Mpc(Xy[, 1:50]), as.numeric(male), Zpop); cat(sprintf("   chrY bins: median |t| = %.0f\n", median(abs(r$t))))
