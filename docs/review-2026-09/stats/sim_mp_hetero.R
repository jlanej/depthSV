# Does the Marchenko-Pastur count in choose_ndim.R report components for PURE NOISE
# when the noise is heteroskedastic (per-bin and per-sample variance differ, as in
# real depth data)? Same gamma = n/p as 3202 x 142070, at n = 800.
suppressPackageStartupMessages(library(data.table))
S <- "/private/tmp/claude-501/-Users-Kitty-git-depthSV/1db72424-c4a1-45f4-b81e-5b466f3ec8a4/scratchpad/review_stats"
# --- fit_one copied verbatim from example/1000G_highcov/R/choose_ndim.R ------
mp_quantile_fn <- function(gamma) {
  a <- (1 - sqrt(gamma))^2; b <- (1 + sqrt(gamma))^2
  x <- seq(a, b, length.out = 20001)
  dens <- sqrt(pmax(0, (b - x) * (x - a))) / (2 * pi * gamma * x)
  cdf <- cumsum(c(0, (dens[-1] + dens[-length(dens)]) / 2 * diff(x)))
  cdf <- cdf / max(cdf)
  function(q) approx(cdf, x, xout = q, ties = "ordered")$y
}
fit_one <- function(sv, n, p, margin, gap) {
  k_rep <- length(sv); gamma <- n / p; qf <- mp_quantile_fn(gamma)
  s_unit <- sqrt(p * qf(1 - (seq_len(k_rep) - 0.5) / n)); edge_unit <- sqrt(p) + sqrt(n)
  step <- function(k, m) {
    ranks <- (k + gap + 1):k_rep
    if (length(ranks) < 10) return(NULL)
    sigma <- exp(mean(log(sv[ranks] / s_unit[ranks]))); edge <- sigma * edge_unit * (1 + m)
    list(sigma = sigma, edge = edge, k = sum(sv > edge))
  }
  solve_k <- function(m) { k <- 0L; for (i in 1:100) { f <- step(k, m); if (is.null(f)) return(list(k = NA_integer_)); if (f$k == k) return(f); k <- f$k }; list(k = NA_integer_) }
  sapply(c(0, 0.01, 0.02, 0.05), function(m) solve_k(m)$k)
}
n <- 800; p <- 35495
cov <- fread(file.path(S, "real_pheno.tsv"))$HQ_MEDIAN_COV
set.seed(5)
cov_s <- sample(cov, n)
sd_sample <- sqrt(median(cov) / cov_s)                     # Poisson: log2 noise sd ~ 1/sqrt(coverage)
top_sv <- function(X) svd(X, nu = 0, nv = 0)$d[1:200]
pa_edge <- function(X) { Xp <- apply(X, 2, sample); svd(Xp, nu = 0, nv = 0)$d[1] }   # parallel analysis: column permutation
report <- function(label, X, extra = "") {
  sv <- top_sv(X); k <- fit_one(sv, n, p, 0.01, 20)
  cat(sprintf("%-58s MP count @0/1/2/5%%: %3s %3s %3s %3s | PA (column-permutation) count: %3d %s\n",
              label, k[1], k[2], k[3], k[4], sum(sv > pa_edge(X)), extra))
}
cat(sprintf("n=%d, p=%d (gamma=%.4f). TRUE number of non-noise components is stated per row.\n", n, p, n / p))
E <- matrix(rnorm(n * p), n, p)
report("iid noise (0 true components)", E)
report("per-sample sd from real coverage 27-70x (0 true)", E * sd_sample)
for (s in c(0.25, 0.5)) {
  sd_bin <- exp(rnorm(p, sd = s))
  report(sprintf("per-bin lognormal sd (sdlog=%.2f) (0 true)", s), sweep(E, 2, sd_bin, "*"))
  report(sprintf("per-bin sdlog=%.2f AND per-sample coverage (0 true)", s), sweep(E * sd_sample, 2, sd_bin, "*"))
}
# with real structure: 30 technical factors (decaying) + 12 CNV blocks; noise per-bin sdlog 0.5 + per-sample
sd_bin <- exp(rnorm(p, sd = 0.5))
tech <- matrix(0, n, p)
for (k in 1:30) tech <- tech + (2.0 * 0.9^k) * tcrossprod(rnorm(n), rnorm(p) / sqrt(p) * sqrt(p) / sqrt(p) * 1)
# scale technical factors so that factor k has singular value ~ (6 - 0.15k) x the noise edge... simpler: set explicit sv
tech <- matrix(0, n, p); U <- qr.Q(qr(matrix(rnorm(n * 30), n, 30))); V <- qr.Q(qr(matrix(rnorm(p * 30), p, 30)))
edge0 <- sqrt(n) + sqrt(p)
tech <- U %*% diag(edge0 * seq(4, 1.3, length.out = 30)) %*% t(V)
cnv <- matrix(0, n, p)
cnv_info <- character(0)
set.seed(6)
for (j in 1:12) {
  af <- c(0.05, 0.1, 0.2, 0.3)[(j - 1) %% 4 + 1]; width <- c(3, 10, 30)[(j - 1) %/% 4 + 1]
  g <- rbinom(n, 2, af); code <- ifelse(g == 1, -1, ifelse(g == 2, -8, 0))
  cols <- (200 * j):(200 * j + width - 1)
  cnv[, cols] <- code
  cnv_info <- c(cnv_info, sprintf("AF%.2f x %d bins", af, width))
}
X <- tech + cnv + sweep(E * sd_sample, 2, sd_bin, "*")
sv <- top_sv(X); k <- fit_one(sv, n, p, 0.01, 20)
cat(sprintf("\n30 technical factors (sv 4x..1.3x the iid edge) + 12 CNV blocks + heteroskedastic noise:\n   MP count @0/1/2/5%%: %s %s %s %s ; PA count: %d ; (30 technical + 12 CNVs = 42 real components)\n",
            k[1], k[2], k[3], k[4], sum(sv > pa_edge(X))))
# which of the counted components are CNVs? correlate the top sample-side singular vectors with each CNV genotype
s <- svd(X, nu = 60, nv = 0)
cnv_geno <- sapply(1:12, function(j) { cols <- (200 * j); cnv[, cols] })
r2max <- apply(cnv_geno, 2, function(g) max(cor(s$u[, 1:60], g)^2))
which_pc <- apply(cnv_geno, 2, function(g) which.max(cor(s$u[, 1:60], g)^2))
cat("   CNV blocks and the sample-side PC that carries each (max r^2 with the genotype):\n")
for (j in 1:12) cat(sprintf("     %-18s -> PC%-3d r^2=%.2f %s\n", cnv_info[j], which_pc[j], r2max[j], if (which_pc[j] <= k[2]) "REMOVED at the 1% MP count" else "survives"))
