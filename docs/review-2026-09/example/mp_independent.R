#!/usr/bin/env Rscript
# Independent Marchenko-Pastur check of choose_ndim.R on the committed spectrum.
suppressPackageStartupMessages(library(data.table))
S <- "/private/tmp/claude-501/-Users-Kitty-git-depthSV/1db72424-c4a1-45f4-b81e-5b466f3ec8a4/scratchpad/review_example/spectrum"
sv <- fread(file.path(S, "svd.singularvalues.txt"))[[2]]
n <- 3202L; p <- 142070L
gamma <- n / p
a <- (1 - sqrt(gamma))^2; b <- (1 + sqrt(gamma))^2
cat(sprintf("n=%d p=%d gamma=%.5f a=%.5f b=%.5f edge_unit=%.3f\n", n, p, gamma, a, b, sqrt(p) + sqrt(n)))

# Independent MP quantile function: numerical inversion of the CDF via integrate().
mp_cdf <- function(x) {
  dens <- function(t) sqrt(pmax(0, (b - t) * (t - a))) / (2 * pi * gamma * t)
  if (x <= a) return(0); if (x >= b) return(1)
  integrate(dens, a, x, rel.tol = 1e-10, subdivisions = 2000L)$value
}
total <- mp_cdf(b - 1e-12)
mp_q <- function(q) uniroot(function(x) mp_cdf(x) / total - q, c(a, b), tol = 1e-12)$root

# Model A (the script's): observed rank j <-> noise rank j.
# Model B (spiked model / interlacing): observed rank j <-> noise rank j - k.
fit <- function(sv, margin, gap = 20L, offset = FALSE) {
  k_rep <- length(sv)
  k <- 0L
  for (it in 1:100) {
    ranks <- (k + gap + 1):k_rep
    if (length(ranks) < 10) return(NA_integer_)
    noise_rank <- if (offset) ranks - k else ranks
    s_unit <- sqrt(p * sapply(1 - (noise_rank - 0.5) / n, mp_q))
    sigma <- exp(mean(log(sv[ranks] / s_unit)))
    edge <- sigma * (sqrt(p) + sqrt(n)) * (1 + margin)
    knew <- sum(sv > edge)
    if (knew == k) return(c(k = k, sigma = sigma, edge = edge))
    k <- knew
  }
  NA
}
for (m in c(0, 0.01, 0.02, 0.05)) {
  fa <- fit(sv, m); fb <- fit(sv, m, offset = TRUE)
  cat(sprintf("margin %.2f: script-model k=%s sigma=%.5f edge=%.3f | rank-offset model k=%s sigma=%.5f edge=%.3f\n",
              m, fa["k"], fa["sigma"], fa["edge"], fb["k"], fb["sigma"], fb["edge"]))
}

# How flat is the observed tail vs MP? print sv ranks 30..60 with the fitted bulk (script model, margin 1%)
fa <- fit(sv, 0.01)
ranks <- 25:60
s_unit <- sqrt(p * sapply(1 - (ranks - 0.5) / n, mp_q))
cat("\nrank  sv  sigma*s_unit  ratio\n")
for (i in seq_along(ranks)) cat(sprintf("%3d %8.3f %8.3f %6.4f\n", ranks[i], sv[ranks[i]], fa["sigma"] * s_unit[i], sv[ranks[i]] / (fa["sigma"] * s_unit[i])))

# Tracy-Widom sd of the largest noise singular value, relative to the edge (Johnstone 2001).
mu_tw <- (sqrt(n - 1) + sqrt(p))^2
sig_tw <- (sqrt(n - 1) + sqrt(p)) * (1 / sqrt(n - 1) + 1 / sqrt(p))^(1 / 3)
sd_tw1 <- 1.2680   # sd of the TW1 distribution
rel_eig <- sig_tw * sd_tw1 / mu_tw
rel_sv <- rel_eig / 2
cat(sprintf("\nTW: eigenvalue centre %.1f scale %.3f -> relative sd of the largest noise eigenvalue %.2e; of the singular value %.2e\n",
            mu_tw, sig_tw, rel_eig, rel_sv))
cat(sprintf("A 1%% margin on the singular value is %.1f TW sd (README/config say 'about four')\n", 0.01 / rel_sv))
cat(sprintf("Ranks 1..%d span a CDF range 1-(0.5/n)=%.5f .. 1-(199.5/n)=%.5f of the bulk\n", 200, 1 - 0.5 / n, 1 - 199.5 / n))

# Edge cases: fewer than gap+10 singular values, and all above the edge
src <- "/Users/Kitty/git/depthSV/.claude/worktrees/agent-aa17d545ebe28b8cc/example/1000G_highcov/R/choose_ndim.R"
cat("\n--- edge case: 15 singular values (fewer than gap+1=21) via the real script\n")
