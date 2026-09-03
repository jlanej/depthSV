# Logistic path: Wald (what analyze.R reports) vs likelihood-ratio under separation.
set.seed(8)
n <- 3199; sex <- rep(c(0, 1), length.out = n)      # 1 = male
fit <- function(x, y) {                              # analyze.R fit_logistic, verbatim logic
  X <- cbind(1, cov_resids = x)
  f <- suppressWarnings(glm.fit(X, y, family = binomial()))
  p <- f$rank; piv <- f$qr$pivot[seq_len(p)]; pos <- match(2L, piv)
  cov_unscaled <- chol2inv(f$qr$qr[seq_len(p), seq_len(p), drop = FALSE])
  se <- sqrt(cov_unscaled[pos, pos]); beta <- unname(f$coefficients[2]); z <- beta / se
  f0 <- suppressWarnings(glm.fit(X[, 1, drop = FALSE], y, family = binomial()))
  lrt <- f0$deviance - f$deviance
  c(beta = unname(beta), se = se, z = unname(z), p_wald = 2 * pnorm(-abs(z)), p_lrt = pchisq(lrt, 1, lower.tail = FALSE), converged = f$converged, iter = f$iter)
}
cat("chrY-like bin (males log2 ~ N(-1.3, 0.1), females ~ N(-4.1, 0.5)), response = sex: complete separation\n")
x <- ifelse(sex == 1, rnorm(n, -1.3, 0.1), rnorm(n, -4.1, 0.5)); print(round(fit(x, sex), 4))
cat("\nchrX-like bin (males ~ N(-1, 0.1), females ~ N(0, 0.1)): complete separation\n")
x <- ifelse(sex == 1, rnorm(n, -1, 0.1), rnorm(n, 0, 0.1)); print(round(fit(x, sex), 4))
cat("\nstrong but non-separating effect: y ~ Bernoulli(logit = -2 + 2.5 x), x = deletion coding at AF 0.2 (het -1, hom floor -12)\n")
g <- rbinom(n, 2, 0.2); x <- ifelse(g == 1, -1, ifelse(g == 2, -12, 0)) + rnorm(n, 0, 0.1)
y <- rbinom(n, 1, plogis(-2 + 2.5 * (-x) / 12 * 3)); print(round(fit(x, y), 4))
cat("\nrare homozygous deletion: 3 carriers at the floor, all 3 cases (5% case rate); a Wald p near 1 for the strongest possible signal\n")
y <- rbinom(n, 1, 0.05); x <- rnorm(n, 0, 0.1); i <- sample(which(y == 1), 3); x[i] <- -12.5; print(round(fit(x, y), 4))
cat("\nnull check: 2000 null bins, 5% cases, one floored sample per bin: Wald vs LRT tail rates\n")
res <- t(sapply(1:2000, function(b) { y <- rbinom(n, 1, 0.05); x <- rnorm(n, 0, 0.1); x[sample.int(n, 1)] <- -12.5; fit(x, y)[c("p_wald", "p_lrt")] }))
cat(sprintf("P(p<1e-4): Wald %.4f  LRT %.4f ; P(p<0.05): Wald %.3f LRT %.3f ; Wald NA (analyze.R would skip the bin): %d of 2000\n",
            mean(res[, 1] < 1e-4, na.rm = TRUE), mean(res[, 2] < 1e-4, na.rm = TRUE), mean(res[, 1] < .05, na.rm = TRUE), mean(res[, 2] < .05, na.rm = TRUE), sum(is.na(res[, 1]))))
