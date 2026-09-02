set.seed(7); n <- 3000; k <- 10; B <- 1500
Q <- qr.Q(qr(cbind(1, matrix(rnorm(n * k), n, k))))          # orthonormal basis incl. intercept, as correct.R builds it
a <- c(0, rnorm(k)); y_pc <- as.numeric(Q %*% a); y_pc <- y_pc / sd(y_pc) * sqrt(0.5)
y <- y_pc + rnorm(n, sd = sqrt(0.5))                            # phenotype ~50% explained by the PCs
cat(sprintf("R2(y ~ PCs) = %.3f\n", summary(lm(y ~ Q[, -1]))$r.squared))
resid_Q <- function(v) as.numeric(v - Q %*% crossprod(Q, v))    # exactly correct.R's projection
lam <- function(t) median(t^2) / qchisq(0.5, 1)
tstat <- function(fit) summary(fit)$coefficients[2, 3]

# (i) technical null bins: x independent of y. Pipeline = y ~ x_res ; full = y ~ x + PCs
tp <- tf <- numeric(B)
for (b in 1:B) { x <- rnorm(n); xr <- resid_Q(x)
  tp[b] <- tstat(lm(y ~ xr)); tf[b] <- tstat(lm(y ~ x + Q[, -1])) }
cat(sprintf("(i)  technical null:        lambda pipeline=%.3f  full=%.3f   P(p<.05) pipeline=%.3f full=%.3f\n",
            lam(tp), lam(tf), mean(2*pt(-abs(tp), n-2) < .05), mean(2*pt(-abs(tf), n-k-2) < .05)))

# (ii) a covariate z correlated with a coverage PC; bins correlated with z but NULL for y given (z, PCs)
z <- 3 * Q[, 2] + rnorm(n)
tp <- tf <- numeric(B)
for (b in 1:B) { x <- 0.8 * z + rnorm(n); xr <- resid_Q(x)
  tp[b] <- tstat(lm(y ~ xr + z)); tf[b] <- tstat(lm(y ~ x + z + Q[, -1])) }
cat(sprintf("(ii) z-correlated null bins: lambda pipeline=%.3f (mean t=%+.2f, %.0f%% same sign)  full=%.3f (mean t=%+.2f)\n",
            lam(tp), mean(tp), 100*max(mean(tp > 0), mean(tp < 0)), lam(tf), mean(tf)))

# (iii) the remedy: PCs appended to the association model's covariates, depth still residualised (FWL)
tp <- numeric(B)
for (b in 1:B) { x <- 0.8 * z + rnorm(n); xr <- resid_Q(x); tp[b] <- tstat(lm(y ~ xr + z + Q[, -1])) }
cat(sprintf("(iii) same bins, PCs in the model: lambda=%.3f (mean t=%+.2f)\n", lam(tp), mean(tp)))
