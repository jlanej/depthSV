# Cox path (R/analyze.R fit_coxph): behaviour on a floored sample and under separation.
suppressPackageStartupMessages(library(survival))
set.seed(41); n <- 3199
fit_coxph <- function(dt) {                                  # analyze.R:268-276 logic
  f <- try(coxph(Surv(time, event) ~ cov_resids + age, data = dt), silent = TRUE)
  if (inherits(f, "try-error")) return("try-error")
  cf <- coef(summary(f)); round(as.numeric(cf["cov_resids", ]), 4)
}
base <- function() { age <- rnorm(n, 50, 10); time <- rexp(n, 0.02); event <- rbinom(n, 1, 0.3); data.frame(time, event, age, cov_resids = rnorm(n, 0, 0.1)) }
cat("columns: coef exp(coef) se(coef) z p\n")
dt <- base(); cat("null bin, no outlier:                              ", fit_coxph(dt), "\n")
dt <- base(); i <- which(dt$event == 1)[1]; dt$cov_resids[i] <- -12.7; dt$time[i] <- 0.01
w <- NULL; r <- withCallingHandlers(fit_coxph(dt), warning = function(x) { w <<- conditionMessage(x); invokeRestart("muffleWarning") })
cat("one floored sample (-12.7) with the earliest event:  ", r, " warning:", if (is.null(w)) "none" else w, "\n")
dt <- base(); j <- sample.int(n, 3); dt$cov_resids[j] <- -12.7; dt$event[j] <- 1; dt$time[j] <- c(0.01, 0.02, 0.03)
w <- NULL; r <- withCallingHandlers(fit_coxph(dt), warning = function(x) { w <<- conditionMessage(x); invokeRestart("muffleWarning") })
cat("3 floored carriers, all with the 3 earliest events:  ", r, " warning:", if (is.null(w)) "none" else w, "\n")
# null tail rates with one floored sample per bin, 2000 bins, exponential (non-normal) survival
p <- replicate(2000, { dt <- base(); dt$cov_resids[sample.int(n, 1)] <- -12.7; suppressWarnings(fit_coxph(dt))[5] })
cat(sprintf("2000 null bins, one floored sample each: P(p<0.05)=%.3f  P(p<1e-3)=%.4f  P(p<1e-4)=%.4f  min p=%.2e\n", mean(p < .05), mean(p < 1e-3), mean(p < 1e-4), min(p)))
