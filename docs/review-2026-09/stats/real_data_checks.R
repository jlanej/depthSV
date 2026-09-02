suppressPackageStartupMessages(library(data.table))
D <- "/Users/Kitty/git/NGS-PCA/example/1000G_highcov/output"
qc  <- fread(file.path(D, "qc_output/sample_qc.tsv"))
pcs <- fread(file.path(D, "ngspca_output/svd.pcs.txt"))
sv  <- fread(file.path(D, "ngspca_output/svd.singularvalues.txt"))[[2]]
pcs[, SAMPLE := sub("\\.by1000\\.$", "", SAMPLE)]
d <- merge(qc, pcs, by.x = "SAMPLE_ID", by.y = "SAMPLE")
cat("merged samples:", nrow(d), "\n")
cat("INFERRED_SEX:", paste(names(table(d$INFERRED_SEX)), table(d$INFERRED_SEX)), "\n")
cat("REPORTED vs INFERRED mismatches:", sum(d$REPORTED_SEX != d$INFERRED_SEX, na.rm = TRUE), "\n")
cat("SUPERPOP:", paste(names(table(d$SUPERPOPULATION)), table(d$SUPERPOPULATION)), "\n")
cat("RELATEDNESS:", paste(names(table(d$RELATEDNESS)), table(d$RELATEDNESS)), "\n")
cat("FAMILY_ROLE:", paste(names(table(d$FAMILY_ROLE)), table(d$FAMILY_ROLE)), "\n")
cat("RELEASE_BATCH:", paste(names(table(d$RELEASE_BATCH)), table(d$RELEASE_BATCH)), "\n")
cat("INSTRUMENT:", paste(names(table(d$INSTRUMENT_MODEL)), table(d$INSTRUMENT_MODEL)), "\n")
cat("n LIBRARY_NAME levels:", length(unique(d$LIBRARY_NAME)), "\n")

# --- phenotype distribution ------------------------------------------------
sk <- function(x) { x <- x[is.finite(x)]; m <- mean(x); s <- sd(x); c(skew = mean((x-m)^3)/s^3, kurt = mean((x-m)^4)/s^4 - 3) }
y <- d$MTDNA_CN
cat(sprintf("\nMTDNA_CN: n=%d median=%.1f mean=%.1f sd=%.1f min=%.1f max=%.1f max/median=%.2f skew=%.2f exkurt=%.2f\n",
  sum(is.finite(y)), median(y), mean(y), sd(y), min(y), max(y), max(y)/median(y), sk(y)[1], sk(y)[2]))
ly <- log2(y)
cat(sprintf("log2(MTDNA_CN): sd=%.3f skew=%.2f exkurt=%.2f\n", sd(ly), sk(ly)[1], sk(ly)[2]))
cat(sprintf("HQ_MEDIAN_COV: median=%.1f min=%.1f max=%.1f ; log2 noise sd for a 1kb bin ~ 1/sqrt(reads): at min %.3f, median %.3f, max %.3f (150bp reads)\n",
  median(d$HQ_MEDIAN_COV), min(d$HQ_MEDIAN_COV), max(d$HQ_MEDIAN_COV),
  1/sqrt(min(d$HQ_MEDIAN_COV)*1000/150)/log(2), 1/sqrt(median(d$HQ_MEDIAN_COV)*1000/150)/log(2), 1/sqrt(max(d$HQ_MEDIAN_COV)*1000/150)/log(2)))
cat(sprintf("X_COV_RATIO by sex: F median %.3f, M median %.3f ; Y_COV_RATIO: F median %.4f (min %.4f max %.4f), M median %.3f\n",
  median(d$X_COV_RATIO[d$INFERRED_SEX=="F"]), median(d$X_COV_RATIO[d$INFERRED_SEX=="M"]),
  median(d$Y_COV_RATIO[d$INFERRED_SEX=="F"]), min(d$Y_COV_RATIO[d$INFERRED_SEX=="F"]), max(d$Y_COV_RATIO[d$INFERRED_SEX=="F"]),
  median(d$Y_COV_RATIO[d$INFERRED_SEX=="M"])))
cat(sprintf("log2 of female chrY ratio: median %.2f, range %.2f .. %.2f (vs floor log2(0.005/30) = %.2f)\n",
  median(log2(d$Y_COV_RATIO[d$INFERRED_SEX=="F"])), min(log2(d$Y_COV_RATIO[d$INFERRED_SEX=="F"])), max(log2(d$Y_COV_RATIO[d$INFERRED_SEX=="F"])), log2(0.005/30)))

# --- how much of each phenotype/covariate the removed PCs explain ---------
pcmat <- as.matrix(d[, grep("^PC[0-9]+$", names(d)), with = FALSE])
r2 <- function(yv, k) { ok <- is.finite(yv); X <- cbind(1, pcmat[ok, seq_len(k), drop = FALSE]); f <- lm.fit(X, yv[ok]); 1 - sum(f$residuals^2)/sum((yv[ok]-mean(yv[ok]))^2) }
ks <- c(4, 10, 20, 32, 42, 52, 100, 200)
cat("\nR^2 of phenotype / covariate on the leading coverage PCs (the ones correct.R removes):\n")
cat(sprintf("%-22s %s\n", "k =", paste(sprintf("%6d", ks), collapse = "")))
sexn <- as.numeric(d$INFERRED_SEX == "M")
sp <- model.matrix(~ SUPERPOPULATION, d)[, -1]
pop <- model.matrix(~ POPULATION, d)[, -1]
rel <- as.numeric(d$RELATEDNESS != "unrelated")
for (nm in c("MTDNA_CN", "log2(MTDNA_CN)", "SEX", "related(0/1)")) {
  v <- switch(nm, "MTDNA_CN" = y, "log2(MTDNA_CN)" = ly, "SEX" = sexn, "related(0/1)" = rel)
  cat(sprintf("%-22s %s\n", nm, paste(sprintf("%6.3f", sapply(ks, function(k) r2(v, k))), collapse = "")))
}
# superpopulation: R^2 of each of the 4 dummies, take the mean and max
for (k in c(10, 42, 100)) {
  rr <- apply(sp, 2, r2, k = k)
  cat(sprintf("superpop dummies on PC1..%d: R^2 %s\n", k, paste(sprintf("%s=%.3f", sub("SUPERPOPULATION", "", names(rr)), rr), collapse = " ")))
}
# reverse: how much of each PC is explained by superpop / sex / batch / relatedness
cat("\nR^2 of each coverage PC on categorical covariates (PC ~ factor):\n")
rf <- function(pc, f) { s <- summary(lm(pc ~ f)); s$r.squared }
cat(sprintf("%-14s %s\n", "PC", paste(sprintf("%5d", 1:20), collapse = "")))
for (nm in c("SUPERPOPULATION", "POPULATION", "INFERRED_SEX", "RELATEDNESS", "RELEASE_BATCH", "INSTRUMENT_MODEL")) {
  f <- factor(d[[nm]])
  cat(sprintf("%-14s %s\n", substr(nm, 1, 14), paste(sprintf("%5.2f", sapply(1:20, function(i) rf(pcmat[, i], f))), collapse = "")))
}
# relatives vs unrelated: is "related" a batch? RELEASE_BATCH x RELATEDNESS
cat("\nRELEASE_BATCH x RELATEDNESS:\n"); print(table(d$RELEASE_BATCH, d$RELATEDNESS))
cat("\nMTDNA_CN by superpop (median):\n"); print(d[, .(median = median(MTDNA_CN), n = .N), by = SUPERPOPULATION][order(SUPERPOPULATION)])
cat("\nMTDNA_CN by relatedness (median):\n"); print(d[, .(median = median(MTDNA_CN), n = .N), by = RELATEDNESS])
cat(sprintf("\nMTDNA_CN ~ sex: t = %.2f ; log2 ~ sex t = %.2f\n", coef(summary(lm(y ~ sexn)))[2, 3], coef(summary(lm(ly ~ sexn)))[2, 3]))
s <- summary(lm(y ~ factor(d$SUPERPOPULATION))); cat(sprintf("MTDNA_CN ~ superpop: R^2 %.3f F p = %.2g\n", s$r.squared, pf(s$fstatistic[1], s$fstatistic[2], s$fstatistic[3], lower.tail = FALSE)))
s <- summary(lm(y ~ factor(d$RELEASE_BATCH))); cat(sprintf("MTDNA_CN ~ release batch: R^2 %.3f\n", s$r.squared))

# --- spectrum shape --------------------------------------------------------
cat(sprintf("\nspectrum: sv1=%.1f sv2=%.1f sv10=%.1f sv42=%.2f sv52=%.2f sv60=%.2f sv100=%.2f sv200=%.2f\n", sv[1], sv[2], sv[10], sv[42], sv[52], sv[60], sv[100], sv[200]))
cat(sprintf("count of sv within 1%%/2%%/5%%/10%% above sv[200]*(sqrt(p)+sqrt(n))/s_unit200 not needed; simple gaps: sv[k]/sv[k+1]-1 for k=40..55: %s\n",
  paste(sprintf("%.4f", sv[40:55]/sv[41:56] - 1), collapse = " ")))
fwrite(d[, .(SAMPLE = SAMPLE_ID, MTDNA_CN, SEX = sexn, SUPERPOPULATION, POPULATION, RELATEDNESS, FAMILY_ROLE, HQ_MEDIAN_COV, X_COV_RATIO, Y_COV_RATIO, MITO_COV_RATIO, MEAN_AUTOSOMAL_COV)],
       "/private/tmp/claude-501/-Users-Kitty-git-depthSV/1db72424-c4a1-45f4-b81e-5b466f3ec8a4/scratchpad/review_stats/real_pheno.tsv", sep = "\t")
