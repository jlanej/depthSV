# Size of the fixed artefact -c * M_Z P_PC(sex) that every chrX/chrY bin carries after
# PC removal, on the real 1000G PCs, and the t-statistic it induces against MTDNA_CN.
suppressPackageStartupMessages(library(data.table))
S <- "/private/tmp/claude-501/-Users-Kitty-git-depthSV/1db72424-c4a1-45f4-b81e-5b466f3ec8a4/scratchpad/review_stats"
d <- fread(file.path(S, "real_pheno.tsv"))
pcs <- fread("/Users/Kitty/git/NGS-PCA/example/1000G_highcov/output/ngspca_output/svd.pcs.txt")
pcs[, SAMPLE := sub("\\.by1000\\.$", "", SAMPLE)]
d <- merge(d, pcs, by = "SAMPLE"); n <- nrow(d)
for (K in c(10, 20, 42, 100)) {
  P <- as.matrix(d[, paste0("PC", 1:K), with = FALSE]); Qpc <- qr.Q(qr(cbind(1, P)))
  Z <- model.matrix(~ SEX + POPULATION, d); Qz <- qr.Q(qr(Z))
  Mz <- function(v) v - Qz %*% crossprod(Qz, v)
  s <- as.numeric(d$SEX)
  y <- d$MTDNA_CN; yt <- as.numeric(Mz(y))
  for (c in c(2.84, 1.0)) {
    art <- as.numeric(Mz(-c * (Qpc %*% crossprod(Qpc, s))))
    noise_ss <- if (c > 2) 1596 * 0.01 + 1603 * 0.16 else n * 0.01     # chrY: male sd .1, female sd .4 ; chrX: sd .1
    t_art <- sum(art * yt) / sqrt(sum(art^2) + noise_ss) / sqrt(sum(yt^2) / (n - ncol(Z) - 1))
    cat(sprintf("K=%3d  %s: ||artefact||^2 = %6.1f  vs noise SS ~%5.0f (artefact share %.0f%%) ; artefact-induced t vs MTDNA_CN = %5.2f  (identical on every bin)\n",
                K, if (c > 2) "chrY (c=2.84)" else "chrX (c=1.00)", sum(art^2), noise_ss, 100 * sum(art^2) / (sum(art^2) + noise_ss), t_art))
  }
}
cat("\nR^2 of SEX on PC1..K (the fraction of the sex vector the projection removes):",
    paste(sapply(c(10, 20, 42, 100), function(K) { P <- as.matrix(d[, paste0("PC", 1:K), with = FALSE]); sprintf("K=%d %.3f", K, summary(lm(d$SEX ~ P))$r.squared) }), collapse = "  "), "\n")
