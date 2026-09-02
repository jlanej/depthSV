#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# depthSV 1000G example — assemble genotype PCs into a covariate table
#
# plink2 has done the heavy lifting (LD-pruned SNPs, KING-unrelated set,
# PCA with allele weights, projection of every sample onto those weights).
# This step turns that into what the association models consume, and shows
# the evidence:
#
#   * calibrates the projected scores against the in-sample PCs of the
#     unrelated set (they are proportional; the constant is fitted per PC
#     and the fit's r^2 is the check that the projection is sound), then
#     applies it to everyone so relatives and unrelated samples are on one
#     scale;
#   * counts the PCs above the Marchenko-Pastur edge of the unrelated
#     samples x pruned SNPs matrix, the same test used for the coverage
#     PCs, as a suggestion beside the fixed default the models use;
#   * plots PC1/2 and PC3/4 by superpopulation (relatives as open points)
#     and the eigenvalue scree against the MP edge.
#
#   Rscript genotype_pcs.R --eigenvec unrel.eigenvec --eigenval unrel.eigenval \
#       --sscore all.sscore --psam hg38_corrected.psam --nsnp 123456 \
#       --unrelated unrel.king.cutoff.in.id --out DIR
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

option_list <- list(
  make_option("--eigenvec",  type = "character", help = "plink2 .eigenvec of the unrelated set"),
  make_option("--eigenval",  type = "character", help = "plink2 .eigenval of the unrelated set"),
  make_option("--sscore",    type = "character", help = "plink2 .sscore projecting every sample"),
  make_option("--psam",      type = "character", help = "sample table with SuperPop / Population columns"),
  make_option("--unrelated", type = "character", help = "plink2 .king.cutoff.in.id (the PCA sample set)"),
  make_option("--nsnp",      type = "integer",   help = "number of pruned SNPs the PCA used"),
  make_option("--margin",    type = "double",    default = 0.01),
  make_option("--out",       type = "character", help = "output directory")
)
opt <- parse_args(OptionParser(option_list = option_list))
for (req in c("eigenvec", "eigenval", "sscore", "psam", "unrelated", "nsnp", "out")) {
  if (is.null(opt[[req]])) stop("--", req, " is required", call. = FALSE)
}
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)

read_plink <- function(f) { d <- fread(f); setnames(d, 1, sub("^#", "", names(d)[1])); d }
iid_col <- function(d) if ("IID" %in% names(d)) "IID" else names(d)[1]

ev  <- read_plink(opt$eigenvec)
val <- fread(opt$eigenval, header = FALSE)[[1]]
sc  <- read_plink(opt$sscore)
ps  <- read_plink(opt$psam)
unrel <- read_plink(opt$unrelated)

pc_cols <- grep("^PC[0-9]+$", names(ev), value = TRUE)
k <- length(pc_cols)
# plink2 names the per-PC score columns after the weight file's header
# (PC1_AVG ...) or SCORE<n>_AVG without one; only score columns carry _AVG.
score_cols <- grep("_AVG$", names(sc), value = TRUE)
if (length(score_cols) < k) stop("sscore has ", length(score_cols), " score columns for ", k, " PCs", call. = FALSE)
score_cols <- score_cols[seq_len(k)]

ev[, SAMPLE := as.character(get(iid_col(ev)))]
sc[, SAMPLE := as.character(get(iid_col(sc)))]
ps[, SAMPLE := as.character(get(iid_col(ps)))]
unrel_ids <- as.character(unrel[[iid_col(unrel)]])

# --- calibrate the projection ----------------------------------------------
# For the unrelated samples the projected score is proportional to their own
# PC; fit the constant per PC through the origin and report r^2.

both <- merge(ev[, c("SAMPLE", pc_cols), with = FALSE],
              sc[, c("SAMPLE", score_cols), with = FALSE], by = "SAMPLE")
scale <- numeric(k); r2 <- numeric(k)
for (i in seq_len(k)) {
  x <- both[[score_cols[i]]]; y <- both[[pc_cols[i]]]
  scale[i] <- sum(x * y) / sum(x * x)
  r2[i] <- suppressWarnings(cor(x, y))^2
}
calib <- data.table(PC = pc_cols, scale = scale, r2 = r2, eigenvalue = val[seq_len(k)])
fwrite(calib, file.path(opt$out, "gpc_calibration.tsv"), sep = "\t")
if (min(r2) < 0.99) {
  message(sprintf("[gpc] WARN: projection reproduces in-sample PCs with r^2 down to %.3f (PC%d)",
                  min(r2), which.min(r2)))
}

proj <- sc[, c("SAMPLE", score_cols), with = FALSE]
for (i in seq_len(k)) proj[[score_cols[i]]] <- proj[[score_cols[i]]] * scale[i]
setnames(proj, score_cols, sub("^PC", "GPC", pc_cols))
proj[, GPC_PROJECTED := as.integer(!(SAMPLE %in% unrel_ids))]

meta_cols <- intersect(c("SuperPop", "Population", "SEX"), names(ps))
cov <- merge(proj, ps[, c("SAMPLE", meta_cols), with = FALSE], by = "SAMPLE", all.x = TRUE)
if ("SEX" %in% names(cov)) setnames(cov, "SEX", "PSAM_SEX")
if ("SuperPop" %in% names(cov)) setnames(cov, "SuperPop", "GPC_SUPERPOP")
if ("Population" %in% names(cov)) setnames(cov, "Population", "GPC_POP")
setcolorder(cov, c("SAMPLE", grep("^GPC[0-9]+$", names(cov), value = TRUE), "GPC_PROJECTED"))
fwrite(cov, file.path(opt$out, "covariates.tsv"), sep = "\t")

# --- how many genotype PCs carry structure (MP edge) -----------------------
# plink2's eigenvalues are those of the variance-standardised relationship
# matrix over m SNPs; the noise bulk for n unrelated samples has edge
# (1 + sqrt(n/m))^2 in units of the noise variance, fitted from the trailing
# reported eigenvalues exactly as for the coverage spectrum.

n_unrel <- length(unrel_ids); m <- opt$nsnp
gamma <- n_unrel / m
k_mp <- NA_integer_; s2 <- NA_real_; edge <- NA_real_; status <- "undetermined"
q_unit <- rep(NA_real_, k)
if (gamma < 1) {
  a <- (1 - sqrt(gamma))^2; b <- (1 + sqrt(gamma))^2
  x <- seq(a, b, length.out = 20001)
  dens <- sqrt(pmax(0, (b - x) * (x - a))) / (2 * pi * gamma * x)
  cdf <- cumsum(c(0, (dens[-1] + dens[-length(dens)]) / 2 * diff(x))); cdf <- cdf / max(cdf)
  q_unit <- approx(cdf, x, xout = 1 - (seq_len(k) - 0.5) / n_unrel, ties = "ordered")$y
  gap <- 5L
  kk <- 0L
  for (it in 1:100) {
    ranks <- (kk + gap + 1):k
    if (length(ranks) < 5) break
    s2 <- exp(mean(log(val[ranks] / q_unit[ranks])))
    edge <- s2 * b * (1 + opt$margin)
    knew <- sum(val[seq_len(k)] > edge)
    if (knew == kk) { k_mp <- kk; status <- "ok"; break }
    kk <- knew
  }
} else {
  # Fewer SNPs than samples (a one-chromosome smoke run): the relationship
  # matrix is rank-deficient and the MP bulk takes a different form; the
  # count is not attempted rather than reported wrong.
  status <- "not applicable (SNPs < samples)"
}
message(sprintf("[gpc] %d PCs computed on %d unrelated samples x %d SNPs; MP-edge count = %s (%s)",
                k, n_unrel, m, k_mp, status))

# --- plots -----------------------------------------------------------------

pal <- c(AFR = "#E69F00", AMR = "#56B4E9", EAS = "#009E73", EUR = "#D55E00", SAS = "#CC79A7")
grp <- if ("GPC_SUPERPOP" %in% names(cov)) cov$GPC_SUPERPOP else rep("all", nrow(cov))
col <- ifelse(grp %in% names(pal), pal[grp], "grey50")
pch <- ifelse(cov$GPC_PROJECTED == 1L, 1, 16)

png(file.path(opt$out, "gpc_plots.png"), width = 2100, height = 700, res = 130)
par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3, 1))
plot(cov$GPC1, cov$GPC2, col = col, pch = pch, cex = 0.7, xlab = "GPC1", ylab = "GPC2",
     main = sprintf("genotype PCs (%d unrelated + %d projected)", n_unrel, sum(cov$GPC_PROJECTED)))
legend("topright", c(names(pal), "projected relative"), col = c(pal, "black"),
       pch = c(rep(16, 5), 1), bty = "n", cex = 0.9)
plot(cov$GPC3, cov$GPC4, col = col, pch = pch, cex = 0.7, xlab = "GPC3", ylab = "GPC4", main = "GPC3 vs GPC4")
plot(seq_len(k), val[seq_len(k)], type = "b", pch = 16, cex = 0.7, log = "y",
     xlab = "PC", ylab = "eigenvalue (log)", main = sprintf("scree; MP-edge count = %s", k_mp))
if (is.finite(edge)) abline(h = edge, col = "firebrick", lty = 2)
if (is.finite(s2)) lines(seq_len(k), s2 * q_unit, col = "steelblue", lwd = 2)
dev.off()

writeLines(c("# Genotype PCs", "",
             sprintf("- PCA on %d KING-unrelated samples x %d LD-pruned SNPs; %d PCs computed", n_unrel, m, k),
             sprintf("- %d related samples projected; calibration r^2 min %.4f (PC%d)",
                     sum(cov$GPC_PROJECTED), min(r2), which.min(r2)),
             sprintf("- PCs above the Marchenko-Pastur edge (+%.0f%%): %s (%s)", 100 * opt$margin, k_mp, status),
             "", "Files: covariates.tsv (SAMPLE, GPC1..GPCk, GPC_PROJECTED, superpopulation/population),",
             "gpc_calibration.tsv, gpc_plots.png."),
           file.path(opt$out, "summary.md"))
cat(sprintf("[gpc] covariates for %d samples -> %s\n", nrow(cov), file.path(opt$out, "covariates.tsv")))
