suppressPackageStartupMessages(library(data.table))
a <- commandArgs(TRUE); dir <- a[1]
rd <- function(pat) rbindlist(lapply(Sys.glob(file.path(dir, pat)), function(f) fread(cmd = paste("gzip -cd", shQuote(f)))))
lam <- function(p) { p <- p[is.finite(p) & p > 0]; median(qchisq(p, 1, lower.tail = FALSE)) / qchisq(0.5, 1) }
cls <- function(ch) { b <- toupper(sub("^chr", "", ch)); ifelse(b %in% c("M","MT"), "chrM", ifelse(b %in% c("X","Y"), "sex", "auto")) }
for (nm in c("mtdna_cn_null.linear.*.txt.gz", "mtdna_cn_null_adj.linear.*.txt.gz")) {
  d <- rd(nm); if (!nrow(d)) next
  setnames(d, 1, "CHROM"); setnames(d, ncol(d), "P"); d[, cl := cls(CHROM)]
  cat(sprintf("%-28s all=%.3f  auto=%.3f (n=%d)  sex=%.3f (n=%d)  chrM=%.3f\n", nm, lam(d$P),
              lam(d[cl=="auto"]$P), sum(d$cl=="auto"), lam(d[cl=="sex"]$P), sum(d$cl=="sex"), lam(d[cl=="chrM"]$P)))
}
