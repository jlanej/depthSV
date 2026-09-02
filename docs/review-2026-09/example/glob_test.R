#!/usr/bin/env Rscript
# Does "<name>.<method>.*.txt.gz" match only its own shards, with the real layout?
S <- "/private/tmp/claude-501/-Users-Kitty-git-depthSV/1db72424-c4a1-45f4-b81e-5b466f3ec8a4/scratchpad/review_example"
src <- file.path(S, "smoke/work/standard/assoc_ndim4")
d <- file.path(S, "glob_test"); unlink(d, recursive = TRUE); dir.create(d)
file.copy(list.files(src, full.names = TRUE), d)
# add the names a covariate-adjusted manifest would produce, plus leftovers a failed unit leaves
for (nm in c("mtdna_cn_adj", "log2_mtdna_cn_adj", "mtdna_cn_null_adj")) {
  for (slug in c("chr20_1-1000000", "chrM_1-16569")) {
    file.create(file.path(d, sprintf("%s.linear.%s.txt.gz", nm, slug)))
    file.create(file.path(d, sprintf("%s.linear.%s.txt.gz.tbi", nm, slug)))
    file.create(file.path(d, sprintf("%s.linear.%s.txt.gz.done", nm, slug)))
    file.create(file.path(d, sprintf("%s.linear.%s.log", nm, slug)))
  }
}
file.create(file.path(d, "mtdna_cn.linear.chrX_1-1000000.txt.gz.tmp.gz"))
cat("files in dir:", length(list.files(d)), "\n")
for (nm in c("mtdna_cn", "mtdna_cn_adj", "mtdna_cn_null", "mtdna_cn_null_adj", "log2_mtdna_cn", "log2_mtdna_cn_adj", "sex_linear")) {
  f <- Sys.glob(file.path(d, sprintf("%s.%s.*.txt.gz", nm, "linear")))
  f2 <- f[!grepl("\\.tmp\\.gz$", f)]
  cat(sprintf("%-20s glob=%d after tmp-filter=%d :: %s\n", nm, length(f), length(f2),
              paste(basename(f2), collapse = " ")))
}
f <- Sys.glob(file.path(d, "inferred_sex.logistic.*.txt.gz")); cat("inferred_sex:", basename(f), "\n")
