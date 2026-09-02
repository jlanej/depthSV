#!/usr/bin/env Rscript
# evaluate.R's read_shards with a header-only shard among normal ones (a unit whose regions were all QC-skipped)
suppressPackageStartupMessages(library(data.table))
S <- "/private/tmp/claude-501/-Users-Kitty-git-depthSV/1db72424-c4a1-45f4-b81e-5b466f3ec8a4/scratchpad/review_example"
src <- file.path(S, "smoke/work/standard/assoc_ndim4")
d <- file.path(S, "hdr_test"); unlink(d, recursive = TRUE); dir.create(d)
fs <- Sys.glob(file.path(src, "mtdna_cn.linear.*.txt.gz"))
file.copy(fs, d)
# make the chrX shard header-only
x <- file.path(d, basename(grep("chrX", fs, value = TRUE)[1]))
hdr <- system(paste("gzip -cd", shQuote(x), "| head -1"), intern = TRUE)
con <- gzfile(x, "w"); writeLines(hdr, con); close(con)
read_shards <- function(dir, name, method) {
  files <- Sys.glob(file.path(dir, sprintf("%s.%s.*.txt.gz", name, method)))
  files <- files[!grepl("\\.tmp\\.gz$", files)]
  if (!length(files)) return(list(files = character(0), dt = NULL))
  dt <- rbindlist(lapply(files, function(f) fread(cmd = paste("gzip -cd", shQuote(f)))))
  setnames(dt, 1:7, c("CHROM", "START", "END", "Region", "N", "NCase", "NControl"))
  stat_idx <- if (method == "coxph") 11L else 10L
  setnames(dt, ncol(dt), "P"); setnames(dt, 8L, "Estimate"); setnames(dt, stat_idx, "STAT")
  dt[, `:=`(P = as.numeric(P), Estimate = as.numeric(Estimate), STAT = as.numeric(STAT))]
  list(files = files, dt = dt)
}
r <- tryCatch(read_shards(d, "mtdna_cn", "linear"), error = function(e) e)
if (inherits(r, "error")) cat("ERROR:", conditionMessage(r), "\n") else cat("ok:", nrow(r$dt), "rows from", length(r$files), "files; classes:", paste(sapply(r$dt, class)[1:8], collapse=","), "\n")
# and a shard where every value is present but the header-only shard comes FIRST in glob order
y <- file.path(d, basename(grep("chr20_1-", fs, value = TRUE)[1]))
con <- gzfile(y, "w"); writeLines(hdr, con); close(con)
r <- tryCatch(read_shards(d, "mtdna_cn", "linear"), error = function(e) e)
if (inherits(r, "error")) cat("ERROR (header-only first):", conditionMessage(r), "\n") else cat("ok (header-only first):", nrow(r$dt), "rows; STAT class", class(r$dt$STAT), "; P NA count", sum(is.na(r$dt$P)), "\n")
