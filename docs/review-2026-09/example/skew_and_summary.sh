#!/usr/bin/env bash
S=/private/tmp/claude-501/-Users-Kitty-git-depthSV/1db72424-c4a1-45f4-b81e-5b466f3ec8a4/scratchpad/review_example
echo "=== 04's summary when a concordance FAILs (stale tree):"
cat "$S/smoke_stale/compare/summary.md"
echo; echo "=== standard_vs_fast summary status column in that run:"
cut -f1,2,5 "$S/smoke_stale/compare/standard_vs_fast/concordance.tsv"
echo; echo "=== MTDNA_CN distribution in the committed QC table (skew matters for REVIEW 1.2's minDepth artefact)"
Rscript -e '
suppressPackageStartupMessages(library(data.table))
q <- fread("'"$S"'/spectrum/std.sample_qc.tsv")
x <- q$MTDNA_CN; x <- x[is.finite(x)]
sk <- function(v) mean((v-mean(v))^3)/sd(v)^3
cat(sprintf("n=%d median=%.0f mean=%.0f sd=%.0f min=%.0f max=%.0f skewness=%.2f  (log2: skewness=%.2f)\n", length(x), median(x), mean(x), sd(x), min(x), max(x), sk(x), sk(log2(x))))
cat(sprintf("HQ_MEDIAN_COV: min=%.1f median=%.1f max=%.1f; samples with HQ_MEDIAN_COV NA: %d\n", min(q$HQ_MEDIAN_COV,na.rm=TRUE), median(q$HQ_MEDIAN_COV,na.rm=TRUE), max(q$HQ_MEDIAN_COV,na.rm=TRUE), sum(is.na(q$HQ_MEDIAN_COV))))
cat("INFERRED_SEX values:", paste(names(table(q$INFERRED_SEX, useNA="ifany")), table(q$INFERRED_SEX, useNA="ifany"), collapse=", "), "\n")
p <- fread("'"$S"'/spectrum/svd.pcs.txt")
cat("PC table samples:", nrow(p), " QC samples:", nrow(q), " PCs in QC but not PC table:", length(setdiff(q$SAMPLE_ID, sub("\\.by1000\\.$","",p$SAMPLE))), " in PC table but not QC:", paste(setdiff(sub("\\.by1000\\.$","",p$SAMPLE), q$SAMPLE_ID), collapse=","), "\n")
cat("QC rows lacking HQ_MEDIAN_COV -> dropped from coverage -> PC/coverage overlap:", sum(sub("\\.by1000\\.$","",p$SAMPLE) %in% q$SAMPLE_ID[is.finite(q$HQ_MEDIAN_COV) & q$HQ_MEDIAN_COV>0]), "\n")
'
echo; echo "=== README claim: 'MTDNA_CN = 2 x chrM mean / HQ autosomal median' -- and log2 slope: does the PC correction shrink the slope (README Notes) or the t-stat?"
grep -h 'chrM_log2_slope\|chrM_top_hit' "$S/smoke/eval/standard/summary.md" | cut -d'|' -f2-5
