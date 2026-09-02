#!/bin/bash
for c in bgzip tabix parallel Rscript gtime shellcheck plink2 curl sbatch shasum; do
  printf '%s: ' "$c"; command -v "$c" || echo MISSING
done
echo "--- GitHub availability (HTTP codes) ---"
B=https://raw.githubusercontent.com/jlanej/NGS-PCA/master/example/1000G_highcov
for p in output/ngspca_output/svd.pcs.txt output_fast/ngspca_output/svd.pcs.txt \
         output/qc_output/sample_qc.tsv output_fast/qc_output/sample_qc.tsv \
         output/ngspca_output/svd.singularvalues.txt output_fast/ngspca_output/svd.singularvalues.txt \
         output_fast/ngspca_output/svd.samples.txt output_fast/ngspca_output/svd.bins.txt \
         output/ngspca_output/autosomal.median.txt; do
  printf '%s: ' "$p"
  curl -s -o /dev/null -w '%{http_code}\n' --max-time 20 -r 0-0 "$B/$p" || echo curl-failed
done
