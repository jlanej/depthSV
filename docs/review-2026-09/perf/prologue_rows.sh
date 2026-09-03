#!/usr/bin/env bash
# Step-by-step prologue breakdown and per-row costs at each width.
cd "$(dirname "$0")"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
for N in 3202 50000 200000 500000; do
  echo "=== prologue N=$N ndim=40"
  Rscript fixed_cost.R tables/svd.pcs_$N.txt tables/autosomal.median_$N.txt 40 tables/raw_$N.txt
done
for N in 3202 50000 200000 500000; do
  echo "=== per-row N=$N ndim=40"
  Rscript row_cost.R tables/raw_$N.txt 40
done
