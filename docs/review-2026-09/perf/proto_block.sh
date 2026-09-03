#!/usr/bin/env bash
cd "$(dirname "$0")"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
head -n 23 tables/raw_500000.txt > chunk_raw_500k.txt      # header + 22 rows = one 64 MB block
head -n 14 tables/corr_500000.txt > chunk_corr_500k.txt    # header + 13 rows
head -n 41 tables/raw_50000.txt > chunk_raw_50k.txt
head -n 41 tables/corr_50000.txt > chunk_corr_50k.txt
rm -f q_500k.rds q_50k.rds
echo "=== 500k"; Rscript proto_block.R chunk_raw_500k.txt chunk_corr_500k.txt tables/svd.pcs_500000.txt tables/autosomal.median_500000.txt 40 q_500k.rds
echo "=== 500k again (Q cached)"; Rscript proto_block.R chunk_raw_500k.txt chunk_corr_500k.txt tables/svd.pcs_500000.txt tables/autosomal.median_500000.txt 40 q_500k.rds
echo "=== 50k"; Rscript proto_block.R chunk_raw_50k.txt chunk_corr_50k.txt tables/svd.pcs_50000.txt tables/autosomal.median_50000.txt 40 q_50k.rds
echo "=== shipped scripts on the same chunks (wall, single process):"
ROOT=/Users/Kitty/git/depthSV/.claude/worktrees/agent-a265e6aa144253f7b
/usr/bin/time -p Rscript "$ROOT/R/correct.R" --inputPCs tables/svd.pcs_500000.txt --inputFile chunk_raw_500k.txt --coverageStats tables/autosomal.median_500000.txt --ndim 40 --skipOutputHeader 2>&1 >/dev/null | grep real | sed 's/^/correct.R 22 rows @500k: /'
/usr/bin/time -p Rscript "$ROOT/R/analyze.R" -f chunk_corr_500k.txt -p tables/phenotypes_500000.tsv -m "y~cov_resids+age+sex+gpc1+gpc2+gpc3+gpc4+gpc5+gpc6+gpc7+gpc8+gpc9+gpc10" -r linear 2>&1 >/dev/null | grep real | sed 's/^/analyze.R 13 rows @500k linear: /'
/usr/bin/time -p Rscript "$ROOT/R/correct.R" --inputPCs tables/svd.pcs_50000.txt --inputFile chunk_raw_50k.txt --coverageStats tables/autosomal.median_50000.txt --ndim 40 --skipOutputHeader 2>&1 >/dev/null | grep real | sed 's/^/correct.R 40 rows @50k: /'
/usr/bin/time -p Rscript "$ROOT/R/analyze.R" -f chunk_corr_50k.txt -p tables/phenotypes_50000.tsv -m "y~cov_resids+age+sex+gpc1+gpc2+gpc3+gpc4+gpc5+gpc6+gpc7+gpc8+gpc9+gpc10" -r linear 2>&1 >/dev/null | grep real | sed 's/^/analyze.R 40 rows @50k linear: /'
rm -f chunk_*.txt
