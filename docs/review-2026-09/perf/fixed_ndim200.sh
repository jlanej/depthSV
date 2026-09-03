#!/usr/bin/env bash
# Fixed cost of correct.R at ndim=200 (PLAN.md: "validated at scale") with a
# 200-column PC table, 1 row, at 200k and 500k samples.
cd "$(dirname "$0")"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
ROOT=/Users/Kitty/git/depthSV/.claude/worktrees/agent-a265e6aa144253f7b
for N in 200000 500000; do
  for nd in 40 200; do
    head -n 2 tables/raw_$N.txt > in.tmp
    /usr/bin/time -l -o time.tmp Rscript "$ROOT/R/correct.R" --inputPCs tables/svd.pcs200_$N.txt --inputFile in.tmp --coverageStats tables/autosomal.median_$N.txt --ndim $nd --skipOutputHeader > out.tmp 2> err.tmp
    printf 'N=%-7s 200-PC table, ndim=%-3s rows=1  wall %6.2fs  cpu %6.2fs  maxRSS %6.0f MB  outrows %s\n' "$N" "$nd" "$(awk '/real/{print $1}' time.tmp)" "$(awk '/user/{u=$1} /sys/{s=$1} END{print u+s}' time.tmp)" "$(awk '/maximum resident/{printf "%.0f", $1/1048576}' time.tmp)" "$(wc -l < out.tmp | tr -d ' ')"
  done
done
rm -f in.tmp out.tmp err.tmp time.tmp
