#!/usr/bin/env bash
# Fixed cost of one R/correct.R process: header + ONE row, exactly the probe
# call in scripts/correct.sh and the floor every parallel chunk pays.
# Prints wall / user+sys / max RSS for each width and ndim.
cd "$(dirname "$0")"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
ROOT=/Users/Kitty/git/depthSV/.claude/worktrees/agent-a265e6aa144253f7b
REAL=/private/tmp/claude-501/-Users-Kitty-git-depthSV/1db72424-c4a1-45f4-b81e-5b466f3ec8a4/scratchpad/runwork/inputs/standard
run() {  # run <label> <pcs> <cov> <rawfile> <ndim> <nrows>
  local label="$1" pcs="$2" cov="$3" raw="$4" ndim="$5" nrows="$6"
  head -n $((nrows + 1)) "$raw" > in.tmp
  /usr/bin/time -l -o time.tmp Rscript "$ROOT/R/correct.R" --inputPCs "$pcs" --inputFile in.tmp --coverageStats "$cov" --ndim "$ndim" --skipOutputHeader > out.tmp 2> err.tmp
  local wall user sys rss
  wall=$(awk '/real/{print $1}' time.tmp); user=$(awk '/user/{print $1}' time.tmp); sys=$(awk '/sys/{print $1}' time.tmp)
  rss=$(awk '/maximum resident/{print $1}' time.tmp)
  printf '%-32s ndim=%-3s rows=%-3s wall %6.2fs  cpu %6.2fs  maxRSS %6.0f MB  outrows %s\n' "$label" "$ndim" "$nrows" "$wall" "$(echo "$user + $sys" | bc)" "$(echo "$rss / 1048576" | bc)" "$(wc -l < out.tmp | tr -d ' ')"
}
run "real 3202 (200-PC table)"   "$REAL/svd.pcs.txt" "$REAL/autosomal.median.txt" tables/raw_3202.txt 16 1
run "real 3202 (200-PC table)"   "$REAL/svd.pcs.txt" "$REAL/autosomal.median.txt" tables/raw_3202.txt 40 1
run "real 3202 (200-PC table)"   "$REAL/svd.pcs.txt" "$REAL/autosomal.median.txt" tables/raw_3202.txt 200 1
for N in 3202 50000 200000 500000; do
  for nd in 16 40; do
    run "synthetic $N (40-PC table)" tables/svd.pcs_$N.txt tables/autosomal.median_$N.txt tables/raw_$N.txt $nd 1
  done
done
# one full chunk's worth of rows at 500k (64 MB / row bytes) to see marginal cost
run "synthetic 500000 (40-PC table)" tables/svd.pcs_500000.txt tables/autosomal.median_500000.txt tables/raw_500000.txt 40 21
run "synthetic 200000 (40-PC table)" tables/svd.pcs_200000.txt tables/autosomal.median_200000.txt tables/raw_200000.txt 40 40
rm -f in.tmp out.tmp err.tmp time.tmp
