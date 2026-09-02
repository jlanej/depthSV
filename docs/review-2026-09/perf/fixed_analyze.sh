#!/usr/bin/env bash
# Fixed cost of one R/analyze.R process: header + ONE corrected row, exactly
# the per-phenotype probe in scripts/analyze.sh and the floor per chunk.
cd "$(dirname "$0")"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
ROOT=/Users/Kitty/git/depthSV/.claude/worktrees/agent-a265e6aa144253f7b
REAL=/private/tmp/claude-501/-Users-Kitty-git-depthSV/1db72424-c4a1-45f4-b81e-5b466f3ec8a4/scratchpad/runwork/inputs/standard
run() {  # run <label> <pheno> <corrfile> <model> <method> <nrows>
  local label="$1" ph="$2" corr="$3" model="$4" meth="$5" nrows="$6"
  head -n $((nrows + 1)) "$corr" > in.tmp
  /usr/bin/time -l -o time.tmp Rscript "$ROOT/R/analyze.R" -f in.tmp -p "$ph" -m "$model" -r "$meth" > out.tmp 2> err.tmp
  local wall user sys rss
  wall=$(awk '/real/{print $1}' time.tmp); user=$(awk '/user/{print $1}' time.tmp); sys=$(awk '/sys/{print $1}' time.tmp)
  rss=$(awk '/maximum resident/{print $1}' time.tmp)
  printf '%-24s %-9s rows=%-3s wall %6.2fs  cpu %6.2fs  maxRSS %6.0f MB  outrows %s %s\n' "$label" "$meth" "$nrows" "$wall" "$(echo "$user + $sys" | bc)" "$(echo "$rss / 1048576" | bc)" "$(wc -l < out.tmp | tr -d ' ')" "$(grep -o 'ERROR.*\|Error.*' err.tmp | head -1)"
}
run "real 3202" "$REAL/phenotypes.tsv" tables/corr_3202.txt "MTDNA_CN~cov_resids+SEX+GPC1+GPC2+GPC3+GPC4+GPC5+GPC6+GPC7+GPC8+GPC9+GPC10" linear 1
run "real 3202" "$REAL/phenotypes.tsv" tables/corr_3202.txt "SEX~cov_resids" logistic 1
for N in 3202 50000 200000 500000; do
  run "synthetic $N" tables/phenotypes_$N.tsv tables/corr_$N.txt "y~cov_resids+age+sex+gpc1+gpc2+gpc3+gpc4+gpc5+gpc6+gpc7+gpc8+gpc9+gpc10" linear 1
  run "synthetic $N" tables/phenotypes_$N.tsv tables/corr_$N.txt "ybin~cov_resids+age+sex+gpc1+gpc2+gpc3+gpc4+gpc5+gpc6+gpc7+gpc8+gpc9+gpc10" logistic 1
done
run "synthetic 500000" tables/phenotypes_500000.tsv tables/corr_500000.txt "y~cov_resids+age+sex+gpc1+gpc2+gpc3+gpc4+gpc5+gpc6+gpc7+gpc8+gpc9+gpc10" linear 11
run "synthetic 500000" tables/phenotypes_500000.tsv tables/corr_500000.txt "ybin~cov_resids+age+sex+gpc1+gpc2+gpc3+gpc4+gpc5+gpc6+gpc7+gpc8+gpc9+gpc10" logistic 11
rm -f in.tmp out.tmp err.tmp time.tmp
