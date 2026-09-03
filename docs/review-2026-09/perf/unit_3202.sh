#!/usr/bin/env bash
# One SLURM array unit as example/1000G_highcov/02_run_depthsv.sh runs it:
# correct.sh then analyze.sh (6-phenotype manifest) over a 25 Mb window at
# 3,202 samples, --jobs 8 --threads 4, ndim 20, real PC/phenotype tables.
cd "$(dirname "$0")"
ROOT=/Users/Kitty/git/depthSV/.claude/worktrees/agent-a265e6aa144253f7b
REAL=/private/tmp/claude-501/-Users-Kitty-git-depthSV/1db72424-c4a1-45f4-b81e-5b466f3ec8a4/scratchpad/runwork/inputs/standard
M=g3202/join/depth.matrix.txt.gz
REGION=chr1:1-25000000
rm -rf g3202/corrected g3202/assoc
t() { local lab="$1"; shift; /usr/bin/time -l -o t.tmp "$@"; printf '%-40s wall %7.1fs  cpu %7.1fs  maxRSS %6.0f MB\n' "$lab" "$(awk '/real/{print $1}' t.tmp)" "$(awk '/user/{u=$1} /sys/{s=$1} END{print u+s}' t.tmp)" "$(awk '/maximum resident/{printf "%.0f", $1/1048576}' t.tmp)"; }
echo "matrix: $(ls -la $M | awk '{print $5}') bytes; region rows: $(tabix $M $REGION | wc -l | tr -d ' ')"
t "tabix region -> /dev/null (one read)" sh -c "tabix $M $REGION > /dev/null"
t "correct.sh ndim=20 jobs=8" bash "$ROOT/scripts/correct.sh" --matrix "$M" --pcs "$REAL/svd.pcs.txt" --coverage "$REAL/autosomal.median.txt" --region "$REGION" --out g3202/corrected --ndim 20 --jobs 8 --threads 4 2> g3202/correct.err
grep -E 'region|parity|wrote|stats' g3202/correct.err | tail -4
C=g3202/corrected/corrected_ndim20.chr1_1-25000000.txt.gz
ls -la g3202/corrected/
echo "corrected chunks (stats files merged): $(bgzip -dc g3202/corrected/stats_ndim20.chr1_1-25000000.txt.gz | wc -l | tr -d ' ') rows"
t "tabix corrected region -> /dev/null" sh -c "tabix $C $REGION > /dev/null"
t "analyze.sh 6-phenotype manifest jobs=8" bash "$ROOT/scripts/analyze.sh" --corrected "$C" --pheno "$REAL/phenotypes.tsv" --pheno-manifest "$REAL/analyses.tsv" --region "$REGION" --out g3202/assoc --jobs 8 --threads 4 --min-obs 100 2> g3202/analyze.err
grep -E 'wrote|complete' g3202/analyze.err
for f in g3202/assoc/*.txt.gz; do printf '%-60s %s rows\n' "$(basename $f)" "$(bgzip -dc $f | grep -vc '^#')"; done
# per-phenotype timing: run the manifest rows one at a time to see linear vs logistic
rm -rf g3202/assoc1
t "analyze.sh mtdna_cn_adj (linear, 12 cov)" bash "$ROOT/scripts/analyze.sh" --corrected "$C" --pheno "$REAL/phenotypes.tsv" --model "MTDNA_CN~cov_resids+SEX+GPC1+GPC2+GPC3+GPC4+GPC5+GPC6+GPC7+GPC8+GPC9+GPC10" --method linear --name a1 --region "$REGION" --out g3202/assoc1 --jobs 8 --threads 4 2> /dev/null
t "analyze.sh inferred_sex (logistic)" bash "$ROOT/scripts/analyze.sh" --corrected "$C" --pheno "$REAL/phenotypes.tsv" --model "SEX~cov_resids" --method logistic --name a2 --region "$REGION" --out g3202/assoc1 --jobs 8 --threads 4 2> /dev/null
t "analyze.sh linear jobs=1" bash "$ROOT/scripts/analyze.sh" --corrected "$C" --pheno "$REAL/phenotypes.tsv" --model "MTDNA_CN~cov_resids+SEX+GPC1+GPC2+GPC3+GPC4+GPC5+GPC6+GPC7+GPC8+GPC9+GPC10" --method linear --name a3 --region "$REGION" --out g3202/assoc1 --jobs 1 --threads 1 2> /dev/null
t "correct.sh ndim=20 jobs=1" bash "$ROOT/scripts/correct.sh" --matrix "$M" --pcs "$REAL/svd.pcs.txt" --coverage "$REAL/autosomal.median.txt" --region "$REGION" --out g3202/corrected1 --ndim 20 --jobs 1 --threads 1 2> /dev/null
rm -f t.tmp
