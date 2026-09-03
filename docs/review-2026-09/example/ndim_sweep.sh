#!/usr/bin/env bash
# README Notes claim: "the chrM_log2_slope shrinks accordingly [with the PCs' absorption of the phenotype]".
# Sweep ndim on the smoke chrM unit and print the slope and t at the best chrM bin.
S=/private/tmp/claude-501/-Users-Kitty-git-depthSV/1db72424-c4a1-45f4-b81e-5b466f3ec8a4/scratchpad/review_example
R=/Users/Kitty/git/depthSV/.claude/worktrees/agent-aa17d545ebe28b8cc
W="$S/smoke"
export DSV_MATRIX="$W/work/standard/join/depth.matrix.txt.gz"
export DSV_PCS="$W/inputs/standard/svd.pcs.txt"
export DSV_COVERAGE="$W/inputs/standard/autosomal.median.txt"
export DSV_JOBS=2 DSV_THREADS=2 DSV_MIN_OBS=30
printf 'ndim\tbest_chrM_bin\tslope\tt\tp\tnull_lambda_auto\n'
for nd in 0 2 4 10 20 40; do
  O="$S/ndim_sweep/nd$nd"; rm -rf "$O"; mkdir -p "$O"
  export DSV_NDIM=$nd
  for reg in chrM:1-16569 chr20:1-1000000; do
    slug=$(printf '%s' "$reg" | tr ':' '_')
    bash "$R/scripts/correct.sh" --region "$reg" --out "$O" >/dev/null 2>&1 || { echo "correct failed nd=$nd $reg"; continue; }
    bash "$R/scripts/analyze.sh" --corrected "$O/corrected_ndim$nd.$slug.txt.gz" --pheno "$W/inputs/standard/phenotypes.tsv" \
         --region "$reg" --out "$O" --pheno-manifest "$W/inputs/standard/analyses.tsv" >/dev/null 2>&1 || { echo "analyze failed nd=$nd $reg"; continue; }
  done
  best=$(gzip -cd "$O/log2_mtdna_cn.linear.chrM_1-16569.txt.gz" | awk -F'\t' 'NR>1{a=$10<0?-$10:$10; if(a>m){m=a; line=$4"\t"$8"\t"$10"\t"$11}} END{print line}')
  lam=$(gzip -cd "$O/mtdna_cn_null.linear.chr20_1-1000000.txt.gz" | awk -F'\t' 'NR>1{print $11}' | Rscript -e 'p<-scan("stdin",quiet=TRUE); cat(sprintf("%.3f", median(qchisq(p,1,lower.tail=FALSE))/qchisq(0.5,1)))')
  printf '%s\t%s\t%s\n' "$nd" "$best" "$lam"
done
