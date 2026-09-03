#!/usr/bin/env bash
# (1) compare_modes.R with only 2 common regions in one analysis -> PASS?
# (2) calibration_summary.R when every ratio is undetermined -> headline verdict?
S=/private/tmp/claude-501/-Users-Kitty-git-depthSV/1db72424-c4a1-45f4-b81e-5b466f3ec8a4/scratchpad/review_example
EX=/Users/Kitty/git/depthSV/.claude/worktrees/agent-aa17d545ebe28b8cc/example/1000G_highcov
W="$S/smoke"
A="$W/work/standard/assoc_ndim4"
B="$S/tiny_b"; rm -rf "$B"; mkdir -p "$B"
cp "$A"/*.txt.gz "$B/"
# keep only 2 regions of mtdna_cn_null in mode b, and only 3 rows for sex_linear
for f in "$B"/mtdna_cn_null.linear.*.txt.gz; do gzip -cd "$f" | head -1 > "$f.tmp"; rm "$f"; done
gzip -cd "$A/mtdna_cn_null.linear.chr20_1-1000000.txt.gz" | head -3 | gzip > "$B/mtdna_cn_null.linear.chr20_1-1000000.txt.gz"
rm -f "$B"/*.tmp
out="$S/tiny_cmp"; rm -rf "$out"; mkdir -p "$out"
Rscript "$EX/R/compare_modes.R" --a-dir "$A" --b-dir "$B" --a-name standard --b-name fast \
  --analyses "$W/inputs/standard/analyses.tsv" --top-k 100 --profile real --out "$out" 2>&1 | tail -2
echo "--- concordance for mtdna_cn_null (n_common, status, r_stat):"
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i} $1=="mtdna_cn_null"{print "n_common="$h["n_common"], "status="$h["status"], "r_stat="$h["r_stat"], "topk="$h["topk_jaccard"]}' "$out/concordance.tsv"
grep verdict "$out/summary.md"
echo
echo "=== calibration_summary with a seed control whose r_stat is exactly 1 (or NA) for every analysis"
ctl="$S/tiny_ctl.tsv"
awk -F'\t' 'BEGIN{OFS="\t"} NR==1{for(i=1;i<=NF;i++)h[$i]=i; print; next} {$h["r_stat"]=1; print}' "$W/compare/standard_vs_fast/concordance.tsv" > "$ctl"
out2="$S/tiny_cal"; rm -rf "$out2"; mkdir -p "$out2"
Rscript "$EX/R/calibration_summary.R" --primary "$W/compare/standard_vs_fast/concordance.tsv" --control "$ctl" --factor 1.5 --out "$out2" 2>&1 | tail -1
grep -E 'Verdict|undetermined' "$out2/summary.md" | head -8
echo
echo "=== same with the control's r_stat = NA (a FAIL row: missing results)"
awk -F'\t' 'BEGIN{OFS="\t"} NR==1{for(i=1;i<=NF;i++)h[$i]=i; print; next} {$h["r_stat"]="NA"; $h["status"]="FAIL"; print}' "$W/compare/standard_vs_fast/concordance.tsv" > "$ctl"
Rscript "$EX/R/calibration_summary.R" --primary "$W/compare/standard_vs_fast/concordance.tsv" --control "$ctl" --factor 1.5 --out "$out2" 2>&1 | tail -1
grep -E 'Verdict' "$out2/summary.md"
