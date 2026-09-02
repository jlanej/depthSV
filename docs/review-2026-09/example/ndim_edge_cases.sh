#!/usr/bin/env bash
S=/private/tmp/claude-501/-Users-Kitty-git-depthSV/1db72424-c4a1-45f4-b81e-5b466f3ec8a4/scratchpad/review_example
EX=/Users/Kitty/git/depthSV/.claude/worktrees/agent-aa17d545ebe28b8cc/example/1000G_highcov
SP="$S/spectrum"
mk() { # mk <dir> <n_values> [transform]
  local d="$1" n="$2"; mkdir -p "$d"
  cp "$SP/svd.samples.txt" "$SP/svd.bins.txt" "$d/"
  head -n $((n + 1)) "$SP/svd.singularvalues.txt" > "$d/svd.singularvalues.txt"
}
echo "=== (a) 15 singular values (k_reported < gap+1)"
mk "$S/ndim_a" 15
Rscript "$EX/R/choose_ndim.R" --runs "std=$S/ndim_a" --out "$S/ndim_a/out" 2>&1 | tail -4; echo "exit=${PIPESTATUS[0]}"
echo "=== (b) 25 singular values (gap+10 > available fit ranks)"
mk "$S/ndim_b" 25
Rscript "$EX/R/choose_ndim.R" --runs "std=$S/ndim_b" --out "$S/ndim_b/out" 2>&1 | tail -3; echo "exit=${PIPESTATUS[0]}"; ls "$S/ndim_b/out"
echo "=== (c) all 200 values scaled 3x above the noise edge except the fit is on the same ranks (spectrum never reaches the bulk)"
mkdir -p "$S/ndim_c"; cp "$SP/svd.samples.txt" "$SP/svd.bins.txt" "$S/ndim_c/"
awk -F'\t' 'NR==1{print;next}{printf "%s\t%.6f\n",$1,200-0.5*$1}' "$SP/svd.singularvalues.txt" > "$S/ndim_c/svd.singularvalues.txt"
Rscript "$EX/R/choose_ndim.R" --runs "std=$S/ndim_c" --out "$S/ndim_c/out" 2>&1 | tail -3; echo "exit=${PIPESTATUS[0]}"; ls "$S/ndim_c/out"
echo "=== (d) columns swapped (SINGULAR_VALUES first, PC second): silently uses the rank index as the spectrum?"
mkdir -p "$S/ndim_d"; cp "$SP/svd.samples.txt" "$SP/svd.bins.txt" "$S/ndim_d/"
awk -F'\t' 'BEGIN{OFS="\t"}{print $2,$1}' "$SP/svd.singularvalues.txt" > "$S/ndim_d/svd.singularvalues.txt"
head -2 "$S/ndim_d/svd.singularvalues.txt"
Rscript "$EX/R/choose_ndim.R" --runs "std=$S/ndim_d" --out "$S/ndim_d/out" 2>&1 | tail -3; echo "exit=${PIPESTATUS[0]}"; cat "$S/ndim_d/out/ndim.txt" 2>/dev/null
echo "=== (e) two modes where only one is determined: what does ndim.txt say and does the summary say 'averaged'?"
Rscript "$EX/R/choose_ndim.R" --runs "standard=$SP,fast=$S/ndim_c" --out "$S/ndim_e" 2>&1 | tail -4; cat "$S/ndim_e/ndim.txt"; grep -E 'ndim =|No mode' "$S/ndim_e/summary.md"
echo "=== (f) two modes, both determined but different (fast = standard spectrum x 1.03 on ranks 30-60 only)"
mkdir -p "$S/ndim_f"; cp "$SP/svd.samples.txt" "$SP/svd.bins.txt" "$S/ndim_f/"
awk -F'\t' 'BEGIN{OFS="\t"}NR==1{print;next}{v=$2; if($1>=30&&$1<=60) v=v*1.03; print $1,v}' "$SP/svd.singularvalues.txt" > "$S/ndim_f/svd.singularvalues.txt"
Rscript "$EX/R/choose_ndim.R" --runs "standard=$SP,fast=$S/ndim_f" --out "$S/ndim_f_out" 2>&1 | grep '\[ndim\]'; cat "$S/ndim_f_out/ndim.txt"
