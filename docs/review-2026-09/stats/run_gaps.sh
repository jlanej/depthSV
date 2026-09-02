#!/usr/bin/env bash
S=/private/tmp/claude-501/-Users-Kitty-git-depthSV/1db72424-c4a1-45f4-b81e-5b466f3ec8a4/scratchpad/review_stats
D=/Users/Kitty/git/NGS-PCA/example/1000G_highcov/output/ngspca_output
R=/Users/Kitty/git/depthSV/.claude/worktrees/agent-a6ee8885847baece0/example/1000G_highcov/R/choose_ndim.R
for g in 5 10 20 40 80 120; do
  echo "== gap $g"
  Rscript "$R" --runs standard=$D --margin 0.01 --gap $g --out "$S/ndim_gap$g" 2>&1 | grep '\[ndim\]'
done
echo
cat "$S/ndim_gap20/summary.md"
