#!/bin/bash
# X2: compare_modes.R FAIL verdict swallowed by 04_compare_modes.sh; also 03 without EX_SMOKE after a smoke run
S=/private/tmp/claude-501/-Users-Kitty-git-depthSV/1db72424-c4a1-45f4-b81e-5b466f3ec8a4/scratchpad/review_shell
REPO=/Users/Kitty/git/depthSV/.claude/worktrees/agent-a08479f00845fba52
export EX_WORK_DIR="$S/exB"
cd "$REPO/example/1000G_highcov" || exit 1
rm -f "$EX_WORK_DIR"/work/fast/assoc_ndim4/inferred_sex.logistic.*
cp "$EX_WORK_DIR/compare/summary.md" "$S/x2_summary_before.md"
EX_SMOKE=1 /bin/bash 04_compare_modes.sh > "$S/x2_04.log" 2>&1; echo "04_compare_modes.sh exit=$?"
cat "$S/x2_04.log"
echo "--- compare/summary.md after:"; cat "$EX_WORK_DIR/compare/summary.md"
echo "--- standard_vs_fast/concordance.tsv status column:"; cut -f1,2 "$EX_WORK_DIR/compare/standard_vs_fast/concordance.tsv"
echo "--- standard_vs_fast/summary.md verdict line:"; grep -m1 verdict "$EX_WORK_DIR/compare/standard_vs_fast/summary.md"
echo
echo "=== X10: stage scripts rerun by hand after a smoke run, without EX_SMOKE=1 ==="
/bin/bash 03_evaluate.sh > "$S/x10.log" 2>&1; echo "03_evaluate.sh exit=$?"; cat "$S/x10.log"
