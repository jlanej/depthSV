#!/usr/bin/env bash
# Real (non-smoke) run with a local mosdepth tree but no local NGS-PCA outputs: is the
# GitHub fallback (local depths + committed PCs/medians) silent?
S=/private/tmp/claude-501/-Users-Kitty-git-depthSV/1db72424-c4a1-45f4-b81e-5b466f3ec8a4/scratchpad/review_example
EX=/Users/Kitty/git/depthSV/.claude/worktrees/agent-aa17d545ebe28b8cc/example/1000G_highcov
cd "$EX" || exit 9
export EX_WORK_DIR="$S/realfallback" EX_MODES=standard
export EX_MOSDEPTH_DIR_STANDARD="$S/smoke/smoke_mosdepth/standard"
export EX_NGSPCA_DIR_STANDARD=/nonexistent/ngspca_output EX_QC_DIR_STANDARD=/nonexistent/qc_output
export EX_CACHE_DIR="$S/smoke/github_cache"   # reuse the already-fetched tables (no network)
rm -rf "$EX_WORK_DIR"
echo "=== 00"; bash 00_fetch_inputs.sh 2>&1 | sed 's/^\[[^]]*\] //'
echo "=== paths.env"; cat "$EX_WORK_DIR/inputs/standard/paths.env"
echo "=== 01"; bash 01_prepare_inputs.sh 2>&1 | sed 's/^\[[^]]*\] //' | grep -vE '^\[(qc|pcs|coverage|phenotypes|align)\]'
echo "=== any WARN about the PCs/medians not coming from the same run as the depths?"
grep -ci 'warn' "$EX_WORK_DIR/inputs/standard/prepare.summary.txt"
echo "=== run.sh --prepare-only summary line:"
bash run.sh --prepare-only 2>&1 | tail -2
