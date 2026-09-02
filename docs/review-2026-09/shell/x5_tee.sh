#!/bin/bash
# X5: does a failing prepare_inputs.R stop 01 despite the `| tee`?
S=/private/tmp/claude-501/-Users-Kitty-git-depthSV/1db72424-c4a1-45f4-b81e-5b466f3ec8a4/scratchpad/review_shell
REPO=/Users/Kitty/git/depthSV/.claude/worktrees/agent-a08479f00845fba52
export EX_WORK_DIR="$S/exE"
cd "$REPO/example/1000G_highcov" || exit 1
sed '1s/MTDNA_CN/MTDNA_XX/' "$EX_WORK_DIR/github_cache/standard/sample_qc.tsv" > "$S/x5_badqc.tsv"
sed -i.bak "s|^EX_M_QC_TABLE=.*|EX_M_QC_TABLE=$S/x5_badqc.tsv|" "$EX_WORK_DIR/inputs/standard/paths.env"
EX_SMOKE=1 EX_MODES=standard /bin/bash 01_prepare_inputs.sh > "$S/x5.log" 2>&1; echo "01 exit=$?"
tail -4 "$S/x5.log"
echo "--- stale tables left from the earlier run (01 does not clear them):"; ls -la "$EX_WORK_DIR/inputs/standard/" | awk '{print $5, $9}' | grep -v '^$'
