#!/bin/bash
# X13: the whole example with EX_WORK_DIR containing a space (fresh work dir; GitHub cache reused)
S=/private/tmp/claude-501/-Users-Kitty-git-depthSV/1db72424-c4a1-45f4-b81e-5b466f3ec8a4/scratchpad/review_shell
REPO=/Users/Kitty/git/depthSV/.claude/worktrees/agent-a08479f00845fba52
export EX_WORK_DIR="$S/ex space"
rm -rf "$EX_WORK_DIR"; mkdir -p "$EX_WORK_DIR"; cp -a "$S/exwork/github_cache" "$EX_WORK_DIR/"
cd "$REPO/example/1000G_highcov" || exit 1
start=$(date +%s)
/bin/bash run.sh --smoke --runner local > "$S/x13.log" 2>&1; echo "run.sh (spaced EX_WORK_DIR, bash 3.2) exit=$? in $(( $(date +%s) - start ))s"
grep -n 'FAILED\|ERROR\|verdict' "$S/x13.log" | head
ls "$EX_WORK_DIR/compare/summary.md" "$EX_WORK_DIR/profile/profile_report.md" 2>&1 | sed "s|$S/||"
