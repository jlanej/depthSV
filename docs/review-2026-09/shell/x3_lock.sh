#!/bin/bash
# X3: stale .finalize.lock
S=/private/tmp/claude-501/-Users-Kitty-git-depthSV/1db72424-c4a1-45f4-b81e-5b466f3ec8a4/scratchpad/review_shell
REPO=/Users/Kitty/git/depthSV/.claude/worktrees/agent-a08479f00845fba52
export EX_WORK_DIR="$S/exC"
cd "$REPO/example/1000G_highcov" || exit 1
mkdir -p "$EX_WORK_DIR/.finalize.lock"
rm -f "$EX_WORK_DIR/compare/summary.md" "$EX_WORK_DIR/profile/profile_report.md"
EX_SMOKE=1 /bin/bash 02_run_depthsv.sh --runner local --mode all > "$S/x3.log" 2>&1; echo "02 (resubmission with a stale lock) exit=$?"
grep -n 'finalize\|verdict' "$S/x3.log"
ls "$EX_WORK_DIR/compare/summary.md" "$EX_WORK_DIR/profile/profile_report.md" 2>&1
echo "--- code: maybe_finalize has no trap; lock removed only on the happy path:"; grep -n 'finalize.lock' 02_run_depthsv.sh
