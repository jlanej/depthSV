#!/bin/bash
# Mirror the CI example job under bash 3.2: preamble --smoke --ndim-only, then run.sh --smoke --runner local.
S=/private/tmp/claude-501/-Users-Kitty-git-depthSV/1db72424-c4a1-45f4-b81e-5b466f3ec8a4/scratchpad/review_shell
REPO=/Users/Kitty/git/depthSV/.claude/worktrees/agent-a08479f00845fba52
export EX_WORK_DIR="$S/exwork"
cd "$REPO/example/1000G_highcov" || exit 1
start=$(date +%s)
/bin/bash preamble.sh --smoke --ndim-only > "$S/ex_preamble.log" 2>&1
echo "preamble exit=$? in $(( $(date +%s) - start ))s"
tail -3 "$S/ex_preamble.log"
start=$(date +%s)
/bin/bash run.sh --smoke --runner local > "$S/ex_run.log" 2>&1
echo "run.sh exit=$? in $(( $(date +%s) - start ))s"
tail -8 "$S/ex_run.log"
echo "--- repo tree litter check (untracked files in the checkout after the run) ---"
cd "$REPO" && git status --short --untracked-files=all | head
