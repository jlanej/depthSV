#!/usr/bin/env bash
S=/private/tmp/claude-501/-Users-Kitty-git-depthSV/1db72424-c4a1-45f4-b81e-5b466f3ec8a4/scratchpad/review_example
EX=/Users/Kitty/git/depthSV/.claude/worktrees/agent-aa17d545ebe28b8cc/example/1000G_highcov
cd "$EX" || exit 9
export EX_WORK_DIR="$S/preamble_work"
bash preamble.sh --smoke --ndim-only > "$S/preamble.ndim.log" 2>&1
echo "exit=$?" >> "$S/preamble.ndim.log"
