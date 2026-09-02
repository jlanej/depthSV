#!/bin/bash
# Run the core smoke test under bash 3.2 in a work dir containing a space.
S=/private/tmp/claude-501/-Users-Kitty-git-depthSV/1db72424-c4a1-45f4-b81e-5b466f3ec8a4/scratchpad/review_shell
REPO=/Users/Kitty/git/depthSV/.claude/worktrees/agent-a08479f00845fba52
cd "$REPO" || exit 1
echo "bash for suite: $(/bin/bash --version | head -1)"
start=$(date +%s)
/bin/bash tests/smoke_test.sh "$S/smoke dir" > "$S/smoke_space.log" 2>&1
echo "smoke_test.sh (space dir, bash 3.2) exit=$? in $(( $(date +%s) - start ))s"
tail -3 "$S/smoke_space.log"
