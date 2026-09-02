#!/bin/bash
# X7: preamble ndim.txt persistence when the MP fit is undetermined; and choose_ndim.R with a large gap
S=/private/tmp/claude-501/-Users-Kitty-git-depthSV/1db72424-c4a1-45f4-b81e-5b466f3ec8a4/scratchpad/review_shell
REPO=/Users/Kitty/git/depthSV/.claude/worktrees/agent-a08479f00845fba52
export EX_WORK_DIR="$S/exG"
cd "$REPO/example/1000G_highcov" || exit 1
echo "ndim.txt before: $(cat "$EX_WORK_DIR/preamble/ndim.txt")"
EX_MP_GAP=195 /bin/bash preamble.sh --smoke --ndim-only > "$S/x7.log" 2>&1; echo "preamble (gap=195) exit=$?"
grep -n 'ndim' "$S/x7.log" | tail -4
echo "ndim.txt after: $(cat "$EX_WORK_DIR/preamble/ndim.txt")   ndim/summary.md says:"; grep -n 'No mode could be determined\|ndim =' "$EX_WORK_DIR/preamble/ndim/summary.md"
grep -n 'coverage PCs' "$EX_WORK_DIR/preamble/preamble_summary.md"
echo
echo "--- choose_ndim.R with --gap 300 (ranks past k+gap+1 run backwards):"
EX_MP_GAP=300 /bin/bash preamble.sh --smoke --ndim-only > "$S/x7b.log" 2>&1; echo "preamble (gap=300) exit=$?"; grep -n 'Error\|error\|missing value\|ndim' "$S/x7b.log" | tail -5
