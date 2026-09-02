#!/bin/bash
# X12: what --force redoes locally; and --stage unit outside SLURM
S=/private/tmp/claude-501/-Users-Kitty-git-depthSV/1db72424-c4a1-45f4-b81e-5b466f3ec8a4/scratchpad/review_shell
REPO=/Users/Kitty/git/depthSV/.claude/worktrees/agent-a08479f00845fba52
export EX_WORK_DIR="$S/exG"
cd "$REPO/example/1000G_highcov" || exit 1
rmdir "$EX_WORK_DIR/.finalize.lock" 2>/dev/null
m_before="$(stat -f %m "$EX_WORK_DIR/work/standard/join/depth.matrix.txt.gz")"
t_before="$(ls "$EX_WORK_DIR/smoke_mosdepth/standard" | wc -l | tr -d ' ')"
sim_before="$(stat -f %m "$EX_WORK_DIR/smoke_mosdepth/standard/smoke.params.txt")"
n_rec_before="$(ls "$EX_WORK_DIR/profile/timings.d" | wc -l | tr -d ' ')"
/bin/bash run.sh --smoke --runner local --mode standard --force > "$S/x12.log" 2>&1; echo "run.sh --force exit=$?"
echo "matrix rebuilt: $([ "$(stat -f %m "$EX_WORK_DIR/work/standard/join/depth.matrix.txt.gz")" != "$m_before" ] && echo yes || echo no); simulated tree re-made: $([ "$(stat -f %m "$EX_WORK_DIR/smoke_mosdepth/standard/smoke.params.txt")" != "$sim_before" ] && echo yes || echo no)"
echo "wrote lines: $(grep -c 'wrote ' "$S/x12.log"); skipped lines: $(grep -c 'already complete' "$S/x12.log")"
echo "timing records: $n_rec_before -> $(ls "$EX_WORK_DIR/profile/timings.d" | wc -l | tr -d ' ')"
echo "--- --stage unit outside SLURM (SLURM_ARRAY_TASK_ID unset -> line 1):"
EX_SMOKE=1 /bin/bash 02_run_depthsv.sh --stage unit --mode standard > "$S/x12b.log" 2>&1; echo "exit=$?"; grep -n 'task 1\|already complete' "$S/x12b.log" | head -3
echo "--- run.sh --prepare-only:"
/bin/bash run.sh --smoke --runner local --prepare-only > "$S/x12c.log" 2>&1; echo "exit=$?"; tail -1 "$S/x12c.log"
