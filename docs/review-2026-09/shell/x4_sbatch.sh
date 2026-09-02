#!/bin/bash
# X4: what the SLURM driver actually passes to sbatch, and what environment a job would inherit
S=/private/tmp/claude-501/-Users-Kitty-git-depthSV/1db72424-c4a1-45f4-b81e-5b466f3ec8a4/scratchpad/review_shell
REPO=/Users/Kitty/git/depthSV/.claude/worktrees/agent-a08479f00845fba52
export EX_WORK_DIR="$S/exD"
export SBATCH_SHIM_DIR="$S/x4_calls"; rm -rf "$SBATCH_SHIM_DIR"; mkdir -p "$SBATCH_SHIM_DIR"
export PATH="$S/bin:$PATH"
cd "$REPO/example/1000G_highcov" || exit 1
# a real (non-smoke) configuration: ndim comes from preamble/ndim.txt (42)
export EX_SBATCH_EXTRA="--partition=short --account=acct"
/bin/bash 02_run_depthsv.sh --runner slurm --mode all > "$S/x4.log" 2>&1; echo "driver exit=$?"
cat "$S/x4.log"
for f in "$SBATCH_SHIM_DIR"/call.*.args; do echo "--- $(basename "$f"):"; tr '\n' ' ' < "$f"; echo; done
echo "--- EX_/DSV_ variables the join job would inherit (call.0.env):"; cat "$SBATCH_SHIM_DIR/call.0.env"
echo "--- driver's own view: EX_NDIM=$(EX_SMOKE=0 /bin/bash -c 'source ./lib.sh; echo $EX_NDIM') ; is EX_NDIM exported to jobs? $(grep -c '^EX_NDIM=' "$SBATCH_SHIM_DIR/call.0.env")"
echo "--- jobs.tsv recorded ids:"; cat "$EX_WORK_DIR/profile/jobs.tsv"
echo
echo "--- same driver, now simulating the dispatch stage inside a job (SLURM_ARRAY_TASK_ID unset): what array/eval get:"
rm -rf "$SBATCH_SHIM_DIR"; mkdir -p "$SBATCH_SHIM_DIR"
EX_SMOKE=1 /bin/bash 02_run_depthsv.sh --runner slurm --stage dispatch --mode standard > "$S/x4b.log" 2>&1; echo "dispatch exit=$?"; tail -3 "$S/x4b.log"
for f in "$SBATCH_SHIM_DIR"/call.*.args; do echo "--- $(basename "$f"):"; tr '\n' ' ' < "$f"; echo; done
