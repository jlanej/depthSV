#!/usr/bin/env bash
# Experiment 3: does DSV_MODULES help when the tools are only on PATH after
# `module load`? Emulate an HPC login: a `module` shell function that adds
# the tool directory to PATH, exported into the stage script's environment,
# and a PATH that lacks bgzip/tabix/parallel/Rscript until it runs.
set -u
here="$(cd "$(dirname "$0")" && pwd)"
ROOT=/Users/Kitty/git/depthSV/.claude/worktrees/agent-a8beeeb9b1a9d1356
TOOLS=/opt/homebrew/bin            # where bgzip/tabix/parallel/Rscript live here

module() { echo "[fake module] load $* -> adding $TOOLS to PATH" >&2; export PATH="$TOOLS:$PATH"; }
export -f module
export DSV_MODULES="htslib parallel R"

mkdir -p "$here/m" && cd "$here/m" || exit 1
printf 'x\n' > f.txt   # dsv_require_file only needs non-empty files

echo "== correct.sh with tools absent from PATH until 'module load' =="
env PATH=/usr/bin:/bin bash "$ROOT/scripts/correct.sh" \
    --matrix f.txt --pcs f.txt --coverage f.txt --region chr1 --out o 2>&1 | head -3
echo "exit=${PIPESTATUS[0]}"

echo
echo "== same, but with the module function invoked BEFORE the script (what a user would have to do by hand) =="
env PATH=/usr/bin:/bin bash -c 'module htslib; bash "$0" --matrix f.txt --pcs f.txt --coverage f.txt --region chr1 --out o' \
    "$ROOT/scripts/correct.sh" 2>&1 | head -3
echo "exit=${PIPESTATUS[0]}"
echo "(the second run gets past the command check and fails later, on the fake matrix — expected)"
