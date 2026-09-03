#!/usr/bin/env bash
# The README's "lagging fast tree" scenario: the fast tree gains a sample after the first join.
S=/private/tmp/claude-501/-Users-Kitty-git-depthSV/1db72424-c4a1-45f4-b81e-5b466f3ec8a4/scratchpad/review_example
EX=/Users/Kitty/git/depthSV/.claude/worktrees/agent-aa17d545ebe28b8cc/example/1000G_highcov
cd "$EX" || exit 9
rm -rf "$S/smoke_grow"; cp -R "$S/smoke" "$S/smoke_grow"
export EX_WORK_DIR="$S/smoke_grow" EX_SMOKE=1 EX_RUNNER=local EX_MODES="standard fast"
export EX_SMOKE_DIR="$S/smoke_grow/smoke_mosdepth"
# the smoke tree lives under the copied work dir; 00 will see 64 >= EX_SMOKE_SAMPLES and keep it.
cp "$EX_SMOKE_DIR/fast/HG00096.by1000.regions.bed.gz" "$EX_SMOKE_DIR/fast/HG99999.by1000.regions.bed.gz"
echo "fast tree now has $(ls "$EX_SMOKE_DIR/fast"/*.regions.bed.gz | wc -l | tr -d ' ') samples"
bash 00_fetch_inputs.sh 2>&1 | grep -E 'simulated tree|SKIP' | sed 's/^\[[^]]*\] //'
bash 01_prepare_inputs.sh 2>&1 | grep -E 'region files|sample set|WARN' | sed 's/^\[[^]]*\] //'
cat "$EX_WORK_DIR/inputs/cross_mode_samples.txt"
bash 02_run_depthsv.sh --mode fast --runner local 2>&1 | grep -E 'join|already complete, skipping: .*depth.matrix' | sed 's/^\[[^]]*\] //' | head -3
echo "matrix columns (samples) in the fast matrix after the rerun: $(( $(gzip -cd "$EX_WORK_DIR/work/fast/join/depth.matrix.txt.gz" | head -1 | tr '\t' '\n' | wc -l | tr -d ' ') - 3 ))"
echo "manifest lists: $(grep -c . "$EX_WORK_DIR/inputs/fast/mosdepth.manifest.txt") samples"
grep -E 'verdict' "$EX_WORK_DIR/eval/fast/summary.md"
