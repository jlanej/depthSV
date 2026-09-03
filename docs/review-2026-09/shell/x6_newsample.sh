#!/bin/bash
# X6: a sample added to the mosdepth tree after the join (REVIEW 1.7 through the example)
S=/private/tmp/claude-501/-Users-Kitty-git-depthSV/1db72424-c4a1-45f4-b81e-5b466f3ec8a4/scratchpad/review_shell
REPO=/Users/Kitty/git/depthSV/.claude/worktrees/agent-a08479f00845fba52
export EX_WORK_DIR="$S/exF"
cd "$REPO/example/1000G_highcov" || exit 1
src="$(ls "$EX_WORK_DIR"/smoke_mosdepth/standard/*.regions.bed.gz | head -1)"
# a real 1000G ID that has a QC row but was not among the first 64 simulated
newid="$(awk -F'\t' 'NR>1{print $1}' "$EX_WORK_DIR/github_cache/standard/sample_qc.tsv" | grep -v -F -f <(sed 's/\.by1000.*//' <(ls "$EX_WORK_DIR/smoke_mosdepth/standard" | grep regions)) | head -1)"
echo "adding $newid to the standard tree (copy of $(basename "$src"))"
cp "$src" "$EX_WORK_DIR/smoke_mosdepth/standard/$newid.by1000.regions.bed.gz"
/bin/bash run.sh --smoke --runner local > "$S/x6.log" 2>&1; echo "run.sh exit=$?"
grep -n 'region files\|present in the coverage\|already complete\|join:\|verdict\|WARN\|sample set' "$S/x6.log" | head
echo "--- manifest lines: $(grep -c . "$EX_WORK_DIR/inputs/standard/mosdepth.manifest.txt"); matrix manifest says samples: $(awk -F'\t' '$1=="samples"{print $2}' "$EX_WORK_DIR/work/standard/join/depth.matrix.manifest"); matrix columns: $(( $(bgzip -dc "$EX_WORK_DIR/work/standard/join/depth.matrix.txt.gz" | head -1 | awk -F'\t' '{print NF}') - 3 ))"
echo "--- cross_mode_samples.txt:"; cat "$EX_WORK_DIR/inputs/cross_mode_samples.txt"
