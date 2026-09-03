#!/bin/bash
S=/private/tmp/claude-501/-Users-Kitty-git-depthSV/1db72424-c4a1-45f4-b81e-5b466f3ec8a4/scratchpad/review_shell
REPO=/Users/Kitty/git/depthSV/.claude/worktrees/agent-a08479f00845fba52
cd "$REPO" || exit 1
FX="$S/fx"
[ -s "$FX/mosdepth.input.txt" ] || Rscript tests/make_fixtures.R "$FX" 60 200 >/dev/null 2>&1

echo "=== E1: join.sh with a manifest of many duplicates (SIGPIPE in the dup check) ==="
awk '{for(i=0;i<300;i++) print}' "$FX/mosdepth.input.txt" > "$S/manydup.txt"
wc -l < "$S/manydup.txt"
/bin/bash scripts/join.sh --manifest "$S/manydup.txt" --out "$S/manydup.out" > "$S/manydup.log" 2>&1; echo "exit=$?"
tail -3 "$S/manydup.log"

echo
echo "=== E2: duplicate SAMPLE in --coverage with different medians (correct.R) ==="
{ cat "$FX/autosomal.median.txt"; printf 'SAMPLE001\t1000\nSAMPLE002\t1000\n'; } > "$S/dupcov.txt"
Rscript R/correct.R -i "$FX/svd.pcs.txt" -f "$S/smoke dir/join/depth.matrix.txt.gz" -c "$S/dupcov.txt" -d 0 2>"$S/dupcov.err" | head -2 | cut -f1-6 > "$S/dupcov.out"; echo "exit=${PIPESTATUS[0]}"
Rscript R/correct.R -i "$FX/svd.pcs.txt" -f "$S/smoke dir/join/depth.matrix.txt.gz" -c "$FX/autosomal.median.txt" -d 0 2>/dev/null | head -2 | cut -f1-6 > "$S/cleancov.out"
echo "stderr: $(cat "$S/dupcov.err")"
echo "with dup rows:"; cat "$S/dupcov.out"
echo "clean:"; cat "$S/cleancov.out"
# Reverse the duplicate order: the 1000 rows first
{ head -1 "$FX/autosomal.median.txt"; printf 'SAMPLE001\t1000\nSAMPLE002\t1000\n'; tail -n +2 "$FX/autosomal.median.txt"; } > "$S/dupcov2.txt"
Rscript R/correct.R -i "$FX/svd.pcs.txt" -f "$S/smoke dir/join/depth.matrix.txt.gz" -c "$S/dupcov2.txt" -d 0 2>/dev/null | head -2 | cut -f1-6
echo "(SAMPLE001 first data column differs between orderings if match() silently takes the first row)"

echo
echo "=== E13: USER unset (containers/cron) ==="
cd example/1000G_highcov
env -u USER EX_WORK_DIR="$S/nouser" /bin/bash 00_fetch_inputs.sh --help > "$S/nouser.log" 2>&1; echo "exit=$?"; head -3 "$S/nouser.log"
env -u USER EX_WORK_DIR="$S/nouser" NGSPCA_WORK_DIR=/nonexistent /bin/bash 00_fetch_inputs.sh --help > "$S/nouser2.log" 2>&1; echo "with NGSPCA_WORK_DIR set too: exit=$?"; head -2 "$S/nouser2.log"

echo
echo "=== E14: --help under bash 3.2 for every example script and stage script ==="
for f in run.sh preamble.sh 00_fetch_inputs.sh 01_prepare_inputs.sh 02_run_depthsv.sh 03_evaluate.sh 04_compare_modes.sh 05_profile.sh; do
  EX_WORK_DIR="$S/helpwork" /bin/bash "$f" --help > "$S/help.$f.log" 2>&1; printf '%s: exit=%s lines=%s first=%s\n' "$f" "$?" "$(wc -l < "$S/help.$f.log" | tr -d ' ')" "$(head -1 "$S/help.$f.log")"
done
cd "$REPO"
for f in scripts/join.sh scripts/correct.sh scripts/analyze.sh scripts/regions.sh; do
  /bin/bash "$f" --help > "$S/help.$(basename "$f").log" 2>&1; printf '%s: exit=%s lines=%s\n' "$f" "$?" "$(wc -l < "$S/help.$(basename "$f").log" | tr -d ' ')"
done
