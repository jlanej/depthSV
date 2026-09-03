#!/bin/bash
S=/private/tmp/claude-501/-Users-Kitty-git-depthSV/1db72424-c4a1-45f4-b81e-5b466f3ec8a4/scratchpad/review_shell
REPO=/Users/Kitty/git/depthSV/.claude/worktrees/agent-a08479f00845fba52
cd "$REPO" || exit 1

echo "=== X1b: forced adjusted sex_linear (why did the earlier forced run fail?) ==="
tail -3 "$S/x1_forced.log"
rm -rf "$S/x1_forced"
/bin/bash scripts/analyze.sh --corrected "$S/exA/work/standard/corrected/corrected_ndim4.chrX_3000001-4000000.txt.gz" \
   --pheno "$S/exA/inputs/standard/phenotypes.tsv" --region chrX:3000001-4000000 --out "$S/x1_forced" \
   --name sex_linear --method linear --model "SEX~cov_resids+GPC1+GPC2+GPC3+GPC4+GPC5+GPC6+GPC7+GPC8+GPC9+GPC10" --jobs 2 -- --minObs 30 > "$S/x1_forced.log" 2>&1; echo "exit=$?"
echo "skipped-as-complete shard (old unadjusted fit)   vs   adjusted fit actually run:"
paste <(bgzip -dc "$S/exA/work/standard/assoc_ndim4/sex_linear.linear.chrX_3000001-4000000.txt.gz" | sed -n '2,4p' | cut -f4,8,10) \
      <(bgzip -dc "$S/x1_forced/sex_linear.linear.chrX_3000001-4000000.txt.gz" | sed -n '2,4p' | cut -f8,10)

echo
echo "=== X9: workflows/slurm_array.sh login-node run (no SLURM), bash 3.2 ==="
FX="$S/fx"; W="$S/x9"; rm -rf "$W"; mkdir -p "$W"
printf 'chr1\nchr2\n' > "$W/regions.txt"
env -i PATH="$PATH" HOME="$HOME" DSV_MATRIX="$S/smoke dir/join/depth.matrix.txt.gz" DSV_PCS="$FX/svd.pcs.txt" DSV_COVERAGE="$FX/autosomal.median.txt" \
    DSV_PHENO="$FX/phenotypes.tsv" DSV_PHENO_MANIFEST="$REPO/conf/phenotypes.example.tsv" DSV_CORRECTED_DIR="$W/corr" DSV_RESULTS_DIR="$W/assoc" DSV_NDIM=4 DSV_MIN_OBS=30 \
    /bin/bash workflows/slurm_array.sh "$W/regions.txt" > "$W/log" 2>&1; echo "exit=$? ; outputs: $(ls "$W/assoc" 2>/dev/null | grep -c 'txt.gz$')"; tail -2 "$W/log"

echo
echo "=== X11: join.sh resume after being killed mid-batch (README: 'resumes from its last finished batch') ==="
J="$S/x11"; rm -rf "$J"; mkdir -p "$J"
perl -e 'setpgrp(0,0); exec @ARGV' /bin/bash scripts/join.sh --manifest "$FX/mosdepth.input.txt" --out "$J/out" --batch-size 10 --jobs 1 --threads 1 > "$J/run1.log" 2>&1 &
pid=$!
for i in $(seq 1 200); do [ -f "$J/out/.join.work/batch.000002.tsv.gz.done" ] && break; sleep 0.05; done
kill -TERM -- -"$pid" 2>/dev/null; wait "$pid" 2>/dev/null; echo "killed after: $(ls "$J/out/.join.work" | grep -c 'tsv.gz.done$') batches done; cols left: $(ls "$J/out/.join.work/cols" | wc -l | tr -d ' ')"
ls "$J/out/.join.work" | tr '\n' ' '; echo
/bin/bash scripts/join.sh --manifest "$FX/mosdepth.input.txt" --out "$J/out" --batch-size 10 --jobs 1 --threads 1 > "$J/run2.log" 2>&1; echo "rerun exit=$?"
grep -n 'resuming\|wrote\|ERROR\|FAILED' "$J/run2.log"
cmp "$J/out/depth.matrix.txt.gz" "$S/smoke dir/join/depth.matrix.txt.gz" && echo "resumed matrix byte-identical to the clean one" || echo "resumed matrix DIFFERS"

echo
echo "=== README / doc line checks ==="
grep -n 'cov_resids' README.md | head -3; echo "(README mentions of cov_resids: $(grep -c cov_resids README.md))"
grep -n 'chr1:0-10000000\|byte-identical\|bash. ≥ 3.2\|DSV_MODULES' README.md
grep -n 'sbatch preamble.sh\|Minutes end to end\|can never silently reuse\|runs both on a laptop\|one submission' example/1000G_highcov/README.md
grep -n 'depthsv:dev' workflows/depthsv.wdl Dockerfile README.md
grep -n 'SBATCH --output' example/1000G_highcov/preamble.sh workflows/slurm_array.sh
