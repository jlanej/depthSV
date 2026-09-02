#!/usr/bin/env bash
# Stale-shard reproduction: same ndim, changed covariates (a column already in phenotypes.tsv)
# and a changed null-permutation seed. Does the rerun reuse the old sex_linear / mtdna_cn_null shards?
S=/private/tmp/claude-501/-Users-Kitty-git-depthSV/1db72424-c4a1-45f4-b81e-5b466f3ec8a4/scratchpad/review_example
EX=/Users/Kitty/git/depthSV/.claude/worktrees/agent-aa17d545ebe28b8cc/example/1000G_highcov
rm -rf "$S/smoke_stale"; cp -R "$S/smoke" "$S/smoke_stale"
cd "$EX" || exit 9
export EX_WORK_DIR="$S/smoke_stale" EX_SMOKE=1 EX_RUNNER=local EX_MODES=standard
export EX_COVARIATES=MEAN_AUTOSOMAL_COV EX_PHENO_SEED=1
A="$EX_WORK_DIR/work/standard/assoc_ndim4"
echo "=== before: md5 of the sex_linear and mtdna_cn_null chr20 shards, and phenotypes MTDNA_CN_NULL head"
md5 -q "$A/sex_linear.linear.chr20_1-1000000.txt.gz" "$A/mtdna_cn_null.linear.chr20_1-1000000.txt.gz"
cut -f1,4 "$EX_WORK_DIR/inputs/standard/phenotypes.tsv" | head -3
echo "=== rerun 01 with EX_COVARIATES=MEAN_AUTOSOMAL_COV EX_PHENO_SEED=1"
bash 01_prepare_inputs.sh 2>&1 | grep -E 'analyses|warn|WARN'
grep -v '^#' "$EX_WORK_DIR/inputs/standard/analyses.tsv"
echo "=== after 01: phenotypes MTDNA_CN_NULL head (changed?)"
cut -f1,4 "$EX_WORK_DIR/inputs/standard/phenotypes.tsv" | head -3
echo "=== rerun 02 (local): which analyses are skipped as already complete?"
bash 02_run_depthsv.sh --mode standard --runner local 2>&1 | grep -E 'already complete|manifest complete|wrote .*assoc|FAIL|evaluate\]' | sort | uniq -c | sort -rn | head -20
echo "=== after: md5 of the same shards (unchanged => stale results reported under a changed model)"
md5 -q "$A/sex_linear.linear.chr20_1-1000000.txt.gz" "$A/mtdna_cn_null.linear.chr20_1-1000000.txt.gz"
echo "=== evaluation summary now claims (analyses listed):"
grep -E 'verdict|sex_linear \| sex_top_hit|mtdna_cn_null \| null_lambda' "$EX_WORK_DIR/eval/standard/summary.md"
echo "=== manifest says sex_linear model is:"; grep '^sex_linear' "$EX_WORK_DIR/inputs/standard/analyses.tsv"
echo "=== but the shard's log (from the ORIGINAL run) says the fit was:"; head -3 "$A/sex_linear.linear.chr20_1-1000000.log"
