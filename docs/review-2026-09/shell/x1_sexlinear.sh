#!/bin/bash
# X1: covariates.tsv appears after an unadjusted run; does the adjusted sex_linear get recomputed?
S=/private/tmp/claude-501/-Users-Kitty-git-depthSV/1db72424-c4a1-45f4-b81e-5b466f3ec8a4/scratchpad/review_shell
REPO=/Users/Kitty/git/depthSV/.claude/worktrees/agent-a08479f00845fba52
export EX_WORK_DIR="$S/exA"
cd "$REPO/example/1000G_highcov" || exit 1
# fabricate preamble/covariates.tsv for every QC sample
Rscript -e '
suppressPackageStartupMessages(library(data.table)); set.seed(7)
ph <- fread(file.path(Sys.getenv("EX_WORK_DIR"), "inputs/standard/phenotypes.tsv"))
cv <- data.table(SAMPLE = ph$SAMPLE)
for (i in 1:10) cv[[paste0("GPC", i)]] <- round(rnorm(nrow(cv)), 5)
cv[, GPC_PROJECTED := 0L]
fwrite(cv, file.path(Sys.getenv("EX_WORK_DIR"), "preamble/covariates.tsv"), sep = "\t")'
before_sex="$(shasum "$EX_WORK_DIR/work/standard/assoc_ndim4/sex_linear.linear.chrX_3000001-4000000.txt.gz" | cut -c1-12)"
before_mt="$(shasum "$EX_WORK_DIR/work/standard/assoc_ndim4/mtdna_cn.linear.chrX_3000001-4000000.txt.gz" | cut -c1-12)"
/bin/bash run.sh --smoke --runner local > "$S/x1_run.log" 2>&1; echo "run.sh exit=$?"
grep -n 'covariates:\|analyses ->\|already complete, skipping: sex_linear\|already complete, skipping: mtdna_cn_adj\|wrote .*sex_linear\|wrote .*mtdna_cn_adj' "$S/x1_run.log" | head -12
echo "--- analyses.tsv now:"; grep -v '^#' "$EX_WORK_DIR/inputs/standard/analyses.tsv"
after_sex="$(shasum "$EX_WORK_DIR/work/standard/assoc_ndim4/sex_linear.linear.chrX_3000001-4000000.txt.gz" | cut -c1-12)"
echo "sex_linear shard before=$before_sex after=$after_sex (identical => the adjusted model was never fitted)"
echo "--- what the adjusted sex_linear would actually give (forced, one region):"
EX_SMOKE=1 /bin/bash -c 'source ./lib.sh; ex_export_dsv_env standard; bash "$DSV_ROOT/scripts/analyze.sh" --corrected "$DSV_CORRECTED_DIR/corrected_ndim4.chrX_3000001-4000000.txt.gz" --region chrX:3000001-4000000 --out "'"$S"'/x1_forced" --name sex_linear --method linear --model "SEX~cov_resids+GPC1+GPC2+GPC3+GPC4+GPC5+GPC6+GPC7+GPC8+GPC9+GPC10" --jobs 2 -- --minObs 30' > "$S/x1_forced.log" 2>&1
paste <(bgzip -dc "$EX_WORK_DIR/work/standard/assoc_ndim4/sex_linear.linear.chrX_3000001-4000000.txt.gz" | sed -n '2,4p' | cut -f4,5,8,10) \
      <(bgzip -dc "$S/x1_forced/sex_linear.linear.chrX_3000001-4000000.txt.gz" | sed -n '2,4p' | cut -f5,8,10)
echo "(left: Region N Estimate t from the skipped 'adjusted' run = old unadjusted; right: N Estimate t when actually fitted with the GPCs)"
echo "--- eval summary now claims:"; grep -c 'sex_linear' "$EX_WORK_DIR/eval/standard/summary.md"; grep 'sex_linear | sex_top_hit' "$EX_WORK_DIR/eval/standard/summary.md"
echo "--- stray shards from the previous (unadjusted) manifest still in assoc dir:"; ls "$EX_WORK_DIR/work/standard/assoc_ndim4" | grep -c 'log2_mtdna_cn\.linear\|mtdna_cn_null\.linear'
