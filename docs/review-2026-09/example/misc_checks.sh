#!/usr/bin/env bash
S=/private/tmp/claude-501/-Users-Kitty-git-depthSV/1db72424-c4a1-45f4-b81e-5b466f3ec8a4/scratchpad/review_example
EX=/Users/Kitty/git/depthSV/.claude/worktrees/agent-aa17d545ebe28b8cc/example/1000G_highcov
R=/Users/Kitty/git/depthSV/.claude/worktrees/agent-aa17d545ebe28b8cc
cd "$EX" || exit 9

echo "=== (1) compare FAIL: exit code of 04_compare_modes.sh on the stale tree, and whether 02 swallows it"
export EX_WORK_DIR="$S/smoke_stale" EX_SMOKE=1 EX_RUNNER=local
bash 04_compare_modes.sh > /dev/null 2>&1; echo "04 exit=$?"
grep -n 'compare_modes.sh\|05_profile.sh' 02_run_depthsv.sh | grep -n 'WARN'
grep -c FAIL "$S/smoke_stale/compare/standard_vs_fast/summary.md"
grep 'Verdict' "$S/smoke_stale/compare/summary.md"

echo; echo "=== (2) stale finalize lock: pre-create the lock, rerun 02 locally (all units skip)"
rm -rf "$S/smoke_lock"; cp -R "$S/smoke" "$S/smoke_lock"
export EX_WORK_DIR="$S/smoke_lock" EX_MODES="standard fast"
mkdir "$EX_WORK_DIR/.finalize.lock"
before=$(stat -f %m "$EX_WORK_DIR/compare/summary.md")
bash 02_run_depthsv.sh --runner local 2>&1 | grep -E 'finalize|evaluate\]' ; echo "02 exit=${PIPESTATUS[0]}"
after=$(stat -f %m "$EX_WORK_DIR/compare/summary.md")
echo "compare/summary.md mtime before=$before after=$after (equal => comparison silently not refreshed)"
ls -d "$EX_WORK_DIR/.finalize.lock" && echo "lock still present"

echo; echo "=== (3) is EX_NDIM exported by config.sh (i.e. frozen for SLURM child jobs), or re-derived per job?"
unset EX_WORK_DIR EX_SMOKE EX_MODES EX_NDIM
bash -c 'export EX_WORK_DIR=/tmp/x; source lib.sh >/dev/null 2>&1; echo "EX_NDIM in this shell=$EX_NDIM"; env | grep -c "^EX_NDIM=" ; env | grep "^EX_" | cut -d= -f1 | tr "\n" " "'; echo
echo "--- what a child job would compute once preamble/ndim.txt appears:"
mkdir -p /tmp/dsvx_ndimtest/preamble; echo 37 > /tmp/dsvx_ndimtest/preamble/ndim.txt
bash -c 'export EX_WORK_DIR=/tmp/dsvx_ndimtest; source lib.sh >/dev/null 2>&1; echo "child EX_NDIM=$EX_NDIM DSV dir=$(ex_assoc_dir standard)"'
rm -rf /tmp/dsvx_ndimtest

echo; echo "=== (4) run.sh's ndim provenance line in smoke mode when ndim.txt exists (CI order: preamble --ndim-only, then run --smoke)"
bash -c 'export EX_WORK_DIR="'"$S"'/preamble_work" EX_SMOKE=1; source lib.sh >/dev/null 2>&1; echo "ndim=$EX_NDIM ($([ -s "$EX_PREAMBLE_DIR/ndim.txt" ] && echo from\ the\ preamble || echo default)) ndim.txt says $(cat $EX_PREAMBLE_DIR/ndim.txt)"'

echo; echo "=== (5) .gitignore vs the sbatch output file preamble.sh writes into the example directory"
grep -n 'out' "$R/.gitignore" || echo "no *.out pattern in .gitignore"
grep -n 'SBATCH --output' preamble.sh

echo; echo "=== (6) analyses.tsv fread with the comment block: any warning?"
Rscript -e 'suppressPackageStartupMessages(library(data.table)); m <- fread("'"$S"'/smoke/inputs/standard/analyses.tsv", header=FALSE, sep="\t", col.names=c("name","method","model")); print(m[!grepl("^#", name) & nzchar(name)])' 2>&1

echo; echo "=== (7) EX_COVARIATES sed edge cases in 01_prepare_inputs.sh"
for cov in "SEX+GPC1+GPC2" "GPC1+SEX" "SEX" "SEX+SEX2" "AGE*SEX" "SEXY+GPC1" "GPC1+SEX+GPC2"; do
  adj_nosex="$(printf '%s' "+$cov" | sed -e 's/+SEX+/+/g' -e 's/+SEX$//')"; [ "$adj_nosex" = "+" ] && adj_nosex=""
  printf '%-16s -> sex_linear model: SEX~cov_resids%s\n' "$cov" "$adj_nosex"
done

echo; echo "=== (8) hg38 unit count at EX_WINDOW=25000000 over chr1-22,X,Y,M"
python3 - <<'EOF'
import math
sizes={'1':248956422,'2':242193529,'3':198295559,'4':190214555,'5':181538259,'6':170805979,'7':159345973,'8':145138636,'9':138394717,'10':133797422,'11':135086622,'12':133275309,'13':114364328,'14':107043718,'15':101991189,'16':90338345,'17':83257441,'18':80373285,'19':58617616,'20':64444167,'21':46709983,'22':50818468,'X':156040895,'Y':57227415,'M':16569}
print(sum(math.ceil(v/25e6) for v in sizes.values()), "units")
EOF
