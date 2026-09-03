#!/usr/bin/env bash
# Two chromosomes -> exercises --pmerge-list (the smoke run only ever takes the single --pfile branch).
S=/private/tmp/claude-501/-Users-Kitty-git-depthSV/1db72424-c4a1-45f4-b81e-5b466f3ec8a4/scratchpad/review_example
EX=/Users/Kitty/git/depthSV/.claude/worktrees/agent-aa17d545ebe28b8cc/example/1000G_highcov
cd "$EX" || exit 9
export EX_WORK_DIR="$S/preamble_work2"
mkdir -p "$EX_WORK_DIR/preamble/genotypes"
cp "$S/preamble_work/preamble/genotypes/chr22.pgen.zst" "$S/preamble_work/preamble/genotypes/chr22.pvar.zst" \
   "$S/preamble_work/preamble/genotypes/samples.psam" "$EX_WORK_DIR/preamble/genotypes/"
export EX_GENO_THREADS=4 EX_GENO_MEMORY_MB=8000 EX_GENO_CHROMS="21 22"
bash preamble.sh --smoke --genotypes-only > "$S/preamble.geno2.log" 2>&1
echo "exit=$?" >> "$S/preamble.geno2.log"
grep -v 'present:' "$S/preamble.geno2.log" | grep -E 'fetched|pruned SNPs|merging|pruned SNP set|unrelated|gpc\]|Error|ERROR|Warning|exit='
echo "--- merge log"; grep -n -i -E 'pmerge|Error|Warning|variants|samples' "$EX_WORK_DIR/preamble/genotypes/pruned.log" | head -12
