#!/usr/bin/env bash
# preamble.sh part B on the real chr22 files (1,066,557 variants x 3,202):
# zst decompress, the QC+prune pass (2 vs 8 threads), the extract pass, and
# memory, exactly as preamble.sh invokes plink2 (16000 MB, --silent).
cd "$(dirname "$0")"
mkdir -p pre && cd pre
SRC=/private/tmp/claude-501/-Users-Kitty-git-depthSV/1db72424-c4a1-45f4-b81e-5b466f3ec8a4/scratchpad/geno22
LD=/Users/Kitty/git/depthSV/.claude/worktrees/agent-a265e6aa144253f7b/example/1000G_highcov/resources/long_range_ld_grch38.txt
cp -f "$SRC/hg38_corrected.psam" samples.psam
ln -sf "$SRC/chr22_hg38.pgen.zst" chr22.pgen.zst
ln -sf "$SRC/chr22_hg38_rs_noannot.pvar.zst" chr22.pvar.zst
t() { /usr/bin/time -l -o t.tmp "$@" > /dev/null 2>&1; printf 'wall %ss  cpu %ss  maxRSS %s MB\n' "$(awk '/real/{print $1}' t.tmp)" "$(awk '/user/{u=$1} /sys/{s=$1} END{print u+s}' t.tmp)" "$(awk '/maximum resident/{printf "%.0f", $1/1048576}' t.tmp)"; }
printf 'zst-decompress pgen (%.0f MB zst): ' "$(echo "$(stat -f %z "$SRC/chr22_hg38.pgen.zst")/1e6" | bc -l)"; t plink2 --zst-decompress chr22.pgen.zst chr22.pgen
printf '  -> %.0f MB pgen\n' "$(echo "$(stat -f %z chr22.pgen)/1e6" | bc -l)"
id_args=(--set-all-var-ids '@:#:$r:$a' --new-id-max-allele-len 10 missing --rm-dup exclude-all)
for th in 2 8; do
  printf 'QC+prune pass, %s threads: ' "$th"
  t plink2 --pgen chr22.pgen --pvar chr22.pvar.zst --psam samples.psam --threads "$th" --memory 16000 --silent \
       --snps-only just-acgt --max-alleles 2 "${id_args[@]}" --maf 0.05 --geno 0.01 --exclude range "$LD" \
       --indep-pairwise 1000kb 1 0.1 --out chr22.prune
done
grep -E 'variants (loaded|remaining)|pruned|--indep-pairwise' chr22.prune.log | head -8
printf 'extract pass (%s SNPs), 8 threads: ' "$(wc -l < chr22.prune.prune.in | tr -d ' ')"
t plink2 --pgen chr22.pgen --pvar chr22.pvar.zst --psam samples.psam --threads 8 --memory 16000 --silent "${id_args[@]}" --extract chr22.prune.prune.in --make-pgen --out chr22.pruned
printf 'QC-only pass (no pruning), 8 threads, --make-pgen of the MAF/geno-filtered set: '
t plink2 --pgen chr22.pgen --pvar chr22.pvar.zst --psam samples.psam --threads 8 --memory 16000 --silent \
     --snps-only just-acgt --max-alleles 2 "${id_args[@]}" --maf 0.05 --geno 0.01 --exclude range "$LD" --make-pgen --out chr22.qc
printf '  QC set: %s variants, %.0f MB pgen\n' "$(grep -vc '^#' chr22.qc.pvar)" "$(echo "$(stat -f %z chr22.qc.pgen)/1e6" | bc -l)"
printf 'prune on the QC set (in memory-sized pgen), 8 threads: '
t plink2 --pfile chr22.qc --threads 8 --memory 16000 --silent --indep-pairwise 1000kb 1 0.1 --out chr22.qc.prune
printf 'variant-count scaling: chr22 = %s variants; whole callset per genotype_sources.tsv is ~%s\n' "$(zstd -dc chr22.pvar.zst | grep -vc '^#')" "see README"
cd .. && rm -rf pre
