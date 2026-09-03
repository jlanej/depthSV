#!/usr/bin/env bash
# Does plink2's --exclude range match "22" in the range file against the pvar's chromosome code?
S=/private/tmp/claude-501/-Users-Kitty-git-depthSV/1db72424-c4a1-45f4-b81e-5b466f3ec8a4/scratchpad/review_example
G="$S/preamble_work/preamble/genotypes"
T="$S/exclude_test"; rm -rf "$T"; mkdir -p "$T"
plink2 --zst-decompress "$G/chr22.pgen.zst" "$T/chr22.pgen" >/dev/null 2>&1
code=$(plink2 --zst-decompress "$G/chr22.pvar.zst" 2>/dev/null | grep -v '^#' | head -1 | cut -f1)
echo "pvar chromosome code: '$code'"
for style in "22" "chr22"; do
  printf '%s\t20000000\t30000000\tTEST\n' "$style" > "$T/range_$style.txt"
  plink2 --pgen "$T/chr22.pgen" --pvar "$G/chr22.pvar.zst" --psam "$G/samples.psam" \
         --threads 2 --memory 4000 --silent --snps-only just-acgt --max-alleles 2 \
         --exclude range "$T/range_$style.txt" --write-snplist --out "$T/ex_$style" 2>&1 | tail -2
  echo "range file '$style': $(grep -E 'exclude range|remaining' "$T/ex_$style.log" | tr '\n' ' ')"
done
rm -f "$T/chr22.pgen"
