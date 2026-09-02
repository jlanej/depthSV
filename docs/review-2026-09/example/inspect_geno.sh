#!/usr/bin/env bash
S=/private/tmp/claude-501/-Users-Kitty-git-depthSV/1db72424-c4a1-45f4-b81e-5b466f3ec8a4/scratchpad/review_example
G="$S/preamble_work/preamble/genotypes"
P="$S/preamble_work/preamble"
echo "=== preamble geno log"; grep -v '^\[.*\] present' "$S/preamble.geno.log" | tail -25
echo; echo "=== files"; ls -la "$G" | awk '{print $5, $9}'
echo; echo "=== pvar chromosome code (first data line) and psam header"
plink2 --zst-decompress "$G/chr22.pvar.zst" 2>/dev/null | grep -v '^##' | head -3 | cut -f1-5
head -2 "$G/samples.psam"
echo; echo "=== prune log: order of operations / exclude range / rm-dup / set-all-var-ids"
grep -n -i -E 'exclude|rm-dup|set-all-var-ids|snps-only|max-alleles|maf|geno|indep|variants (loaded|remaining)|removed|pruned|Note|Warning|Error' "$G/chr22.prune.log" | head -40
echo; echo "=== extract log"; grep -n -i -E 'rm-dup|variants (loaded|remaining)|Warning|Error' "$G/chr22.pruned.log" | head
echo; echo "=== king log"; grep -n -i -E 'king|excluded|remain|Warning' "$G/king.log" | head
echo; echo "=== eigenvec.allele header"; head -2 "$G/pca_unrel.eigenvec.allele" | cut -f1-9
echo "columns: $(head -1 "$G/pca_unrel.eigenvec.allele" | tr '\t' '\n' | wc -l | tr -d ' ')"
echo; echo "=== sscore header"; head -2 "$G/pca_proj.sscore" | cut -f1-8
echo "columns: $(head -1 "$G/pca_proj.sscore" | tr '\t' '\n' | wc -l | tr -d ' ')"
echo; echo "=== score log"; grep -n -i -E 'score|Warning|Error|variants' "$G/pca_proj.log" | head
echo; echo "=== calibration"; head -12 "$P/gpc_calibration.tsv"
echo; echo "=== covariates head"; head -3 "$P/covariates.tsv" | cut -f1-4,41-46
echo; echo "=== gpc summary"; cat "$P/summary.md"
echo; echo "=== eigenval"; head -12 "$G/pca_unrel.eigenval" | tr '\n' ' '; echo
