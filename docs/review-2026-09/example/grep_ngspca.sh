#!/usr/bin/env bash
D=/private/tmp/claude-501/-Users-Kitty-git-depthSV/1db72424-c4a1-45f4-b81e-5b466f3ec8a4/scratchpad/review_example/ngspca_docs
cd "$D" || exit 9
echo "=== variables the example README quotes: RANDOM_SEED COMPARE_FAST_MODE NUM_PC MOSDEPTH_DIR NGSPCA_OUTPUT QC_OUTPUT"
for v in RANDOM_SEED COMPARE_FAST_MODE NUM_PC MOSDEPTH_DIR NGSPCA_OUTPUT QC_OUTPUT MOSDEPTH_BIN_SIZE; do
  printf '%-18s README:%3d  02:%2d  03a:%2d  03:%2d\n' "$v" "$(grep -c "$v" README.md)" "$(grep -c "$v" 02_run_ngspca.sh)" "$(grep -c "$v" 03a_mosdepth_coverage_summary.sh)" "$(grep -c "$v" 03_collect_qc.sh)"
done
echo; echo "=== dir names: mosdepth_output_fast ngspca_output_fast qc_output_fast ngspca_output_seed"
grep -n -E 'mosdepth_output_fast|ngspca_output_fast|qc_output_fast|_seed[0-9]|seed control|seed-control|reseed' README.md | head -20
echo; echo "=== autosomal.median.txt / N_BINS / HQ_MEDIAN_COV sourcing"
grep -n -E 'autosomal\.median|N_BINS|HQ_MEDIAN_COV' README.md 03a_mosdepth_coverage_summary.sh 03_collect_qc.sh | head -30
echo; echo "=== randomized / randomSeed / singular values / NUM_PC in README"
grep -n -i -E 'random|singular|svd\.|numPC|NUM_PC|-numPC' README.md | head -30
echo; echo "=== fast-mode section head"
grep -n -i -E 'fast' README.md | head -40
echo; echo "=== MTDNA_CN definition"
grep -n -E 'MTDNA_CN|mtDNA|MITO_COV_RATIO' README.md 03a_mosdepth_coverage_summary.sh 03_collect_qc.sh | head -20
echo; echo "=== 02_run_ngspca.sh key lines"
grep -n -E 'NUM_PC|RANDOM_SEED|randomSeed|numPC|sampleSuffix|ngspca_output|median' 02_run_ngspca.sh | head -20
