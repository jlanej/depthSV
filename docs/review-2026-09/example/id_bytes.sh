#!/usr/bin/env bash
S=/private/tmp/claude-501/-Users-Kitty-git-depthSV/1db72424-c4a1-45f4-b81e-5b466f3ec8a4/scratchpad/review_example
echo "=== raw bytes of the three PC-table IDs the QC table does not match (committed svd.pcs.txt / svd.samples.txt):"
grep -E '^HG02635|^HG03025|^HG03366' "$S/spectrum/svd.pcs.txt" | cut -f1 | cat -A
grep -E '^HG02635|^HG03025|^HG03366' "$S/spectrum/svd.samples.txt" | cat -A
echo "=== and in the QC table:"; grep -E '^HG02635|^HG03025|^HG03366' "$S/spectrum/std.sample_qc.tsv" | cut -f1 | cat -A
echo "=== how the example's prepared PC table carries them (after suffix stripping):"
grep -E '^HG02635|^HG03025|^HG03366' "$S/smoke/inputs/standard/svd.pcs.txt" | cut -f1 | cat -A
echo "=== does the corrected-stage log in a real run report the drop? (smoke: these 3 are not in the 64-sample tree, so nothing shows)"
grep -h align "$S/smoke/work/standard/corrected/"*.log | sort -u
