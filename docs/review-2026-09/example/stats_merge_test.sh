#!/usr/bin/env bash
# correct.sh with >= 10 parallel chunks: is the stats merge (lexicographic glob) still fatal in this tree?
S=/private/tmp/claude-501/-Users-Kitty-git-depthSV/1db72424-c4a1-45f4-b81e-5b466f3ec8a4/scratchpad/review_example
R=/Users/Kitty/git/depthSV/.claude/worktrees/agent-aa17d545ebe28b8cc
W="$S/smoke"
export DSV_MATRIX="$W/work/standard/join/depth.matrix.txt.gz"
export DSV_PCS="$W/inputs/standard/svd.pcs.txt"
export DSV_COVERAGE="$W/inputs/standard/autosomal.median.txt"
export DSV_NDIM=4 DSV_JOBS=4 DSV_THREADS=2
rowbytes=$(gzip -cd "$DSV_MATRIX" | sed -n 2p | wc -c | tr -d ' ')
rows=$(tabix "$DSV_MATRIX" chr20:1-1000000 | wc -l | tr -d ' ')
echo "smoke matrix: $rows rows in chr20:1-1000000, $rowbytes bytes/row -> $((rows * rowbytes)) bytes"
export DSV_BLOCK_BYTES=$(( rows * rowbytes / 12 ))
echo "DSV_BLOCK_BYTES=$DSV_BLOCK_BYTES => ~12 chunks"
rm -rf "$S/stats_test"; mkdir -p "$S/stats_test"
bash "$R/scripts/correct.sh" --region chr20:1-1000000 --out "$S/stats_test" 2>&1 | tail -6
echo "exit=${PIPESTATUS[0]}"
ls "$S/stats_test"
echo
echo "=== the real configuration: bytes per raw-matrix row at 3,202 samples"
echo "mosdepth prints depth with 2 decimals; at 30x most values are 5-6 chars + tab."
python3 - <<'EOF'
for width in (6.5, 7.5):
    row = 3202 * width + 20
    unit_bytes = 25000 * row      # 25 Mb window at 1 kb bins
    print(f"row ~{row:.0f} B -> 25 Mb unit ~{unit_bytes/1e6:.0f} MB -> {unit_bytes/67108864:.1f} x 64 MB blocks")
EOF
