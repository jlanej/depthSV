#!/usr/bin/env bash
# Experiment 1: how many matrix rows does one R worker process get at
# biobank width, given DSV_BLOCK_BYTES=64MB and dsv_chunk_lines?
# Also: does GNU parallel --pipe need a writable TMPDIR (cloud /tmp sizing)?
set -u
here="$(cd "$(dirname "$0")" && pwd)"
ROOT=/Users/Kitty/git/depthSV/.claude/worktrees/agent-a8beeeb9b1a9d1356
source "$ROOT/lib/common.sh"
set +e

echo "== rows per R process (DSV_BLOCK_BYTES=$DSV_BLOCK_BYTES, cap 2000) =="
for n in 3202 50000 100000 200000 500000; do
    awk -v n="$n" 'BEGIN{printf "chr1\t1000\t2000"; for(i=0;i<n;i++) printf "\t31.42"; printf "\n"}' > "$here/row_$n.txt"
    bytes=$(wc -c < "$here/row_$n.txt" | tr -d ' ')
    lines=$(dsv_chunk_lines "$here/row_$n.txt" 2000)
    # corrected rows are wider (signed, 6 sig. digits): emulate "-0.012345"
    awk -v n="$n" 'BEGIN{printf "chr1\t1000\t2000\tchr1:1000-2000"; for(i=0;i<n;i++) printf "\t-0.012345"; printf "\n"}' > "$here/crow_$n.txt"
    cbytes=$(wc -c < "$here/crow_$n.txt" | tr -d ' ')
    clines=$(dsv_chunk_lines "$here/crow_$n.txt" 2000)
    printf 'samples=%-7s raw_row=%9s B -> %4s rows/R-process (correct)   corrected_row=%9s B -> %4s rows/R-process (analyze)\n' \
        "$n" "$bytes" "$lines" "$cbytes" "$clines"
    rm -f "$here/row_$n.txt" "$here/crow_$n.txt"
done

echo
echo "== GNU parallel --pipe with an unwritable TMPDIR =="
printf 'a\nb\nc\n' | TMPDIR=/nonexistent/dir parallel --pipe -k -L 1 cat 2>&1 | head -3
echo "exit=${PIPESTATUS[1]}"
echo "(parallel version: $(parallel --version | head -1))"
