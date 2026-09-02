#!/usr/bin/env bash
# What GNU parallel --pipe actually does with the pipeline's flags at 500k
# width: lines per job, and the cost of moving bytes through parallel itself
# (including -k output buffering) with a `cat` worker.
cd "$(dirname "$0")"
export TMPDIR=/private/tmp/claude-501/-Users-Kitty-git-depthSV/1db72424-c4a1-45f4-b81e-5b466f3ec8a4/scratchpad/review_perf/ptmp
mkdir -p "$TMPDIR"
if [ ! -s stream_500k.txt ]; then
  for i in 1 2 3 4 5 6 7 8 9 10; do tail -n +2 tables/raw_500000.txt; done > stream_500k.txt
fi
bytes=$(wc -c < stream_500k.txt | tr -d ' '); lines=$(wc -l < stream_500k.txt | tr -d ' ')
echo "stream: $lines rows, $bytes bytes ($((bytes / lines)) B/row)"
chunk=$(( 67108864 / (bytes / lines) ))
echo "dsv_chunk_lines would give: $chunk"
echo "--- lines per job with --block 64M -L $chunk -j4 (worker: wc -l):"
cat stream_500k.txt | parallel --pipe -k --halt now,fail=1 --block 67108864 -L "$chunk" -j 4 'wc -l' | tr '\n' ' '; echo
echo "--- lines per job with -L 2000 (the 3,202-width default cap) on the same stream:"
cat stream_500k.txt | parallel --pipe -k --halt now,fail=1 --block 67108864 -L 2000 -j 4 'wc -l' | tr '\n' ' '; echo
echo "--- baseline cat > /dev/null:"
/usr/bin/time -p sh -c 'cat stream_500k.txt > /dev/null' 2>&1 | tr '\n' ' '; echo
echo "--- through parallel --pipe -k (cat worker), -j4:"
/usr/bin/time -p sh -c "cat stream_500k.txt | parallel --pipe -k --halt now,fail=1 --block 67108864 -L $chunk -j 4 cat > /dev/null" 2>&1 | tr '\n' ' '; echo
echo "--- through parallel --pipe WITHOUT -k (cat worker), -j4:"
/usr/bin/time -p sh -c "cat stream_500k.txt | parallel --pipe --halt now,fail=1 --block 67108864 -L $chunk -j 4 cat > /dev/null" 2>&1 | tr '\n' ' '; echo
echo "--- through parallel --pipe -k --block 512M:"
/usr/bin/time -p sh -c "cat stream_500k.txt | parallel --pipe -k --halt now,fail=1 --block 536870912 -L $chunk -j 4 cat > /dev/null" 2>&1 | tr '\n' ' '; echo
echo "--- bgzip -dc throughput on the same bytes (single thread) and tabix on an indexed copy:"
if [ ! -s stream_500k.bgz ]; then
  { head -n 1 tables/raw_500000.txt; awk -F'\t' 'BEGIN{OFS="\t"} {$2=(NR-1)*1000; $3=NR*1000; print}' stream_500k.txt; } | bgzip -@ 8 > stream_500k.bgz
  tabix -f -p bed stream_500k.bgz
fi
ls -la stream_500k.bgz stream_500k.bgz.tbi
/usr/bin/time -p sh -c 'bgzip -dc stream_500k.bgz > /dev/null' 2>&1 | tr '\n' ' '; echo " (bgzip -dc, 1 thread)"
/usr/bin/time -p sh -c 'bgzip -dc -@ 4 stream_500k.bgz > /dev/null' 2>&1 | tr '\n' ' '; echo " (bgzip -dc, 4 threads)"
/usr/bin/time -p sh -c 'tabix stream_500k.bgz chr1 > /dev/null' 2>&1 | tr '\n' ' '; echo " (tabix whole contig)"
/usr/bin/time -p sh -c 'tabix stream_500k.bgz chr1:200001-201000 | wc -c' 2>&1 | tr '\n' ' '; echo " (tabix one bin from the middle)"
/usr/bin/time -p sh -c 'tabix stream_500k.bgz chr1:200001-220000 | wc -l' 2>&1 | tr '\n' ' '; echo " (tabix 20 bins)"
rm -rf "$TMPDIR"
