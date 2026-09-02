#!/usr/bin/env bash
# paste throughput vs number of input columns, bgzip compress throughput, and
# the per-sample extraction cost on a genome-sized (3.1M-line) mosdepth file.
cd "$(dirname "$0")"
mkdir -p pastebench && cd pastebench
mk() {  # mk <ncols> <nlines> : column files of mosdepth-like values
  local n="$1" L="$2" i
  rm -f col_*;
  awk -v n="$n" -v L="$L" 'BEGIN{srand(3); for(i=1;i<=n;i++){f=sprintf("col_%04d",i); for(j=1;j<=L;j++) printf "%.2f\n", 20+rand()*30 > f; close(f)}}'
}
run_paste() {  # run_paste <ncols> <nlines>
  mk "$1" "$2"
  local t0 t1 bytes
  t0=$(python3 -c 'import time;print(time.time())')
  bytes=$(paste col_* | wc -c | tr -d ' ')
  t1=$(python3 -c 'import time;print(time.time())')
  printf 'paste %4d cols x %8d lines: %7.1f MB out in %6.2f s -> %6.0f MB/s\n' "$1" "$2" "$(echo "$bytes/1e6" | bc -l)" "$(echo "$t1-$t0" | bc -l)" "$(echo "$bytes/1e6/($t1-$t0)" | bc -l)"
}
run_paste 2 10000000
run_paste 57 1000000
run_paste 708 100000
run_paste 1024 60000
# bgzip compress throughput on paste output (708 cols x 100k lines ~ 425 MB)
mk 708 100000
paste col_* > big.txt
sz=$(stat -f %z big.txt)
for opt in "-l 1 -@ 1" "-l 1 -@ 4" "-l 1 -@ 8" "-@ 1" "-@ 4" "-@ 8"; do
  t0=$(python3 -c 'import time;print(time.time())'); bgzip $opt -c big.txt > big.bgz; t1=$(python3 -c 'import time;print(time.time())')
  printf 'bgzip %-10s: %6.1f MB -> %6.1f MB (ratio %.2f) in %5.2f s -> %5.0f MB/s in\n' "$opt" "$(echo "$sz/1e6" | bc -l)" "$(echo "$(stat -f %z big.bgz)/1e6" | bc -l)" "$(echo "$sz/$(stat -f %z big.bgz)" | bc -l)" "$(echo "$t1-$t0" | bc -l)" "$(echo "$sz/1e6/($t1-$t0)" | bc -l)"
done
# FIFO fan-in as join.sh does it: 8 bgzip -dc feeding 8 fifos into one paste
rm -f col_* big.bgz
split -l 12500 -d big.txt part_
for p in part_*; do bgzip -l 1 -@ 2 "$p"; done
mkdir -p fifos; rm -f fifos/*
i=0; pids=""
for p in part_*.gz; do i=$((i+1)); mkfifo "fifos/f$i"; bgzip -dc "$p" > "fifos/f$i" & pids="$pids $!"; done
t0=$(python3 -c 'import time;print(time.time())')
bytes=$(paste fifos/* | wc -c | tr -d ' ')
t1=$(python3 -c 'import time;print(time.time())')
wait $pids
printf 'paste over %d FIFOs fed by bgzip -dc: %.1f MB in %.2f s -> %.0f MB/s\n' "$i" "$(echo "$bytes/1e6" | bc -l)" "$(echo "$t1-$t0" | bc -l)" "$(echo "$bytes/1e6/($t1-$t0)" | bc -l)"
# per-sample extraction on a genome-sized mosdepth file (3.1M bins)
awk 'BEGIN{srand(5); for(i=0;i<3100000;i++) printf "chr1\t%d\t%d\t%.2f\n", i*1000, (i+1)*1000, 20+rand()*30}' > genome.bed
bgzip -f -@ 4 genome.bed
printf 'genome-sized mosdepth file: %s bytes bgz (%.1f MB)\n' "$(stat -f %z genome.bed.gz)" "$(echo "$(stat -f %z genome.bed.gz)/1e6" | bc -l)"
/usr/bin/time -p sh -c 'gzip -cd genome.bed.gz | cut -f4 > col.tmp; wc -l < col.tmp > /dev/null' 2>&1 | tr '\n' ' '; echo " (join_extract.sh: gzip -cd | cut -f4, then wc -l)"
/usr/bin/time -p sh -c 'gzip -cd genome.bed.gz | cut -f1-3 | cksum > /dev/null' 2>&1 | tr '\n' ' '; echo " (--strict-coords second pass)"
/usr/bin/time -p sh -c 'gzip -cd genome.bed.gz > /dev/null' 2>&1 | tr '\n' ' '; echo " (gzip -cd alone)"
/usr/bin/time -p sh -c 'bgzip -dc genome.bed.gz > /dev/null' 2>&1 | tr '\n' ' '; echo " (bgzip -dc alone)"
/usr/bin/time -p sh -c 'bgzip -dc genome.bed.gz | cut -f4 > /dev/null' 2>&1 | tr '\n' ' '; echo " (bgzip -dc | cut -f4)"
tabix -f -p bed genome.bed.gz
/usr/bin/time -p sh -c 'tabix genome.bed.gz chr1:1000000001-1025000000 | wc -l' 2>&1 | tr '\n' ' '; echo " (tabix one 25 Mb window from a per-sample file)"
cd .. && rm -rf pastebench
