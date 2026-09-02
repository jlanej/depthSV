#!/usr/bin/env bash
# GNU paste (coreutils gpaste) vs BSD paste, same inputs; and the FIFO fan-in
# with gpaste. What a Linux cluster would actually run.
cd "$(dirname "$0")"
mkdir -p pastegnu && cd pastegnu
mk() { local n="$1" L="$2"; rm -f col_*; awk -v n="$n" -v L="$L" 'BEGIN{srand(3); for(i=1;i<=n;i++){f=sprintf("col_%04d",i); for(j=1;j<=L;j++) printf "%.2f\n", 20+rand()*30 > f; close(f)}}'; }
run_paste() {  # run_paste <ncols> <nlines> <pastecmd>
  mk "$1" "$2"
  local t0 t1 bytes
  t0=$(python3 -c 'import time;print(time.time())')
  bytes=$("$3" col_* | wc -c | tr -d ' ')
  t1=$(python3 -c 'import time;print(time.time())')
  printf '%-6s %4d cols x %8d lines: %7.1f MB out in %6.2f s -> %6.0f MB/s\n' "$3" "$1" "$2" "$(echo "$bytes/1e6" | bc -l)" "$(echo "$t1-$t0" | bc -l)" "$(echo "$bytes/1e6/($t1-$t0)" | bc -l)"
}
run_paste 2 10000000 gpaste
run_paste 57 1000000 gpaste
run_paste 708 100000 gpaste
run_paste 1024 60000 gpaste
# 2-stage as join.sh: 57 batch pastes of 57 cols (3249 cols) then final paste of 57 batch files
mk 708 100000
gpaste col_* > big.txt
rm -f col_*
gsplit -l 12500 -d big.txt part_
for p in part_*; do bgzip -l 1 -@ 2 "$p"; done
mkdir -p fifos; rm -f fifos/*
i=0; pids=""
for p in part_*.gz; do i=$((i+1)); mkfifo "fifos/f$i"; bgzip -dc "$p" > "fifos/f$i" & pids="$pids $!"; done
t0=$(python3 -c 'import time;print(time.time())')
bytes=$(gpaste fifos/* | wc -c | tr -d ' ')
t1=$(python3 -c 'import time;print(time.time())')
wait $pids
printf 'gpaste over %d FIFOs fed by bgzip -dc: %.1f MB in %.2f s -> %.0f MB/s\n' "$i" "$(echo "$bytes/1e6" | bc -l)" "$(echo "$t1-$t0" | bc -l)" "$(echo "$bytes/1e6/($t1-$t0)" | bc -l)"
# gpaste | bgzip -@4 -l1 and default: the batch and final stages as written
t0=$(python3 -c 'import time;print(time.time())'); gpaste fifos/* > /dev/null 2>&1 &
for p in part_*.gz; do :; done
rm -rf fifos
sz=$(stat -f %z big.txt)
for opt in "-l 1 -@ 4" "-@ 4" "-@ 8"; do
  t0=$(python3 -c 'import time;print(time.time())'); bgzip $opt -c big.txt > big.bgz; t1=$(python3 -c 'import time;print(time.time())')
  printf 'bgzip %-10s on paste output: %6.1f MB -> %6.1f MB (ratio %.2f) in %5.2f s -> %5.0f MB/s in\n' "$opt" "$(echo "$sz/1e6" | bc -l)" "$(echo "$(stat -f %z big.bgz)/1e6" | bc -l)" "$(echo "$sz/$(stat -f %z big.bgz)" | bc -l)" "$(echo "$t1-$t0" | bc -l)" "$(echo "$sz/1e6/($t1-$t0)" | bc -l)"
done
cd .. && rm -rf pastegnu
