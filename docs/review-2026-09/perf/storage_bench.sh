#!/usr/bin/env bash
# Compression ratio and throughput: text (bgzip -l1 / default, zstd) vs
# float32 / int16 binaries (zstd), on the same 40 x 500k values.
cd "$(dirname "$0")"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
Rscript storage_bench.R
sz() { stat -f %z "$1"; }
comp() {  # comp <label> <file> <cmd...>
  local label="$1" f="$2"; shift 2
  local t0 t1 out
  t0=$(python3 -c 'import time;print(time.time())')
  "$@" < "$f" > cmp.tmp
  t1=$(python3 -c 'import time;print(time.time())')
  printf '%-34s %-22s in %9.1f MB -> %9.1f MB  ratio %5.2f  %6.1f s  %6.0f MB/s in\n' "$label" "$*" "$(echo "$(sz "$f")/1e6" | bc -l)" "$(echo "$(sz cmp.tmp)/1e6" | bc -l)" "$(echo "$(sz "$f")/$(sz cmp.tmp)" | bc -l)" "$(echo "$t1-$t0" | bc -l)" "$(echo "$(sz "$f")/1e6/($t1-$t0)" | bc -l)"
}
for f in tables/raw_500000.txt tables/corr_500000.txt; do
  comp "$f" "$f" bgzip -l 1 -@ 1 -c
  comp "$f" "$f" bgzip -l 1 -@ 4 -c
  comp "$f" "$f" bgzip -@ 1 -c
  comp "$f" "$f" bgzip -@ 4 -c
  comp "$f" "$f" zstd -1 -q -c
  comp "$f" "$f" zstd -3 -T4 -q -c
done
for f in tables/raw_500000.f32 tables/raw_500000.u16 tables/corr_500000.f32 tables/corr_500000.i16 tables/corr_500000.f64; do
  comp "$f" "$f" zstd -1 -q -c
  comp "$f" "$f" zstd -3 -q -c
  comp "$f" "$f" bgzip -l 1 -@ 1 -c
done
echo "--- decompression throughput (output MB/s), single thread:"
bgzip -l 1 -@ 4 -c tables/corr_500000.txt > corr.l1.bgz; bgzip -@ 4 -c tables/corr_500000.txt > corr.l6.bgz; zstd -3 -q -c tables/corr_500000.f32 > corr.f32.zst
for x in corr.l1.bgz corr.l6.bgz; do
  t0=$(python3 -c 'import time;print(time.time())'); bgzip -dc -@ 1 "$x" > /dev/null; t1=$(python3 -c 'import time;print(time.time())')
  printf '%-14s bgzip -dc: %.2f s -> %.0f MB/s of text\n' "$x" "$(echo "$t1-$t0" | bc -l)" "$(echo "$(sz tables/corr_500000.txt)/1e6/($t1-$t0)" | bc -l)"
done
t0=$(python3 -c 'import time;print(time.time())'); zstd -dc -q corr.f32.zst > /dev/null; t1=$(python3 -c 'import time;print(time.time())')
printf '%-14s zstd -dc:  %.2f s -> %.0f MB/s of float32 (= %.0f M values/s)\n' corr.f32.zst "$(echo "$t1-$t0" | bc -l)" "$(echo "$(sz tables/corr_500000.f32)/1e6/($t1-$t0)" | bc -l)" "$(echo "$(sz tables/corr_500000.f32)/4e6/($t1-$t0)" | bc -l)"
rm -f cmp.tmp corr.l1.bgz corr.l6.bgz corr.f32.zst
