#!/usr/bin/env bash
cd "$(dirname "$0")"
cat gen_tables.log
ls -la tables/
for N in 3202 50000 200000 500000; do
  printf 'N=%s raw row bytes: %s   corr row bytes: %s\n' "$N" "$(sed -n 2p tables/raw_$N.txt | wc -c | tr -d ' ')" "$(sed -n 2p tables/corr_$N.txt | wc -c | tr -d ' ')"
done
