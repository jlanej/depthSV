#!/usr/bin/env bash
# How many jobs does GNU parallel --pipe make when one -L record is 60% of --block?
# (the real example: -L 2000 rows x ~19 KB = 38 MB records against a 64 MiB block)
echo "parallel version: $(parallel --version | head -1)"
# rows of 6 bytes ("12345\n" .. ), record = 2000 rows = 12000 B, block = 1.67 x record = 20040 B
n=$(seq -w 10000 34999 | parallel --pipe -k --block 20040 -L 2000 'wc -l' | wc -l | tr -d ' ')
echo "25,000 rows, record 12000 B, block 20040 B (1.67 records): jobs = $n  (13 => one record per job; 7 => two)"
seq -w 10000 34999 | parallel --pipe -k --block 20040 -L 2000 'wc -l' | tr '\n' ' '; echo
# and the real ratio: 64 MiB block vs 2000 rows x 19,232 B (3,202 x 6 B + 20)
python3 -c "print('records per 64 MiB block at 19.2 KB/row: %.2f -> floor 1' % (67108864/(2000*19232)))"
