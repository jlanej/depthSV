#!/usr/bin/env bash
# Per-region join alternative: cost of one 25 Mb range query against a
# genome-sized per-sample mosdepth file through a csi index (mosdepth writes
# .csi), and GNU cut for the extraction step.
cd "$(dirname "$0")"
mkdir -p csi && cd csi
awk 'BEGIN{srand(5); n=0; for(c=1;c<=22;c++){len=(c<=8?200000000:100000000); for(i=0;i<len/1000;i++) printf "chr%d\t%d\t%d\t%.2f\n", c, i*1000, (i+1)*1000, 20+rand()*30}}' > genome.bed
echo "bins: $(wc -l < genome.bed | tr -d ' ')"
bgzip -f -@ 4 genome.bed
tabix -f -C -p bed genome.bed.gz
printf 'per-sample file: %.1f MB bgz, csi index %s bytes\n' "$(echo "$(stat -f %z genome.bed.gz)/1e6" | bc -l)" "$(stat -f %z genome.bed.gz.csi)"
/usr/bin/time -p sh -c 'tabix genome.bed.gz chr5:100000001-125000000 | wc -l' 2>&1 | tr '\n' ' '; echo " (tabix csi: one 25 Mb window, 25k bins)"
/usr/bin/time -p sh -c 'for r in chr1:1-25000000 chr3:50000001-75000000 chr7:1-25000000 chr12:25000001-50000000 chr20:1-25000000; do tabix genome.bed.gz $r; done | wc -l' 2>&1 | tr '\n' ' '; echo " (5 windows)"
/usr/bin/time -p sh -c 'gzip -cd genome.bed.gz | gcut -f4 | wc -l' 2>&1 | tr '\n' ' '; echo " (gzip -cd | gcut -f4  [GNU cut])"
/usr/bin/time -p sh -c 'gzip -cd genome.bed.gz | cut -f4 | wc -l' 2>&1 | tr '\n' ' '; echo " (gzip -cd | cut -f4  [BSD cut])"
/usr/bin/time -p sh -c 'gzip -cd genome.bed.gz | gawk -F"\t" "{print \$4}" | wc -l' 2>&1 | tr '\n' ' '; echo " (gzip -cd | gawk)"
/usr/bin/time -p sh -c 'gzip -cd genome.bed.gz | gawk -F"\t" "{print \$4 > \"/dev/null\"; print \$1\"\t\"\$2\"\t\"\$3}" | cksum' 2>&1 | tr '\n' ' '; echo " (one gawk pass: column + coordinate checksum)"
cd .. && rm -rf csi
