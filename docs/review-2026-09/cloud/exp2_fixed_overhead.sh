#!/usr/bin/env bash
# Experiment 2: per-process fixed cost of the R drivers at cohort width.
# correct.R: reads the PC table + coverage, QR of an n x (ndim+1) design, qr.Q.
# analyze.R: reads the phenotype table, builds the design, QR.
# Each GNU-parallel --pipe block starts a fresh Rscript, so this cost is paid
# once per block (see exp1 for rows per block).
set -u
here="$(cd "$(dirname "$0")" && pwd)"
ROOT=/Users/Kitty/git/depthSV/.claude/worktrees/agent-a8beeeb9b1a9d1356
N="${1:-100000}"       # samples
NPC="${2:-200}"        # PCs in the table
ROWS="${3:-20}"        # matrix rows per process (exp1: ~21 at 500k)
cd "$here" || exit 1

echo "== generating synthetic tables: N=$N samples, $NPC PCs, $ROWS rows =="
Rscript - "$N" "$NPC" "$ROWS" <<'RS'
a <- commandArgs(TRUE); n <- as.integer(a[1]); npc <- as.integer(a[2]); rows <- as.integer(a[3])
set.seed(1)
ids <- sprintf("S%07d", seq_len(n))
pcs <- matrix(rnorm(n * npc), n, npc); colnames(pcs) <- paste0("PC", seq_len(npc))
pc_dt <- data.frame(SAMPLE = ids, pcs, check.names = FALSE)
write.table(pc_dt, "pcs.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(SAMPLE = ids, AUTO_HQ_median = round(runif(n, 25, 40), 2)),
            "cov.txt", sep = "\t", quote = FALSE, row.names = FALSE)
ph <- data.frame(SAMPLE = ids, y = rnorm(n), age = round(runif(n, 40, 70)), sex = sample(0:1, n, TRUE))
write.table(ph, "pheno.txt", sep = "\t", quote = FALSE, row.names = FALSE)
# raw matrix chunk (correct input) and corrected chunk (analyze input)
con <- file("raw_chunk.txt", "w")
writeLines(paste(c("#CHR", "START", "STOP", ids), collapse = "\t"), con)
for (i in seq_len(rows)) {
  s <- (i - 1) * 1000
  writeLines(paste(c("chr1", s, s + 1000, sprintf("%.2f", runif(n, 20, 45))), collapse = "\t"), con)
}
close(con)
con <- file("corr_chunk.txt", "w")
writeLines(paste(c("#CHROM", "START", "END", "Region", ids), collapse = "\t"), con)
for (i in seq_len(rows)) {
  s <- (i - 1) * 1000
  writeLines(paste(c("chr1", s, s + 1000, sprintf("chr1:%d-%d", s, s + 1000), sprintf("%.6g", rnorm(n, 0, 0.2))), collapse = "\t"), con)
}
close(con)
RS
ls -la pcs.txt cov.txt pheno.txt raw_chunk.txt corr_chunk.txt | awk '{print $5, $9}'

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
tm() { local t0 t1; t0=$(date +%s.%N); "$@" > /dev/null 2> "$here/err.log"; local rc=$?; t1=$(date +%s.%N); printf '%.1f s (exit %d)\n' "$(echo "$t1 - $t0" | bc)" "$rc"; }

echo
echo "== correct.R, $ROWS rows, one process =="
for nd in 0 42 "$NPC"; do
    printf 'ndim=%-4s ' "$nd"
    tm Rscript "$ROOT/R/correct.R" --inputPCs pcs.txt --inputFile raw_chunk.txt --coverageStats cov.txt --ndim "$nd"
done
/usr/bin/time -l Rscript "$ROOT/R/correct.R" --inputPCs pcs.txt --inputFile raw_chunk.txt --coverageStats cov.txt --ndim "$NPC" > /dev/null 2> time.log
echo "peak RSS at ndim=$NPC: $(awk '/maximum resident set size/ {printf "%.2f GB", $1/1024/1024/1024}' time.log)"

echo
echo "== analyze.R, $ROWS rows, one process =="
printf 'linear   ' ; tm Rscript "$ROOT/R/analyze.R" -f corr_chunk.txt -p pheno.txt -m 'y~cov_resids+age+sex' -r linear
printf 'logistic ' ; tm Rscript "$ROOT/R/analyze.R" -f corr_chunk.txt -p pheno.txt -m 'sex~cov_resids+age' -r logistic
echo
echo "== R startup + package load only =="
printf 'baseline ' ; tm Rscript -e 'suppressPackageStartupMessages({library(optparse); library(data.table)})'
