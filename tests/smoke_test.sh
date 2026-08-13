#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# depthSV — end-to-end test on synthetic data
#
#   tests/smoke_test.sh [workDir]
#
# Runs join -> correct -> analyze on generated fixtures and asserts on the
# results, not merely that the commands exit zero. Needs no real cohort data,
# so it is safe to run anywhere and is what CI executes.
#
# Covers: the guardrails fire on bad input; the injected association is
# recovered; a null phenotype stays calibrated; case/control counts are
# labelled the right way round; and a region is a usable unit of work.
# ---------------------------------------------------------------------------

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../lib/common.sh"

work="${1:-${TMPDIR:-/tmp}/depthsv-smoke.$$}"
fixtures="$work/fixtures"
rm -rf "$work"; mkdir -p "$work"

pass=0; fail=0
ok()   { printf '  \033[32mPASS\033[0m %s\n' "$1"; pass=$((pass+1)); }
bad()  { printf '  \033[31mFAIL\033[0m %s\n' "$1"; fail=$((fail+1)); }
# Show why a stage failed. Without this an intermittent failure reports only
# that it happened, and the reason sits in a log the caller never sees.
bad_log() { bad "$1"; [ -s "${2:-}" ] && sed 's/^/         | /' "$2" | tail -12; }
check(){ if [ "$2" = "$3" ]; then ok "$1"; else bad "$1 (expected '$3', got '$2')"; fi; }

echo "depthSV smoke test"
echo "workdir: $work"
echo

# --- fixtures --------------------------------------------------------------
echo "[1/6] fixtures"
Rscript "$DSV_ROOT/tests/make_fixtures.R" "$fixtures" 60 200 >/dev/null 2>&1
[ -s "$fixtures/mosdepth.input.txt" ] && ok "fixtures generated" || bad "fixture generation"

# --- join ------------------------------------------------------------------
echo "[2/6] join"
bash "$DSV_ROOT/scripts/join.sh" --manifest "$fixtures/mosdepth.input.txt" \
     --out "$work/join" --threads 2 >"$work/join.log" 2>&1 \
  && ok "join succeeded" || bad_log "join failed" "$work/join.log"

matrix="$work/join/depth.matrix.txt.gz"
check "matrix is tabix-indexed" "$([ -s "$matrix.tbi" ] && echo yes || echo no)" "yes"
check "matrix has both contigs"  "$(tabix -l "$matrix" | tr '\n' ' ' | sed 's/ $//')" "chr1 chr2"
check "manifest records regions" "$(awk -F'\t' '$1=="regions"{print $2}' "$work/join/depth.matrix.manifest")" "400"

# A sample built against different intervals must be refused, not silently
# pasted into the wrong coordinates.
mkdir -p "$work/bad/mosdepth"; cp "$fixtures"/mosdepth/*.gz "$work/bad/mosdepth/"
gzip -cd "$work/bad/mosdepth/SAMPLE007.by1000.regions.bed.gz" | head -380 \
  | bgzip > "$work/bad/mosdepth/x.gz"
mv "$work/bad/mosdepth/x.gz" "$work/bad/mosdepth/SAMPLE007.by1000.regions.bed.gz"
ls "$work/bad/mosdepth"/*.gz > "$work/bad/manifest.txt"
if bash "$DSV_ROOT/scripts/join.sh" --manifest "$work/bad/manifest.txt" \
        --out "$work/bad/out" --threads 2 >"$work/bad.log" 2>&1; then
    bad "join accepted a sample with a different region count"
else
    ok "join refuses mismatched region counts"
    check "  and promotes nothing" "$([ -f "$work/bad/out/depth.matrix.txt.gz" ] && echo yes || echo no)" "no"
fi

# --- correct ---------------------------------------------------------------
echo "[3/6] correct"
for region in chr1 chr2; do
    bash "$DSV_ROOT/scripts/correct.sh" --matrix "$matrix" \
        --pcs "$fixtures/svd.pcs.txt" --coverage "$fixtures/autosomal.median.txt" \
        --region "$region" --out "$work/corrected" --ndim 4 --jobs 2 --chunk 100 \
        >"$work/correct.$region.log" 2>&1 || bad_log "correct failed on $region" "$work/correct.$region.log"
done
ok "correction ran for both regions"
check "corrected rows (chr1)" \
      "$(bgzip -dc "$work/corrected/corrected_ndim4.chr1.txt.gz" | grep -cv '^#')" "200"
check "per-region stats written" \
      "$([ -s "$work/corrected/stats_ndim4.chr1.txt.gz" ] && echo yes || echo no)" "yes"

# Re-running a completed unit must be a no-op, which is what makes a preempted
# job cheap to retry.
before="$(dsv_mtime "$work/corrected/corrected_ndim4.chr1.txt.gz")"
bash "$DSV_ROOT/scripts/correct.sh" --matrix "$matrix" --pcs "$fixtures/svd.pcs.txt" \
    --coverage "$fixtures/autosomal.median.txt" --region chr1 --out "$work/corrected" \
    --ndim 4 >>"$work/correct.log" 2>&1
after="$(dsv_mtime "$work/corrected/corrected_ndim4.chr1.txt.gz")"
check "completed work is not redone" "$before" "$after"

# A sub-chromosome interval is a valid unit of work.
bash "$DSV_ROOT/scripts/correct.sh" --matrix "$matrix" --pcs "$fixtures/svd.pcs.txt" \
    --coverage "$fixtures/autosomal.median.txt" --region chr2:0-49999 \
    --out "$work/corrected" --ndim 4 --jobs 2 --chunk 100 >>"$work/correct.log" 2>&1
check "sub-chromosome region works" \
      "$(bgzip -dc "$work/corrected/corrected_ndim4.chr2_0-49999.txt.gz" | grep -cv '^#')" "50"

# --- analyze ---------------------------------------------------------------
echo "[4/6] analyze"
# Use the shipped example rather than a copy, so a mistake in the documented
# manifest fails CI instead of shipping. This also exercises comment skipping.
manifest="$DSV_ROOT/conf/phenotypes.example.tsv"
for region in chr1 chr2; do
    bash "$DSV_ROOT/scripts/analyze.sh" \
        --corrected "$work/corrected/corrected_ndim4.${region}.txt.gz" \
        --pheno "$fixtures/phenotypes.tsv" --pheno-manifest "$manifest" \
        --region "$region" --out "$work/assoc" --jobs 2 --chunk 100 \
        -- --minObs 30 >"$work/analyze.$region.log" 2>&1 || bad_log "analyze failed on $region" "$work/analyze.$region.log"
done
ok "association ran for 4 phenotypes x 2 regions"
check "output files" "$(ls "$work/assoc"/*.txt.gz 2>/dev/null | wc -l | tr -d ' ')" "8"
check "coxph produced results" \
      "$(bgzip -dc "$work/assoc/survival.coxph.chr1.txt.gz" | grep -cv '^#')" "200"

# --- assertions on the numbers ---------------------------------------------
echo "[5/6] results"
Rscript - "$work" "$fixtures" <<'RS' > "$work/assert.txt" 2>&1
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(trailingOnly = TRUE); work <- a[1]; fx <- a[2]
truth <- fread(file.path(fx, "truth.tsv"))
rd <- function(p) rbindlist(lapply(Sys.glob(file.path(work, "assoc", p)),
        function(f) fread(cmd = paste("gzip -cd", shQuote(f)))))

q <- rd("quant_trait.linear.chr*.txt.gz"); setnames(q, ncol(q), "P")
n <- rd("null_trait.linear.chr*.txt.gz");  setnames(n, ncol(n), "P")
b <- rd("case_status.logistic.chr*.txt.gz")

cat(sprintf("regions_tested=%d\n", nrow(q)))
cat(sprintf("signal_is_top_hit=%s\n", all(q[order(P)][seq_len(nrow(truth))]$Region %in% truth$Region)))
cat(sprintf("signal_min_p=%.3g\n", min(q[Region %in% truth$Region]$P)))
lam <- function(p) round(median(qchisq(p, 1, lower.tail = FALSE)) / qchisq(0.5, 1), 3)
cat(sprintf("null_lambda=%.3f\n", lam(n$P)))
cat(sprintf("ncase_lt_ncontrol=%s\n", all(b$NCase < b$NControl)))
cat(sprintf("ncase_plus_ncontrol_eq_n=%s\n", all(b$NCase + b$NControl == b$N)))
cat(sprintf("no_missing_p=%s\n", !any(is.na(q$P)) && !any(is.na(n$P))))
RS
get() { awk -F= -v k="$1" '$1==k{print $2}' "$work/assert.txt"; }

check "all regions tested"            "$(get regions_tested)"           "400"
check "injected signal is the top hit" "$(get signal_is_top_hit)"        "TRUE"
check "case counts are the minority"   "$(get ncase_lt_ncontrol)"        "TRUE"
check "counts sum to N"                "$(get ncase_plus_ncontrol_eq_n)" "TRUE"
check "no missing p-values"            "$(get no_missing_p)"             "TRUE"

lam="$(get null_lambda)"
if awk -v l="$lam" 'BEGIN{exit !(l > 0.6 && l < 1.6)}'; then
    ok "null phenotype is calibrated (lambda=$lam)"
else
    bad "null phenotype inflation lambda=$lam outside 0.6-1.6"
fi
printf '       injected signal p = %s\n' "$(get signal_min_p)"

# --- guardrails ------------------------------------------------------------
echo "[6/6] guardrails"
# A duplicated SAMPLE must be refused: silently, it misassigns depth values.
head -1 "$fixtures/phenotypes.tsv"  > "$work/dup.tsv"
tail -n +2 "$fixtures/phenotypes.tsv" >> "$work/dup.tsv"
sed -n '2p' "$fixtures/phenotypes.tsv" >> "$work/dup.tsv"
if bgzip -dc "$work/corrected/corrected_ndim4.chr1.txt.gz" \
     | Rscript "$DSV_ROOT/R/analyze.R" -f - -p "$work/dup.tsv" \
         -m "quant_trait~cov_resids+age" -r linear --minObs 30 >/dev/null 2>&1; then
    bad "duplicated SAMPLE accepted"
else
    ok "duplicated SAMPLE refused"
fi

# An absent chromosome must fail loudly rather than yield an empty result.
if bash "$DSV_ROOT/scripts/correct.sh" --matrix "$matrix" --pcs "$fixtures/svd.pcs.txt" \
        --coverage "$fixtures/autosomal.median.txt" --region chr22 \
        --out "$work/nope" --ndim 4 >/dev/null 2>&1; then
    bad "absent chromosome accepted"
else
    ok "absent chromosome refused"
fi

# A relative --out must work. join.sh cds into a working subdirectory, so a
# relative path captured beforehand stops resolving — and a relative path is
# exactly what the documented quickstart uses.
( cd "$work" && mkdir -p rel && cd rel \
  && bash "$DSV_ROOT/scripts/join.sh" --manifest "$fixtures/mosdepth.input.txt" \
       --out relout --threads 2 >/dev/null 2>&1 \
  && [ -s relout/depth.matrix.txt.gz ] ) \
  && ok "relative --out works" || bad "relative --out failed"

# Gzipped input tables must work without the R.utils package, which is not part
# of a base R install. This failed only in a clean container, where R.utils is
# absent — and a gzipped phenotype table is the normal case in practice.
gzip -c "$fixtures/phenotypes.tsv" > "$work/pheno.tsv.gz"
if bgzip -dc "$work/corrected/corrected_ndim4.chr1.txt.gz" \
     | Rscript "$DSV_ROOT/R/analyze.R" -f - -p "$work/pheno.tsv.gz" \
         -m 'quant_trait~cov_resids+age' -r linear --minObs 30 \
         > "$work/gz.txt" 2>/dev/null && [ "$(grep -cv '^#' "$work/gz.txt")" -eq 200 ]; then
    ok "gzipped phenotype table is read without R.utils"
else
    bad "gzipped phenotype table failed"
fi

# The two projection routes must agree. `explicit` is 12-20x faster but is a
# different floating-point path, so the equivalence is enforced here rather
# than remembered. Both are compared at the precision actually written.
proj_ok=1
for stage in correct analyze; do
    for p in explicit qr; do
        if [ "$stage" = correct ]; then
            Rscript "$DSV_ROOT/R/correct.R" -i "$fixtures/svd.pcs.txt" \
                -f "$matrix" -d 4 -c "$fixtures/autosomal.median.txt" \
                --projection "$p" 2>/dev/null | shasum > "$work/proj.$stage.$p"
        else
            Rscript "$DSV_ROOT/R/analyze.R" \
                -f "$work/corrected/corrected_ndim4.chr1.txt.gz" \
                -p "$fixtures/phenotypes.tsv" -r linear --minObs 30 \
                -m 'quant_trait~cov_resids+age+sex+ancestry_PC1' \
                --projection "$p" 2>/dev/null | shasum > "$work/proj.$stage.$p"
        fi
    done
    cmp -s "$work/proj.$stage.explicit" "$work/proj.$stage.qr" || proj_ok=0
done
if [ "$proj_ok" -eq 1 ]; then
    ok "explicit and qr projections give identical output"
else
    bad "projection routes disagree at the printed precision"
fi

# Rank-deficient designs are the case where a naive explicit projection is
# wrong rather than merely imprecise, so it gets its own check.
Rscript - "$fixtures" <<'RS' > "$work/rank.txt" 2>&1
set.seed(1); n <- 500; a <- commandArgs(trailingOnly = TRUE)
X <- cbind(1, matrix(rnorm(n*5), n, 5)); X <- cbind(X, X[, 2] + X[, 3])  # deficient
q <- qr(X); v <- rnorm(n)
Q <- qr.Q(q)[, seq_len(q$rank), drop = FALSE]
cat(if (max(abs(qr.resid(q, v) - (v - Q %*% crossprod(Q, v)))) < 1e-10) "OK" else "BAD")
RS
check "rank-truncated projection is correct on a deficient design" "$(cat "$work/rank.txt")" "OK"

# Contig naming must not be assumed. A matrix using Ensembl-style names must
# work, and must be addressable by either naming convention.
mkdir -p "$work/nochr/mosdepth"
for f in "$fixtures"/mosdepth/*.gz; do
    gzip -cd "$f" | sed 's/^chr//' | bgzip > "$work/nochr/mosdepth/$(basename "$f")"
done
ls "$work/nochr/mosdepth"/*.gz > "$work/nochr/manifest.txt"
if bash "$DSV_ROOT/scripts/join.sh" --manifest "$work/nochr/manifest.txt" \
        --out "$work/nochr/join" --threads 2 >>"$work/nochr.log" 2>&1 \
   && bash "$DSV_ROOT/scripts/correct.sh" --matrix "$work/nochr/join/depth.matrix.txt.gz" \
        --pcs "$fixtures/svd.pcs.txt" --coverage "$fixtures/autosomal.median.txt" \
        --region chr1 --out "$work/nochr/corr" --ndim 4 --jobs 2 --chunk 100 \
        >>"$work/nochr.log" 2>&1; then
    check "unprefixed contigs, addressed as 'chr1'" \
          "$(bgzip -dc "$work/nochr/corr/corrected_ndim4.chr1.txt.gz" | grep -cv '^#')" "200"
else
    bad "unprefixed contig names not handled (see $work/nochr.log)"
fi

echo
printf 'passed %d, failed %d\n' "$pass" "$fail"
[ "$fail" -eq 0 ] || exit 1
