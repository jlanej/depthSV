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
# recovered with the estimate lm() and coxph() give; a null phenotype stays
# calibrated, including one driven by the PC the correction removed; a
# binary phenotype keeps its direction under either coding; a region is a
# usable unit of work, including one wide enough to need ten chunks.
# ---------------------------------------------------------------------------

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../lib/common.sh"

work="${1:-${TMPDIR:-/tmp}/depthsv-smoke.$$}"
case "$work" in /|"$HOME"|"") echo "refusing to use '$work' as a work directory" >&2; exit 2 ;; esac
fixtures="$work/fixtures"
rm -rf "$work"; mkdir -p "$work"

pass=0; fail=0
ok()   { printf '  \033[32mPASS\033[0m %s\n' "$1"; pass=$((pass+1)); }
bad()  { printf '  \033[31mFAIL\033[0m %s\n' "$1"; fail=$((fail+1)); }
# Show why a stage failed. Without this an intermittent failure reports only
# that it happened, and the reason sits in a log the caller never sees.
# Returns 0 so it can end a `||` list under set -e.
bad_log() { bad "$1"; [ -s "${2:-}" ] && sed 's/^/         | /' "$2" | tail -12; return 0; }
check(){ if [ "$2" = "$3" ]; then ok "$1"; else bad "$1 (expected '$3', got '$2')"; fi; }

echo "depthSV smoke test"
echo "workdir: $work"
echo

# --- fixtures --------------------------------------------------------------
echo "[1/6] fixtures"
Rscript "$DSV_ROOT/tests/make_fixtures.R" "$fixtures" 60 200 >/dev/null 2>&1
[ -s "$fixtures/mosdepth.input.txt" ] && ok "fixtures generated" || bad "fixture generation"
pcs="$fixtures/svd.pcs.txt"; cov="$fixtures/autosomal.median.txt"; pheno="$fixtures/phenotypes.tsv"

# --- join ------------------------------------------------------------------
echo "[2/6] join"
bash "$DSV_ROOT/scripts/join.sh" --manifest "$fixtures/mosdepth.input.txt" \
     --out "$work/join" --threads 2 >"$work/join.log" 2>&1 \
  && ok "join succeeded" || bad_log "join failed" "$work/join.log"

matrix="$work/join/depth.matrix.txt.gz"
check "matrix is tabix-indexed" "$([ -s "$matrix.tbi" ] && echo yes || echo no)" "yes"
check "matrix has both contigs"  "$(tabix -l "$matrix" | tr '\n' ' ' | sed 's/ $//')" "chr1 chr2"
check "manifest records regions" "$(awk -F'\t' '$1=="regions"{print $2}' "$work/join/depth.matrix.manifest")" "400"

# Re-running a completed join must be a no-op.
before="$(dsv_mtime "$matrix")"
bash "$DSV_ROOT/scripts/join.sh" --manifest "$fixtures/mosdepth.input.txt" \
     --out "$work/join" --threads 2 >>"$work/join.log" 2>&1
after="$(dsv_mtime "$matrix")"
check "completed join is not redone" "$before" "$after"

# The matrix must not depend on how the work was batched or parallelised.
bash "$DSV_ROOT/scripts/join.sh" --manifest "$fixtures/mosdepth.input.txt" \
     --out "$work/join.alt" --batch-size 7 --jobs 3 --threads 2 >"$work/join.alt.log" 2>&1 \
  && cmp -s "$work/join.alt/depth.matrix.txt.gz" "$matrix" \
  && ok "matrix is byte-identical across batching and job counts" \
  || bad_log "matrix depends on batching" "$work/join.alt.log"

# A windowed region list must partition the matrix — every bin in exactly one
# window — including when the requested window is not a multiple of the bin
# size (it is rounded up to a bin edge).
printf 'chr1\t200000\nchr2\t200000\n' > "$work/chrom.sizes"
bash "$DSV_ROOT/scripts/regions.sh" --matrix "$matrix" --window 2500 \
     --sizes "$work/chrom.sizes" > "$work/regions.win.txt" 2>/dev/null
win_rows=0
: > "$work/win.keys"
while IFS= read -r r; do
    win_rows=$(( win_rows + $(tabix "$matrix" "$r" | wc -l) ))
    tabix "$matrix" "$r" | cut -f1,2 >> "$work/win.keys"
done < "$work/regions.win.txt"
check "windowed regions cover every bin exactly once" "$win_rows" "400"
check "  and contain no duplicate bins" "$(sort "$work/win.keys" | uniq -d | wc -l | tr -d ' ')" "0"

# The same sample listed twice must be refused before any pasting: it would
# become two columns under one name, and later stages match columns by name.
{ cat "$fixtures/mosdepth.input.txt"; head -1 "$fixtures/mosdepth.input.txt"; } > "$work/dupman.txt"
if bash "$DSV_ROOT/scripts/join.sh" --manifest "$work/dupman.txt" \
        --out "$work/dupjoin" >/dev/null 2>&1; then
    bad "join accepted a manifest listing the same sample twice"
else
    ok "join refuses a duplicated sample in the manifest"
fi

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
        --pcs "$pcs" --coverage "$cov" \
        --region "$region" --out "$work/corrected" --ndim 4 --jobs 2 --chunk 100 \
        >"$work/correct.$region.log" 2>&1 || bad_log "correct failed on $region" "$work/correct.$region.log"
done
ok "correction ran for both regions"
check "corrected rows (chr1)" \
      "$(bgzip -dc "$work/corrected/corrected_ndim4.chr1.txt.gz" | grep -cv '^#')" "200"
check "per-region stats written" \
      "$([ -s "$work/corrected/stats_ndim4.chr1.txt.gz" ] && echo yes || echo no)" "yes"
check "worker log written" \
      "$([ -f "$work/corrected/corrected_ndim4.chr1.log" ] && echo yes || echo no)" "yes"

# Re-running a completed unit must be a no-op, which is what makes a preempted
# job cheap to retry.
before="$(dsv_mtime "$work/corrected/corrected_ndim4.chr1.txt.gz")"
bash "$DSV_ROOT/scripts/correct.sh" --matrix "$matrix" --pcs "$pcs" \
    --coverage "$cov" --region chr1 --out "$work/corrected" \
    --ndim 4 >>"$work/correct.log" 2>&1
after="$(dsv_mtime "$work/corrected/corrected_ndim4.chr1.txt.gz")"
check "completed work is not redone" "$before" "$after"

# A sub-chromosome interval is a valid unit of work (1-based, inclusive).
bash "$DSV_ROOT/scripts/correct.sh" --matrix "$matrix" --pcs "$pcs" \
    --coverage "$cov" --region chr2:1-50000 \
    --out "$work/corrected" --ndim 4 --jobs 2 --chunk 100 >>"$work/correct.log" 2>&1
check "sub-chromosome region works" \
      "$(bgzip -dc "$work/corrected/corrected_ndim4.chr2_1-50000.txt.gz" | grep -cv '^#')" "50"

# A unit wide enough to need ten or more chunks must still index its
# statistics: the per-chunk files are merged in coordinate order, not in the
# lexicographic order of their job numbers.
bash "$DSV_ROOT/scripts/correct.sh" --matrix "$matrix" --pcs "$pcs" --coverage "$cov" \
    --region chr1 --out "$work/manychunks" --ndim 4 --jobs 2 --chunk 20 \
    >"$work/manychunks.log" 2>&1 \
  && [ -s "$work/manychunks/stats_ndim4.chr1.txt.gz.tbi" ] \
  && cmp -s <(bgzip -dc "$work/manychunks/corrected_ndim4.chr1.txt.gz") \
            <(bgzip -dc "$work/corrected/corrected_ndim4.chr1.txt.gz") \
  && ok "a ten-chunk unit indexes its stats and matches the two-chunk result" \
  || bad_log "ten-chunk unit failed" "$work/manychunks.log"

# A duplicated or missing coverage median must be refused, not resolved by
# row order or propagated as NA through the projection.
{ cat "$cov"; printf 'SAMPLE001\t1000\n'; } > "$work/cov.dup.txt"
if bash "$DSV_ROOT/scripts/correct.sh" --matrix "$matrix" --pcs "$pcs" --coverage "$work/cov.dup.txt" \
        --region chr1 --out "$work/covdup" --ndim 4 >/dev/null 2>&1; then
    bad "duplicated coverage row accepted"
else
    ok "duplicated coverage row refused"
fi
awk -F'\t' 'BEGIN{OFS="\t"} NR==2 {$2="NA"} {print}' "$cov" > "$work/cov.na.txt"
if bash "$DSV_ROOT/scripts/correct.sh" --matrix "$matrix" --pcs "$pcs" --coverage "$work/cov.na.txt" \
        --region chr1 --out "$work/covna" --ndim 4 >/dev/null 2>&1; then
    bad "NA coverage median accepted"
else
    ok "NA coverage median refused"
fi

# --- analyze ---------------------------------------------------------------
echo "[4/6] analyze"
# Use the shipped example rather than a copy, so a mistake in the documented
# manifest fails CI instead of shipping. This also exercises comment skipping.
manifest="$DSV_ROOT/conf/phenotypes.example.tsv"
for region in chr1 chr2; do
    bash "$DSV_ROOT/scripts/analyze.sh" \
        --corrected "$work/corrected/corrected_ndim4.${region}.txt.gz" \
        --pheno "$pheno" --pheno-manifest "$manifest" --pcs "$pcs" \
        --case-level case --min-cases 5 \
        --region "$region" --out "$work/assoc" --jobs 2 --chunk 100 \
        -- --minObs 30 >"$work/analyze.$region.log" 2>&1 || bad_log "analyze failed on $region" "$work/analyze.$region.log"
done
ok "association ran for 6 phenotypes x 2 regions"
check "output files" "$(ls "$work/assoc"/*.txt.gz 2>/dev/null | wc -l | tr -d ' ')" "12"
check "coxph produced results" \
      "$(bgzip -dc "$work/assoc/survival.coxph.chr1.txt.gz" | grep -cv '^#')" "200"

# --- assertions on the numbers ---------------------------------------------
echo "[5/6] results"
Rscript - "$work" "$fixtures" <<'RS' > "$work/assert.txt" 2>&1
suppressPackageStartupMessages({ library(data.table); library(survival) })
a <- commandArgs(trailingOnly = TRUE); work <- a[1]; fx <- a[2]
truth <- fread(file.path(fx, "truth.tsv"))
rd <- function(p) rbindlist(lapply(Sys.glob(file.path(work, "assoc", p)),
        function(f) fread(cmd = paste("gzip -cd", shQuote(f)))))

q  <- rd("quant_trait.linear.chr*.txt.gz")
n  <- rd("null_trait.linear.chr*.txt.gz")
pn <- rd("pc_null.linear.chr*.txt.gz")
b  <- rd("case_status.logistic.chr*.txt.gz")
bl <- rd("case_label.logistic.chr*.txt.gz")
cx <- rd("survival.coxph.chr*.txt.gz")

cat(sprintf("regions_tested=%d\n", nrow(q)))
cat(sprintf("signal_is_top_hit=%s\n", all(q[order(P)][seq_len(nrow(truth))]$Region %in% truth$Region)))
cat(sprintf("signal_min_p=%.3g\n", min(q[Region %in% truth$Region]$P)))
cat(sprintf("signal_direction_negative=%s\n", all(q[Region %in% truth$Region]$BETA < 0)))
lam <- function(p) round(median(qchisq(p, 1, lower.tail = FALSE)) / qchisq(0.5, 1), 3)
cat(sprintf("null_lambda=%.3f\n", lam(n$P)))
cat(sprintf("pc_null_lambda=%.3f\n", lam(pn$P)))
cat(sprintf("ncase_lt_ncontrol=%s\n", all(b$NCase < b$NControl)))
cat(sprintf("ncase_plus_ncontrol_eq_n=%s\n", all(b$NCase + b$NControl == b$N)))
cat(sprintf("no_missing_p=%s\n", !any(is.na(q$P)) && !any(is.na(n$P))))
cat(sprintf("log10p_matches_p=%s\n", all(abs(q$LOG10P + log10(q$P)) < 1e-6)))

# The text-coded binary phenotype must give the same test as the 0/1 one:
# same case count, same sign, same statistic.
setkey(b, Region); setkey(bl, Region)
cat(sprintf("case_label_same_ncase=%s\n", all(bl$NCase == b[bl$Region]$NCase)))
cat(sprintf("case_label_same_stat=%s\n", max(abs(bl$STAT - b[bl$Region]$STAT)) < 1e-6))
cat(sprintf("logistic_converged_flag_present=%s\n", all(bl$CONVERGED %in% c(0, 1))))

# Reference values: the pipeline's linear and Cox estimates at the injected
# bin must be what lm() and coxph() give on the same corrected depth with the
# PCs in the model.
ph  <- fread(file.path(fx, "phenotypes.tsv"))
pcs <- fread(file.path(fx, "svd.pcs.txt"))
d   <- merge(ph, pcs[, .(SAMPLE, PC1, PC2, PC3, PC4)], by = "SAMPLE")
cm  <- fread(cmd = paste("gzip -cd", shQuote(file.path(work, "corrected", "corrected_ndim4.chr1.txt.gz"))))
row <- cm[Region == truth$Region[1]]
x   <- as.numeric(row[, -(1:4)]); names(x) <- names(cm)[-(1:4)]
d[, cov_resids := x[SAMPLE]]
fit <- lm(quant_trait ~ cov_resids + age + sex + ancestry_PC1 + ancestry_PC2 + PC1 + PC2 + PC3 + PC4, data = d)
ref <- summary(fit)$coefficients["cov_resids", ]
got <- q[Region == truth$Region[1]]
cat(sprintf("linear_beta_matches_lm=%s\n", abs(got$BETA - ref[1]) < 1e-6 * abs(ref[1])))
cat(sprintf("linear_t_matches_lm=%s\n", abs(got$STAT - ref[3]) < 1e-6 * abs(ref[3])))
cfit <- coxph(Surv(time, event) ~ cov_resids + age + sex + PC1 + PC2 + PC3 + PC4, data = d)
cref <- coef(summary(cfit))["cov_resids", ]
cgot <- cx[Region == truth$Region[1]]
cat(sprintf("cox_beta_matches_coxph=%s\n", abs(cgot$BETA - cref[1]) < 1e-5 * abs(cref[1])))
cat(sprintf("cox_z_matches_coxph=%s\n", abs(cgot$STAT - cref[4]) < 1e-5 * abs(cref[4])))
RS
get() { awk -F= -v k="$1" '$1==k{print $2}' "$work/assert.txt"; }

check "all regions tested"            "$(get regions_tested)"           "400"
check "injected signal is the top hit" "$(get signal_is_top_hit)"        "TRUE"
check "injected deletion has a negative effect" "$(get signal_direction_negative)" "TRUE"
check "linear estimate equals lm() with the PCs in the model" "$(get linear_beta_matches_lm)" "TRUE"
check "linear t equals lm()"           "$(get linear_t_matches_lm)"      "TRUE"
check "Cox estimate equals coxph()"    "$(get cox_beta_matches_coxph)"   "TRUE"
check "Cox z equals coxph()"           "$(get cox_z_matches_coxph)"      "TRUE"
check "case counts are the minority"   "$(get ncase_lt_ncontrol)"        "TRUE"
check "counts sum to N"                "$(get ncase_plus_ncontrol_eq_n)" "TRUE"
check "text-coded cases counted the same" "$(get case_label_same_ncase)" "TRUE"
check "text-coded cases give the same statistic" "$(get case_label_same_stat)" "TRUE"
check "logistic convergence flag reported" "$(get logistic_converged_flag_present)" "TRUE"
check "no missing p-values"            "$(get no_missing_p)"             "TRUE"
check "LOG10P matches P"               "$(get log10p_matches_p)"         "TRUE"

lam="$(get null_lambda)"
if awk -v l="$lam" 'BEGIN{exit !(l > 0.75 && l < 1.35)}'; then
    ok "null phenotype is calibrated (lambda=$lam)"
else
    bad "null phenotype inflation lambda=$lam outside 0.75-1.35"
fi
# A phenotype driven by PC1 is deflated to lambda ~ 0.2 unless the PCs are in
# the association model; with them it is a plain null.
lam="$(get pc_null_lambda)"
if awk -v l="$lam" 'BEGIN{exit !(l > 0.7 && l < 1.4)}'; then
    ok "PC-driven null phenotype is calibrated with the PCs in the model (lambda=$lam)"
else
    bad "PC-driven null phenotype lambda=$lam outside 0.7-1.4: the PCs are not being conditioned on"
fi
printf '       injected signal p = %s\n' "$(get signal_min_p)"

# --- guardrails ------------------------------------------------------------
echo "[6/6] guardrails"
# A duplicated SAMPLE must be refused: silently, it misassigns depth values.
head -1 "$pheno"  > "$work/dup.tsv"
tail -n +2 "$pheno" >> "$work/dup.tsv"
sed -n '2p' "$pheno" >> "$work/dup.tsv"
if bgzip -dc "$work/corrected/corrected_ndim4.chr1.txt.gz" \
     | Rscript "$DSV_ROOT/R/analyze.R" -f - -p "$work/dup.tsv" \
         -m "quant_trait~cov_resids+age" -r linear --minObs 30 >/dev/null 2>&1; then
    bad "duplicated SAMPLE accepted"
else
    ok "duplicated SAMPLE refused"
fi

# A text-coded binary phenotype without --case-level must be refused rather
# than recoded by sort order.
if bgzip -dc "$work/corrected/corrected_ndim4.chr1.txt.gz" \
     | Rscript "$DSV_ROOT/R/analyze.R" -f - -p "$pheno" \
         -m "case_label~cov_resids+age" -r logistic --minObs 30 --minCases 5 >/dev/null 2>&1; then
    bad "text-coded binary phenotype accepted without --case-level"
else
    ok "text-coded binary phenotype refused without --case-level"
fi

# Too few cases must be refused up front, not reported per bin.
if bgzip -dc "$work/corrected/corrected_ndim4.chr1.txt.gz" \
     | Rscript "$DSV_ROOT/R/analyze.R" -f - -p "$pheno" \
         -m "case_status~cov_resids+age" -r logistic --minObs 30 --minCases 1000 >/dev/null 2>&1; then
    bad "phenotype with too few cases accepted"
else
    ok "phenotype with too few cases refused"
fi

# The analysis stage must not run a correction with ndim > 0 without the PCs.
if bash "$DSV_ROOT/scripts/analyze.sh" --corrected "$work/corrected/corrected_ndim4.chr1.txt.gz" \
        --pheno "$pheno" --model "quant_trait~cov_resids+age" --region chr1 \
        --out "$work/nopcs" >/dev/null 2>&1; then
    bad "analysis of a PC-corrected matrix ran without the PCs"
else
    ok "analysis of a PC-corrected matrix refuses to run without the PCs"
fi

# An absent chromosome must fail loudly rather than yield an empty result.
if bash "$DSV_ROOT/scripts/correct.sh" --matrix "$matrix" --pcs "$pcs" \
        --coverage "$cov" --region chr22 \
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
gzip -c "$pheno" > "$work/pheno.tsv.gz"
if bgzip -dc "$work/corrected/corrected_ndim4.chr1.txt.gz" \
     | Rscript "$DSV_ROOT/R/analyze.R" -f - -p "$work/pheno.tsv.gz" \
         -m 'quant_trait~cov_resids+age' -r linear --minObs 30 --pcs "$pcs" --ndim 4 \
         > "$work/gz.txt" 2>/dev/null && [ "$(grep -cv '^#' "$work/gz.txt")" -eq 200 ]; then
    ok "gzipped phenotype table is read without R.utils"
else
    bad "gzipped phenotype table failed"
fi

# Passthrough arguments full of regex metacharacters must reach the R workers
# intact — they are interpolated into a worker shell command, and an unquoted
# interpolation once made any real --sampleIdPattern a syntax error there.
bash "$DSV_ROOT/scripts/correct.sh" --matrix "$matrix" --pcs "$pcs" \
    --coverage "$cov" --region chr1 --out "$work/pattern" \
    --ndim 4 --jobs 2 --chunk 100 -- --sampleIdPattern '(SAMPLE[0-9]+)' \
    >"$work/pattern.log" 2>&1 \
  && cmp -s <(bgzip -dc "$work/pattern/corrected_ndim4.chr1.txt.gz") \
            <(bgzip -dc "$work/corrected/corrected_ndim4.chr1.txt.gz") \
  && ok "regex --sampleIdPattern passes through to the workers" \
  || bad_log "regex --sampleIdPattern passthrough failed" "$work/pattern.log"

# The two projection routes must agree. `explicit` is 12-20x faster but is a
# different floating-point path, so the equivalence is enforced here rather
# than remembered. Both are compared at the precision actually written.
proj_ok=1
for stage in correct analyze; do
    for p in explicit qr; do
        if [ "$stage" = correct ]; then
            Rscript "$DSV_ROOT/R/correct.R" -i "$pcs" \
                -f "$matrix" -d 4 -c "$cov" \
                --projection "$p" 2>/dev/null | shasum > "$work/proj.$stage.$p"
        else
            Rscript "$DSV_ROOT/R/analyze.R" \
                -f "$work/corrected/corrected_ndim4.chr1.txt.gz" \
                -p "$pheno" -r linear --minObs 30 --pcs "$pcs" --ndim 4 \
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
        --pcs "$pcs" --coverage "$cov" \
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
