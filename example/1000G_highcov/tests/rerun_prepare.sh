#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# depthSV 1000G example — rerunning stages 0-1 resolves the run afresh
#
#   bash tests/rerun_prepare.sh <smokeWorkDir> [scratchDir]
#
# Stage 1 freezes ndim, the covariates and the thresholds into
# inputs/run.env, and every later job loads that freeze, so a preamble
# finishing mid-run cannot change ndim under the array. Stages 0 and 1 are
# the ones that MAKE the freeze, so a rerun of them must resolve from the
# environment and the preamble's files, not from the previous freeze. This
# replays the sequence that once went wrong: a prepare before the preamble
# (ndim 20, no covariates), the preamble's files landing, a second prepare —
# which must freeze the preamble's ndim and write the adjusted manifest.
# Then: an unchanged rerun rewrites no table; a later stage keeps the frozen
# values when the preamble changes again; an explicit environment still
# wins; and a stage-0 resolution that switches an upstream table says so.
#
# A real (non-smoke) configuration on a mock upstream tree assembled, read
# only, from a finished smoke run: its cached NGS-PCA tables and its
# simulated standard mosdepth tree. No network. Exit status = the number of
# failed checks.
# ---------------------------------------------------------------------------

set -uo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
smoke="${1:?usage: rerun_prepare.sh <smokeWorkDir> [scratchDir]}"
T="${2:-${TMPDIR:-/tmp}/depthsv-rerun-prepare.$$}"
case "$T" in /|"$HOME"|"") echo "refusing to use '$T' as a scratch directory" >&2; exit 2 ;; esac

pcs="$smoke/github_cache/standard/svd.pcs.txt"
qc="$smoke/github_cache/standard/sample_qc.tsv"
tree="$smoke/smoke_mosdepth/standard"
for f in "$pcs" "$qc"; do
    [ -s "$f" ] || { echo "no $f: run the smoke example into $smoke first" >&2; exit 2; }
done
ls "$tree"/*.regions.bed.gz >/dev/null 2>&1 || { echo "no simulated tree under $tree" >&2; exit 2; }

rm -rf "$T"; mkdir -p "$T"
T="$(cd "$T" && pwd)"

# The mock upstream, in the NGS-PCA example's layout.
up="$T/ngspca"
mkdir -p "$up/ngspca_output" "$up/qc_output" "$up/mosdepth_output"
ln -s "$pcs" "$up/ngspca_output/svd.pcs.txt"
ln -s "$qc"  "$up/qc_output/sample_qc.tsv"
for f in "$tree"/*.regions.bed.gz; do ln -s "$f" "$up/mosdepth_output/"; done

# Nothing from the caller's configuration: one mode, this upstream, a fresh
# work directory, and config.sh's defaults for everything else.
for v in $(compgen -v | grep '^EX_' || true); do unset "$v"; done
export NGSPCA_WORK_DIR="$up" EX_WORK_DIR="$T/work" EX_MODES=standard EX_RUNNER=local
W="$EX_WORK_DIR"; in_dir="$W/inputs/standard"

pass=0; fail=0
ok()    { pass=$((pass + 1)); printf '  PASS %s\n' "$1"; }
bad()   { fail=$((fail + 1)); printf '  FAIL %s\n' "$1"
          [ -z "${2:-}" ] || [ ! -s "$2" ] || tail -n 12 "$2" | sed 's/^/       | /'; }
check() { if [ "$2" = "$3" ]; then ok "$1"; else bad "$1 (expected '$3', got '$2')" "${4:-}"; fi; }

# prepare <log> [VAR=value ...]: stages 0 and 1, as run.sh runs them.
prepare() { local log="$1"; shift; env "$@" bash "$here/run.sh" --prepare-only > "$log" 2>&1; }
# The value inputs/run.env froze for a variable.
frozen()  { ( unset "$1"; . "$W/inputs/run.env"; printf '%s' "${!1-}" ); }
# What a stage after the freeze (02-06) resolves: lib.sh sourced as they
# source it. later <VAR> [VAR=value ...]
later()   { env "${@:2}" bash -c 'source "$1/lib.sh" >/dev/null 2>&1 && printf "%s" "${!2-}"' _ "$here" "$1"; }
inode()   { ls -i "$1" 2>/dev/null | awk '{print $1}'; }
adj="SEX+GPC1+GPC2+GPC3+GPC4+GPC5+GPC6+GPC7+GPC8+GPC9+GPC10"

echo "rerunning stages 0-1 -> $T"

# --- 1. a prepare before the preamble ----------------------------------------
prepare "$T/1.log" || bad "first prepare exited non-zero" "$T/1.log"
check "before the preamble: ndim 20 frozen" "$(frozen EX_NDIM)" "20" "$T/1.log"
check "  and no covariates"                  "$(frozen EX_COVARIATES)" "none"
check "  and an unadjusted manifest" "$(grep -c '_adj' "$in_dir/analyses.tsv" 2>/dev/null)" "0"

# --- 2. the preamble lands; the rerun must take its ndim and covariates ------
mkdir -p "$W/preamble"
echo 7 > "$W/preamble/ndim.txt"
awk -F'\t' '
  NR == 1 { for (i = 1; i <= NF; i++) if ($i == "SAMPLE_ID") c = i
            printf "SAMPLE"; for (k = 1; k <= 10; k++) printf "\tGPC%d", k; print ""; next }
  { id = $c; gsub(/^ +| +$/, "", id); printf "%s", id
    for (k = 1; k <= 10; k++) printf "\t%.4f", sin(NR * k); print "" }' "$qc" > "$W/preamble/covariates.tsv"
prepare "$T/2.log" || bad "prepare after the preamble exited non-zero" "$T/2.log"
check "after the preamble: its ndim is frozen" "$(frozen EX_NDIM)" "7" "$T/2.log"
check "  and its covariates" "$(frozen EX_COVARIATES)" "$adj"
check "  the manifest is adjusted" \
      "$(awk -F'\t' '$1 == "mtdna_cn_adj" {print $3}' "$in_dir/analyses.tsv")" "MTDNA_CN~cov_resids+$adj"
check "  the genotype PCs are in the phenotype table" \
      "$(head -n 1 "$in_dir/phenotypes.tsv" | tr '\t' '\n' | grep -c '^GPC[0-9]*$')" "10"
check "  run.sh reports where ndim came from" \
      "$(grep -c 'ndim=7 (from the preamble)' "$T/2.log")" "1"
check "  and stage 1 says what the rerun changed" \
      "$(grep -c 're-frozen EX_NDIM: 20 -> 7' "$T/2.log")" "1"

# --- 3. an unchanged rerun rewrites no table ----------------------------------
tables="analyses.tsv phenotypes.tsv svd.pcs.txt autosomal.median.txt mosdepth.manifest.txt samples.txt chrom.sizes"
before="$(for t in $tables; do inode "$in_dir/$t"; done | tr '\n' ' ')"
prepare "$T/3.log" || bad "unchanged rerun exited non-zero" "$T/3.log"
check "an unchanged rerun rewrites no input table" \
      "$(for t in $tables; do inode "$in_dir/$t"; done | tr '\n' ' ')" "$before"
check "  and changes nothing in the freeze" "$(grep -c 're-frozen' "$T/3.log")" "0"

# --- 4. later stages keep the freeze when the preamble changes again ---------
echo 9 > "$W/preamble/ndim.txt"
mv "$W/preamble/covariates.tsv" "$W/preamble/covariates.tsv.away"
check "a later stage keeps the frozen ndim when ndim.txt changes" "$(later EX_NDIM)" "7"
check "  and the frozen covariates when covariates.tsv goes" "$(later EX_COVARIATES)" "$adj"
check "  an explicit environment still wins there" "$(later EX_NDIM EX_NDIM=6)" "6"
echo 7 > "$W/preamble/ndim.txt"
mv "$W/preamble/covariates.tsv.away" "$W/preamble/covariates.tsv"

# --- 5. an explicit environment wins over the preamble at the freeze ---------
prepare "$T/5.log" EX_NDIM=5 || bad "prepare with EX_NDIM=5 exited non-zero" "$T/5.log"
check "an explicit EX_NDIM wins over the preamble at the freeze" "$(frozen EX_NDIM)" "5" "$T/5.log"
check "  and run.sh says so" "$(grep -c 'ndim=5 (set in the environment)' "$T/5.log")" "1"

# --- 6. stage 0 says when a rerun switches an upstream table ------------------
mkdir -p "$T/qc_alt"; cp "$qc" "$T/qc_alt/sample_qc.tsv"
prepare "$T/6.log" EX_QC_DIR_STANDARD="$T/qc_alt" || bad "prepare with another QC table exited non-zero" "$T/6.log"
check "a switched phenotype table is warned about" \
      "$(grep -c 'WARN standard: the upstream inputs differ' "$T/6.log")" "1" "$T/6.log"
check "  naming the new table" "$(grep -c "now: EX_M_QC_TABLE=$T/qc_alt/sample_qc.tsv" "$T/6.log")" "1"
prepare "$T/7.log" EX_QC_DIR_STANDARD="$T/qc_alt" || bad "repeat prepare exited non-zero" "$T/7.log"
check "  and only once" "$(grep -c 'WARN standard: the upstream inputs differ' "$T/7.log")" "0"

echo
echo "passed $pass, failed $fail   ($T)"
exit "$fail"
