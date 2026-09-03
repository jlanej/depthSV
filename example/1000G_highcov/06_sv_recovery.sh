#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# depthSV 1000G example — stage 6: recover known deletions as a function of ndim
#
#   bash 06_sv_recovery.sh [--mode standard|fast] [--force]
#
# The Marchenko-Pastur count says how many coverage components stand out
# from noise, not how many should be removed from the depth before an
# association test; the outcome that matters is whether real copy-number
# variation survives the correction. This stage measures that directly: a
# set of known autosomal deletions (the NYGC 3,202-sample SV callset; in
# smoke mode the deletions the simulated tree carries) is corrected at
# several ndims, and for each deletion the corrected depth of carriers is
# compared with non-carriers — AUC, shift in log2 units (one lost copy is
# about -1) and Welch t. The summary per ndim shows where recovery plateaus
# and where it starts to erode, with the MP count and the run's EX_NDIM
# marked on the same axis.
#
# Outputs under EX_WORK_DIR/sv_recovery/<mode>/: selected.tsv (the
# deletions used), sv_recovery.tsv (per deletion x ndim),
# sv_recovery_summary.tsv (per ndim), sv_recovery.png, summary.md and
# recommended_ndim.txt — the smallest ndim within 0.01 of the best median
# AUC. Informational: it does not change EX_NDIM.
#
# A real run downloads the callset VCF once into EX_WORK_DIR/sv_callset
# and extracts the deletions with awk (EX_SV_CALLSET_URL; EX_SV_CALLS
# points at a prepared table instead). Corrections run one small unit per
# deletion and ndim, skipped when already done.
# ---------------------------------------------------------------------------

EX_EXAMPLE_DIR="${EX_EXAMPLE_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}}"
source "$EX_EXAMPLE_DIR/lib.sh"
dsv_enable_error_trace

mode="standard"; force=0; fetch_only=0
while [ $# -gt 0 ]; do
    case "$1" in
        --mode)       mode="$2"; shift 2 ;;
        --force)      force=1; shift ;;
        --fetch-only) fetch_only=1; shift ;;   # download and reduce the callset, then stop
        -h|--help)    dsv_usage ;;
        *)            dsv_die "unknown argument: $1" ;;
    esac
done
[ "$mode" != all ] || mode=standard
ex_check_mode "$mode"
[ "$mode" != seedctl ] || dsv_die "the seed control shares the standard matrix; use --mode standard"

# Modules first: on a cluster the tools may only exist after `module load`.
dsv_load_modules
dsv_require_cmd Rscript gzip awk
[ "$fetch_only" -eq 1 ] || dsv_require_cmd bgzip tabix parallel
ex_export_dsv_env "$mode"
out="$(ex_sv_dir "$mode")"
if [ "$fetch_only" -eq 0 ]; then
    dsv_output_complete "$DSV_MATRIX" || dsv_die "$mode: no finished matrix at $DSV_MATRIX (run 02_run_depthsv.sh first)"
    in_dir="$(ex_inputs_dir "$mode")"
    [ -f "$in_dir/prepared.ok" ] || dsv_die "$mode: inputs are not marked ready; run 01_prepare_inputs.sh first"
    dsv_require_file "$in_dir/samples.txt"
    mkdir -p "$out"
fi
force_flag=()
[ "$force" -eq 0 ] || force_flag=(--force)

# --- the deletions ---------------------------------------------------------

calls=""
if [ -n "$EX_SV_CALLS" ]; then
    dsv_require_file "$EX_SV_CALLS"; calls="$EX_SV_CALLS"
    dsv_log "deletions from EX_SV_CALLS: $calls"
elif [ "$EX_SMOKE" = "1" ]; then
    calls="$(ex_mosdepth_dir "$mode")/sv_calls.tsv"
    dsv_require_file "$calls"
    dsv_log "smoke: deletions the simulated tree carries ($calls)"
else
    cs="$EX_WORK_DIR/sv_callset"; mkdir -p "$cs"
    vcf="$cs/$(basename "$EX_SV_CALLSET_URL")"
    calls="$cs/sv_calls.tsv"
    if [ "$force" -eq 0 ] && [ -s "$calls" ] && [ -f "$calls.done" ]; then
        dsv_log "deletions present: $calls"
    else
        if [ ! -s "$vcf" ]; then
            dsv_require_cmd curl
            dsv_log "downloading the SV callset: $EX_SV_CALLSET_URL"
            curl -fL --retry 5 --retry-delay 10 -C - -o "$vcf.part" "$EX_SV_CALLSET_URL" \
                || dsv_die "download failed; set EX_SV_CALLS to a prepared calls table, or EX_SV_RECOVERY=0"
            mv "$vcf.part" "$vcf"
        fi
        dsv_log "extracting PASS deletions >= $EX_SV_MIN_LEN bp with carrier frequency $EX_SV_MIN_AF-$EX_SV_MAX_AF"
        # One row per deletion: 0-based start, end, id, carrier count and the
        # carriers (any ALT allele in the GT, the first FORMAT field).
        gzip -cd "$vcf" | awk -F'\t' -v min_len="$EX_SV_MIN_LEN" -v min_af="$EX_SV_MIN_AF" -v max_af="$EX_SV_MAX_AF" '
            BEGIN { OFS = "\t"; print "CHROM", "START", "END", "ID", "N_CARRIERS", "CARRIERS" }
            /^##/ { next }
            /^#CHROM/ { for (i = 10; i <= NF; i++) sid[i] = $i; ns = NF - 9; next }
            {
              if ($7 != "PASS" && $7 != ".") next
              if (index($8, "SVTYPE=DEL") == 0) next
              n = split($8, kv, ";"); end = 0; svlen = 0
              for (i = 1; i <= n; i++) {
                  if (kv[i] ~ /^END=/) end = substr(kv[i], 5) + 0
                  else if (kv[i] ~ /^SVLEN=/) svlen = substr(kv[i], 7) + 0
              }
              if (svlen < 0) svlen = -svlen
              if (end == 0 && svlen > 0) end = $2 + svlen
              if (end == 0 || end - $2 < min_len) next
              nc = 0; car = ""
              for (i = 10; i <= NF; i++) {
                  gt = $i; p = index(gt, ":"); if (p) gt = substr(gt, 1, p - 1)
                  if (gt ~ /[1-9]/) { nc++; car = (car == "" ? sid[i] : car "," sid[i]) }
              }
              if (ns == 0 || nc / ns < min_af || nc / ns > max_af) next
              print $1, $2 - 1, end, $3, nc, car
            }' > "$calls.tmp"
        mv "$calls.tmp" "$calls"
        touch "$calls.done"
        dsv_log "$(( $(grep -c . "$calls") - 1 )) deletions -> $calls"
    fi
fi
if [ "$fetch_only" -eq 1 ]; then
    dsv_log "--fetch-only: the deletions are ready at $calls; the recovery runs at the end of 02_run_depthsv.sh"
    exit 0
fi

# --- which deletions, which ndims -------------------------------------------

mp=""
[ ! -s "$EX_PREAMBLE_DIR/ndim.txt" ] || mp="$(tr -cd '0-9' < "$EX_PREAMBLE_DIR/ndim.txt")"
# shellcheck disable=SC2086
ndims="$(printf '%s\n' $EX_SV_NDIMS "$EX_NDIM" ${mp:+"$mp"} | grep -E '^[0-9]+$' | sort -nu | tr '\n' ' ' | sed 's/ $//')"

Rscript "$EX_EXAMPLE_DIR/R/sv_recovery.R" --stage select \
    --calls "$calls" --contigs "$(tabix -l "$DSV_MATRIX" | tr '\n' ',')" \
    --samples "$in_dir/samples.txt" --max-dels "$EX_SV_MAX_DELS" --min-len "$EX_SV_MIN_LEN" \
    --min-af "$EX_SV_MIN_AF" --max-af "$EX_SV_MAX_AF" --out "$out"
n_sel="$(( $(grep -c . "$out/selected.tsv") - 1 ))"
[ "$n_sel" -gt 0 ] || dsv_die "no usable deletion: none on the matrix's autosomes within the length and frequency bands"
dsv_log "$n_sel deletions x ndims [$ndims] -> $out/corrected"

# --- corrections, one small unit per deletion and ndim -------------------------

: > "$out/units.list"
for k in $ndims; do
    while IFS= read -r region; do
        [ -z "$region" ] || printf '%s %s\n' "$k" "$region" >> "$out/units.list"
    done < "$out/units.txt"
done
# The unit's own log carries its diagnostics; the driver output here is the
# progress lines.
ex_timed "$mode" sv-correct all -- \
    parallel --colsep ' ' -j "$EX_JOBS" --halt soon,fail=1 \
        bash "$DSV_ROOT/scripts/correct.sh" --matrix "$DSV_MATRIX" --pcs "$DSV_PCS" --coverage "$DSV_COVERAGE" \
             --region {2} --out "$out/corrected" --ndim {1} --jobs 1 --threads 1 \
             ${force_flag[@]+"${force_flag[@]}"} \
        :::: "$out/units.list" > "$out/correct.log" 2>&1 \
    || dsv_die "a correction unit failed (see $out/correct.log)"

# --- recovery per deletion and ndim ---------------------------------------------

ex_timed "$mode" sv-evaluate all -- \
    Rscript "$EX_EXAMPLE_DIR/R/sv_recovery.R" --stage evaluate \
        --selected "$out/selected.tsv" --corrected "$out/corrected" --ndims "$ndims" \
        --mp "${mp:-NA}" --ndim "$EX_NDIM" --mode "$mode" --smoke "$EX_SMOKE" --out "$out"
dsv_log "SV recovery: $out/summary.md (recommended ndim: $(cat "$out/recommended_ndim.txt"))"
