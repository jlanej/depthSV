#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# depthSV 1000G example — stage 4: compare the modes end to end
#
#   bash 04_compare_modes.sh
#
# Two pairs, each joined per analysis on the regions tested in both:
#
#   standard vs fast      the question — does mosdepth --fast-mode upstream
#                         change what the association sweep concludes?
#   standard vs seedctl   the yardstick — the same standard matrix corrected
#                         with PCs from a reseeded NGS-PCA run, so every
#                         difference is the randomized estimator's own noise
#
# Each pair lands in compare/<a>_vs_<b>/{concordance.tsv, summary.md}, and
# compare/summary.md states the calibration verdict per analysis: whether
# fast mode moved the statistics more than a reseed did. Without a seed
# control the fast comparison stands alone and says so; with only one mode
# on disk the stage exits cleanly after writing a note. A pair whose
# concordance FAILs (an analysis missing on one side) fails this stage.
# ---------------------------------------------------------------------------

EX_EXAMPLE_DIR="${EX_EXAMPLE_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}}"
source "$EX_EXAMPLE_DIR/lib.sh"
dsv_enable_error_trace

while [ $# -gt 0 ]; do
    case "$1" in
        -h|--help) dsv_usage ;;
        *)         dsv_die "unknown argument: $1" ;;
    esac
done

dsv_require_cmd Rscript gzip
mkdir -p "$EX_COMPARE_DIR"

# Aggregate whatever timing records exist so the reports can show runtimes.
timings="$EX_PROFILE_DIR/timings.tsv"
t_opt=()
if [ -d "$EX_PROFILE_DIR/timings.d" ]; then
    {
        printf 'mode\tstage\tunit\tstart\telapsed_s\tmax_rss_kb\texit\thost\tjobid\n'
        cat "$EX_PROFILE_DIR"/timings.d/*.tsv 2>/dev/null || true
    } > "$timings"
    t_opt=(--timings "$timings")
fi

mode_source() {                    # mode_source <mode> -> EX_M_SOURCE from paths.env
    local EX_M_SOURCE=""
    # shellcheck disable=SC1090
    [ ! -s "$(ex_paths_env "$1")" ] || source "$(ex_paths_env "$1")"
    printf '%s\n' "${EX_M_SOURCE:-unknown}"
}

# compare_pair <a> <b>: 2 = a mode has no output (not an error), 1 = the
# comparison ran and FAILed, 0 = ran.
compare_pair() {
    local a="$1" b="$2" a_dir b_dir out
    a_dir="$(ex_assoc_dir "$a")"; b_dir="$(ex_assoc_dir "$b")"
    [ -d "$a_dir" ] && [ -d "$b_dir" ] || return 2
    out="$EX_COMPARE_DIR/${a}_vs_${b}"
    mkdir -p "$out"
    Rscript "$EX_EXAMPLE_DIR/R/compare_modes.R" \
        --a-dir "$a_dir" --b-dir "$b_dir" --a-name "$a" --b-name "$b" \
        --a-source "$(mode_source "$a")" --b-source "$(mode_source "$b")" \
        --analyses "$(ex_inputs_dir "$a")/analyses.tsv" \
        --top-k "$EX_TOP_K" --profile "$EX_EVAL_PROFILE" \
        ${t_opt[@]+"${t_opt[@]}"} \
        --out "$out" || return 1
}

rc=0; primary=""; control=""
# The function's status is read inside an if: a bare non-zero return would
# trip errexit before the case below could interpret it.
if compare_pair standard fast; then st=0; else st=$?; fi
case "$st" in
    0) primary="$EX_COMPARE_DIR/standard_vs_fast/concordance.tsv" ;;
    1) primary="$EX_COMPARE_DIR/standard_vs_fast/concordance.tsv"; rc=1
       dsv_log "standard vs fast: concordance FAILED (see $EX_COMPARE_DIR/standard_vs_fast/summary.md)" ;;
    2) dsv_log "standard vs fast: one of the modes has no association output yet" ;;
esac
if compare_pair standard seedctl; then st=0; else st=$?; fi
case "$st" in
    0) control="$EX_COMPARE_DIR/standard_vs_seedctl/concordance.tsv" ;;
    1) control="$EX_COMPARE_DIR/standard_vs_seedctl/concordance.tsv"; rc=1
       dsv_log "standard vs seedctl: concordance FAILED (see $EX_COMPARE_DIR/standard_vs_seedctl/summary.md)" ;;
    2) dsv_log "no seed-control results; the fast comparison will stand without its yardstick" ;;
esac

if [ -z "$primary" ]; then
    {
        echo "# depthSV 1000G example — mode comparison"
        echo
        echo "The standard-vs-fast comparison needs association output from both modes:"
        echo
        for m in standard fast; do
            if [ -d "$(ex_assoc_dir "$m")" ]; then echo "- $m: $(ex_assoc_dir "$m")"; else echo "- $m: missing"; fi
        done
        echo
        echo "Run the missing mode (00 -> 01 -> 02) and rerun this stage."
    } > "$EX_COMPARE_DIR/summary.md"
    dsv_log "wrote $EX_COMPARE_DIR/summary.md"
    exit 0
fi

note="Inputs — standard: \`$(mode_source standard)\`; fast: \`$(mode_source fast)\`"
[ -z "$control" ] || note="$note; seedctl: \`$(mode_source seedctl)\`"
case "$(mode_source fast)" in
    *synthetic*) note="$note. **SYNTHETIC fast mode** (smoke): the fast depths were simulated from the standard tables, so this verdict tests the machinery, not mosdepth." ;;
esac
c_opt=()
[ -z "$control" ] || c_opt=(--control "$control")
Rscript "$EX_EXAMPLE_DIR/R/calibration_summary.R" \
    --primary "$primary" ${c_opt[@]+"${c_opt[@]}"} \
    --factor "$EX_CALIBRATION_FACTOR" --exclude "inferred_sex" --note "$note" \
    --out "$EX_COMPARE_DIR" || rc=1
exit "$rc"
