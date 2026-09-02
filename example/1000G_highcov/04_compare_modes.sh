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
# on disk the stage exits cleanly after writing a note.
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

# compare_pair <a> <b> — runs when both modes have association output.
compare_pair() {
    local a="$1" b="$2" a_dir b_dir out
    a_dir="$(ex_assoc_dir "$a")"; b_dir="$(ex_assoc_dir "$b")"
    [ -d "$a_dir" ] && [ -d "$b_dir" ] || return 1
    out="$EX_COMPARE_DIR/${a}_vs_${b}"
    mkdir -p "$out"
    Rscript "$EX_EXAMPLE_DIR/R/compare_modes.R" \
        --a-dir "$a_dir" --b-dir "$b_dir" --a-name "$a" --b-name "$b" \
        --analyses "$(ex_inputs_dir "$a")/analyses.tsv" \
        --top-k "$EX_TOP_K" --profile "$EX_EVAL_PROFILE" \
        ${t_opt[@]+"${t_opt[@]}"} \
        --out "$out"
}

primary=""; control=""
if compare_pair standard fast; then
    primary="$EX_COMPARE_DIR/standard_vs_fast/concordance.tsv"
else
    dsv_log "standard vs fast: one of the modes has no association output yet"
fi
if compare_pair standard seedctl; then
    control="$EX_COMPARE_DIR/standard_vs_seedctl/concordance.tsv"
else
    dsv_log "no seed-control results; the fast comparison will stand without its yardstick"
fi

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

c_opt=()
[ -z "$control" ] || c_opt=(--control "$control")
Rscript "$EX_EXAMPLE_DIR/R/calibration_summary.R" \
    --primary "$primary" ${c_opt[@]+"${c_opt[@]}"} \
    --factor "$EX_CALIBRATION_FACTOR" --out "$EX_COMPARE_DIR"
