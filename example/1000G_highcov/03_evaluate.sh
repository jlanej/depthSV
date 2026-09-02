#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# depthSV 1000G example — stage 3: evaluate one or both modes
#
#   bash 03_evaluate.sh [--mode all|standard|fast]
#
# Runs the truth checks over a mode's association output (see R/evaluate.R
# for what is asserted and why) and writes eval/<mode>/{summary.md,
# checks.tsv, top_hits.*.tsv}. Exits non-zero if any machinery check FAILs;
# statistical WARNs do not fail the run.
# ---------------------------------------------------------------------------

EX_EXAMPLE_DIR="${EX_EXAMPLE_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}}"
source "$EX_EXAMPLE_DIR/lib.sh"
dsv_enable_error_trace

mode_arg="all"
while [ $# -gt 0 ]; do
    case "$1" in
        --mode)    mode_arg="$2"; shift 2 ;;
        -h|--help) dsv_usage ;;
        *)         dsv_die "unknown argument: $1" ;;
    esac
done

dsv_require_cmd Rscript gzip

modes="$(ex_ready_modes)"
[ "$mode_arg" = "all" ] || { ex_check_mode "$mode_arg"; modes="$mode_arg"; }

rc=0
for mode in $modes; do
    assoc="$(ex_assoc_dir "$mode")"
    in_dir="$(ex_inputs_dir "$mode")"
    if [ ! -d "$assoc" ]; then
        dsv_log "SKIP $mode: no association output at $assoc"
        continue
    fi
    out="$EX_EVAL_DIR/$mode"
    mkdir -p "$out"
    regions_opt=()
    [ ! -s "$(ex_regions_file "$mode")" ] || regions_opt=(--regions "$(ex_regions_file "$mode")")
    ex_timed "$mode" evaluate evaluate -- \
        Rscript "$EX_EXAMPLE_DIR/R/evaluate.R" \
            --assoc "$assoc" --analyses "$in_dir/analyses.tsv" \
            --mode "$mode" --profile "$EX_EVAL_PROFILE" --out "$out" \
            ${regions_opt[@]+"${regions_opt[@]}"} \
        || rc=1
done

exit "$rc"
