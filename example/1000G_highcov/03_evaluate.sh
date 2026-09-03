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
    # Provenance travels with the verdict: a synthetic smoke tree is labelled
    # as such in the summary, not only in paths.env.
    EX_M_SOURCE=""
    # shellcheck disable=SC1090
    [ ! -s "$(ex_paths_env "$mode")" ] || source "$(ex_paths_env "$mode")"

    # The export first: every analysis over the whole region list, count
    # suppression, and the empirical threshold from the permutation maxima.
    # A failed export is reported but does not stop the truth checks.
    export_opt=()
    if [ -s "$(ex_regions_file "$mode")" ]; then
        ex_export_dsv_env "$mode"
        if ex_timed "$mode" export export -- \
               bash "$DSV_ROOT/scripts/export.sh" --results "$assoc" --regions "$(ex_regions_file "$mode")" \
                    --pheno-manifest "$in_dir/analyses.tsv" --out "$DSV_EXPORT_DIR"; then
            export_opt=(--export "$DSV_EXPORT_DIR")
        else
            dsv_log "WARN: export failed for $mode (see above); evaluating the shards without it"
            rc=1
        fi
    fi

    ex_timed "$mode" evaluate evaluate -- \
        Rscript "$EX_EXAMPLE_DIR/R/evaluate.R" \
            --assoc "$assoc" --analyses "$in_dir/analyses.tsv" \
            --mode "$mode" --profile "$EX_EVAL_PROFILE" --out "$out" \
            --source "${EX_M_SOURCE:-unknown}" --pheno "$in_dir/phenotypes.tsv" \
            --samples "$in_dir/samples.txt" --ploidy "$EX_PLOIDY" \
            ${regions_opt[@]+"${regions_opt[@]}"} ${export_opt[@]+"${export_opt[@]}"} \
        || rc=1
done

exit "$rc"
