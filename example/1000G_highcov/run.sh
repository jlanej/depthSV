#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# depthSV 1000G example — one command for the whole example
#
#   bash run.sh [--smoke] [--prepare-only] [--mode all|standard|fast|seedctl]
#               [--runner auto|slurm|local] [--force]
#
# Runs the numbered stages in order — 00 (resolve the upstream NGS-PCA
# inputs), 01 (input tables and the upstream checks), 02 (join, correct,
# analyze; evaluation, comparison and profile follow on their own) — with
# the site configuration taken from the environment, exactly as the stages
# take it themselves:
#
#   export NGSPCA_WORK_DIR=/scratch/$USER/1000G_highcov
#   export EX_WORK_DIR=/scratch/$USER/depthsv_1000G_highcov
#   export EX_SBATCH_EXTRA="--partition=... --account=..."   # if your site needs it
#   bash run.sh
#
#   --smoke          EX_SMOKE=1: the committed upstream results plus simulated
#                    depth trees, so the whole example runs anywhere in minutes
#   --prepare-only   stop after the upstream checks and the input tables: the
#                    preflight to run while the NGS-PCA comparison finishes
#   --mode           one mode instead of every prepared one
#   --runner         slurm or local (default: slurm wherever sbatch exists)
#   --force          redo completed pipeline units (fetched inputs are kept)
#
# Under SLURM this returns as soon as the chains are submitted; results land
# under $EX_WORK_DIR when the last evaluate job finishes. Every stage can
# still be run on its own — this only strings them together.
# ---------------------------------------------------------------------------

smoke=0; prepare_only=0; want_help=0; run_args=()
while [ $# -gt 0 ]; do
    case "$1" in
        --smoke)        smoke=1; shift ;;
        --prepare-only) prepare_only=1; shift ;;
        --runner)       export EX_RUNNER="$2"; run_args+=("$1" "$2"); shift 2 ;;
        --mode)         run_args+=("$1" "$2"); shift 2 ;;
        --force)        run_args+=(--force); shift ;;
        -h|--help)      want_help=1; shift ;;
        *)              echo "unknown argument: $1 (see --help)" >&2; exit 2 ;;
    esac
done
# Exported before config.sh is read: the smoke defaults (ndim, window,
# thresholds) are decided there from EX_SMOKE.
[ "$smoke" -eq 0 ] || export EX_SMOKE=1

EX_EXAMPLE_DIR="${EX_EXAMPLE_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}}"
source "$EX_EXAMPLE_DIR/lib.sh"
dsv_enable_error_trace
[ "$want_help" -eq 0 ] || dsv_usage

dsv_log "depthSV 1000G example -> $EX_WORK_DIR (smoke=$EX_SMOKE, runner=$(ex_runner), modes: $EX_MODES)"
if [ "$EX_SMOKE" = "1" ]; then ndim_src="smoke default"
elif [ -s "$EX_INPUTS_DIR/run.env" ] && grep -q 'EX_NDIM' "$EX_INPUTS_DIR/run.env"; then ndim_src="frozen in inputs/run.env"
elif [ -s "$EX_PREAMBLE_DIR/ndim.txt" ]; then ndim_src="from the preamble"
else ndim_src="default; preamble.sh not run"; fi
dsv_log "ndim=$EX_NDIM ($ndim_src), covariates: $EX_COVARIATES"
[ -s "$EX_PREAMBLE_DIR/covariates.tsv" ] \
    || dsv_log "no preamble covariates: the mtDNA-CN models run unadjusted (submit preamble.sh for genotype PCs and the MP-derived ndim)"

bash "$EX_EXAMPLE_DIR/00_fetch_inputs.sh"
bash "$EX_EXAMPLE_DIR/01_prepare_inputs.sh"

ready="$(ex_ready_modes | tr '\n' ' ')"
[ -n "$ready" ] || dsv_die "nothing to run: no mode has both its upstream tables and a mosdepth tree (see the messages above)"
dsv_log "prepared modes: $ready"

if [ "$prepare_only" -eq 1 ]; then
    dsv_log "--prepare-only: stopping here. Any WARN above is the preflight finding; inputs are under $EX_INPUTS_DIR"
    exit 0
fi

rc=0
bash "$EX_EXAMPLE_DIR/02_run_depthsv.sh" ${run_args[@]+"${run_args[@]}"} || rc=$?

echo
if [ "$(ex_runner)" = slurm ]; then
    dsv_log "chains submitted. When the last evaluate job finishes, read:"
    dsv_log "  $EX_EVAL_DIR/<mode>/summary.md   $EX_COMPARE_DIR/summary.md   $EX_PROFILE_DIR/profile_report.md"
else
    for f in "$EX_EVAL_DIR"/*/summary.md "$EX_COMPARE_DIR/summary.md"; do
        [ -s "$f" ] || continue
        printf '%s: %s\n' "${f#"$EX_WORK_DIR"/}" "$(grep -m1 -i 'verdict' "$f" | sed 's/^- //')"
    done
    dsv_log "profile: $EX_PROFILE_DIR/profile_report.md"
fi
exit "$rc"
