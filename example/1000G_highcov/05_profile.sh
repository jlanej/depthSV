#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# depthSV 1000G example — stage 5: profile the run
#
#   bash 05_profile.sh
#
# Folds the per-invocation timing records (written by every stage via
# ex_timed) into profile/timings.tsv, augments them with the scheduler's own
# accounting where sacct exists, and writes profile/profile_report.md: time
# by stage and mode, the slowest work units, and peak memory where it was
# recordable. This is the artefact to read when deciding what to speed up —
# and what to hand to a bigger allocation.
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

dsv_require_cmd Rscript
mkdir -p "$EX_PROFILE_DIR"

timings="$EX_PROFILE_DIR/timings.tsv"
[ -d "$EX_PROFILE_DIR/timings.d" ] || dsv_die "no timing records yet (nothing has run under 02_run_depthsv.sh)"
{
    printf 'mode\tstage\tunit\tstart\telapsed_s\tmax_rss_kb\texit\thost\tjobid\n'
    cat "$EX_PROFILE_DIR"/timings.d/*.tsv 2>/dev/null || true
} > "$timings"
dsv_log "$(( $(grep -c . "$timings") - 1 )) timing records -> $timings"

# Scheduler accounting for every job the driver submitted. Array tasks show
# up as <jobid>_<index>, which is how per-unit MaxRSS becomes visible.
sacct_tsv="$EX_PROFILE_DIR/sacct.tsv"
sacct_opt=()
if command -v sacct >/dev/null 2>&1 && [ -s "$EX_PROFILE_DIR/jobs.tsv" ]; then
    ids="$(cut -f3 "$EX_PROFILE_DIR/jobs.tsv" | sort -u | paste -sd, -)"
    if sacct -j "$ids" -P --noconvert \
             --format=JobID,JobName,State,Elapsed,TotalCPU,MaxRSS,ReqMem,AllocCPUS,NodeList \
             > "$sacct_tsv" 2>/dev/null && [ -s "$sacct_tsv" ]; then
        sacct_opt=(--sacct "$sacct_tsv")
        dsv_log "scheduler accounting -> $sacct_tsv"
    else
        dsv_log "sacct gave nothing usable; the report uses the recorder's timings only"
    fi
fi

Rscript "$EX_EXAMPLE_DIR/R/profile_report.R" \
    --timings "$timings" \
    ${sacct_opt[@]+"${sacct_opt[@]}"} \
    --out "$EX_PROFILE_DIR"

dsv_log "report: $EX_PROFILE_DIR/profile_report.md"
