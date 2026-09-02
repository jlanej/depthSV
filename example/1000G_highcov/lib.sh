#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# depthSV 1000G example — shared helpers
#
# Source this at the top of every stage script in this directory:
#
#     EX_EXAMPLE_DIR="${EX_EXAMPLE_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}}"
#     source "$EX_EXAMPLE_DIR/lib.sh"
#
# The EX_EXAMPLE_DIR dance matters under SLURM: sbatch copies the submitted
# script into a spool directory, so dirname of the running script no longer
# points at this directory. The driver exports EX_EXAMPLE_DIR before every
# submission; SLURM_SUBMIT_DIR covers a script submitted by hand from here.
#
# Loads the pipeline's own lib/common.sh (logging, guards, dsv_sample_name)
# and config.sh, then adds what only the example needs: per-mode path
# resolution, the DSV_* environment for one mode, and a timing recorder that
# every stage runs under so the run can be profiled afterwards.
#
# Modes:
#   standard  mosdepth's default mode: its tree, its NGS-PCA run, its QC table
#   fast      mosdepth --fast-mode: its tree, its NGS-PCA run, its QC table
#   seedctl   the STANDARD tree and QC table with the seed-control PCs
#             (NGS-PCA rerun under another -randomSeed). Shares the standard
#             depth matrix; only correction and analysis are redone. The
#             standard-vs-seedctl differences are the estimator's own noise,
#             which is what the fast-vs-standard differences are judged
#             against.
# ---------------------------------------------------------------------------

EX_EXAMPLE_DIR="${EX_EXAMPLE_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
[ -s "$EX_EXAMPLE_DIR/config.sh" ] \
    || { echo "ERROR: $EX_EXAMPLE_DIR does not look like example/1000G_highcov (no config.sh)." >&2
         echo "       Run from that directory, or export EX_EXAMPLE_DIR." >&2; exit 2; }
export EX_EXAMPLE_DIR

source "$EX_EXAMPLE_DIR/../../lib/common.sh"   # sets -euo pipefail, DSV_ROOT

# The parameters stage 1 froze for this run (inputs/run.env) load BEFORE the
# configuration, as defaults that only an explicit environment overrides:
# every job of a submission then agrees on ndim and the covariates however
# the preamble's files change underneath. The work-dir default here must
# match config.sh's.
EX_WORK_DIR="${EX_WORK_DIR:-/scratch/${USER:-$(id -un)}/depthsv_1000G_highcov}"
# shellcheck disable=SC1090
[ ! -s "${EX_INPUTS_DIR:-$EX_WORK_DIR/inputs}/run.env" ] || source "${EX_INPUTS_DIR:-$EX_WORK_DIR/inputs}/run.env"

source "$EX_EXAMPLE_DIR/config.sh"

# --- modes -----------------------------------------------------------------

ex_check_mode() {                  # ex_check_mode <mode>
    case "$1" in standard|fast|seedctl) ;;
        *) dsv_die "unknown mode '$1' (expected: standard | fast | seedctl)" ;;
    esac
}

# The mosdepth tree a mode is built on. seedctl reuses the standard tree;
# smoke runs swap in the simulated trees so later stages are oblivious.
ex_mosdepth_dir() {
    ex_check_mode "$1"
    local tree="$1"
    [ "$tree" != seedctl ] || tree=standard
    if [ "$EX_SMOKE" = "1" ]; then printf '%s/%s\n' "$EX_SMOKE_DIR" "$tree"; return; fi
    case "$tree" in
        standard) printf '%s\n' "$EX_MOSDEPTH_DIR_STANDARD" ;;
        fast)     printf '%s\n' "$EX_MOSDEPTH_DIR_FAST" ;;
    esac
}

ex_ngspca_dir() {
    ex_check_mode "$1"
    case "$1" in
        standard) printf '%s\n' "$EX_NGSPCA_DIR_STANDARD" ;;
        fast)     printf '%s\n' "$EX_NGSPCA_DIR_FAST" ;;
        seedctl)  printf '%s\n' "$EX_NGSPCA_DIR_SEEDCTL" ;;
    esac
}

ex_qc_dir() {
    ex_check_mode "$1"
    case "$1" in
        standard|seedctl) printf '%s\n' "$EX_QC_DIR_STANDARD" ;;
        fast)             printf '%s\n' "$EX_QC_DIR_FAST" ;;
    esac
}

ex_github_path() {
    case "$1" in
        standard) printf '%s\n' "$EX_GITHUB_PATH_STANDARD" ;;
        fast)     printf '%s\n' "$EX_GITHUB_PATH_FAST" ;;
        *)        printf '\n' ;;
    esac
}

# Where 00_fetch_inputs.sh recorded what it resolved for a mode.
ex_paths_env() { printf '%s/%s/paths.env\n' "$EX_INPUTS_DIR" "$1"; }

ex_inputs_dir() { printf '%s/%s\n' "$EX_INPUTS_DIR" "$1"; }

# The joined matrix a mode reads. seedctl shares the standard one — same
# manifest, same rows — so the join is never repeated for it.
ex_join_dir() {
    ex_check_mode "$1"
    local tree="$1"
    [ "$tree" != seedctl ] || tree=standard
    printf '%s/%s/join\n' "$EX_RUN_DIR" "$tree"
}

# Association output directory for one mode. Carries the ndim because the
# analysis filenames do not: without it, rerunning after an EX_NDIM change
# would skip every analysis as already complete and hand back stale results.
ex_assoc_dir() { printf '%s/%s/assoc_ndim%s\n' "$EX_RUN_DIR" "$1" "$EX_NDIM"; }
# The export of that directory: one table, summary and hits per analysis.
ex_export_dir() { printf '%s/%s/export_ndim%s\n' "$EX_RUN_DIR" "$1" "$EX_NDIM"; }
# SV-callset recovery output for one mode.
ex_sv_dir() { printf '%s/sv_recovery/%s\n' "$EX_WORK_DIR" "$1"; }

# A mode is ready once 00 and 01 have produced its inputs and 01 finished
# the mode (prepared.ok); a prepare that died part-way leaves no marker.
ex_mode_ready() {                  # ex_mode_ready <mode>
    local in; in="$(ex_inputs_dir "$1")"
    [ -f "$in/prepared.ok" ] && [ -s "$in/mosdepth.manifest.txt" ] \
        && [ -s "$in/svd.pcs.txt" ] && [ -s "$in/phenotypes.tsv" ]
}

ex_ready_modes() {
    local m
    for m in $EX_MODES; do
        if ex_mode_ready "$m"; then printf '%s\n' "$m"; fi
    done
    return 0
}

# --- DSV_* environment for one mode ----------------------------------------
#
# The stage scripts read their fixed inputs from DSV_* (see conf/example.env);
# exporting them per mode is what lets one driver, one sbatch script and one
# region list serve every mode without repeating flags anywhere.

ex_export_dsv_env() {              # ex_export_dsv_env <mode>
    local mode="$1" in run
    ex_check_mode "$mode"
    in="$(ex_inputs_dir "$mode")"
    run="$EX_RUN_DIR/$mode"

    export DSV_MANIFEST="$in/mosdepth.manifest.txt"
    export DSV_PCS="$in/svd.pcs.txt"
    export DSV_COVERAGE="$in/autosomal.median.txt"
    export DSV_PHENO="$in/phenotypes.tsv"
    export DSV_PHENO_MANIFEST="$in/analyses.tsv"

    DSV_JOIN_DIR="$(ex_join_dir "$mode")"; export DSV_JOIN_DIR
    export DSV_MATRIX="$DSV_JOIN_DIR/depth.matrix.txt.gz"
    export DSV_CORRECTED_DIR="$run/corrected"
    DSV_RESULTS_DIR="$(ex_assoc_dir "$mode")"; export DSV_RESULTS_DIR

    export DSV_NDIM="$EX_NDIM"
    export DSV_MIN_OBS="$EX_MIN_OBS"
    export DSV_MAX_SHARE="$EX_MAX_SHARE"
    export DSV_WINSOR_LOG2="$EX_WINSOR_LOG2"
    export DSV_PERMS="$EX_PERMS"
    export DSV_PERM_SEED="$EX_PERM_SEED"
    export DSV_MIN_COUNT="$EX_MIN_COUNT"
    DSV_EXPORT_DIR="$(ex_export_dir "$mode")"; export DSV_EXPORT_DIR
    export DSV_JOBS="$EX_JOBS"
    export DSV_THREADS="$EX_THREADS"

    # The ploidy model reads the inferred sex the phenotype table carries as
    # M/F; the PAR BED comes from the depthSV checkout.
    if [ "$EX_PLOIDY" = "1" ]; then
        export DSV_SEX="$in/phenotypes.tsv"
        export DSV_SEX_COL="SEX_MF"
        export DSV_PAR="$EX_PAR"
    else
        unset DSV_SEX DSV_SEX_COL DSV_PAR
    fi
}

# --- timing ----------------------------------------------------------------
#
# Every stage invocation runs under ex_timed, which appends one record to a
# per-process file under profile/timings.d/ (one file per record: concurrent
# array tasks on a shared filesystem never contend). 05_profile.sh folds
# them into one table. Wall time is portable; peak RSS is recorded where GNU
# time exists (Linux clusters), and NA elsewhere.

if [ -z "${EX_GTIME+x}" ]; then
    EX_GTIME=""
    if /usr/bin/time --version >/dev/null 2>&1; then EX_GTIME="/usr/bin/time"; fi
fi

ex_timed() {                       # ex_timed <mode> <stage> <unit> -- cmd args...
    local mode="$1" stage="$2" unit="$3"; shift 3
    [ "${1:-}" = "--" ] && shift
    local dir="$EX_PROFILE_DIR/timings.d"
    mkdir -p "$dir"
    local rec rss="NA" rc=0 t0 t1 start_iso use_gtime="$EX_GTIME"
    # GNU time can only exec a binary; shell functions run bare (wall clock only).
    [ "$(type -t "$1" 2>/dev/null)" = "function" ] && use_gtime=""
    rec="$dir/$(date -u +%Y%m%dT%H%M%SZ).$$.${RANDOM}.tsv"
    start_iso="$(dsv_now)"
    t0=$(date +%s)
    if [ -n "$use_gtime" ]; then
        "$use_gtime" -f '%M' -o "$rec.rss" "$@" || rc=$?
        [ -s "$rec.rss" ] && rss="$(tail -n 1 "$rec.rss")"
        rm -f "$rec.rss"
    else
        "$@" || rc=$?
    fi
    t1=$(date +%s)
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$mode" "$stage" "$unit" "$start_iso" "$((t1 - t0))" "$rss" "$rc" \
        "$(hostname -s 2>/dev/null || echo unknown)" "${SLURM_JOB_ID:-NA}" > "$rec"
    return "$rc"
}

# Record a submitted SLURM job so 05_profile.sh can ask sacct about it.
ex_record_job() {                  # ex_record_job <mode> <stage> <jobid>
    mkdir -p "$EX_PROFILE_DIR"
    printf '%s\t%s\t%s\n' "$1" "$2" "$3" >> "$EX_PROFILE_DIR/jobs.tsv"
}

# --- misc ------------------------------------------------------------------

ex_runner() {
    case "$EX_RUNNER" in
        slurm|local) printf '%s\n' "$EX_RUNNER" ;;
        auto) if command -v sbatch >/dev/null 2>&1; then echo slurm; else echo local; fi ;;
        *) dsv_die "EX_RUNNER must be slurm, local or auto (got '$EX_RUNNER')" ;;
    esac
}

ex_regions_file() { printf '%s/%s.regions.txt\n' "$EX_REGIONS_DIR" "$1"; }
