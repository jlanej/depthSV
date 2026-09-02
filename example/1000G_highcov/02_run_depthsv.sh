#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# depthSV 1000G example — stage 2: run the pipeline per mode
#
#   bash 02_run_depthsv.sh [--mode all|standard|fast|seedctl] [--runner auto|slurm|local] [--force]
#
# Per mode this runs join -> regions -> (correct + analyze per region), with
# every stage under the timing recorder so 05_profile.sh can show where the
# time went. `all` (the default) runs every mode 00 and 01 prepared, in
# EX_MODES order; a mode without prepared inputs is skipped with a note.
#
# With the SLURM runner each mode becomes a self-scheduling chain:
#
#   join job ──afterok──> dispatch job ──> array over the region list
#                                          └─afterok──> evaluate job
#
# The dispatch job exists because the region list can only be produced from
# the joined matrix's index, and the array cannot be sized before the list
# exists. The seed-control mode shares the standard matrix, so its chain has
# no join of its own: its dispatch waits on the standard join submitted
# alongside it (or requires the matrix to exist already). Whichever mode's
# evaluate job finishes last also runs the cross-mode comparison and the
# profile report, so one submission takes the whole example to its end
# state. Every job receives the parameters this submission resolved
# (inputs/run.env, frozen by stage 1), so a preamble landing mid-run cannot
# change ndim under the array. Completed units are skipped on resubmission
# (the stage scripts' parameter-aware .done markers), so a partial failure
# is rerun by resubmitting the same command.
#
# The local runner executes the same stages in-process, which is what smoke
# mode uses.
#
# Internal entry points (used by this script's own submissions):
#   --stage join-exec|dispatch|unit|eval-exec --mode <mode>
# ---------------------------------------------------------------------------

EX_EXAMPLE_DIR="${EX_EXAMPLE_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}}"
source "$EX_EXAMPLE_DIR/lib.sh"
dsv_enable_error_trace

mode_arg="all"; stage="driver"; runner=""
while [ $# -gt 0 ]; do
    case "$1" in
        --mode)   mode_arg="$2"; shift 2 ;;
        --runner) runner="$2"; shift 2 ;;
        --stage)  stage="$2"; shift 2 ;;
        --force)  export EX_FORCE=1; shift ;;
        -h|--help) dsv_usage ;;
        *)        dsv_die "unknown argument: $1" ;;
    esac
done
[ -z "$runner" ] || EX_RUNNER="$runner"

force_flag=()
[ "${EX_FORCE:-0}" != "1" ] || force_flag=(--force)

selected_file="$EX_WORK_DIR/.selected_modes"
lock_dir="$EX_WORK_DIR/.finalize.lock"

modes_selected() {
    local m
    if [ "$mode_arg" = all ]; then
        for m in $EX_MODES; do
            if ex_mode_ready "$m"; then
                printf '%s\n' "$m"
            else
                dsv_log "skipping mode '$m': no prepared inputs (00_fetch_inputs.sh / 01_prepare_inputs.sh did not produce them)"
            fi
        done
    else
        ex_check_mode "$mode_arg"
        printf '%s\n' "$mode_arg"
    fi
}

require_inputs() {                 # require_inputs <mode>
    local in; in="$(ex_inputs_dir "$1")"
    dsv_require_file "$in/svd.pcs.txt" "$in/autosomal.median.txt" \
                     "$in/phenotypes.tsv" "$in/analyses.tsv"
    [ -s "$in/mosdepth.manifest.txt" ] \
        || dsv_die "$1: no mosdepth manifest in $in - run 00_fetch_inputs.sh and 01_prepare_inputs.sh first (EX_SMOKE=1 for a simulated tree)"
    [ -f "$in/prepared.ok" ] || dsv_die "$1: inputs are not marked ready; rerun 01_prepare_inputs.sh"
}

# --- stage bodies ----------------------------------------------------------

do_join() {                        # do_join <mode>
    local mode="$1" ff=()
    ex_export_dsv_env "$mode"
    # seedctl joins into the standard matrix directory: complete there means
    # complete here, and --force must not rebuild a matrix it merely shares.
    [ "$mode" = seedctl ] || ff=(${force_flag[@]+"${force_flag[@]}"})
    dsv_log "[$mode] join: $(grep -c . "$DSV_MANIFEST") samples -> $DSV_MATRIX"
    ex_timed "$mode" join join -- \
        bash "$DSV_ROOT/scripts/join.sh" ${ff[@]+"${ff[@]}"}
}

do_regions() {                     # do_regions <mode>
    local mode="$1" out tmp sizes
    ex_export_dsv_env "$mode"
    out="$(ex_regions_file "$mode")"
    mkdir -p "$EX_REGIONS_DIR"
    tmp="$out.tmp"
    sizes="$(ex_inputs_dir "$mode")/chrom.sizes"
    if [ "$EX_WINDOW" -gt 0 ]; then
        dsv_require_file "$sizes"
        bash "$DSV_ROOT/scripts/regions.sh" --matrix "$DSV_MATRIX" \
             --window "$EX_WINDOW" --sizes "$sizes" > "$tmp"
    else
        bash "$DSV_ROOT/scripts/regions.sh" --matrix "$DSV_MATRIX" > "$tmp"
    fi
    grep -E "$EX_CONTIG_REGEX" "$tmp" > "$out" \
        || dsv_die "$mode: no region matched EX_CONTIG_REGEX ('$EX_CONTIG_REGEX')"
    rm -f "$tmp"
    dsv_log "[$mode] $(grep -c . "$out") work units -> $out"
}

do_unit() {                        # do_unit <mode> <region>
    local mode="$1" region="$2" slug
    ex_export_dsv_env "$mode"
    slug="$(printf '%s' "$region" | tr ':' '_' | tr -d ' ')"
    ex_timed "$mode" correct "$region" -- \
        bash "$DSV_ROOT/scripts/correct.sh" --region "$region" ${force_flag[@]+"${force_flag[@]}"}
    ex_timed "$mode" analyze "$region" -- \
        bash "$DSV_ROOT/scripts/analyze.sh" \
            --corrected "$DSV_CORRECTED_DIR/corrected_ndim${EX_NDIM}.${slug}.txt.gz" \
            --region "$region" ${force_flag[@]+"${force_flag[@]}"}
}

# The comparison and the profile run once every mode THIS submission
# selected has been evaluated. Under SLURM the chains run concurrently, so
# whichever evaluate job finishes last is the one that sees a complete set.
# The lock is owned and dated: a finalize that was killed leaves a lock the
# next attempt can recognise as stale, and a new submission clears it.
finalize_modes() {
    if [ -s "$selected_file" ]; then cat "$selected_file"; else ex_ready_modes | tr '\n' ' '; fi
}

lock_is_stale() {
    local stamp age
    [ -d "$lock_dir" ] || return 1
    stamp="$(cat "$lock_dir/epoch" 2>/dev/null || echo 0)"
    age=$(( $(date +%s) - stamp ))
    [ "$age" -gt $(( ${EX_FINALIZE_LOCK_HOURS:-4} * 3600 )) ]
}

maybe_finalize() {
    local m rc=0
    for m in $(finalize_modes); do
        [ -s "$EX_EVAL_DIR/$m/summary.md" ] || { dsv_log "finalize: waiting on mode '$m'"; return 0; }
    done
    if lock_is_stale; then
        dsv_log "finalize: removing a stale lock held by $(cat "$lock_dir/owner" 2>/dev/null || echo '?')"
        rm -rf "$lock_dir"
    fi
    if ! mkdir "$lock_dir" 2>/dev/null; then
        dsv_log "finalize: another task holds the lock ($(cat "$lock_dir/owner" 2>/dev/null || echo '?'))"
        return 0
    fi
    printf '%s\n' "${SLURM_JOB_ID:-pid$$}@$(hostname -s 2>/dev/null || echo unknown)" > "$lock_dir/owner"
    date +%s > "$lock_dir/epoch"
    bash "$EX_EXAMPLE_DIR/04_compare_modes.sh" || { rc=1; dsv_log "WARN: 04_compare_modes.sh exited non-zero"; }
    bash "$EX_EXAMPLE_DIR/05_profile.sh"       || { rc=1; dsv_log "WARN: 05_profile.sh exited non-zero"; }
    rm -rf "$lock_dir"
    return "$rc"
}

# --- SLURM plumbing --------------------------------------------------------

submit() {                         # submit <mode> <stage> <sbatchArgs...> -- <scriptArgs...>
    local mode="$1" st="$2" jid; shift 2
    local sb=()
    while [ $# -gt 0 ] && [ "$1" != "--" ]; do sb+=("$1"); shift; done
    shift
    mkdir -p "$EX_LOG_DIR"
    # Resource strings are word-split on purpose (globbing off while they
    # are); EX_SBATCH_EXTRA may be empty. The environment travels
    # explicitly so a site default of SBATCH_EXPORT=NONE cannot strip it.
    set -f
    # shellcheck disable=SC2086
    jid="$(sbatch --parsable --export=ALL --kill-on-invalid-dep=yes "${sb[@]}" $EX_SBATCH_EXTRA \
                  "$EX_EXAMPLE_DIR/02_run_depthsv.sh" "$@")" \
        || { set +f; dsv_die "sbatch failed for $mode/$st"; }
    set +f
    jid="${jid%%;*}"               # federated clusters append ;clustername
    ex_record_job "$mode" "$st" "$jid"
    dsv_log "[$mode] submitted $st as job $jid"
    printf '%s\n' "$jid"
}

# The standard join submitted by this invocation, if any: the seed-control
# chain depends on it rather than joining the same matrix a second time.
join_jid_standard=""

slurm_chain() {                    # slurm_chain <mode>
    local mode="$1" j_join="" j_disp dep=()
    require_inputs "$mode"
    if [ "$mode" = seedctl ]; then
        if [ -n "$join_jid_standard" ]; then
            dep=(--dependency "afterok:$join_jid_standard")
        elif dsv_output_complete "$(ex_join_dir standard)/depth.matrix.txt.gz"; then
            dsv_log "[seedctl] reusing the completed standard matrix"
        else
            dsv_die "seedctl shares the standard matrix: run --mode standard (or all) first"
        fi
    else
        # shellcheck disable=SC2086
        j_join="$(submit "$mode" join $EX_SBATCH_JOIN \
            --job-name "dsvx-join-$mode" --output "$EX_LOG_DIR/%x.%j.out" \
            -- --stage join-exec --mode "$mode")"
        [ "$mode" != standard ] || join_jid_standard="$j_join"
        dep=(--dependency "afterok:$j_join")
    fi
    # shellcheck disable=SC2086
    j_disp="$(submit "$mode" dispatch $EX_SBATCH_LIGHT ${dep[@]+"${dep[@]}"} \
        --job-name "dsvx-dispatch-$mode" --output "$EX_LOG_DIR/%x.%j.out" \
        -- --stage dispatch --mode "$mode")"
    dsv_log "[$mode] chain submitted: join=${j_join:-shared} dispatch=$j_disp (array + evaluate follow automatically)"
}

# --- entry points ----------------------------------------------------------

case "$stage" in

driver)
    [ -s "$EX_INPUTS_DIR/run.env" ] || dsv_die "no $EX_INPUTS_DIR/run.env: run 01_prepare_inputs.sh first (it freezes ndim and the covariates for this run)"
    # A preamble still running would rewrite ndim.txt under the array; the
    # frozen run.env protects this run, but the user almost certainly wants
    # its result, so say so and stop.
    if [ "$(ex_runner)" = slurm ] && command -v squeue >/dev/null 2>&1 \
       && squeue --me --name=dsvx-preamble -h 2>/dev/null | grep -q .; then
        dsv_die "a dsvx-preamble job is queued or running; let it finish, rerun 01_prepare_inputs.sh to pick up its ndim and covariates, then submit"
    fi
    modes_selected | tr '\n' ' ' > "$selected_file"
    rm -rf "$lock_dir"
    case "$(ex_runner)" in
    local)
        rc=0
        for mode in $(cat "$selected_file"); do
            require_inputs "$mode"
            do_join "$mode"
            ex_timed "$mode" regions regions -- do_regions "$mode"
            regions_file="$(ex_regions_file "$mode")"
            n=0; total="$(grep -c . "$regions_file")"
            while IFS= read -r region; do
                n=$((n + 1))
                dsv_log "[$mode] unit $n/$total: $region"
                do_unit "$mode" "$region"
            done < "$regions_file"
            # An evaluation FAIL must not stop the other modes or the
            # comparison — the failure is easiest to read beside them.
            bash "$EX_EXAMPLE_DIR/03_evaluate.sh" --mode "$mode" || rc=1
        done
        maybe_finalize || rc=1
        exit "$rc"
        ;;
    slurm)
        for mode in $(cat "$selected_file"); do slurm_chain "$mode"; done
        dsv_log "monitor with:  squeue --me | grep dsvx-"
        dsv_log "results land in $EX_EVAL_DIR, $EX_COMPARE_DIR and $EX_PROFILE_DIR"
        ;;
    esac
    ;;

join-exec)
    do_join "$mode_arg"
    ;;

dispatch)
    ex_timed "$mode_arg" regions regions -- do_regions "$mode_arg"
    regions_file="$(ex_regions_file "$mode_arg")"
    n="$(grep -c . "$regions_file")"
    # shellcheck disable=SC2086
    j_arr="$(submit "$mode_arg" unit-array $EX_SBATCH_UNIT \
        --array "1-${n}%${EX_ARRAY_THROTTLE}" \
        --job-name "dsvx-unit-$mode_arg" --output "$EX_LOG_DIR/%x.%A_%a.out" \
        -- --stage unit --mode "$mode_arg")"
    # shellcheck disable=SC2086
    submit "$mode_arg" evaluate $EX_SBATCH_LIGHT \
        --dependency "afterok:$j_arr" \
        --job-name "dsvx-eval-$mode_arg" --output "$EX_LOG_DIR/%x.%j.out" \
        -- --stage eval-exec --mode "$mode_arg" > /dev/null
    ;;

unit)
    regions_file="$(ex_regions_file "$mode_arg")"
    dsv_require_file "$regions_file"
    idx="${SLURM_ARRAY_TASK_ID:-1}"
    region="$(sed -n "${idx}p" "$regions_file")"
    [ -n "$region" ] || dsv_die "no region on line $idx of $regions_file"
    dsv_log "[$mode_arg] task $idx -> $region"
    do_unit "$mode_arg" "$region"
    ;;

eval-exec)
    rc=0
    bash "$EX_EXAMPLE_DIR/03_evaluate.sh" --mode "$mode_arg" || rc=$?
    maybe_finalize || rc=1
    exit "$rc"
    ;;

*)
    dsv_die "unknown --stage '$stage'"
    ;;
esac
