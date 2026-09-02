#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# depthSV — association test one region against one phenotype
#
#   scripts/analyze.sh --corrected FILE --pheno FILE --model FORMULA \
#                      --method linear|logistic|coxph \
#                      --region chr1[:start-end] --out DIR
#
# This replaces the three near-identical wrappers that differed only by their
# regression method and their defaults. Everything site-specific belongs in an
# env file (see conf/example.env), not here.
#
# Phenotype sweeps: pass --pheno-manifest instead of --model/--method to run a
# tab-separated table of `name<TAB>method<TAB>model` rows against the same
# region. Each row writes its own output and is skipped if already complete, so
# adding a phenotype is a line in a file and a re-run costs only the new work.
# ---------------------------------------------------------------------------

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../lib/common.sh"
dsv_enable_error_trace

corrected=""; pheno="${DSV_PHENO:-}"; model=""; method="linear"; region=""
out_dir="${DSV_RESULTS_DIR:-}"; manifest="${DSV_PHENO_MANIFEST:-}"; name=""
jobs="${DSV_JOBS:-4}"; chunk="${DSV_CHUNK:-2000}"; threads="${DSV_THREADS:-2}"
min_obs="${DSV_MIN_OBS:-}"; min_var="${DSV_MIN_VARIANCE:-}"
pcs="${DSV_PCS:-}"; ndim=""; case_level=""; min_cases="${DSV_MIN_CASES:-}"
force=0; extra=(); model_flag=0; manifest_flag=0

while [ $# -gt 0 ]; do
    case "$1" in
        --corrected)      corrected="$2"; shift 2 ;;
        --pheno)          pheno="$2";     shift 2 ;;
        --model)          model="$2"; model_flag=1; shift 2 ;;
        --method)         method="$2";    shift 2 ;;
        --name)           name="$2";      shift 2 ;;
        --region)         region="$2";    shift 2 ;;
        --out)            out_dir="$2";   shift 2 ;;
        --pheno-manifest) manifest="$2"; manifest_flag=1; shift 2 ;;
        --pcs)            pcs="$2";       shift 2 ;;
        --ndim)           ndim="$2";      shift 2 ;;
        --case-level)     case_level="$2"; shift 2 ;;
        --min-cases)      min_cases="$2"; shift 2 ;;
        --jobs)           jobs="$2";      shift 2 ;;
        --chunk)          chunk="$2";     shift 2 ;;
        --threads)        threads="$2";   shift 2 ;;
        --min-obs)        min_obs="$2";   shift 2 ;;
        --min-variance)   min_var="$2";   shift 2 ;;
        --force)          force=1;        shift ;;
        --)               shift; extra=("$@"); break ;;
        -h|--help)        dsv_usage ;;
        *)                dsv_die "unknown argument: $1" ;;
    esac
done

# A --model on the command line is a single-analysis run even when the
# environment names a manifest: flags win over the environment.
if [ "$model_flag" -eq 1 ] && [ "$manifest_flag" -eq 0 ]; then manifest=""; fi
dsv_require_opt corrected pheno region out_dir
[ -z "$manifest" ] || [ -z "$model" ] \
    || dsv_die "--model and --pheno-manifest are mutually exclusive"
dsv_load_modules
dsv_require_file "$corrected" "$pheno"
dsv_require_cmd bgzip tabix parallel Rscript

# The PCs removed by the correction must be in the association model, or the
# test is deflated and — with covariates — biased. The count is read from the
# corrected file's name (corrected_ndim<N>.<region>) unless given, and the
# table from DSV_PCS or --pcs. --ndim 0 is the explicit way to opt out.
if [ -z "$ndim" ]; then
    ndim="$(printf '%s' "$(basename "$corrected")" | sed -nE 's/^corrected_ndim([0-9]+)\..*/\1/p')"
    [ -n "$ndim" ] || ndim=0
fi
if [ "$ndim" -gt 0 ]; then
    [ -n "$pcs" ] || dsv_die "the corrected matrix removed $ndim PCs; pass --pcs (or set DSV_PCS) so the test conditions on them, or --ndim 0 to opt out"
    dsv_require_file "$pcs"
fi

# Fold the thresholds and the PC arguments into the pass-through, BEFORE
# whatever the user put after `--`: optparse keeps the last occurrence, so
# the user's flag wins over the environment.
front=(--minObs "${min_obs:-100}")
[ -z "$min_var" ]    || front+=(--minVariance "$min_var")
[ -z "$min_cases" ]  || front+=(--minCases "$min_cases")
[ -z "$case_level" ] || front+=(--caseLevel "$case_level")
[ "$ndim" -eq 0 ]    || front+=(--pcs "$pcs" --ndim "$ndim")
extra=(${front[@]+"${front[@]}"} ${extra[@]+"${extra[@]}"})

rscript="${DSV_ANALYZE_R:-$DSV_ROOT/R/analyze.R}"
dsv_require_file "$rscript"
mkdir -p "$out_dir"
# Absolute from here on: join.sh cds into a working subdirectory, and a
# relative path captured beforehand would stop resolving after that.
out_dir="$(cd "$out_dir" && pwd)"

slug="$(printf '%s' "$region" | tr ':' '_' | tr -d ' ')"

header_file="$(mktemp "${TMPDIR:-/tmp}/dsv.hdr.XXXXXX")"
probe_row="$(mktemp "${TMPDIR:-/tmp}/dsv.row.XXXXXX")"
out_header="$(mktemp "${TMPDIR:-/tmp}/dsv.outhdr.XXXXXX")"
trap 'rm -f "$header_file" "$probe_row" "$out_header"' EXIT

dsv_header "$corrected" > "$header_file"

# One real data row: proves the region is non-empty, sizes the parallel chunk,
# and probes each analysis's output header. Reading one row rather than counting
# them all matters — the count was only ever a log denominator, and a full pass
# over a region is minutes of decompression at biobank width.
set +o pipefail
dsv_read_region "$corrected" "$region" | head -n 1 > "$probe_row"
set -o pipefail
[ -s "$probe_row" ] || dsv_die "region $region contains no rows"
chunk="$(dsv_chunk_lines "$probe_row" "$chunk")"

# --- one phenotype ---------------------------------------------------------

# Passthrough args are re-quoted with %q because the worker command below is a
# string a worker shell parses; interpolating them raw broke any argument
# containing shell metacharacters (every realistic --sampleIdPattern).
extra_q=""
[ ${#extra[@]} -eq 0 ] || extra_q="$(printf '%q ' "${extra[@]}")"

run_one() {                        # run_one <name> <method> <model>
    local name="$1" meth="$2" form="$3"
    local final="$out_dir/${name}.${meth}.${slug}.txt.gz"
    local log="$out_dir/${name}.${meth}.${slug}.log"

    if [ "$force" -eq 0 ] && dsv_output_complete "$final"; then
        dsv_log "already complete, skipping: $(basename "$final")"
        return 0
    fi
    dsv_output_reset "$final"

    # R worker diagnostics — the [align] sample-drop counts and any real
    # error — go to a per-analysis log rather than /dev/null.
    : > "$log"

    set +o pipefail
    cat "$header_file" "$probe_row" \
      | Rscript "$rscript" -f - -p "$pheno" -m "$form" -r "$meth" \
          ${extra[@]+"${extra[@]}"} 2>>"$log" | head -n 1 > "$out_header"
    set -o pipefail
    [ -s "$out_header" ] || dsv_die "could not derive the output header for $name; check the model formula against the phenotype table (see $log)"

    # The probe succeeded, so what it logged is only the SIGPIPE noise from
    # having its output truncated by `head`. Clear it: the log should hold the
    # workers' real diagnostics, not an expected artefact.
    : > "$log"

    # Regions are legitimately dropped by the QC filters, so the row count is
    # not a parity check here. Rscript is the LAST command of the worker
    # pipeline on purpose: with a trailing `| tail` GNU parallel only ever saw
    # tail's exit status and an R crash part-way through a chunk went
    # unnoticed. The header is suppressed in R instead.
    local worker
    worker="cat $(dsv_q "$header_file") - | Rscript $(dsv_q "$rscript") \
      -f - -p $(dsv_q "$pheno") -m $(dsv_q "$form") -r $(dsv_q "$meth") --skipOutputHeader \
      ${extra_q}2>>$(dsv_q "$log")"
    (
        cat "$out_header"
        dsv_read_region "$corrected" "$region" \
          | parallel "${DSV_PARALLEL_FLAGS[@]}" --block "$DSV_BLOCK_BYTES" -L "$chunk" -j "$jobs" "$worker"
    ) | bgzip -@ "$threads" > "$(dsv_output_tmp "$final")"

    # A result with no rows is not a result. The QC filters can drop every
    # region of a small window legitimately, but at that point the caller
    # must say so explicitly rather than inherit a green run.
    if [ "$(bgzip -dc "$(dsv_output_tmp "$final")" | grep -cv '^#')" -eq 0 ] && [ "${DSV_ALLOW_EMPTY:-0}" != "1" ]; then
        dsv_die "$name produced no result rows for $region (every region failed QC?). Set DSV_ALLOW_EMPTY=1 to accept an empty shard."
    fi

    dsv_output_commit "$final"
}

# --- dispatch --------------------------------------------------------------
if [ -n "$manifest" ]; then
    dsv_require_file "$manifest"
    n=0
    # Fields are cut positionally rather than read with IFS=$'\t'. Tab is an
    # IFS *whitespace* character, so `read` collapses runs of tabs and an empty
    # field silently shifts every later field left — which would hand the wrong
    # string to --model. The manifest is also read on fd 3 so nothing inside
    # the loop can consume it.
    while IFS= read -r line <&3; do
        case "$line" in ''|'#'*) continue ;; esac
        name="$(printf '%s' "$line" | cut -f1)"
        meth="$(printf '%s' "$line" | cut -f2)"
        form="$(printf '%s' "$line" | cut -f3)"
        [ -n "$name" ] || dsv_die "manifest row has an empty name: $line"
        [ -n "$meth" ] || dsv_die "manifest row '$name' is missing a method (expected name<TAB>method<TAB>model)"
        [ -n "$form" ] || dsv_die "manifest row '$name' is missing a model formula"
        n=$((n+1))
        run_one "$name" "$meth" "$form"
    done 3< "$manifest"
    dsv_log "manifest complete: $n phenotypes over $region"
else
    [ -n "$model" ] || dsv_die "--model is required (or use --pheno-manifest)"
    # Name the analysis after its response variable, not the phenotype file:
    # two models against the same file would otherwise resolve to the same
    # output path and the second would be skipped as already complete.
    # Non-alphanumerics are stripped so Surv(time,event) yields a usable name.
    run_one "${name:-$(printf '%s' "${model%%~*}" | tr -cd '[:alnum:]_.')}" "$method" "$model"
fi
