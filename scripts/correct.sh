#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# depthSV — normalise and correct one region of the depth matrix
#
#   scripts/correct.sh --matrix FILE --pcs FILE --coverage FILE \
#                      --region chr1[:start-end] --out DIR [--ndim N]
#
# The unit of work is a region, not a chromosome. That is what lets the same
# script be driven by a SLURM array, a WDL scatter or a CSV fan-out without any
# of them knowing about the others, and it is what makes a preempted job cheap
# to retry.
#
# Input regions are read through the tabix index rather than by scanning, so
# this reads only what it needs.
# ---------------------------------------------------------------------------

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../lib/common.sh"
dsv_enable_error_trace

matrix="${DSV_MATRIX:-}"; pcs="${DSV_PCS:-}"; coverage="${DSV_COVERAGE:-}"
region=""; out_dir="${DSV_CORRECTED_DIR:-}"
ndim="${DSV_NDIM:-16}"; jobs="${DSV_JOBS:-4}"; chunk="${DSV_CHUNK:-2000}"
threads="${DSV_THREADS:-2}"; force=0; extra=()

while [ $# -gt 0 ]; do
    case "$1" in
        --matrix)   matrix="$2";   shift 2 ;;
        --pcs)      pcs="$2";      shift 2 ;;
        --coverage) coverage="$2"; shift 2 ;;
        --region)   region="$2";   shift 2 ;;
        --out)      out_dir="$2";  shift 2 ;;
        --ndim)     ndim="$2";     shift 2 ;;
        --jobs)     jobs="$2";     shift 2 ;;
        --chunk)    chunk="$2";    shift 2 ;;
        --threads)  threads="$2";  shift 2 ;;
        --force)    force=1;       shift ;;
        --)         shift; extra=("$@"); break ;;
        -h|--help)  dsv_usage ;;
        *)          dsv_die "unknown argument: $1" ;;
    esac
done

dsv_require_opt matrix pcs coverage region out_dir
dsv_require_file "$matrix" "$pcs" "$coverage"
dsv_require_cmd bgzip tabix parallel Rscript
dsv_load_modules

rscript="${DSV_CORRECT_R:-$DSV_ROOT/R/correct.R}"
dsv_require_file "$rscript"

mkdir -p "$out_dir"
# Absolute from here on: join.sh cds into a working subdirectory, and a
# relative path captured beforehand would stop resolving after that.
out_dir="$(cd "$out_dir" && pwd)"
slug="$(printf '%s' "$region" | tr ':' '_' | tr -d ' ')"
final="$out_dir/corrected_ndim${ndim}.${slug}.txt.gz"
stats="$out_dir/stats_ndim${ndim}.${slug}.txt"

if [ "$force" -eq 0 ] && dsv_output_complete "$final"; then
    dsv_log "already complete, skipping: $final"
    exit 0
fi
dsv_output_reset "$final"

# R worker diagnostics — the [align] sample-drop counts, rank warnings and any
# real error — go to a per-unit log rather than /dev/null. A wrong
# --sampleIdPattern that silently drops most of the cohort is visible here.
log="$out_dir/corrected_ndim${ndim}.${slug}.log"
: > "$log"

# Column header of the matrix, needed by every parallel chunk.
header_file="$(mktemp "${TMPDIR:-/tmp}/dsv.header.XXXXXX")"
stats_dir="$(mktemp -d "${TMPDIR:-/tmp}/dsv.stats.XXXXXX")"
trap 'rm -rf "$header_file" "$stats_dir"' EXIT

dsv_header "$matrix" > "$header_file"

n_in="$(dsv_read_region "$matrix" "$region" | wc -l | tr -d ' ')"
dsv_log "region $region: $n_in input rows, ndim=$ndim"
[ "$n_in" -gt 0 ] || dsv_die "region $region contains no rows"

# One real data row, used both to probe the output header and to size the
# parallel chunk. Read once rather than once per purpose.
probe_row="$(mktemp "${TMPDIR:-/tmp}/dsv.row.XXXXXX")"
out_header="$(mktemp "${TMPDIR:-/tmp}/dsv.outhdr.XXXXXX")"
trap 'rm -rf "$header_file" "$stats_dir" "$probe_row" "$out_header"' EXIT
set +o pipefail
dsv_read_region "$matrix" "$region" | head -n 1 > "$probe_row"
# Note the ${arr[@]+"${arr[@]}"} idiom: under `set -u`, expanding an empty
# array is an error on bash 3.2, which is what macOS still ships.
{ cat "$header_file" "$probe_row"; } \
  | Rscript "$rscript" --inputPCs "$pcs" --inputFile - --coverageStats "$coverage" \
      --ndim "$ndim" ${extra[@]+"${extra[@]}"} 2>>"$log" | head -n 1 > "$out_header"
set -o pipefail
[ -s "$out_header" ] || dsv_die "could not derive the output header from $rscript (see $log)"

# The probe succeeded, so what it logged is only the SIGPIPE noise from having
# its output truncated by `head`. Clear it: the log should hold the workers'
# real diagnostics, not an expected artefact that reads like an error.
: > "$log"

chunk="$(dsv_chunk_lines "$probe_row" "$chunk")"

# The worker command is a string a worker shell parses, so every interpolated
# value is quoted with dsv_q. {#} is GNU parallel's job number and must stay
# bare for the per-chunk stats files.
extra_q=""
[ ${#extra[@]} -eq 0 ] || extra_q="$(printf '%q ' "${extra[@]}")"
worker="cat $(dsv_q "$header_file") - | Rscript $(dsv_q "$rscript") \
  --inputPCs $(dsv_q "$pcs") --inputFile - --coverageStats $(dsv_q "$coverage") \
  --ndim $(dsv_q "$ndim") --skipOutputHeader \
  --statsFile $(dsv_q "$stats_dir")/stats.{#}.txt ${extra_q}2>>$(dsv_q "$log")"

(
    cat "$out_header"
    dsv_read_region "$matrix" "$region" \
      | parallel "${DSV_PARALLEL_FLAGS[@]}" --block "$DSV_BLOCK_BYTES" -L "$chunk" -j "$jobs" "$worker"
) | bgzip -@ "$threads" > "$(dsv_output_tmp "$final")"

# Merge the per-chunk statistics BEFORE committing the main output. The .done
# marker written by dsv_output_commit is what makes a re-run skip this unit, so
# anything produced after it would be permanently missing if the job were
# interrupted in between.
if compgen -G "$stats_dir/stats.*.txt" > /dev/null; then
    head -n 1 "$(ls "$stats_dir"/stats.*.txt | head -1)" > "$stats"
    for f in "$stats_dir"/stats.*.txt; do tail -n +2 "$f"; done >> "$stats"
    bgzip -f "$stats" && tabix -f -p bed "${stats}.gz"
    dsv_log "stats: ${stats}.gz"
fi

dsv_output_commit "$final" "$n_in"
