#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# depthSV — export one analysis from its shards
#
#   scripts/export.sh --results DIR --regions FILE --out DIR \
#                     (--name NAME --method METHOD | --pheno-manifest FILE) \
#                     [--min-count N] [--alpha X] [--allow-missing] [--force]
#
# The association stage writes one shard per region. This is the step that
# turns them into a result:
#
#   coverage     every region in the list must have a finished shard (a
#                `.done` marker); a missing one is an error unless
#                --allow-missing, so a lost array task cannot vanish quietly.
#   order        shards are concatenated in region-list order, which is
#                coordinate order, then bgzipped and indexed.
#   suppression  rows whose N, NCase or NControl is below --min-count are
#                dropped and counted. All of Us forbids disseminating counts
#                of 1-20 or anything they can be derived from; 20 is the
#                default and 0 disables it.
#   threshold    the per-shard permutation maxima (analyze.sh --perms) are
#                folded into one genome-wide max-|t| distribution; from it
#                R/export.R derives the empirical family-wise threshold at
#                --alpha, the effective number of tests it implies, and the
#                adjusted p of every hit — beside Bonferroni and lambda.
#
# Outputs, per analysis, under --out:
#   <name>.<method>.txt.gz (+ index)   the exported table
#   <name>.<method>.summary.tsv        counts, lambda, thresholds, M_eff
#   <name>.<method>.hits.tsv           rows past either threshold, with P_ADJ
#   <name>.<method>.permmax.txt        the folded permutation maxima
# ---------------------------------------------------------------------------

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../lib/common.sh"
dsv_enable_error_trace

results="${DSV_RESULTS_DIR:-}"; regions=""; out_dir="${DSV_EXPORT_DIR:-}"
name=""; method="linear"; manifest="${DSV_PHENO_MANIFEST:-}"; name_flag=0
min_count="${DSV_MIN_COUNT:-20}"; alpha="${DSV_ALPHA:-0.05}"
allow_missing=0; threads="${DSV_THREADS:-2}"; force=0

while [ $# -gt 0 ]; do
    case "$1" in
        --results)        results="$2";   shift 2 ;;
        --regions)        regions="$2";   shift 2 ;;
        --out)            out_dir="$2";   shift 2 ;;
        --name)           name="$2"; name_flag=1; shift 2 ;;
        --method)         method="$2";    shift 2 ;;
        --pheno-manifest) manifest="$2";  shift 2 ;;
        --min-count)      min_count="$2"; shift 2 ;;
        --alpha)          alpha="$2";     shift 2 ;;
        --allow-missing)  allow_missing=1; shift ;;
        --threads)        threads="$2";   shift 2 ;;
        --force)          force=1;        shift ;;
        -h|--help)        dsv_usage ;;
        *)                dsv_die "unknown argument: $1" ;;
    esac
done

# A --name on the command line exports that one analysis even when the
# environment names a manifest: flags win over the environment.
if [ "$name_flag" -eq 1 ]; then manifest=""; fi
dsv_require_opt results regions out_dir
[ -n "$manifest" ] || [ -n "$name" ] || dsv_die "pass --name/--method or --pheno-manifest"
dsv_load_modules
dsv_require_cmd bgzip tabix Rscript
dsv_require_file "$regions"
[ -d "$results" ] || dsv_die "no results directory: $results"

rscript="${DSV_EXPORT_R:-$DSV_ROOT/R/export.R}"
dsv_require_file "$rscript"
mkdir -p "$out_dir"
out_dir="$(cd "$out_dir" && pwd)"
results="$(cd "$results" && pwd)"
n_regions="$(grep -c . "$regions" || true)"
[ "$n_regions" -gt 0 ] || dsv_die "the region list is empty: $regions"

export_one() {                     # export_one <name> <method>
    local name="$1" meth="$2"
    local final="$out_dir/${name}.${meth}.txt.gz"
    local summary="$out_dir/${name}.${meth}.summary.tsv"
    local hits="$out_dir/${name}.${meth}.hits.tsv"
    local permmax="$out_dir/${name}.${meth}.permmax.txt"
    local shards=() perms=() missing=() region slug shard

    while IFS= read -r region; do
        [ -n "$region" ] || continue
        slug="$(printf '%s' "$region" | tr ':' '_' | tr -d ' ')"
        shard="$results/${name}.${meth}.${slug}.txt.gz"
        if [ -s "$shard" ] && [ -f "$(dsv_output_done "$shard")" ]; then
            shards+=("$shard")
            [ ! -s "${shard%.txt.gz}.permmax.txt" ] || perms+=("${shard%.txt.gz}.permmax.txt")
        else
            missing+=("$region")
        fi
    done < "$regions"

    if [ ${#missing[@]} -gt 0 ]; then
        if [ "$allow_missing" -eq 1 ]; then
            dsv_log "WARNING: $name.$meth: ${#missing[@]} of $n_regions regions have no finished shard (e.g. ${missing[0]}); exporting the rest"
        else
            dsv_die "$name.$meth: ${#missing[@]} of $n_regions regions have no finished shard, e.g. $(printf '%s ' "${missing[@]:0:5}")-- rerun those units, or pass --allow-missing"
        fi
    fi
    [ ${#shards[@]} -gt 0 ] || dsv_die "$name.$meth: no finished shards under $results"

    # The shards' own signatures identify the inputs; the rest is this step's.
    local sig
    sig="$(printf 'stage=export\nshards=%s\nshard_params=%s\nmin_count=%s\nalpha=%s\nallow_missing=%s\nscript=%s\nrscript=%s' \
           "${#shards[@]}" \
           "$(for s in "${shards[@]}"; do cat "$(dsv_output_params "$s")" 2>/dev/null || true; echo; done | cksum | awk '{print $1}')" \
           "$min_count" "$alpha" "$allow_missing" "$(dsv_script_sig "$0")" "$(dsv_script_sig "$rscript")")"
    if [ "$force" -eq 0 ] && dsv_output_complete "$final" "$sig"; then
        dsv_log "already complete, skipping: $(basename "$final")"
        return 0
    fi
    dsv_output_reset "$final" "$force"
    rm -f "$summary" "$hits" "$permmax"

    local counts tmp
    counts="$(mktemp "${TMPDIR:-/tmp}/dsv.counts.XXXXXX")"
    tmp="$(dsv_output_tmp "$final")"

    # Header from the first shard, bodies in region-list order. Columns 5-7
    # are N, NCase, NControl for every method.
    {
        dsv_header "${shards[0]}"
        for shard in "${shards[@]}"; do bgzip -dc "$shard" | grep -v '^#' || true; done
    } | awk -F'\t' -v mc="$min_count" -v cnt="$counts" '
        BEGIN { OFS = "\t" }
        /^#/  { print; next }
        { total++
          if (mc > 0 && ($5+0 < mc || $6+0 < mc || $7+0 < mc)) { sup++; next }
          print }
        END { printf "rows_in\t%d\nrows_suppressed\t%d\n", total+0, sup+0 > cnt }' \
      | bgzip -@ "$threads" > "$tmp"

    # Permutation maxima: fold the shards, but only a complete set is a
    # genome-wide distribution.
    local n_perm=0 perm_opt=()
    if [ ${#perms[@]} -gt 0 ]; then
        if [ ${#perms[@]} -ne ${#shards[@]} ]; then
            dsv_log "WARNING: $name.$meth: permutation maxima exist for ${#perms[@]} of ${#shards[@]} shards; no empirical threshold"
        else
            {
                grep '^#' "${perms[0]}" | grep -v '^#regions='
                printf '#shards=%d\n' "${#perms[@]}"
                printf 'perm\tmax_abs_stat\n'
                dsv_perm_fold "${perms[@]}"
            } > "$permmax"
            n_perm="$(grep -v -e '^#' -e '^perm' "$permmax" | grep -c . || true)"
            [ "$n_perm" -gt 0 ] && perm_opt=(--permmax "$permmax") || rm -f "$permmax"
        fi
    fi

    # Provenance travels with the summary: the pipeline version and commit,
    # R and htslib. The exported table itself carries no meta-lines, so any
    # reader that takes the first line as the header keeps working.
    Rscript "$rscript" --input "$tmp" --counts "$counts" --alpha "$alpha" --minCount "$min_count" \
        --name "$name" --method "$meth" --regionsListed "$n_regions" --shards "${#shards[@]}" \
        --version "$(dsv_version)" --htslib "$(tabix --version 2>/dev/null | head -n 1 || echo unknown)" \
        ${perm_opt[@]+"${perm_opt[@]}"} --out "$summary" --hits "$hits" \
        || dsv_die "$name.$meth: the summary step failed"
    rm -f "$counts"

    dsv_output_commit "$final" "" "$sig"
    dsv_log "summary: $summary"
}

if [ -n "$manifest" ]; then
    dsv_require_file "$manifest"
    n=0
    while IFS= read -r line <&3; do
        case "$line" in ''|'#'*) continue ;; esac
        name="$(printf '%s' "$line" | cut -f1)"
        meth="$(printf '%s' "$line" | cut -f2)"
        [ -n "$name" ] && [ -n "$meth" ] || dsv_die "manifest row is missing a name or method: $line"
        n=$((n+1))
        export_one "$name" "$meth"
    done 3< "$manifest"
    dsv_log "exported $n analyses -> $out_dir"
else
    export_one "$name" "$method"
fi
