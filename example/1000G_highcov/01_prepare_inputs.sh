#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# depthSV 1000G example — stage 1: prepare per-mode depthSV inputs
#
#   bash 01_prepare_inputs.sh
#
# For each resolved mode this writes, under $EX_WORK_DIR/inputs/<mode>/:
#
#   svd.pcs.txt            PC table, sample suffix stripped to bare IDs
#   autosomal.median.txt   SAMPLE, AUTO_HQ_median — from NGS-PCA's own
#                          median table when the run wrote one, else from
#                          the QC table's HQ_MEDIAN_COV
#   phenotypes.tsv         MTDNA_CN, LOG2_MTDNA_CN, MTDNA_CN_NULL, SEX, ...
#   analyses.tsv           the phenotype manifest the sweep runs
#   mosdepth.manifest.txt  one absolute path per sample region file
#   samples.txt            the matrix column names that manifest implies
#   chrom.sizes            contig lengths read from the first region file
#
# The phenotypes carry their own truth: MTDNA_CN is 2 x chrM / HQ-median by
# construction upstream, so chrM must dominate its association; INFERRED_SEX
# must light up the sex chromosomes; the permuted null must stay flat. That
# is what 03_evaluate.sh asserts. The prepare step checks that the
# phenotype's denominator is the median this pipeline normalises against,
# and that the modes' sample sets agree, before any joining starts.
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

dsv_require_cmd Rscript awk gzip find sort comm

# Replace a generated file only when its content changed, so finished
# pipeline units keyed on these inputs are not invalidated by a rerun that
# produced the same bytes.
install_if_changed() {             # install_if_changed <tmp> <dest>
    if [ -f "$2" ] && cmp -s "$1" "$2"; then rm -f "$1"; else mv "$1" "$2"; fi
}

prepared=""
for mode in $EX_MODES; do
    ex_check_mode "$mode"
    paths_env="$(ex_paths_env "$mode")"
    if [ ! -s "$paths_env" ]; then
        dsv_log "SKIP $mode: no $paths_env (00_fetch_inputs.sh did not resolve it)"
        continue
    fi
    EX_M_MEDIAN_TABLE=""
    # shellcheck disable=SC1090
    source "$paths_env"
    in_dir="$(ex_inputs_dir "$mode")"
    mkdir -p "$in_dir"
    # Not ready until this block finishes: a failed prepare must not leave
    # the previous run's tables looking usable.
    rm -f "$in_dir/prepared.ok"
    dsv_log "=== mode: $mode (source: $EX_M_SOURCE) ==="

    # --- tables ------------------------------------------------------------
    median_opt=()
    [ -z "$EX_M_MEDIAN_TABLE" ] || median_opt=(--median "$EX_M_MEDIAN_TABLE")
    cov_opt=()
    [ ! -s "$EX_PREAMBLE_DIR/covariates.tsv" ] || cov_opt=(--covariates "$EX_PREAMBLE_DIR/covariates.tsv")
    # One permuted null for every mode: the permutation is drawn once (for the
    # first mode prepared, normally standard) and copied by sample ID, so a
    # tree with a different sample set gets the same null values, not a
    # different permutation.
    null_opt=()
    if [ "$mode" != standard ] && [ -s "$(ex_inputs_dir standard)/phenotypes.tsv" ]; then
        null_opt=(--null-from "$(ex_inputs_dir standard)/phenotypes.tsv")
    fi
    Rscript "$EX_EXAMPLE_DIR/R/prepare_inputs.R" \
        --qc "$EX_M_QC_TABLE" --pcs "$EX_M_PCS_FILE" \
        ${median_opt[@]+"${median_opt[@]}"} ${cov_opt[@]+"${cov_opt[@]}"} ${null_opt[@]+"${null_opt[@]}"} \
        --suffix "$EX_SAMPLE_SUFFIX" --seed "$EX_PHENO_SEED" \
        --out "$in_dir" 2>&1 | tee "$in_dir/prepare.summary.txt" >&2

    # --- phenotype manifest -------------------------------------------------
    # With covariates (preamble.sh ran): the unadjusted MTDNA_CN run stays as
    # the pure truth check and the adjusted runs are the science models.
    # The sex analyses take the covariates without SEX itself.
    adj=""; adj_nosex=""
    if [ "$EX_COVARIATES" != "none" ] && [ -n "$EX_COVARIATES" ]; then
        adj="+$EX_COVARIATES"
        adj_nosex="$(printf '%s' "+$EX_COVARIATES" | sed -e 's/+SEX+/+/g' -e 's/+SEX$//')"
        [ "$adj_nosex" = "+" ] && adj_nosex=""
    fi
    {
        printf '# depthSV 1000G example analyses. Format: name<TAB>method<TAB>model\n'
        printf '# The depth term must be named cov_resids (see conf/phenotypes.example.tsv).\n'
        printf '#\n'
        printf '# SEX appears twice on purpose. The linear run carries the truth check:\n'
        printf '# sex separates depth on chrX/Y so completely that the logistic Wald z\n'
        printf '# collapses there (Hauck-Donner), so the logistic run exercises that\n'
        printf '# engine and documents the collapse rather than asserting rank.\n'
        if [ -n "$adj" ]; then
            printf '# Covariates from the preamble: %s (EX_COVARIATES).\n' "$EX_COVARIATES"
            printf 'mtdna_cn\tlinear\tMTDNA_CN~cov_resids\n'
            printf 'mtdna_cn_adj\tlinear\tMTDNA_CN~cov_resids%s\n' "$adj"
            printf 'log2_mtdna_cn_adj\tlinear\tLOG2_MTDNA_CN~cov_resids%s\n' "$adj"
            printf 'mtdna_cn_null_adj\tlinear\tMTDNA_CN_NULL~cov_resids%s\n' "$adj"
            printf 'sex_linear_adj\tlinear\tSEX~cov_resids%s\n' "$adj_nosex"
            printf 'inferred_sex\tlogistic\tSEX~cov_resids\n'
        else
            printf 'mtdna_cn\tlinear\tMTDNA_CN~cov_resids\n'
            printf 'log2_mtdna_cn\tlinear\tLOG2_MTDNA_CN~cov_resids\n'
            printf 'mtdna_cn_null\tlinear\tMTDNA_CN_NULL~cov_resids\n'
            printf 'sex_linear\tlinear\tSEX~cov_resids\n'
            printf 'inferred_sex\tlogistic\tSEX~cov_resids\n'
        fi
    } > "$in_dir/analyses.tsv.tmp"
    install_if_changed "$in_dir/analyses.tsv.tmp" "$in_dir/analyses.tsv"
    dsv_log "$mode: analyses -> $(grep -vc '^#' "$in_dir/analyses.tsv") models (covariates: $EX_COVARIATES; ndim $EX_NDIM)"

    # --- mosdepth manifest --------------------------------------------------
    manifest="$in_dir/mosdepth.manifest.txt"
    find "$EX_M_MOSDEPTH_DIR" -maxdepth 1 -name '*.regions.bed.gz' 2>/dev/null | sort > "$manifest.tmp" || true
    install_if_changed "$manifest.tmp" "$manifest"
    n_md="$(grep -c . "$manifest" || true)"

    if [ "$n_md" -eq 0 ]; then
        dsv_log "NOTE $mode: no mosdepth files under $EX_M_MOSDEPTH_DIR; tables are ready but the join stage cannot run yet"
        continue
    fi
    dsv_log "$mode: $n_md region files -> mosdepth.manifest.txt"

    # The matrix columns will be dsv_sample_name of each path; verify they
    # actually meet the coverage table before hours of joining, because an
    # ID-suffix mismatch is the one way this wiring silently shrinks to an
    # empty cohort.
    while IFS= read -r f; do dsv_sample_name "$f"; done < "$manifest" | sort > "$in_dir/samples.txt.tmp"
    install_if_changed "$in_dir/samples.txt.tmp" "$in_dir/samples.txt"
    matched="$(
        join "$in_dir/samples.txt" <(awk -F'\t' 'NR>1{print $1}' "$in_dir/autosomal.median.txt" | sort) | grep -c . || true
    )"
    dsv_log "$mode: $matched of $n_md mosdepth samples present in the coverage table"
    [ "$matched" -gt 0 ] || dsv_die "no mosdepth sample matches the coverage table; check EX_SAMPLE_SUFFIX ('$EX_SAMPLE_SUFFIX') and the mosdepth file naming"
    if [ $((matched * 10)) -lt $((n_md * 9)) ]; then
        dsv_log "WARN $mode: below 90% overlap; the correction stage will drop the unmatched samples"
    fi

    # --- contig lengths ------------------------------------------------------
    # mosdepth's final bin is truncated at the contig end, so the largest STOP
    # per contig is the contig length. regions.sh needs these to window.
    first="$(awk 'NF {print; exit}' "$manifest")"
    gzip -cd "$first" | awk -F'\t' '
        !($1 in max) { order[++n] = $1 }
        $3 > max[$1] { max[$1] = $3 }
        END { for (i = 1; i <= n; i++) print order[i] "\t" max[order[i]] }
    ' > "$in_dir/chrom.sizes.tmp"
    install_if_changed "$in_dir/chrom.sizes.tmp" "$in_dir/chrom.sizes"
    dsv_log "$mode: $(grep -c . "$in_dir/chrom.sizes") contigs -> chrom.sizes"

    : > "$in_dir/prepared.ok"
    prepared="$prepared $mode"
done

[ -n "$prepared" ] || dsv_die "no mode was prepared"

# --- freeze this run's parameters ---------------------------------------------
# Every later job re-sources config.sh, which re-derives EX_NDIM from
# preamble/ndim.txt and EX_COVARIATES from covariates.tsv at that moment.
# The values resolved HERE are what the tables above were built with, so
# they are written as defaults that lib.sh loads before config.sh: an
# explicit environment still wins, config defaults no longer do. Rerun this
# stage after the preamble to re-freeze.
{
    printf '# frozen by 01_prepare_inputs.sh %s - explicit environment still wins\n' "$(dsv_now)"
    for v in EX_SMOKE EX_MODES EX_NDIM EX_COVARIATES EX_N_GPCS EX_WINDOW EX_MIN_OBS \
             EX_MEDIAN_SOURCE EX_CONTIG_REGEX EX_PHENO_SEED EX_EVAL_PROFILE EX_SAMPLE_SUFFIX; do
        # A guarded plain assignment: %q quoting is exact there, whereas inside
        # a double-quoted ${VAR:=word} the backslashes would survive.
        printf '[ -n "${%s:+x}" ] || %s=%q\n' "$v" "$v" "${!v}"
    done
} > "$EX_INPUTS_DIR/run.env.tmp"
install_if_changed "$EX_INPUTS_DIR/run.env.tmp" "$EX_INPUTS_DIR/run.env"
dsv_log "frozen for this run: ndim=$EX_NDIM covariates=$EX_COVARIATES window=$EX_WINDOW min-obs=$EX_MIN_OBS -> $EX_INPUTS_DIR/run.env"

# --- cross-mode sample sets --------------------------------------------------
# The comparison joins on region, so differing sample sets do not break it —
# but they change what "the same association" means, and a lagging fast tree
# (samples processed before COMPARE_FAST_MODE was on are re-queued upstream)
# is the usual cause. Say so here, once, with numbers.
set -- $prepared
ref="$1"
report="$EX_INPUTS_DIR/cross_mode_samples.txt"
{
    printf '# samples per prepared mode, and differences from %s\n' "$ref"
    for m in "$@"; do
        printf '%s\t%s\n' "$m" "$(grep -c . "$(ex_inputs_dir "$m")/samples.txt")"
    done
    for m in "$@"; do
        [ "$m" != "$ref" ] || continue
        only_ref="$(comm -23 "$(ex_inputs_dir "$ref")/samples.txt" "$(ex_inputs_dir "$m")/samples.txt" | grep -c . || true)"
        only_m="$(comm -13 "$(ex_inputs_dir "$ref")/samples.txt" "$(ex_inputs_dir "$m")/samples.txt" | grep -c . || true)"
        printf 'only_in_%s\t%s\nonly_in_%s\t%s\n' "$ref" "$only_ref" "$m" "$only_m"
        if [ "$only_ref" -gt 0 ] || [ "$only_m" -gt 0 ]; then
            dsv_log "WARN sample sets differ: $only_ref only in $ref, $only_m only in $m (see $report)"
        else
            dsv_log "$m: sample set identical to $ref"
        fi
    done
} > "$report"

dsv_log "done. Next: bash 02_run_depthsv.sh"
