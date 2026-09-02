#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# depthSV 1000G example — stage 0: resolve the upstream NGS-PCA inputs
#
#   bash 00_fetch_inputs.sh [--force]
#
# For each mode in EX_MODES this locates the NGS-PCA outputs depthSV needs —
# svd.pcs.txt, sample_qc.tsv and, when the run wrote one, autosomal.median.txt
# — preferring the local trees under NGSPCA_WORK_DIR and falling back to the
# results committed in the NGS-PCA repository on GitHub. What was resolved
# is recorded in inputs/<mode>/paths.env for the later stages. A mode whose
# inputs cannot be found is skipped, not fatal: the seed control
# (standard tree, PCs from NGS-PCA's step-2b reseed) simply may not exist yet.
#
# Mosdepth region files are never fetched: they are not committed upstream.
# They come from the local NGS-PCA run — or, with EX_SMOKE=1, this stage
# simulates a small per-mode tree consistent with each sample's real QC row
# (HQ median, X/Y ratios, chrM ratio), which is what makes the whole example
# runnable, positive controls included, on a laptop or a login node.
# ---------------------------------------------------------------------------

EX_EXAMPLE_DIR="${EX_EXAMPLE_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}}"
source "$EX_EXAMPLE_DIR/lib.sh"
dsv_enable_error_trace

force=0
while [ $# -gt 0 ]; do
    case "$1" in
        --force)   force=1; shift ;;
        -h|--help) dsv_usage ;;
        *)         dsv_die "unknown argument: $1" ;;
    esac
done

dsv_require_cmd curl awk gzip find
mkdir -p "$EX_WORK_DIR" "$EX_CACHE_DIR" "$EX_INPUTS_DIR"

# --- helpers ---------------------------------------------------------------

# ex_fetch <url> <dest> <mustContain>
# Download one file unless already cached; verify the header mentions the
# expected column so an HTML error page can never impersonate a table.
ex_fetch() {
    local url="$1" dest="$2" must="$3"
    if [ "$force" -eq 0 ] && [ -s "$dest" ]; then
        dsv_log "cached: $dest"
        return 0
    fi
    mkdir -p "$(dirname "$dest")"
    if ! curl -fsSL --retry 3 --retry-delay 2 -o "$dest.tmp" "$url" 2>/dev/null; then
        rm -f "$dest.tmp"
        return 1
    fi
    if ! head -n 1 "$dest.tmp" | grep -q "$must"; then
        dsv_log "unexpected content from $url (no '$must' in the header)"
        rm -f "$dest.tmp"
        return 1
    fi
    mv "$dest.tmp" "$dest"
    dsv_log "fetched: $dest"
}

# github_fetch_mode <mode> — pull svd.pcs.txt + sample_qc.tsv into the cache.
# Returns non-zero when that mode's results are not in the repository (the
# fast tree may simply not be committed yet; the seed control never is).
github_fetch_mode() {
    local mode="$1" path base cache
    path="$(ex_github_path "$mode")"
    [ -n "$path" ] || return 1
    base="$EX_GITHUB_RAW_BASE/$EX_GITHUB_REPO/$EX_GITHUB_REF/$path"
    cache="$EX_CACHE_DIR/$mode"
    ex_fetch "$base/ngspca_output/svd.pcs.txt" "$cache/svd.pcs.txt"   "SAMPLE" || return 1
    ex_fetch "$base/qc_output/sample_qc.tsv"   "$cache/sample_qc.tsv" "SAMPLE_ID" || return 1
}

count_mosdepth() { find "$1" -maxdepth 1 -name '*.regions.bed.gz' 2>/dev/null | grep -c . || true; }

# The coverage-median table for a run, honouring EX_MEDIAN_SOURCE. Empty
# means "use the QC table's HQ_MEDIAN_COV".
resolve_median() {                 # resolve_median <ngspcaDir>
    local t="$1/autosomal.median.txt"
    case "$EX_MEDIAN_SOURCE" in
        qc)     printf '\n' ;;
        ngspca) [ -s "$t" ] || dsv_die "EX_MEDIAN_SOURCE=ngspca but $t is missing (a run reusing a cached matrix does not write it)"
                printf '%s\n' "$t" ;;
        auto)   if [ -s "$t" ]; then printf '%s\n' "$t"; else printf '\n'; fi ;;
        *)      dsv_die "EX_MEDIAN_SOURCE must be auto, ngspca or qc" ;;
    esac
}

write_paths_env() {                # write_paths_env <mode> <pcs> <qc> <mosdepthDir> <median> <source>
    local envf; envf="$(ex_paths_env "$1")"
    mkdir -p "$(dirname "$envf")"
    {
        printf '# written by 00_fetch_inputs.sh %s\n' "$(dsv_now)"
        printf 'EX_M_PCS_FILE=%q\n'     "$2"
        printf 'EX_M_QC_TABLE=%q\n'     "$3"
        printf 'EX_M_MOSDEPTH_DIR=%q\n' "$4"
        printf 'EX_M_MEDIAN_TABLE=%q\n' "$5"
        printf 'EX_M_SOURCE=%q\n'       "$6"
    } > "$envf"
    dsv_log "$1: pcs+qc from $6; medians from ${5:-the QC table}; mosdepth dir: $4"
}

# --- per-mode resolution ---------------------------------------------------

resolved=0
for mode in $EX_MODES; do
    ex_check_mode "$mode"
    dsv_log "=== mode: $mode ==="

    if [ "$EX_SMOKE" = "1" ]; then
        if [ "$mode" = seedctl ]; then
            # Never fetched: a seed control is a local calibration run or
            # nothing. Resolved when one exists beside the smoke inputs.
            seed_dir="$(ex_ngspca_dir seedctl)"
            if [ -s "$seed_dir/svd.pcs.txt" ] && [ -s "$EX_CACHE_DIR/standard/sample_qc.tsv" ] \
               && [ "$(count_mosdepth "$EX_SMOKE_DIR/standard")" -gt 0 ]; then
                write_paths_env seedctl "$seed_dir/svd.pcs.txt" "$EX_CACHE_DIR/standard/sample_qc.tsv" \
                                "$EX_SMOKE_DIR/standard" "$(resolve_median "$seed_dir")" "smoke:local-seed-control"
                resolved=$((resolved + 1))
            else
                dsv_log "SKIP seedctl: no local seed-control PCs at $seed_dir (set EX_NGSPCA_DIR_SEEDCTL to include one)"
            fi
            continue
        fi

        # Smoke: committed outputs + a simulated mosdepth tree.
        qc_for_sim=""
        source_tag="github"
        if github_fetch_mode "$mode"; then
            qc_for_sim="$EX_CACHE_DIR/$mode/sample_qc.tsv"
            pcs="$EX_CACHE_DIR/$mode/svd.pcs.txt"
            qc="$qc_for_sim"
        elif [ "$mode" = "fast" ] && [ -s "$EX_CACHE_DIR/standard/sample_qc.tsv" ]; then
            # The fast tree is not committed upstream (yet). Reuse the
            # standard tables so the two-mode machinery can still be
            # exercised; the simulated fast depths get a small deterministic
            # perturbation. Clearly synthetic, and labelled as such.
            dsv_log "fast outputs not on GitHub; smoke falls back to the standard tables (synthetic fast mode)"
            qc_for_sim="$EX_CACHE_DIR/standard/sample_qc.tsv"
            pcs="$EX_CACHE_DIR/standard/svd.pcs.txt"
            qc="$qc_for_sim"
            source_tag="github-standard-tables-synthetic-fast"
        else
            dsv_die "could not fetch $mode outputs from GitHub ($EX_GITHUB_REPO@$EX_GITHUB_REF); run modes in the order listed in EX_MODES"
        fi

        smoke_tree="$EX_SMOKE_DIR/$mode"
        have="$(count_mosdepth "$smoke_tree")"
        if [ "$force" -eq 0 ] && [ "$have" -ge "$EX_SMOKE_SAMPLES" ]; then
            dsv_log "simulated tree present ($have samples): $smoke_tree"
        else
            dsv_require_cmd Rscript bgzip
            jitter=0
            [ "$mode" = "fast" ] && jitter=0.01
            dsv_log "simulating $EX_SMOKE_SAMPLES samples into $smoke_tree (jitter=$jitter)"
            Rscript "$EX_EXAMPLE_DIR/R/make_smoke_inputs.R" \
                --qc "$qc_for_sim" --out "$smoke_tree" \
                --samples "$EX_SMOKE_SAMPLES" --seed "$EX_SMOKE_SEED" --jitter "$jitter"
        fi
        write_paths_env "$mode" "$pcs" "$qc" "$smoke_tree" "" "smoke:$source_tag"
        resolved=$((resolved + 1))
        continue
    fi

    # Real run: prefer the local NGS-PCA trees.
    ngspca="$(ex_ngspca_dir "$mode")"
    qcdir="$(ex_qc_dir "$mode")"
    mosdepth_dir="$(ex_mosdepth_dir "$mode")"
    n_md="$(count_mosdepth "$mosdepth_dir")"

    if [ -s "$ngspca/svd.pcs.txt" ] && [ -s "$qcdir/sample_qc.tsv" ]; then
        write_paths_env "$mode" "$ngspca/svd.pcs.txt" "$qcdir/sample_qc.tsv" "$mosdepth_dir" \
                        "$(resolve_median "$ngspca")" "local"
    elif [ "$mode" != seedctl ] && github_fetch_mode "$mode"; then
        dsv_log "$mode: no local NGS-PCA outputs under $NGSPCA_WORK_DIR; using the committed GitHub results"
        write_paths_env "$mode" "$EX_CACHE_DIR/$mode/svd.pcs.txt" "$EX_CACHE_DIR/$mode/sample_qc.tsv" \
                        "$mosdepth_dir" "" "github"
    elif [ "$mode" = seedctl ]; then
        dsv_log "SKIP seedctl: no seed-control PCs at $ngspca/svd.pcs.txt (needs $qcdir/sample_qc.tsv too)."
        dsv_log "       Produce them with NGS-PCA's step 2b, e.g."
        dsv_log "       NGSPCA_OUTPUT=\$WORK_DIR/ngspca_output_seed$EX_SEED_CONTROL_SEED RANDOM_SEED=$EX_SEED_CONTROL_SEED sbatch 02_run_ngspca.sh"
        continue
    else
        dsv_log "SKIP $mode: no local outputs (looked in $ngspca and $qcdir) and none committed on GitHub."
        dsv_log "       Run the NGS-PCA example for this mode, or drop it from EX_MODES."
        continue
    fi

    if [ "$n_md" -eq 0 ]; then
        dsv_log "NOTE $mode: no mosdepth files in $mosdepth_dir yet. The join stage needs them:"
        dsv_log "       run NGS-PCA's 01_download_and_mosdepth.sh (COMPARE_FAST_MODE=1), or set EX_SMOKE=1 for a simulated tree."
    else
        dsv_log "$mode: $n_md mosdepth region files in $mosdepth_dir"
    fi
    resolved=$((resolved + 1))
done

[ "$resolved" -gt 0 ] || dsv_die "no mode could be resolved"
dsv_log "done. Next: bash 01_prepare_inputs.sh"
