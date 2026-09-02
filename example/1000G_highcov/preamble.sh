#!/usr/bin/env bash
#SBATCH --job-name=dsvx-preamble
#SBATCH --output=dsvx-preamble.%j.out
#SBATCH --cpus-per-task=24
#SBATCH --mem=248G
#SBATCH --time=24:00:00
# ---------------------------------------------------------------------------
# depthSV 1000G example — preamble: the number of PCs, and the covariates
#
#   sbatch preamble.sh                      # on the cluster, from this directory
#   bash preamble.sh [--ndim-only | --genotypes-only] [--smoke] [--force]
#
# Two decisions the association models depend on, made once, before the
# pipeline runs, from data that is already on disk or a few gigabytes away:
#
#   A. How many coverage PCs to remove. Marchenko-Pastur on each mode's
#      NGS-PCA spectrum (svd.singularvalues.txt), averaged across modes so
#      both trees are corrected identically -> preamble/ndim.txt, which
#      config.sh picks up as the EX_NDIM default. See R/choose_ndim.R.
#
#   B. Genotype PCs as covariates. The NYGC 30x GRCh38 callset for exactly
#      these 3,202 samples is published in PLINK 2 format (resources/
#      genotype_sources.tsv, ~4 GB in total, no CRAMs or VCFs), so this is
#      the textbook recipe end to end: biallelic autosomal SNPs, MAF and
#      missingness filters, long-range LD regions excluded, LD pruning,
#      KING-robust relatedness, PCA on the unrelated set, projection of the
#      relatives, calibration, plots -> preamble/covariates.tsv, which the
#      prepare stage merges into the phenotype table. See R/genotype_pcs.R.
#
# Both parts are idempotent: downloads resume, finished chromosomes are
# skipped, and --force redoes the genotype pipeline from the pruning step.
# Part A needs only R; part B needs plink2 (>= 2.00a5) and curl.
#
# Sized for a 24-core / 248 GB / 24 h allocation to match the upstream
# NGS-PCA stages; the actual need is far smaller (the PCA is seconds at
# 3,202 samples; the wall time is downloading and filtering ~70M variants).
# ---------------------------------------------------------------------------

do_ndim=1; do_geno=1; smoke=0; force=0; want_help=0
while [ $# -gt 0 ]; do
    case "$1" in
        --ndim-only)      do_geno=0; shift ;;
        --genotypes-only) do_ndim=0; shift ;;
        --smoke)          smoke=1; shift ;;
        --force)          force=1; shift ;;
        -h|--help)        want_help=1; shift ;;
        *)                echo "unknown argument: $1 (see --help)" >&2; exit 2 ;;
    esac
done
[ "$smoke" -eq 0 ] || export EX_SMOKE=1

EX_EXAMPLE_DIR="${EX_EXAMPLE_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}}"
source "$EX_EXAMPLE_DIR/lib.sh"
dsv_enable_error_trace
[ "$want_help" -eq 0 ] || dsv_usage

mkdir -p "$EX_PREAMBLE_DIR"
dsv_log "preamble -> $EX_PREAMBLE_DIR (smoke=$EX_SMOKE)"

# Download with resume and retries; refuse an empty or HTML result.
fetch_url() {                      # fetch_url <url> <dest>
    local url="$1" dest="$2"
    if [ -s "$dest" ]; then dsv_log "present: $(basename "$dest")"; return 0; fi
    mkdir -p "$(dirname "$dest")"
    curl -fsSL --retry 5 --retry-delay 5 -C - -o "$dest.part" "$url" \
        || dsv_die "download failed: $url"
    if head -c 512 "$dest.part" | grep -qi '<html'; then
        rm -f "$dest.part"; dsv_die "got an HTML page instead of a file from $url"
    fi
    mv "$dest.part" "$dest"
    dsv_log "fetched: $(basename "$dest") ($(du -h "$dest" | cut -f1))"
}

# =============================================================================
# A. ndim from the Marchenko-Pastur edge
# =============================================================================
if [ "$do_ndim" -eq 1 ]; then
    dsv_require_cmd Rscript
    runs=""
    for mode in standard fast; do
        dir="$(ex_ngspca_dir "$mode")"
        if [ ! -s "$dir/svd.singularvalues.txt" ] || [ "$EX_SMOKE" = "1" ]; then
            # Not on disk (or a smoke run): the committed GitHub results carry
            # the spectrum too. The seed control shares the standard spectrum
            # up to its reseed and is not consulted.
            path="$(ex_github_path "$mode")"
            cache="$EX_CACHE_DIR/$mode"
            if [ -n "$path" ]; then
                base="$EX_GITHUB_RAW_BASE/$EX_GITHUB_REPO/$EX_GITHUB_REF/$path/ngspca_output"
                ok=1
                for f in svd.singularvalues.txt svd.samples.txt svd.bins.txt; do
                    [ -s "$cache/$f" ] && continue
                    mkdir -p "$cache"
                    curl -fsSL --retry 3 -o "$cache/$f.tmp" "$base/$f" 2>/dev/null \
                        && mv "$cache/$f.tmp" "$cache/$f" || { rm -f "$cache/$f.tmp"; ok=0; break; }
                done
                [ "$ok" -eq 1 ] && dir="$cache"
            fi
        fi
        if [ -s "$dir/svd.singularvalues.txt" ] && [ -s "$dir/svd.samples.txt" ] && [ -s "$dir/svd.bins.txt" ]; then
            runs="${runs:+$runs,}$mode=$dir"
        else
            dsv_log "ndim: no spectrum for $mode (looked in $dir)"
        fi
    done
    [ -n "$runs" ] || dsv_die "no NGS-PCA spectrum found for any mode"
    Rscript "$EX_EXAMPLE_DIR/R/choose_ndim.R" --runs "$runs" \
        --margin "$EX_MP_MARGIN" --gap "$EX_MP_GAP" --out "$EX_PREAMBLE_DIR/ndim"
    if [ -s "$EX_PREAMBLE_DIR/ndim/ndim.txt" ]; then
        cp "$EX_PREAMBLE_DIR/ndim/ndim.txt" "$EX_PREAMBLE_DIR/ndim.txt"
        dsv_log "ndim = $(cat "$EX_PREAMBLE_DIR/ndim.txt") -> $EX_PREAMBLE_DIR/ndim.txt (EX_NDIM overrides it)"
    fi
fi

# =============================================================================
# B. genotype PCs
# =============================================================================
if [ "$do_geno" -eq 1 ]; then
    if command -v module >/dev/null 2>&1; then
        for m in ${EX_PREAMBLE_MODULES:-}; do module load "$m" 2>/dev/null || true; done
    fi
    dsv_require_cmd plink2 curl Rscript awk
    dsv_require_file "$EX_GENO_SOURCES" "$EX_GENO_LD_REGIONS"
    G="$EX_PREAMBLE_DIR/genotypes"
    mkdir -p "$G"
    threads="$EX_GENO_THREADS"; mem="$EX_GENO_MEMORY_MB"
    dsv_log "genotypes: chromosomes [$EX_GENO_CHROMS], MAF>=$EX_GENO_MAF, missing<=$EX_GENO_GENO, LD $EX_GENO_LD_WINDOW r2<$EX_GENO_LD_R2, KING cutoff $EX_GENO_KING_CUTOFF, $EX_GENO_NPC PCs; plink2 $threads threads / ${mem} MB"

    source_url() {                 # source_url <key> <kind>
        awk -F'\t' -v k="$1" -v t="$2" '!/^#/ && $1 == k && $2 == t {print $3; exit}' "$EX_GENO_SOURCES"
    }

    psam="$G/samples.psam"
    fetch_url "$(source_url psam psam)" "$psam"
    [ "$(grep -vc '^#' "$psam")" -ge 3000 ] || dsv_die "$psam has fewer than 3,000 samples"

    # IDs are set from position and alleles so both pvar flavours upstream
    # (rsID-annotated or not) yield one ID space; only SNPs survive, so no
    # allele is long enough to hit the length cap.
    id_args=(--set-all-var-ids '@:#:$r:$a' --new-id-max-allele-len 10 missing --rm-dup exclude-all)

    merge_list="$G/merge.list"; : > "$merge_list"
    for c in $EX_GENO_CHROMS; do
        key="chr$c"; out="$G/$key"
        if [ "$force" -eq 0 ] && [ -s "$out.pruned.pgen" ]; then
            dsv_log "$key: pruned set present"
            printf '%s.pruned\n' "$out" >> "$merge_list"
            continue
        fi
        pgen_url="$(source_url "$key" pgen)"; pvar_url="$(source_url "$key" pvar)"
        [ -n "$pgen_url" ] && [ -n "$pvar_url" ] || dsv_die "no download source for $key in $EX_GENO_SOURCES"
        fetch_url "$pgen_url" "$out.pgen.zst"
        fetch_url "$pvar_url" "$out.pvar.zst"
        [ -s "$out.pgen" ] || plink2 --zst-decompress "$out.pgen.zst" "$out.pgen" >/dev/null
        dsv_log "$key: QC + LD pruning"
        ex_timed preamble geno-prune "$key" -- \
            plink2 --pgen "$out.pgen" --pvar "$out.pvar.zst" --psam "$psam" \
                   --threads "$threads" --memory "$mem" --silent \
                   --snps-only just-acgt --max-alleles 2 "${id_args[@]}" \
                   --maf "$EX_GENO_MAF" --geno "$EX_GENO_GENO" \
                   --exclude range "$EX_GENO_LD_REGIONS" \
                   --indep-pairwise "$EX_GENO_LD_WINDOW" 1 "$EX_GENO_LD_R2" \
                   --out "$out.prune"
        ex_timed preamble geno-extract "$key" -- \
            plink2 --pgen "$out.pgen" --pvar "$out.pvar.zst" --psam "$psam" \
                   --threads "$threads" --memory "$mem" --silent "${id_args[@]}" \
                   --extract "$out.prune.prune.in" --make-pgen --out "$out.pruned"
        rm -f "$out.pgen"      # the decompressed copy; the .zst download is kept
        dsv_log "$key: $(grep -vc '^#' "$out.pruned.pvar") pruned SNPs"
        printf '%s.pruned\n' "$out" >> "$merge_list"
    done

    dsv_log "merging $(grep -c . "$merge_list") chromosome sets"
    if [ "$(grep -c . "$merge_list")" -eq 1 ]; then
        plink2 --pfile "$(cat "$merge_list")" --threads "$threads" --memory "$mem" --silent \
               --make-pgen --out "$G/pruned"
    else
        ex_timed preamble geno-merge all -- \
            plink2 --pmerge-list "$merge_list" pfile --threads "$threads" --memory "$mem" --silent \
                   --make-pgen --out "$G/pruned"
    fi
    nsnp="$(grep -vc '^#' "$G/pruned.pvar")"
    dsv_log "pruned SNP set: $nsnp variants"

    dsv_log "KING-robust relatedness (cutoff $EX_GENO_KING_CUTOFF)"
    ex_timed preamble geno-king all -- \
        plink2 --pfile "$G/pruned" --threads "$threads" --memory "$mem" --silent \
               --king-cutoff "$EX_GENO_KING_CUTOFF" --out "$G/king"
    n_unrel="$(grep -vc '^#' "$G/king.king.cutoff.in.id")"
    dsv_log "unrelated set: $n_unrel samples ($(grep -vc '^#' "$G/king.king.cutoff.out.id") removed)"

    dsv_log "PCA on the unrelated set ($EX_GENO_NPC PCs, allele weights)"
    ex_timed preamble geno-pca all -- \
        plink2 --pfile "$G/pruned" --keep "$G/king.king.cutoff.in.id" \
               --threads "$threads" --memory "$mem" --silent \
               --freq --pca "$EX_GENO_NPC" allele-wts vcols=chrom,pos,ref,alt --out "$G/pca_unrel"

    # .eigenvec.allele with vcols=chrom,pos,ref,alt: #CHROM POS ID REF ALT A1 PC1..
    dsv_log "projecting every sample onto the allele weights"
    ex_timed preamble geno-project all -- \
        plink2 --pfile "$G/pruned" --read-freq "$G/pca_unrel.afreq" \
               --threads "$threads" --memory "$mem" --silent \
               --score "$G/pca_unrel.eigenvec.allele" 3 6 header-read no-mean-imputation variance-standardize \
               --score-col-nums "7-$((6 + EX_GENO_NPC))" --out "$G/pca_proj"

    Rscript "$EX_EXAMPLE_DIR/R/genotype_pcs.R" \
        --eigenvec "$G/pca_unrel.eigenvec" --eigenval "$G/pca_unrel.eigenval" \
        --sscore "$G/pca_proj.sscore" --psam "$psam" --unrelated "$G/king.king.cutoff.in.id" \
        --nsnp "$nsnp" --margin "$EX_MP_MARGIN" --out "$EX_PREAMBLE_DIR"
fi

# =============================================================================
# summary
# =============================================================================
{
    echo "# depthSV 1000G example — preamble"
    echo
    echo "- written: $(dsv_now)"
    if [ -s "$EX_PREAMBLE_DIR/ndim.txt" ]; then
        echo "- coverage PCs to remove (MP, averaged over modes): **$(cat "$EX_PREAMBLE_DIR/ndim.txt")** — see ndim/summary.md"
    else
        echo "- coverage ndim: not determined (see ndim/summary.md if present); EX_NDIM default applies"
    fi
    if [ -s "$EX_PREAMBLE_DIR/covariates.tsv" ]; then
        echo "- genotype PCs: covariates.tsv for $(( $(grep -c . "$EX_PREAMBLE_DIR/covariates.tsv") - 1 )) samples; the models use SEX + GPC1..GPC${EX_N_GPCS} (EX_N_GPCS / EX_COVARIATES)"
        echo
        cat "$EX_PREAMBLE_DIR/summary.md" 2>/dev/null | sed -n '2,$p' || true
    else
        echo "- genotype PCs: not produced; the models run unadjusted"
    fi
} > "$EX_PREAMBLE_DIR/preamble_summary.md"
dsv_log "done: $EX_PREAMBLE_DIR/preamble_summary.md"
