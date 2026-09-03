#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# depthSV — paste one region out of block matrices into a window matrix
#
#   scripts/join_paste.sh --blocks FILE --region chr1:1-10000000 --out DIR \
#                         [--threads N] [--force]
#
# The cloud shape of the join. scripts/join.sh run on a block of ~1,000
# samples writes an indexed matrix of every bin x that block; this script
# reads one region from each block through its index and pastes the slices
# side by side, so the whole-cohort matrix for a window exists only for the
# task that corrects and analyses it and the multi-terabyte whole-genome
# matrix is never materialised. --blocks lists the block matrices, one path
# per line, in the column order wanted; each must be tabix-indexed. The
# coordinates of every block must agree row for row, and a sample name may
# appear in one block only.
#
# Output: DIR/depth.matrix.txt.gz (+ index, .done, .params) in exactly the
# format scripts/join.sh writes, so correct.sh reads it unchanged.
# ---------------------------------------------------------------------------

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../lib/common.sh"
dsv_enable_error_trace

blocks=""; region=""; out_dir=""; threads="${DSV_THREADS:-2}"; force=0
while [ $# -gt 0 ]; do
    case "$1" in
        --blocks)  blocks="$2";  shift 2 ;;
        --region)  region="$2";  shift 2 ;;
        --out)     out_dir="$2"; shift 2 ;;
        --threads) threads="$2"; shift 2 ;;
        --force)   force=1;      shift ;;
        -h|--help) dsv_usage ;;
        *)         dsv_die "unknown argument: $1" ;;
    esac
done
dsv_require_opt blocks region out_dir
dsv_load_modules
dsv_require_cmd bgzip tabix paste cut cksum
dsv_require_file "$blocks"

blk=()
while IFS= read -r line; do
    [ -n "$line" ] || continue
    blk+=("$line")
done < "$blocks"
[ ${#blk[@]} -gt 0 ] || dsv_die "no block matrices listed in $blocks"
for b in "${blk[@]}"; do
    dsv_require_file "$b"
    [ -s "$b.tbi" ] || [ -s "$b.csi" ] || dsv_die "block matrix is not indexed: $b"
done

mkdir -p "$out_dir"
out_dir="$(cd "$out_dir" && pwd)"
final="$out_dir/depth.matrix.txt.gz"

sig="$(printf 'stage=join_paste\nregion=%s\nblocks=%s\nscript=%s' "$region" \
       "$(for b in "${blk[@]}"; do dsv_file_sig "$b"; echo; done | cksum | awk '{print $1}')" \
       "$(dsv_script_sig "$0")")"
if [ "$force" -eq 0 ] && dsv_output_complete "$final" "$sig"; then
    dsv_log "already complete, skipping: $final"
    exit 0
fi
dsv_output_reset "$final" "$force"

work="$(mktemp -d "${TMPDIR:-/tmp}/dsv.paste.XXXXXX")"
trap 'rm -rf "$work"' EXIT

# Header: the first block's three coordinate columns and every block's samples.
{
    printf '%s' "$(dsv_header "${blk[0]}")"
    i=0
    for b in "${blk[@]}"; do
        i=$((i + 1)); [ "$i" -gt 1 ] || continue
        printf '\t%s' "$(dsv_header "$b" | cut -f4-)"
    done
    printf '\n'
} > "$work/header"
dups="$(tr '\t' '\n' < "$work/header" | tail -n +4 | sort | uniq -d | head -n 3)"
[ -z "$dups" ] || dsv_die "sample name(s) appear in more than one block: $(printf '%s ' $dups)"

# Slices: the first block keeps its coordinates, the others contribute their
# sample columns; every block must agree on the bins.
i=0; n_rows=""; ref_coords=""
for b in "${blk[@]}"; do
    i=$((i + 1))
    dsv_read_region "$b" "$region" > "$work/full.$i"
    rows="$(grep -c . "$work/full.$i" || true)"
    coords="$(cut -f1-3 "$work/full.$i" | cksum | awk '{print $1}')"
    if [ "$i" -eq 1 ]; then
        n_rows="$rows"; ref_coords="$coords"
        [ "$n_rows" -gt 0 ] || dsv_die "region $region has no rows in ${blk[0]}"
        mv "$work/full.$i" "$work/body.$i"
    else
        [ "$rows" -eq "$n_rows" ] || dsv_die "block $b has $rows rows for $region, the first block $n_rows"
        [ "$coords" = "$ref_coords" ] || dsv_die "block $b: bin coordinates differ from the first block over $region"
        cut -f4- "$work/full.$i" > "$work/body.$i"
        rm -f "$work/full.$i"
    fi
done

# paste in groups, so hundreds of blocks stay inside the open-file limit.
group=100
files=()
for ((j = 1; j <= i; j++)); do files+=("$work/body.$j"); done
level=0
while [ ${#files[@]} -gt "$group" ]; do
    level=$((level + 1)); next=(); k=0
    while [ "$k" -lt ${#files[@]} ]; do
        chunk=("${files[@]:$k:$group}")
        paste "${chunk[@]}" > "$work/g.$level.$k"
        next+=("$work/g.$level.$k")
        k=$((k + group))
    done
    files=("${next[@]}")
done
{
    cat "$work/header"
    paste "${files[@]}"
} | bgzip -@ "$threads" > "$(dsv_output_tmp "$final")"

dsv_output_commit "$final" "$n_rows" "$sig"
printf 'region\t%s\nblocks\t%d\nsamples\t%d\nregions\t%d\n' "$region" "${#blk[@]}" \
    "$(( $(tr '\t' '\n' < "$work/header" | grep -c .) - 3 ))" "$n_rows" > "$out_dir/depth.matrix.manifest"
