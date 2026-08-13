#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# depthSV — build the depth matrix from per-sample mosdepth region files
#
#   scripts/join.sh --manifest FILE --out DIR [--threads N] [--batch-size N]
#
# Output:
#   <out>/depth.matrix.txt.gz      bgzip-compressed, tabix-indexed
#   <out>/depth.matrix.manifest    sample list, row count, coordinate signature
#
# Two deliberate choices:
#
#   bgzip, not pigz.  The matrix is indexed so every later stage can read a
#   region instead of scanning the file. That single change removes the
#   repeated whole-genome decompression from the correction stage and makes an
#   arbitrary interval a valid unit of work on any scheduler.
#
#   Row counts are verified.  paste() aligns by position and cannot detect that
#   one sample was processed against a different contig set; the offending
#   sample's depths would simply be attributed to the wrong coordinates. The
#   count is captured from the decompression already being performed, so it
#   costs nothing measurable.
# ---------------------------------------------------------------------------

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../lib/common.sh"
dsv_enable_error_trace

manifest="${DSV_MANIFEST:-}"; out_dir="${DSV_JOIN_DIR:-}"
threads="${DSV_THREADS:-4}"; batch_size="${DSV_BATCH_SIZE:-200}"; strict_coords=0

while [ $# -gt 0 ]; do
    case "$1" in
        --manifest)      manifest="$2"; shift 2 ;;
        --out)           out_dir="$2"; shift 2 ;;
        --threads)       threads="$2"; shift 2 ;;
        --batch-size)    batch_size="$2"; shift 2 ;;
        --strict-coords) strict_coords=1; shift ;;   # full coordinate checksum
        -h|--help)       dsv_usage ;;
        *)               dsv_die "unknown argument: $1" ;;
    esac
done

dsv_require_opt manifest out_dir
dsv_require_file "$manifest"
dsv_require_cmd bgzip tabix paste awk

dsv_load_modules

mkdir -p "$out_dir"
# Absolute from here on: join.sh cds into a working subdirectory, and a
# relative path captured beforehand would stop resolving after that.
out_dir="$(cd "$out_dir" && pwd)"
work="$out_dir/.join.work"
rm -rf "$work"; mkdir -p "$work/cols"

n_samples=$(grep -c . "$manifest")
dsv_log "joining $n_samples samples into $out_dir"

# --- extract one depth column per sample, capturing the row count ----------
# Sample name is the file's basename with mosdepth's suffixes removed.
sample_name() {
    local b="${1##*/}"
    b="${b%.gz}"; b="${b%.bed}"; b="${b%.regions}"
    b="${b%.by1000}"; b="${b%.src}"
    printf '%s\n' "$b"
}

first_file="$(head -n1 "$manifest")"
dsv_require_file "$first_file"

# Coordinate columns come from the first sample; every other sample is then
# checked against them.
# The header is '#'-prefixed so tabix treats it as a comment rather than trying
# to parse it as an interval. Downstream readers address samples positionally,
# so the name of the first column is immaterial to them.
{ printf '#CHR\tSTART\tSTOP\n'; gzip -cd "$first_file" | cut -f1-3; } > "$work/cols/000000_coords"
ref_rows=$(( $(wc -l < "$work/cols/000000_coords") - 1 ))
if [ "$strict_coords" -eq 1 ]; then
    ref_sig=$(gzip -cd "$first_file" | cut -f1-3 | cksum | awk '{print $1"-"$2}')
fi

i=0
: > "$work/rowcounts.txt"
while IFS= read -r f; do
    [ -n "$f" ] || continue
    i=$((i+1))
    s="$(sample_name "$f")"
    [ -s "$f" ] || dsv_die "sample file missing or empty: $f"
    col="$(printf '%s/cols/%06d_%s' "$work" "$i" "$s")"

    # The row count is free here: decompression is the bottleneck and wc runs
    # concurrently in the pipeline. Doing the same count in awk costs ~2.3x.
    { printf '%s\n' "$s"; gzip -cd "$f" | cut -f4; } > "$col"
    rows=$(( $(wc -l < "$col") - 1 ))
    printf '%s\t%s\t%s\n' "$s" "$rows" "$f" >> "$work/rowcounts.txt"

    if [ "$rows" -ne "$ref_rows" ]; then
        dsv_die "$s has $rows regions but the first sample has $ref_rows. Inputs were not produced against the same intervals; joining them would attribute depths to the wrong coordinates."
    fi
    if [ "$strict_coords" -eq 1 ]; then
        sig=$(gzip -cd "$f" | cut -f1-3 | cksum | awk '{print $1"-"$2}')
        [ "$sig" = "$ref_sig" ] || dsv_die "$s has the same row count but different coordinates (signature $sig vs $ref_sig)"
    fi
done < "$manifest"

dsv_log "all $i samples agree on $ref_rows regions"

# --- paste in batches, then paste the batches ------------------------------
# paste has an open-file limit; batching keeps it well clear of it.
ls "$work/cols" | sort > "$work/col.order"
split -l "$batch_size" "$work/col.order" "$work/batch."

for b in "$work"/batch.*; do
    ( cd "$work/cols" && paste $(awk '{printf "%s ", $0}' "$b") ) > "${b}.tsv"
done

final="$out_dir/depth.matrix.txt.gz"
dsv_output_reset "$final"
paste "$work"/batch.*.tsv | bgzip -@ "$threads" > "$(dsv_output_tmp "$final")"

# Every row must carry a coordinate triple plus one value per sample. A short
# row here means paste ran out of input for some column.
expected_cols=$(( 3 + i ))
bad=$(bgzip -dc "$(dsv_output_tmp "$final")" | awk -F'\t' -v e="$expected_cols" 'NF != e {print NR; exit}')
[ -z "$bad" ] || dsv_die "row $bad does not have $expected_cols columns"

dsv_output_commit "$final" "$ref_rows"

# --- manifest --------------------------------------------------------------
{
    printf '# depthSV depth matrix\n'
    printf 'created\t%s\n' "$(dsv_now)"
    printf 'samples\t%s\n' "$i"
    printf 'regions\t%s\n' "$ref_rows"
    printf 'coordinate_source\t%s\n' "$(sample_name "$first_file")"
    if [ "$strict_coords" -eq 1 ]; then printf 'coordinate_signature\t%s\n' "$ref_sig"; fi
    printf '#\n# sample\tregions\tsource\n'
    cat "$work/rowcounts.txt"
} > "$out_dir/depth.matrix.manifest"

rm -rf "$work"
dsv_log "manifest: $out_dir/depth.matrix.manifest"
