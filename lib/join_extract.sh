#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# depthSV — extract one sample's depth column (worker for scripts/join.sh)
#
#   join_extract.sh <index> <file> <workDir> <strict01> <expectedRows>
#
# Writes, under <workDir>:
#   cols/<index>_<sample>       the sample name, then one depth value per region
#   meta/<index>_<sample>.meta  sample<TAB>rows<TAB>source
#   meta/<index>_<sample>.sig   coordinate checksum (only with strict)
#
# The column is written to a temporary name and moved, so a killed worker never
# leaves a file the batch step would mistake for complete. The row count is
# verified here, inside the worker, so one wrong sample stops the whole
# extraction immediately (the driver runs workers under --halt now,fail=1)
# instead of after days of further decompression.
# ---------------------------------------------------------------------------

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/common.sh"

idx="$1"; f="$2"; work="$3"; strict="$4"; expected="$5"

s="$(dsv_sample_name "$f")"

# Column files are pasted by name from a worker shell, so a name that word-
# splits or globs would corrupt the batch. Refuse it with the sample named.
case "$s" in
    *[\ \	\*\?\[]* ) dsv_die "sample name contains characters unsafe for the join: '$s' (from $f)" ;;
esac

[ -s "$f" ] || dsv_die "sample file missing or empty: $f"

tag="$(printf '%06d_%s' "$idx" "$s")"
col="$work/cols/$tag"

{ printf '%s\n' "$s"; gzip -cd "$f" | cut -f4; } > "$col.tmp"
rows=$(( $(wc -l < "$col.tmp") - 1 ))

if [ "$rows" -ne "$expected" ]; then
    rm -f "$col.tmp"
    dsv_die "$s has $rows regions but the first sample has $expected. Inputs were not produced against the same intervals; joining them would attribute depths to the wrong coordinates."
fi

if [ "$strict" -eq 1 ]; then
    gzip -cd "$f" | cut -f1-3 | cksum | awk '{print $1"-"$2}' > "$work/meta/$tag.sig"
fi

printf '%s\t%s\t%s\n' "$s" "$rows" "$f" > "$work/meta/$tag.meta"
mv "$col.tmp" "$col"
