#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# depthSV — extract one sample's depth column (worker for scripts/join.sh)
#
#   join_extract.sh <index> <file> <workDir> <expectedRows>
#
# Writes, under <workDir>:
#   cols/<index>_<sample>       the sample name, then one depth value per region
#   meta/<index>_<sample>.meta  sample<TAB>rows<TAB>source
#   meta/<index>_<sample>.sig   checksum of the coordinate columns
#
# The column is written to a temporary name and moved, so a killed worker never
# leaves a file the batch step would mistake for complete. The row count is
# verified here, inside the worker, so one wrong sample stops the whole
# extraction immediately (the driver runs workers under --halt now,fail=1)
# instead of after days of further decompression.
#
# The coordinate checksum is computed in the same pass as the column, so it
# costs nothing extra and is always on: two inputs with the same bin set in a
# different order have the same row count and would otherwise paste cleanly
# into misattributed coordinates.
# ---------------------------------------------------------------------------

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/common.sh"

idx="$1"; f="$2"; work="$3"; expected="$4"

s="$(dsv_sample_name "$f")"

# Column files are pasted by name from a worker shell, so a name that word-
# splits or globs would corrupt the batch. Refuse it with the sample named.
case "$s" in
    *[\ \	\*\?\[]* ) dsv_die "sample name contains characters unsafe for the join: '$s' (from $f)" ;;
esac

[ -s "$f" ] || dsv_die "sample file missing or empty: $f"

tag="$(printf '%06d_%s' "$idx" "$s")"
col="$work/cols/$tag"
sig="$work/meta/$tag.sig"

# One decompression feeds both the depth column and the coordinate checksum:
# tee copies the stream into a FIFO read by a background checksum job that is
# waited for explicitly (bash 3.2 cannot wait on a process substitution).
# The first depth value is checked to be numeric so a five-column BED (a
# named-region mosdepth run) fails here, in the first second, rather than at
# the correction stage after the whole join.
mkdir -p "$work/fifos"
fifo="$work/fifos/$tag"
rm -f "$fifo"; mkfifo "$fifo"
( cut -f1-3 < "$fifo" | cksum | awk '{print $1"-"$2}' > "$sig.tmp" ) &
sig_pid=$!
{ printf '%s\n' "$s"; gzip -cd "$f" | tee "$fifo" | cut -f4; } > "$col.tmp"
wait "$sig_pid" || { rm -f "$fifo" "$col.tmp" "$sig.tmp"; dsv_die "$s: coordinate checksum failed"; }
rm -f "$fifo"
rows=$(( $(wc -l < "$col.tmp") - 1 ))

first="$(sed -n '2p' "$col.tmp")"
case "$first" in
    ''|*[!0-9.eE+-]*) rm -f "$col.tmp" "$sig.tmp"
        dsv_die "$s: column 4 of $f is not numeric ('$first'); depthSV expects mosdepth --by <bin size> output" ;;
esac

if [ "$rows" -ne "$expected" ]; then
    rm -f "$col.tmp" "$sig.tmp"
    dsv_die "$s has $rows regions but the first sample has $expected. Inputs were not produced against the same intervals; joining them would attribute depths to the wrong coordinates."
fi

mv "$sig.tmp" "$sig"
printf '%s\t%s\t%s\n' "$s" "$rows" "$f" > "$work/meta/$tag.meta"
mv "$col.tmp" "$col"
