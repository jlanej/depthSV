#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# depthSV — build the depth matrix from per-sample mosdepth region files
#
#   scripts/join.sh --manifest FILE --out DIR [--jobs N] [--threads N]
#                   [--batch-size N|auto] [--force]
#
# Output:
#   <out>/depth.matrix.txt.gz      bgzip-compressed, tabix-indexed
#   <out>/depth.matrix.manifest    sample list, row count, coordinate signature
#
# Deliberate choices:
#
#   bgzip, not pigz.  The matrix is indexed so every later stage can read a
#   region instead of scanning the file. That single change removes the
#   repeated whole-genome decompression from the correction stage and makes an
#   arbitrary interval a valid unit of work on any scheduler.
#
#   Row counts AND coordinates are verified, inside each extraction worker.
#   paste() aligns by position and cannot detect that one sample was processed
#   against a different contig set or bin order; the offending sample's depths
#   would simply be attributed to the wrong coordinates. The coordinate
#   checksum costs nothing (same decompression pass), so it is always on, and
#   a mismatch stops the run at the offending sample, not at the end.
#
#   A finished matrix is reused only for the same manifest CONTENT: adding a
#   sample and re-running rebuilds it rather than serving the old columns.
#
#   The work is batched, parallel and resumable.  Columns are extracted
#   --jobs at a time and pasted in batches; each finished batch is compressed
#   and its columns deleted, so scratch stays near two batches of columns plus
#   the compressed batches rather than a full uncompressed copy of the matrix.
#   Batch size defaults to about sqrt(#samples), which also keeps every paste
#   under the open-file limit; the limit is checked up front with an
#   actionable error. A re-run resumes from the last finished batch as long
#   as the manifest and batching are unchanged.
# ---------------------------------------------------------------------------

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../lib/common.sh"
dsv_enable_error_trace

manifest="${DSV_MANIFEST:-}"; out_dir="${DSV_JOIN_DIR:-}"
threads="${DSV_THREADS:-4}"; jobs="${DSV_JOBS:-4}"
batch_size="${DSV_BATCH_SIZE:-auto}"; force=0

while [ $# -gt 0 ]; do
    case "$1" in
        --manifest)      manifest="$2"; shift 2 ;;
        --out)           out_dir="$2"; shift 2 ;;
        --jobs)          jobs="$2"; shift 2 ;;
        --threads)       threads="$2"; shift 2 ;;
        --batch-size)    batch_size="$2"; shift 2 ;;
        --strict-coords) shift ;;                    # now always on; accepted for compatibility
        --force)         force=1; shift ;;
        -h|--help)       dsv_usage ;;
        *)               dsv_die "unknown argument: $1" ;;
    esac
done

dsv_require_opt manifest out_dir
dsv_require_file "$manifest"
# Modules first: on a cluster the tools may only exist after `module load`.
dsv_load_modules
dsv_require_cmd bgzip tabix paste awk parallel mkfifo

mkdir -p "$out_dir"
# Absolute from here on: the batch step cds into the column directory, and a
# relative path captured beforehand would stop resolving after that.
out_dir="$(cd "$out_dir" && pwd)"
final="$out_dir/depth.matrix.txt.gz"
work="$out_dir/.join.work"

# --- plan the batches -------------------------------------------------------

n_samples="$(grep -c . "$manifest" || true)"
[ "$n_samples" -gt 0 ] || dsv_die "manifest contains no paths"

# The signature of this join: the manifest's content (not just its path — a
# grown cohort re-using the old output is exactly the failure this prevents),
# the batching, and the scripts. Checked before any skip.
manifest_sig="$(cksum < "$manifest" | awk '{print $1"-"$2}')"
sig="$(printf 'stage=join\nmanifest=%s\nsamples=%s\nscript=%s\nworker=%s' \
       "$manifest_sig" "$n_samples" "$(dsv_script_sig "$0")" "$(dsv_script_sig "$DSV_ROOT/lib/join_extract.sh")")"

if [ "$force" -eq 0 ] && dsv_output_complete "$final" "$sig"; then
    dsv_log "already complete, skipping: $final"
    exit 0
fi
if [ "$force" -eq 1 ]; then rm -rf "$work"; fi
dsv_output_reset "$final" "$force"

# Duplicate names would become duplicate matrix columns and every later stage
# matches samples by name, so refuse them before any heavy work. The list is
# written to a file rather than piped into `head`, which would kill `uniq`
# with SIGPIPE (exit 141, no message) on a long duplicate list.
mkdir -p "$out_dir"
dup_file="$out_dir/.dsv.dups.$$"
while IFS= read -r f; do
    [ -n "$f" ] || continue
    dsv_sample_name "$f"
done < "$manifest" | sort | uniq -d > "$dup_file"
if [ -s "$dup_file" ]; then
    dup="$(head -3 "$dup_file" | tr '\n' ' ')"
    rm -f "$dup_file"
    dsv_die "duplicate sample name(s) in the manifest after suffix stripping: $dup"
fi
rm -f "$dup_file"

if [ "$batch_size" = "auto" ]; then
    # sqrt(N) balances the two file-count ceilings: each batch paste opens
    # batch_size columns, the final paste opens one pipe per batch.
    batch_size="$(awk -v n="$n_samples" 'BEGIN{b=int(sqrt(n)); if(b*b<n) b++; if(b<8) b=8; if(b>1024) b=1024; print b}')"
fi
n_batches=$(( (n_samples + batch_size - 1) / batch_size ))

fd_need=$(( (batch_size > n_batches ? batch_size : n_batches) + 16 ))
fd_limit="$(ulimit -n)"
if [ "$fd_limit" != "unlimited" ] && [ "$fd_need" -gt "$fd_limit" ]; then
    dsv_die "joining $n_samples samples with batch size $batch_size needs ~$fd_need open files but the limit is $fd_limit. Raise it (ulimit -n $fd_need) or pass a --batch-size nearer sqrt($n_samples)."
fi

dsv_log "joining $n_samples samples into $out_dir ($n_batches batches of <=$batch_size, $jobs extraction jobs)"

# --- resume or fresh start --------------------------------------------------
# A previous partial run can be resumed only if it was cut from the same
# manifest into the same batches; otherwise the leftovers are meaningless.

resume_meta="$work/resume.meta"
wanted_meta="$(printf 'manifest\t%s\nbatch_size\t%s\n' "$manifest_sig" "$batch_size")"

if [ -f "$resume_meta" ] && [ "$(cat "$resume_meta")" = "$wanted_meta" ]; then
    dsv_log "resuming: $(ls "$work" | grep -c '\.done$' || true) of $n_batches batches already finished"
else
    rm -rf "$work"
    mkdir -p "$work/cols" "$work/meta"
    printf '%s' "$wanted_meta" > "$resume_meta"
fi
mkdir -p "$work/cols" "$work/meta"

# Background jobs (a batch paste, the final fan-in) must not linger if the
# script dies between starting and reaping them. The work directory itself is
# kept on failure so the run can resume.
bg_pids=""
cleanup() {
    [ -z "$bg_pids" ] || kill $bg_pids 2>/dev/null || true
    rm -rf "$work/fifos"
}
trap cleanup EXIT

# --- coordinates ------------------------------------------------------------
# Coordinate columns come from the first sample; every worker then verifies
# its own row count against them (and its coordinate checksum, with strict).
# The header is '#'-prefixed so tabix treats it as a comment rather than
# trying to parse it as an interval. Downstream readers address samples
# positionally, so the name of the first column is immaterial to them.

first_file="$(awk 'NF {print; exit}' "$manifest")"
dsv_require_file "$first_file"

coords="$work/coords.tsv"
{ printf '#CHR\tSTART\tSTOP\n'; gzip -cd "$first_file" | cut -f1-3; } > "$coords"
ref_rows=$(( $(wc -l < "$coords") - 1 ))
[ "$ref_rows" -gt 0 ] || dsv_die "first sample has no regions: $first_file"

# --- extract and paste, one batch at a time ---------------------------------
# The paste of batch k runs in the background while batch k+1 extracts, so the
# single-threaded paste overlaps the parallel decompression instead of
# serialising after it.

awk -v bs="$batch_size" -v dir="$work" 'NF {
    f = sprintf("%s/batch.%06d.list", dir, int(i / bs) + 1); i++
    if (f != prev) { if (prev) close(prev); prev = f }
    printf "%d\t%s\n", i, $0 > f
}' "$manifest"

batch_cols() {                     # batch_cols <lo> <hi> — column files, in order
    ls "$work/cols" | awk -F_ -v lo="$1" -v hi="$2" '$1+0 >= lo && $1+0 <= hi' | sort
}

paste_batch() {                    # paste_batch <lo> <hi> <outGz>
    local lo="$1" hi="$2" out="$3" cols
    cols="$(batch_cols "$lo" "$hi")"
    ( cd "$work/cols" && paste $cols ) | bgzip -l 1 -@ "$threads" > "$out.tmp"
    mv "$out.tmp" "$out"
    : > "$out.done"
    ( cd "$work/cols" && rm -f $cols )
}

paste_pid=""
lo=1
ref_sig=""
for list in "$work"/batch.*.list; do
    bn="$(basename "$list" .list)"
    out="$work/$bn.tsv.gz"
    n_in_batch="$(grep -c . "$list")"
    hi=$(( lo + n_in_batch - 1 ))

    if [ -f "$out.done" ]; then
        lo=$(( hi + 1 ))
        continue
    fi

    # -q makes parallel re-quote the command tokens it composes, so paths
    # survive intact; replacement strings are still substituted.
    parallel -q --colsep '\t' -j "$jobs" --halt now,fail=1 \
        bash "$DSV_ROOT/lib/join_extract.sh" '{1}' '{2}' "$work" "$ref_rows" \
        :::: "$list"

    # Every sample's coordinate columns must checksum like the first sample's:
    # the same bins in a different order have the same row count and would
    # paste cleanly into misattributed coordinates.
    [ -n "$ref_sig" ] || ref_sig="$(cat "$work"/meta/000001_*.sig)"
    for sig_file in $(ls "$work/meta" | awk -F_ -v lo="$lo" -v hi="$hi" '/\.sig$/ && $1+0 >= lo && $1+0 <= hi'); do
        this_sig="$(cat "$work/meta/$sig_file")"
        [ "$this_sig" = "$ref_sig" ] || dsv_die "$(basename "$sig_file" .sig | cut -d_ -f2-) has the same row count but different coordinates (signature $this_sig vs $ref_sig)"
    done

    if [ -n "$paste_pid" ]; then
        wait "$paste_pid" || dsv_die "batch paste failed"
    fi
    paste_batch "$lo" "$hi" "$out" &
    paste_pid=$!
    bg_pids="$paste_pid"

    dsv_log "batch ${bn#batch.}/$n_batches extracted (samples $lo-$hi of $n_samples)"
    lo=$(( hi + 1 ))
done
if [ -n "$paste_pid" ]; then
    wait "$paste_pid" || dsv_die "batch paste failed"
fi
bg_pids=""

# --- fan the batches into the final matrix ----------------------------------
# One decompressor per batch feeds a FIFO and a single paste assembles them
# beside the coordinates. FIFO count equals batch count, which the fd check
# above already bounded.

rm -rf "$work/fifos"; mkdir "$work/fifos"
i=0
for bgz in "$work"/batch.*.tsv.gz; do
    i=$(( i + 1 ))
    p="$(printf '%s/fifos/f%06d' "$work" "$i")"
    mkfifo "$p"
    bgzip -dc "$bgz" > "$p" &
    bg_pids="$bg_pids $!"
done

paste "$coords" "$work"/fifos/* | bgzip -@ "$threads" > "$(dsv_output_tmp "$final")"

for pid in $bg_pids; do
    wait "$pid" || dsv_die "batch decompression failed during the final paste"
done
bg_pids=""

# Every row must carry a coordinate triple plus one value per sample. A short
# row here means paste ran out of input for some column. awk stops at the
# first bad row, which kills the decompressor with SIGPIPE — so pipefail is
# lifted for this one assignment, or the check itself would abort the script
# before it could say which row was wrong.
expected_cols=$(( 3 + n_samples ))
set +o pipefail
bad=$(bgzip -dc "$(dsv_output_tmp "$final")" | awk -F'\t' -v e="$expected_cols" 'NF != e {print NR; exit}')
set -o pipefail
[ -z "$bad" ] || dsv_die "row $bad does not have $expected_cols columns"

dsv_output_commit "$final" "$ref_rows" "$sig"

# --- manifest --------------------------------------------------------------

( cd "$work/meta" && ls | grep '\.meta$' | sort | xargs cat ) > "$work/rowcounts.txt"
[ "$(grep -c . "$work/rowcounts.txt")" -eq "$n_samples" ] \
    || dsv_die "internal error: $(grep -c . "$work/rowcounts.txt") extraction records for $n_samples samples"

{
    printf '# depthSV depth matrix\n'
    printf 'created\t%s\n' "$(dsv_now)"
    printf 'samples\t%s\n' "$n_samples"
    printf 'regions\t%s\n' "$ref_rows"
    printf 'coordinate_source\t%s\n' "$(dsv_sample_name "$first_file")"
    printf 'coordinate_signature\t%s\n' "$ref_sig"
    printf '#\n# sample\tregions\tsource\n'
    cat "$work/rowcounts.txt"
} > "$out_dir/depth.matrix.manifest"

rm -rf "$work"
dsv_log "manifest: $out_dir/depth.matrix.manifest"
