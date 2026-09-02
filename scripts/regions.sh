#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# depthSV — emit the list of work units for an indexed matrix
#
#   scripts/regions.sh --matrix FILE [--window BP --sizes FILE] > regions.txt
#
# One region per line. This is the only thing a scheduler needs: a SLURM array,
# a WDL scatter and a CSV fan-out all iterate the same list, which is what keeps
# them interchangeable. Handing someone the list also makes a run reproducible.
#
# By default one region per contig, taken from the tabix index — so the list
# reflects what the matrix actually contains rather than an assumed karyotype,
# and contig naming is whatever the file uses.
#
# --window splits each contig into fixed-size intervals, which makes each unit
# small enough to be cheap to lose on a preemptible instance. It needs a
# chrom-sizes file (contig<TAB>length) because a tabix index does not record
# contig lengths; empty intervals are skipped so the list never contains a unit
# that would fail as "contains no rows".
#
# The emitted intervals partition the matrix rows: consecutive windows never
# both contain the same bin, so concatenating per-window outputs yields each
# region exactly once.
# ---------------------------------------------------------------------------

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../lib/common.sh"
dsv_enable_error_trace

matrix="${DSV_MATRIX:-}"; window=0; sizes=""

while [ $# -gt 0 ]; do
    case "$1" in
        --matrix)  matrix="$2"; shift 2 ;;
        --window)  window="$2"; shift 2 ;;
        --sizes)   sizes="$2";  shift 2 ;;
        -h|--help) dsv_usage ;;
        *)         dsv_die "unknown argument: $1" ;;
    esac
done

dsv_require_opt matrix
dsv_require_file "$matrix"
dsv_load_modules
dsv_require_cmd tabix awk

[ -s "${matrix}.tbi" ] || [ -s "${matrix}.csi" ] \
    || dsv_die "no tabix index beside $matrix"

if [ "$window" -le 0 ]; then
    tabix -l "$matrix"
    exit 0
fi

[ -n "$sizes" ] || dsv_die "--window needs --sizes (a contig<TAB>length file)"
dsv_require_file "$sizes"

# Window boundaries must land on bin edges: regions are matched by overlap, so
# a boundary falling inside a bin would hand that bin to both neighbouring
# windows. Bins are fixed-size (the premise of the whole pipeline), so read the
# bin size from the first data row and round the window up to a multiple of it.
bin=$( set +o pipefail; gzip -cd "$matrix" | awk -F'\t' 'NR==2 {print $3-$2; exit}' )
if [ -n "$bin" ] && [ "$bin" -gt 0 ] && [ $((window % bin)) -ne 0 ]; then
    window=$(( (window / bin + 1) * bin ))
    dsv_log "window rounded up to $window to align with the ${bin} bp bins"
fi

# Only contigs present in both the matrix and the sizes file, so the list can
# never name a unit the matrix cannot serve. The sizes file may use the other
# naming convention (`1` vs `chr1`); contigs are matched either way and the
# list is written in the matrix's own names. A contig the matrix has but the
# sizes file lacks is reported, because silently dropping it would analyse a
# smaller genome than the one that was joined; no match at all is an error.
contigs="$(mktemp "${TMPDIR:-/tmp}/dsv.contigs.XXXXXX")"
trap 'rm -f "$contigs"' EXIT
tabix -l "$matrix" > "$contigs"

missing="$(awk '
  NR == FNR { s[$1] = 1; sub(/^chr/, "", $1); s[$1] = 1; next }
  { c = $1; b = c; sub(/^chr/, "", b); if (!(c in s) && !(b in s) && !("chr" b in s)) print c }
' "$sizes" "$contigs" | tr '\n' ' ')"
[ -z "$missing" ] || dsv_log "WARNING: $(printf '%s' "$missing" | wc -w | tr -d ' ') contig(s) in the matrix have no entry in $sizes and get no work units: $missing"

# Windows are emitted as 1-based inclusive tabix queries. The matrix rows are
# BED (0-based, half-open), and a tabix query [B,E] matches any feature that
# overlaps it — so a window written as its 0-based start would also pick up
# the bin that ENDS exactly there, duplicating one bin at every boundary.
# Emitting [s+1, s+w] for a 0-based window start s makes consecutive windows
# a partition: a bin ending on a boundary belongs to the earlier window only.
n_units=0
while IFS= read -r region; do
    # Skip intervals with no rows: an empty unit is a job that fails rather
    # than a job with nothing to do.
    if [ -n "$(tabix "$matrix" "$region" | head -n 1)" ]; then
        printf '%s\n' "$region"
        n_units=$(( n_units + 1 ))
    fi
done < <(awk -v w="$window" '
  NR == FNR { name[$1] = $1; b = $1; sub(/^chr/, "", b); name[b] = $1; name["chr" b] = $1; next }
  {
    if (!($1 in name)) next
    contig = name[$1]; len = $2 + 0
    if (len <= 0) next
    for (s = 0; s < len; s += w) {
      e = s + w
      if (e > len) e = len
      print contig ":" (s + 1) "-" e
    }
  }
' "$contigs" "$sizes")
[ "$n_units" -gt 0 ] || dsv_die "no contig of $matrix matched $sizes (contigs: $(tr '\n' ' ' < "$contigs"))"
