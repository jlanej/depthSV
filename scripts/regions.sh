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
dsv_require_cmd tabix awk

[ -s "${matrix}.tbi" ] || [ -s "${matrix}.csi" ] \
    || dsv_die "no tabix index beside $matrix"

if [ "$window" -le 0 ]; then
    tabix -l "$matrix"
    exit 0
fi

[ -n "$sizes" ] || dsv_die "--window needs --sizes (a contig<TAB>length file)"
dsv_require_file "$sizes"

# Only contigs present in both the matrix and the sizes file, so the list can
# never name a unit the matrix cannot serve.
tabix -l "$matrix" > "$(dirname "$matrix")/.dsv.contigs.$$"
trap 'rm -f "$(dirname "$matrix")/.dsv.contigs.$$"' EXIT

awk -v w="$window" -v matrix="$matrix" '
  NR == FNR { present[$1] = 1; next }
  {
    contig = $1; len = $2 + 0
    if (!(contig in present) || len <= 0) next
    for (s = 0; s < len; s += w) {
      e = s + w - 1
      if (e >= len) e = len - 1
      print contig ":" s "-" e
    }
  }
' "$(dirname "$matrix")/.dsv.contigs.$$" "$sizes" \
| while IFS= read -r region; do
    # Skip intervals with no rows: an empty unit is a job that fails rather
    # than a job with nothing to do.
    if [ -n "$(tabix "$matrix" "$region" | head -n 1)" ]; then
        printf '%s\n' "$region"
    fi
done
