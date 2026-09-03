#!/usr/bin/env bash
#SBATCH --job-name=depthsv
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=8:00:00
# ---------------------------------------------------------------------------
# depthSV — SLURM array dispatcher
#
#   scripts/regions.sh --matrix "$DSV_MATRIX" > regions.txt
#   mkdir -p logs
#   sbatch --array=1-$(wc -l < regions.txt)%100 --export=ALL \
#          --output=logs/depthsv.%A_%a.out workflows/slurm_array.sh regions.txt
#
# One array task per region, correcting then analysing it. This is the same
# unit of work the WDL scatters over and the same scripts it calls; the
# dispatchers agree on the unit of work, which is what makes results from
# either comparable.
#
# %100 throttles concurrency; MaxArraySize (often 1,000) caps N — split a
# longer list into several arrays. No #SBATCH --output here: that would
# write logs into the submit directory, which may be a read-only checkout.
# --export=ALL carries the DSV_* environment through sites whose default is
# SBATCH_EXPORT=NONE.
#
# Completed units are skipped, so a partially failed array can be resubmitted
# whole and only redoes what is missing.
#
# Configuration comes from the environment; source your site's env file first.
# ---------------------------------------------------------------------------
set -euo pipefail

region_list="${1:?usage: sbatch --array=1-N slurm_array.sh <regions.txt>}"

# sbatch runs a spooled COPY of this script, so its own directory is not the
# checkout. Resolve the stage scripts through DSV_ROOT (exported by any
# sourced stage script, or set it yourself), then the submit directory.
DSV_ROOT="${DSV_ROOT:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}}"
[ -x "$DSV_ROOT/scripts/correct.sh" ] \
    || { echo "cannot find scripts/correct.sh under DSV_ROOT=$DSV_ROOT; export DSV_ROOT to the depthSV checkout" >&2; exit 2; }
here="$DSV_ROOT/workflows"

# Outside SLURM (a login-node smoke check), do the first region.
idx="${SLURM_ARRAY_TASK_ID:-1}"
region="$(sed -n "${idx}p" "$region_list")"
[ -n "$region" ] || { echo "no region on line $idx of $region_list" >&2; exit 2; }

: "${DSV_MATRIX:?set DSV_MATRIX (source your conf/*.env)}"
: "${DSV_CORRECTED_DIR:?set DSV_CORRECTED_DIR}"
: "${DSV_RESULTS_DIR:?set DSV_RESULTS_DIR}"

cpus="${SLURM_CPUS_PER_TASK:-4}"
# Workers and compression threads share the allocation; asking for both at
# the full CPU count oversubscribes the node twofold.
threads=$(( cpus / 2 )); [ "$threads" -ge 1 ] || threads=1
# GNU parallel's --pipe buffers and the R workers' temporaries go under
# TMPDIR; on a node with a small /tmp point it at the job's scratch.
export TMPDIR="${TMPDIR:-${SLURM_TMPDIR:-/tmp}}"

echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] task $idx -> region $region"

"$here/../scripts/correct.sh" \
    --region "$region" \
    --jobs "$cpus" --threads "$threads"

"$here/../scripts/analyze.sh" \
    --corrected "${DSV_CORRECTED_DIR}/corrected_ndim${DSV_NDIM:-16}.$(printf '%s' "$region" | tr ':' '_' | tr -d ' ').txt.gz" \
    --region "$region" \
    --jobs "$cpus" --threads "$threads"

echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] task $idx done"
