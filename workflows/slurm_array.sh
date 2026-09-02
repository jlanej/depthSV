#!/usr/bin/env bash
#SBATCH --job-name=depthsv
#SBATCH --output=depthsv.%A_%a.out
#SBATCH --error=depthsv.%A_%a.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=8:00:00
# ---------------------------------------------------------------------------
# depthSV — SLURM array dispatcher
#
#   scripts/regions.sh --matrix "$DSV_MATRIX" > regions.txt
#   sbatch --array=1-$(wc -l < regions.txt) workflows/slurm_array.sh regions.txt
#
# One array task per region, correcting then analysing it. This is the same
# unit of work the WDL scatters over and the same scripts it calls, so the two
# dispatchers are interchangeable and neither is the "real" one.
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

echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] task $idx -> region $region"

"$here/../scripts/correct.sh" \
    --region "$region" \
    --jobs "$cpus" --threads "$threads"

"$here/../scripts/analyze.sh" \
    --corrected "${DSV_CORRECTED_DIR}/corrected_ndim${DSV_NDIM:-16}.$(printf '%s' "$region" | tr ':' '_').txt.gz" \
    --region "$region" \
    --jobs "$cpus" --threads "$threads"

echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] task $idx done"
