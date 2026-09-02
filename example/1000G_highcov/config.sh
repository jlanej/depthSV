# ---------------------------------------------------------------------------
# depthSV 1000G example — shared configuration
#
# Sourced by every stage script (via lib.sh). Override any variable by
# exporting it before running a script:
#
#   export EX_WORK_DIR=/scratch/$USER/depthsv_1000G_highcov
#   export NGSPCA_WORK_DIR=/scratch/$USER/1000G_highcov
#   bash 00_fetch_inputs.sh
#
# Every variable uses ${VAR:-default}, so exported values always win.
# ---------------------------------------------------------------------------

# --- upstream: the NGS-PCA 1000G_highcov example ---------------------------
#
# This example consumes the output of github.com/jlanej/NGS-PCA's
# example/1000G_highcov pipeline, run with COMPARE_FAST_MODE=1 so both
# mosdepth modes exist. NGSPCA_WORK_DIR is that pipeline's WORK_DIR; the
# per-mode trees below default to its layout and can be overridden
# individually.
NGSPCA_WORK_DIR="${NGSPCA_WORK_DIR:-/scratch/${USER}/1000G_highcov}"

EX_MOSDEPTH_DIR_STANDARD="${EX_MOSDEPTH_DIR_STANDARD:-${NGSPCA_WORK_DIR}/mosdepth_output}"
EX_NGSPCA_DIR_STANDARD="${EX_NGSPCA_DIR_STANDARD:-${NGSPCA_WORK_DIR}/ngspca_output}"
EX_QC_DIR_STANDARD="${EX_QC_DIR_STANDARD:-${NGSPCA_WORK_DIR}/qc_output}"

EX_MOSDEPTH_DIR_FAST="${EX_MOSDEPTH_DIR_FAST:-${NGSPCA_WORK_DIR}/mosdepth_output_fast}"
EX_NGSPCA_DIR_FAST="${EX_NGSPCA_DIR_FAST:-${NGSPCA_WORK_DIR}/ngspca_output_fast}"
EX_QC_DIR_FAST="${EX_QC_DIR_FAST:-${NGSPCA_WORK_DIR}/qc_output_fast}"

# Seed-control calibration (NGS-PCA's step 2b): the STANDARD tree's PCs
# recomputed under another -randomSeed. Same data, different draw of the
# randomized estimator - so how far the association statistics move between
# standard and seedctl is the estimator's own noise, the yardstick the
# fast-vs-standard differences are judged against. Resolved only when the
# run exists; skipped with a note otherwise.
EX_SEED_CONTROL_SEED="${EX_SEED_CONTROL_SEED:-43}"
EX_NGSPCA_DIR_SEEDCTL="${EX_NGSPCA_DIR_SEEDCTL:-${NGSPCA_WORK_DIR}/ngspca_output_seed${EX_SEED_CONTROL_SEED}}"

# Sample-ID suffix carried by the NGS-PCA outputs (mosdepth prefix `.by1000`
# plus the trailing dot NGS-PCA leaves when it removes `regions.bed.gz`).
# The depth matrix strips the same mosdepth suffixes itself, so after the
# prepare stage everything joins on bare IDs like HG00096. If you changed
# MOSDEPTH_BIN_SIZE upstream, change this to match (and see the README note
# on --sampleIdPattern).
EX_SAMPLE_SUFFIX="${EX_SAMPLE_SUFFIX:-.by1000.}"

# Where the per-sample coverage median (the correction stage's denominator)
# comes from:
#   auto    NGS-PCA's own autosomal.median.txt when the run wrote one - the
#           median over exactly the bins PCA used - else the QC table
#   ngspca  insist on autosomal.median.txt (error if absent)
#   qc      HQ_MEDIAN_COV from sample_qc.tsv, as the earlier upstream did
EX_MEDIAN_SOURCE="${EX_MEDIAN_SOURCE:-auto}"

# --- GitHub fallback -------------------------------------------------------
#
# When a mode's NGS-PCA outputs are not on disk, 00_fetch_inputs.sh can pull
# the committed results from the NGS-PCA repository instead (PCs and the QC
# table only — mosdepth files are not committed there). This is what makes
# the example testable without rerunning the upstream pipeline. The seed
# control never comes from GitHub: it is a local calibration run or nothing.
EX_GITHUB_RAW_BASE="${EX_GITHUB_RAW_BASE:-https://raw.githubusercontent.com}"
EX_GITHUB_REPO="${EX_GITHUB_REPO:-jlanej/NGS-PCA}"
EX_GITHUB_REF="${EX_GITHUB_REF:-master}"
EX_GITHUB_PATH_STANDARD="${EX_GITHUB_PATH_STANDARD:-example/1000G_highcov/output}"
EX_GITHUB_PATH_FAST="${EX_GITHUB_PATH_FAST:-example/1000G_highcov/output_fast}"

# --- this example's working area -------------------------------------------

EX_WORK_DIR="${EX_WORK_DIR:-/scratch/${USER}/depthsv_1000G_highcov}"

# Which modes to run, in order. standard first: seedctl shares its matrix.
# A mode whose upstream inputs are missing is skipped, so the default is
# safe before the seed-control run exists; EX_MODES=standard runs one tree.
EX_MODES="${EX_MODES:-standard fast seedctl}"

# --- smoke mode ------------------------------------------------------------
#
# EX_SMOKE=1 makes the example runnable anywhere, with no cluster and no
# CRAMs: 00_fetch_inputs.sh downloads the committed NGS-PCA outputs and then
# SIMULATES a small mosdepth tree per mode that is numerically consistent
# with each sample's real QC row (HQ median, X/Y ratios, chrM ratio), so
# every positive control still fires. The numbers it produces exercise the
# machinery; they are not results.
EX_SMOKE="${EX_SMOKE:-0}"
EX_SMOKE_SAMPLES="${EX_SMOKE_SAMPLES:-64}"
EX_SMOKE_SEED="${EX_SMOKE_SEED:-20260818}"

# --- pipeline parameters ---------------------------------------------------

# PCs removed by the correction stage. The upstream run computes 200; how
# many to remove is an empirical choice (see the per-region stats the
# correct stage writes). 20 is a starting point, not a recommendation.
#
# Smoke runs default lower: across a 64-sample subset the real leading PCs
# explain a large share of MTDNA_CN itself (measured R^2 ~ 0.67 at ndim=20,
# ~0.29 at ndim=4), so removing 20 of them from 64 samples blunts the very
# signal the checks assert on. The same absorption exists at full scale —
# that is what makes the ndim choice measurable — but there it dilutes over
# 3,202 samples instead of overwhelming the test.
if [ "${EX_SMOKE}" = "1" ]; then
    EX_NDIM="${EX_NDIM:-4}"
else
    EX_NDIM="${EX_NDIM:-20}"
fi

# Work-unit size for the SLURM array / local loop, in bp. 0 = one unit per
# contig. The default gives ~140 units over chr1-22,X,Y,M.
if [ "${EX_SMOKE}" = "1" ]; then
    EX_WINDOW="${EX_WINDOW:-1000000}"
else
    EX_WINDOW="${EX_WINDOW:-25000000}"
fi

# Only regions matching this pattern are corrected and analysed. The matrix
# still contains every contig mosdepth emitted (decoys, alts, EBV); the
# association sweep is restricted to the primary assembly plus chrM.
# (Assigned via a guard rather than ${:-...}: a `}` inside a brace-expansion
# default — as in a regex quantifier — would end the expansion early.)
if [ -z "${EX_CONTIG_REGEX:-}" ]; then
    EX_CONTIG_REGEX='^(chr)?([0-9][0-9]?|X|Y|M|MT)(:|$)'
fi

# Per-region QC applied at the analysis stage.
if [ "${EX_SMOKE}" = "1" ]; then
    EX_MIN_OBS="${EX_MIN_OBS:-30}"
else
    EX_MIN_OBS="${EX_MIN_OBS:-100}"
fi

# Seed for the permuted null phenotype (MTDNA_CN_NULL).
EX_PHENO_SEED="${EX_PHENO_SEED:-20260818}"

# --- parallelism (within one unit of work) ---------------------------------

EX_JOBS="${EX_JOBS:-8}"
EX_THREADS="${EX_THREADS:-4}"

# --- scheduler -------------------------------------------------------------

# slurm | local | auto (auto = slurm when sbatch is on PATH)
EX_RUNNER="${EX_RUNNER:-auto}"

# Resource requests. Site specifics (partition, account, QoS) belong in
# EX_SBATCH_EXTRA so the defaults stay portable.
EX_SBATCH_JOIN="${EX_SBATCH_JOIN:---cpus-per-task=12 --mem=48G --time=12:00:00}"
EX_SBATCH_UNIT="${EX_SBATCH_UNIT:---cpus-per-task=8 --mem=16G --time=4:00:00}"
EX_SBATCH_LIGHT="${EX_SBATCH_LIGHT:---cpus-per-task=2 --mem=8G --time=2:00:00}"
EX_SBATCH_EXTRA="${EX_SBATCH_EXTRA:-}"
EX_ARRAY_THROTTLE="${EX_ARRAY_THROTTLE:-100}"

# --- evaluation ------------------------------------------------------------

# Threshold profile for the evaluation checks: real | smoke. Smoke runs are
# tiny, so the calibration bands are wider there.
if [ "${EX_SMOKE}" = "1" ]; then
    EX_EVAL_PROFILE="${EX_EVAL_PROFILE:-smoke}"
else
    EX_EVAL_PROFILE="${EX_EVAL_PROFILE:-real}"
fi

# Top-K (by |t|) used for the cross-mode overlap metric.
EX_TOP_K="${EX_TOP_K:-100}"

# Calibration verdict: fast mode is "within seed noise" for an analysis when
# its distance from standard, 1 - r(stat), is at most this many times the
# seed-control's distance. 1 would demand fast mode be no worse than a
# reseed; the default leaves room for the reseed being an unusually close
# draw.
EX_CALIBRATION_FACTOR="${EX_CALIBRATION_FACTOR:-1.5}"

# --- derived layout (not usually overridden) -------------------------------

EX_INPUTS_DIR="${EX_INPUTS_DIR:-${EX_WORK_DIR}/inputs}"
EX_RUN_DIR="${EX_RUN_DIR:-${EX_WORK_DIR}/work}"
EX_REGIONS_DIR="${EX_REGIONS_DIR:-${EX_WORK_DIR}/regions}"
EX_EVAL_DIR="${EX_EVAL_DIR:-${EX_WORK_DIR}/eval}"
EX_COMPARE_DIR="${EX_COMPARE_DIR:-${EX_WORK_DIR}/compare}"
EX_PROFILE_DIR="${EX_PROFILE_DIR:-${EX_WORK_DIR}/profile}"
EX_LOG_DIR="${EX_LOG_DIR:-${EX_WORK_DIR}/logs}"
EX_CACHE_DIR="${EX_CACHE_DIR:-${EX_WORK_DIR}/github_cache}"
EX_SMOKE_DIR="${EX_SMOKE_DIR:-${EX_WORK_DIR}/smoke_mosdepth}"
