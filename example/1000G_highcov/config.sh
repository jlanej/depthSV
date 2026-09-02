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
NGSPCA_WORK_DIR="${NGSPCA_WORK_DIR:-/scratch/${USER:-$(id -un)}/1000G_highcov}"

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
#
# In a real (non-smoke) run that fallback pairs YOUR mosdepth depths with
# PCs and medians computed from SOMEONE ELSE'S run of the same samples — the
# PCs are only valid for the matrix they came from — so it is opt-in.
EX_ALLOW_GITHUB_FALLBACK="${EX_ALLOW_GITHUB_FALLBACK:-0}"
EX_GITHUB_RAW_BASE="${EX_GITHUB_RAW_BASE:-https://raw.githubusercontent.com}"
EX_GITHUB_REPO="${EX_GITHUB_REPO:-jlanej/NGS-PCA}"
EX_GITHUB_REF="${EX_GITHUB_REF:-master}"
EX_GITHUB_PATH_STANDARD="${EX_GITHUB_PATH_STANDARD:-example/1000G_highcov/output}"
EX_GITHUB_PATH_FAST="${EX_GITHUB_PATH_FAST:-example/1000G_highcov/output_fast}"

# --- this example's working area -------------------------------------------

# (lib.sh resolves the same default before loading inputs/run.env.)
EX_WORK_DIR="${EX_WORK_DIR:-/scratch/${USER:-$(id -un)}/depthsv_1000G_highcov}"

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

# --- preamble (preamble.sh): ndim and covariates ---------------------------

EX_PREAMBLE_DIR="${EX_PREAMBLE_DIR:-${EX_WORK_DIR}/preamble}"

# Marchenko-Pastur cut for the coverage PCs: relative margin above the noise
# edge (1% is about four Tracy-Widom sd at 3,202 x 142,070), and the ranks
# skipped between the signal estimate and the noise fit.
EX_MP_MARGIN="${EX_MP_MARGIN:-0.01}"
EX_MP_GAP="${EX_MP_GAP:-20}"

# Genotype PCs from the published PLINK-format NYGC callset.
EX_GENO_SOURCES="${EX_GENO_SOURCES:-${EX_EXAMPLE_DIR}/resources/genotype_sources.tsv}"
EX_GENO_LD_REGIONS="${EX_GENO_LD_REGIONS:-${EX_EXAMPLE_DIR}/resources/long_range_ld_grch38.txt}"
if [ "${EX_SMOKE}" = "1" ]; then
    EX_GENO_CHROMS="${EX_GENO_CHROMS:-22}"                  # one chromosome: minutes, ~60 MB
else
    EX_GENO_CHROMS="${EX_GENO_CHROMS:-1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22}"
fi
EX_GENO_MAF="${EX_GENO_MAF:-0.05}"
EX_GENO_GENO="${EX_GENO_GENO:-0.01}"
EX_GENO_LD_WINDOW="${EX_GENO_LD_WINDOW:-1000kb}"
EX_GENO_LD_R2="${EX_GENO_LD_R2:-0.1}"
EX_GENO_KING_CUTOFF="${EX_GENO_KING_CUTOFF:-0.0884}"      # removes 2nd-degree and closer
EX_GENO_NPC="${EX_GENO_NPC:-40}"                            # computed, so the MP tail is visible
EX_N_GPCS="${EX_N_GPCS:-10}"                                # used in the models
EX_GENO_THREADS="${EX_GENO_THREADS:-${SLURM_CPUS_PER_TASK:-8}}"
if [ -z "${EX_GENO_MEMORY_MB:-}" ]; then
    if [ -n "${SLURM_MEM_PER_NODE:-}" ]; then
        EX_GENO_MEMORY_MB=$(( SLURM_MEM_PER_NODE * 85 / 100 ))
    else
        EX_GENO_MEMORY_MB=16000
    fi
fi
# Environment modules to load before plink2, if the host provides `module`.
EX_PREAMBLE_MODULES="${EX_PREAMBLE_MODULES:-plink2}"

# Covariate terms appended to the mtDNA-CN models (SEX is dropped from the
# sex analyses automatically). Defaults to SEX plus the first EX_N_GPCS
# genotype PCs once the preamble has produced covariates.tsv; set to "none"
# for unadjusted models, or to any '+'-joined list of phenotype columns.
if [ -z "${EX_COVARIATES:-}" ]; then
    if [ -s "${EX_PREAMBLE_DIR}/covariates.tsv" ]; then
        EX_COVARIATES="SEX"
        i=1
        while [ "$i" -le "$EX_N_GPCS" ]; do EX_COVARIATES="${EX_COVARIATES}+GPC${i}"; i=$((i + 1)); done
        unset i
    else
        EX_COVARIATES="none"
    fi
fi

# --- pipeline parameters ---------------------------------------------------

# PCs removed by the correction stage. The upstream run computes 200; how
# many to remove is the preamble's Marchenko-Pastur decision when it has
# run (preamble/ndim.txt, averaged over modes so every tree is corrected
# identically); 20 otherwise, as a starting point rather than a
# recommendation. EX_NDIM in the environment overrides both.
#
# Smoke runs default lower: across a 64-sample subset the real leading PCs
# explain a large share of MTDNA_CN itself (measured R^2 ~ 0.67 at ndim=20,
# ~0.29 at ndim=4), so removing 20 of them from 64 samples blunts the very
# signal the checks assert on. The same absorption exists at full scale —
# that is what makes the ndim choice measurable — but there it dilutes over
# 3,202 samples instead of overwhelming the test.
if [ -z "${EX_NDIM:-}" ]; then
    if [ "${EX_SMOKE}" = "1" ]; then
        EX_NDIM=4
    elif [ -s "${EX_PREAMBLE_DIR}/ndim.txt" ]; then
        EX_NDIM="$(tr -cd '0-9' < "${EX_PREAMBLE_DIR}/ndim.txt")"
    else
        EX_NDIM=20
    fi
fi

# Work-unit size for the SLURM array / local loop, in bp. 0 = one unit per
# contig. 10 Mb gives ~310 units over chr1-22,X,Y,M — about 150 s and under
# 1 GB each at 3,202 samples, small enough to be cheap to lose.
if [ "${EX_SMOKE}" = "1" ]; then
    EX_WINDOW="${EX_WINDOW:-1000000}"
else
    EX_WINDOW="${EX_WINDOW:-10000000}"
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
# A region where one sample carries more than this share of the residual
# depth is a test of one participant and is skipped.
EX_MAX_SHARE="${EX_MAX_SHARE:-0.5}"

# The ploidy model for chrX/chrY (EX_PLOIDY=0 turns it off, and every
# sex-chromosome bin becomes a sex indicator again), the pseudo-autosomal
# BED for the build, and the floor on the log2 ratio.
EX_PLOIDY="${EX_PLOIDY:-1}"
EX_PAR="${EX_PAR:-$EX_EXAMPLE_DIR/../../conf/par.grch38.bed}"
EX_WINSOR_LOG2="${EX_WINSOR_LOG2:--3}"

# Seed for the permuted null phenotype (MTDNA_CN_NULL); the structured null
# and the coverage-PC null derive their seeds from it.
EX_PHENO_SEED="${EX_PHENO_SEED:-20260818}"
# Heritability of the structured null y ~ MVN(0, h2 * 2K + (1 - h2) I)
# drawn from the preamble's KING kinship.
EX_STRUCTURED_H2="${EX_STRUCTURED_H2:-0.5}"

# Permutations of the response per linear analysis (scripts/analyze.sh
# --perms), folded by the export step into the empirical genome-wide
# threshold; the seed is shared by every shard.
if [ "${EX_SMOKE}" = "1" ]; then
    EX_PERMS="${EX_PERMS:-50}"
else
    EX_PERMS="${EX_PERMS:-100}"
fi
EX_PERM_SEED="${EX_PERM_SEED:-1}"
# Export: rows with N, NCase or NControl below this are suppressed.
EX_MIN_COUNT="${EX_MIN_COUNT:-20}"

# --- SV-callset recovery (06_sv_recovery.sh) ------------------------------
#
# How well the corrected depth recovers known deletions as a function of
# ndim: carriers vs non-carriers of NYGC-callset deletions at several
# correction depths, so the Marchenko-Pastur count can be judged on the
# outcome that matters. Real runs download the callset VCF (~1 GB) once
# into EX_WORK_DIR/sv_callset; smoke runs use the deletions the simulated
# tree carries. EX_SV_CALLS points at a prepared calls table instead.
EX_SV_RECOVERY="${EX_SV_RECOVERY:-1}"
EX_SV_CALLSET_URL="${EX_SV_CALLSET_URL:-https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20210124.SV_Illumina_Integration/1KGP_3202.gatksv_svtools_novelins.freeze_V3.wAF.vcf.gz}"
EX_SV_CALLS="${EX_SV_CALLS:-}"
EX_SV_MIN_LEN="${EX_SV_MIN_LEN:-5000}"          # bp; bins are 1 kb
EX_SV_MIN_AF="${EX_SV_MIN_AF:-0.02}"            # carrier frequency band
EX_SV_MAX_AF="${EX_SV_MAX_AF:-0.5}"
if [ "${EX_SMOKE}" = "1" ]; then
    EX_SV_MAX_DELS="${EX_SV_MAX_DELS:-12}"
    EX_SV_NDIMS="${EX_SV_NDIMS:-0 2 4 8}"
else
    EX_SV_MAX_DELS="${EX_SV_MAX_DELS:-200}"
    EX_SV_NDIMS="${EX_SV_NDIMS:-0 5 10 20 40 60}"  # plus the MP count and EX_NDIM
fi

# --- parallelism (within one unit of work) ---------------------------------
#
# Defaults follow the allocation when there is one: workers take every CPU
# and the compression threads half of them (they overlap with the workers).
EX_JOBS="${EX_JOBS:-${SLURM_CPUS_PER_TASK:-8}}"
EX_THREADS="${EX_THREADS:-$(( EX_JOBS / 2 > 0 ? EX_JOBS / 2 : 1 ))}"

# A finalize lock older than this is treated as abandoned (an evaluate job
# killed by its time limit mid-comparison would otherwise block every later
# comparison).
EX_FINALIZE_LOCK_HOURS="${EX_FINALIZE_LOCK_HOURS:-4}"

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
