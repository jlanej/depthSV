#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# depthSV — shared shell helpers
#
# Source this at the top of every stage script:
#
#     source "$(dirname "${BASH_SOURCE[0]}")/../lib/common.sh"
#
# It sets the failure-handling policy for the whole pipeline and provides the
# few operations every stage needs: resolving a contig name from a tabix index,
# reading a region, and promoting an output atomically only after it has been
# verified.
#
# The policy these helpers exist to enforce: a stage that does not produce a
# complete, valid output must exit non-zero and must not leave anything behind
# that a later stage would mistake for a finished result.
# ---------------------------------------------------------------------------

set -euo pipefail

DSV_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export DSV_ROOT

# Pin linear-algebra threading to one thread per worker. Parallelism here comes
# from GNU parallel fanning out workers, so a threaded BLAS (OpenBLAS is common
# on clusters) multiplies: jobs x cores threads thrash the node. Pinning also
# keeps floating-point summation order — and therefore the output bytes —
# independent of how many cores the host happens to have. Already-set values
# are respected so a site can override deliberately.
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export VECLIB_MAXIMUM_THREADS="${VECLIB_MAXIMUM_THREADS:-1}"

# --- logging ---------------------------------------------------------------

# ISO-8601 UTC. `date -Is` is GNU-only and fails on BSD/macOS, which matters
# because the same scripts are expected to run on a workstation and a cluster.
dsv_now() { date -u +%Y-%m-%dT%H:%M:%SZ; }

# Modification time, portably. Order matters: on GNU coreutils `stat -f` means
# "filesystem status" and SUCCEEDS, printing free-block counts rather than a
# time — so a BSD-first fallback silently returns the wrong thing on Linux.
dsv_mtime() { stat -c %Y "$1" 2>/dev/null || stat -f %m "$1"; }

dsv_log()  { printf '[%s] %s\n' "$(dsv_now)" "$*" >&2; }
dsv_die()  { printf '[%s] ERROR: %s\n' "$(dsv_now)" "$*" >&2; exit 2; }

# Print the leading comment block as usage text. Deriving it from the comment
# itself means editing the header cannot silently truncate --help, which is
# what hand-tuned `sed` line ranges did.
dsv_usage() {
    awk 'NR>1 && /^#/ { sub(/^# ?/, ""); print; next } NR>1 { exit }' "$0" \
      | sed -e '/^-\{3,\}$/d' -e '/^$/{ $d; }'
    exit 0
}

# Report where we failed. Without this a mid-script failure under `set -e` is
# silent about which command died.
dsv_enable_error_trace() {
    trap 'rc=$?; [ $rc -ne 0 ] && printf "[%s] FAILED (exit %d) at line %d: %s\n" \
        "$(dsv_now)" "$rc" "$LINENO" "$BASH_COMMAND" >&2; exit $rc' ERR
}

# --- preconditions ---------------------------------------------------------

dsv_require_cmd() {
    local c
    for c in "$@"; do
        command -v "$c" >/dev/null 2>&1 || dsv_die "required command not found: $c"
    done
}

dsv_require_file() {
    local f
    for f in "$@"; do
        [ -s "$f" ] || dsv_die "required file missing or empty: $f"
    done
}

# Each name is a shell variable holding a CLI option's value; the error names
# the flag, not the variable, so it matches what the user typed.
dsv_require_opt() {
    local v
    for v in "$@"; do
        [ -n "${!v:-}" ] || dsv_die "--${v//_/-} is required"
    done
}

# Load environment modules if the host provides them. A machine without
# `module` — a container, a laptop — is not an error.
#
# Module names are site-specific, so they come from DSV_MODULES rather than
# being guessed here. This is advisory only: dsv_require_cmd is what actually
# decides whether the tools are present.
dsv_load_modules() {
    command -v module >/dev/null 2>&1 || return 0
    local m
    for m in ${DSV_MODULES:-}; do module load "$m" 2>/dev/null || true; done
}

# --- SIGPIPE-safe header read ----------------------------------------------
#
# `zcat big.gz | head -n 1` makes zcat die of SIGPIPE (exit 141). Under
# `pipefail` that aborts the script even though the read succeeded. Scope the
# exception narrowly rather than disabling pipefail globally.

dsv_header() {                     # dsv_header <file.gz>
    set +o pipefail
    gzip -cd "$1" | head -n 1
    set -o pipefail
}

# --- contig resolution -----------------------------------------------------
#
# Never assume a `chr` prefix. Ask the index what the file actually contains.
# Accepts either style on input and returns the name as it appears in the file.

dsv_resolve_contig() {             # dsv_resolve_contig <indexed.gz> <chrom>
    local f="$1" c="$2" hit
    [ -s "${f}.tbi" ] || [ -s "${f}.csi" ] || dsv_die "no tabix index beside $f (run: tabix -p bed $f)"
    hit="$(tabix -l "$f" | awk -v c="$c" '$0 == c || $0 == "chr" c || "chr" $0 == c {print; exit}')"
    [ -n "$hit" ] || dsv_die "chromosome '$c' is not present in $f (contigs: $(tabix -l "$f" | tr '\n' ' '))"
    printf '%s\n' "$hit"
}

# Emit the body rows for a region. `region` is a contig, optionally with a
# range: "1", "chr1", "chr1:1-10000000". This replaces full-file scans and is
# what makes an arbitrary interval a valid unit of work.
dsv_read_region() {                # dsv_read_region <indexed.gz> <region>
    local f="$1" region="$2" chrom range contig
    chrom="${region%%:*}"
    range=""
    [ "$region" != "$chrom" ] && range=":${region#*:}"
    contig="$(dsv_resolve_contig "$f" "$chrom")"
    tabix "$f" "${contig}${range}"
}

# --- atomic output ---------------------------------------------------------
#
# Nothing is written under its final name until it has been verified. A stage
# that dies part-way leaves a .tmp.gz for inspection and no .done marker, so a
# resume-aware dispatcher re-runs it rather than treating it as finished.

dsv_output_tmp()  { printf '%s.tmp.gz\n' "$1"; }
dsv_output_done() { printf '%s.done\n'   "$1"; }

# Clear anything a previous partial run left behind, unless it completed.
dsv_output_reset() {               # dsv_output_reset <finalPath>
    local final="$1"
    if [ ! -f "$(dsv_output_done "$final")" ]; then
        rm -f "$(dsv_output_tmp "$final")" "$final" "${final}.tbi" "${final}.csi"
    fi
}

# True if this unit of work is already finished.
dsv_output_complete() {            # dsv_output_complete <finalPath>
    [ -f "$(dsv_output_done "$1")" ] && [ -s "$1" ]
}

# Verify the temp output, index it, and only then promote it.
#
# Indexing happens before the rename on purpose: if tabix rejects the file the
# run must not leave a plausible-looking output behind under its final name.
dsv_output_commit() {              # dsv_output_commit <finalPath> [expected_body_rows]
    local final="$1" expected="${2:-}" tmp actual
    tmp="$(dsv_output_tmp "$final")"
    [ -s "$tmp" ] || dsv_die "no output produced at $tmp"

    bgzip -t "$tmp" || dsv_die "output failed bgzip integrity check: $tmp"

    if [ -n "$expected" ]; then
        actual="$(bgzip -dc "$tmp" | grep -cv '^#' || true)"
        if [ "$actual" -ne "$expected" ]; then
            dsv_die "incomplete output: $actual body rows, expected $expected. Leaving $tmp for inspection."
        fi
        dsv_log "row parity ok: $actual rows"
    fi

    tabix -f -p bed "$tmp" || dsv_die "tabix rejected $tmp (is the header line prefixed with '#'?)"

    mv "$tmp" "$final"
    mv "${tmp}.tbi" "${final}.tbi"
    echo "done" > "$(dsv_output_done "$final")"
    dsv_log "wrote $final"
}

# Shell-quote one value for embedding in a command string that a worker shell
# will parse (the GNU parallel templates below). Interpolating values raw broke
# any passthrough argument containing shell metacharacters — which is every
# realistic --sampleIdPattern, since PCREs are made of them.
dsv_q() { printf '%q' "$1"; }

# Sample name from a mosdepth file path: the basename with mosdepth's suffixes
# removed. Lives here because the join driver and its extraction worker must
# agree on it exactly.
dsv_sample_name() {
    local b="${1##*/}"
    b="${b%.gz}"; b="${b%.bed}"; b="${b%.regions}"
    b="${b%.by1000}"; b="${b%.src}"
    printf '%s\n' "$b"
}

# --- parallel --------------------------------------------------------------
#
# GNU parallel buffers a whole --block in memory and will grow the block at
# runtime if a record does not fit. Stating the block explicitly avoids that
# regrow entirely; it is also reported to drop records at a boundary on some
# versions, which is not something we have reproduced but is cheap to avoid.
#
# It must NOT be derived by multiplying a fixed line count by the row width:
# a row is one value per sample, so 2000 rows of a 500,000-sample matrix is
# ~15 GB of buffer. Bound the buffer instead and derive how many lines fit,
# which keeps memory flat as the cohort grows. Chunk boundaries do not affect
# results — every region is independent — so this is purely a memory decision.

DSV_BLOCK_BYTES="${DSV_BLOCK_BYTES:-67108864}"   # 64 MB

# Lines per chunk: as many as fit in one block, capped by the caller's request.
dsv_chunk_lines() {                # dsv_chunk_lines <sampleRowFile> <maxLines>
    local row_bytes n max="$2"
    row_bytes="$(wc -c < "$1" | tr -d ' ')"
    [ "${row_bytes:-0}" -gt 0 ] || row_bytes=1
    n=$(( DSV_BLOCK_BYTES / row_bytes ))
    [ "$n" -ge 1 ] || n=1
    [ "$n" -le "$max" ] || n="$max"
    printf '%d\n' "$n"
}

# --halt now,fail=1 means one dead worker tears the whole run down non-zero
# instead of silently producing a short file.
DSV_PARALLEL_FLAGS=(--pipe -k --halt now,fail=1)
