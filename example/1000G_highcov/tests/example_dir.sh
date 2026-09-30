#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# depthSV 1000G example — every entry script finds this directory
#
#   bash tests/example_dir.sh [scratchDir]
#
# Each stage script, run.sh and the preamble must find this directory before
# they can source lib.sh: an exported EX_EXAMPLE_DIR, else the script's own
# directory, else - for the copy sbatch runs from its spool - the submit
# directory. SLURM_SUBMIT_DIR once came before the script's own directory,
# and an interactive job (srun --pty, salloc) sets it too, to wherever the
# session started: a stage run by hand inside one looked for lib.sh there
# and fell through a cascade of "command not found" to exit 0. For every
# entry script this checks that
#
#   - run from here, or from anywhere by path, with SLURM_SUBMIT_DIR
#     elsewhere, it finds its own directory;
#   - as sbatch's spooled copy submitted from here, it finds this directory;
#   - as a spooled copy submitted from elsewhere, it stops with one clear
#     error (exit 2);
#   - an exported EX_EXAMPLE_DIR wins over both, and a wrong one is refused
#     by name rather than replaced.
#
# Every run is the script's --help, which loads lib.sh and config.sh and
# exits: no network, no R, nothing written. Exit status = the number of
# failed checks.
# ---------------------------------------------------------------------------

# The case_* functions run by name, through every(), and the stand-in for
# dsv_enable_error_trace in the scripts under test.
# shellcheck disable=SC2329
set -uo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
T="${1:-${TMPDIR:-/tmp}/depthsv-example-dir.$$}"
case "$T" in /|"$HOME"|"") echo "refusing to use '$T' as a scratch directory" >&2; exit 2 ;; esac
rm -rf "$T"; mkdir -p "$T"
T="$(cd "$T" && pwd)"

scripts="00_fetch_inputs.sh 01_prepare_inputs.sh 02_run_depthsv.sh 03_evaluate.sh
         04_compare_modes.sh 05_profile.sh 06_sv_recovery.sh preamble.sh run.sh"

# Nothing from the caller's shell: an EX_EXAMPLE_DIR exported there (a guide
# may set one) would mask exactly what is tested, and a real SLURM_SUBMIT_DIR
# would stand in for the ones set below. lib.sh only reads under EX_WORK_DIR.
for v in $(compgen -v | grep '^EX_' || true); do unset "$v"; done
unset SLURM_SUBMIT_DIR
export EX_WORK_DIR="$T/work"

# A script that failed to load lib.sh carries on without its functions, and
# its --help spins on "dsv_usage: command not found" for ever. Every entry
# script calls this right after loading lib.sh, which replaces this stand-in;
# without lib.sh the stand-in stops the script there, and says why.
dsv_enable_error_trace() { echo "lib.sh was not loaded"; exit 3; }
export -f dsv_enable_error_trace

elsewhere="$T/elsewhere"           # where an interactive job started: no lib.sh
decoy="$T/decoy"                   # an exported EX_EXAMPLE_DIR other than this one
mkdir -p "$elsewhere" "$decoy"
printf 'echo "decoy lib.sh"; exit 0\n' > "$decoy/lib.sh"
# What sbatch runs: a copy of the script, alone in a spool directory.
for s in $scripts; do mkdir -p "$T/spool/$s"; cp "$here/$s" "$T/spool/$s/slurm_script"; done

pass=0; fail=0
ok()  { pass=$((pass + 1)); printf '  PASS %s\n' "$1"; }
bad() { fail=$((fail + 1)); printf '  FAIL %s\n' "$1"
        [ -z "${2:-}" ] || [ ! -s "$2" ] || tail -n 12 "$2" | sed 's/^/       | /'; }

# run <log> <cwd> <script> [VAR=value ...]: the script's --help from <cwd>,
# under the bash running this test; its status is returned, its output logged.
run()   { local log="$1" cwd="$2" script="$3"; shift 3
          ( cd "$cwd" && env "$@" "$BASH" "$script" --help ) > "$log" 2>&1; }
# The usage title --help prints once lib.sh has loaded.
title() { sed -n 's/^# \(depthSV 1000G example .*\)/\1/p' "$here/$1" | head -n 1; }
usage() { [ "$1" -eq 0 ] && grep -qF "$(title "$2")" "$3"; }          # usage <rc> <script> <log>
# One clear line naming the directory, not a cascade of "command not found".
refused() { [ "$1" -eq 2 ] && [ "$(grep -c . "$3")" -eq 1 ] \
                && grep -qF "no lib.sh in $2;" "$3" && ! grep -q 'command not found' "$3"; }
decoyed() { [ "$1" -eq 0 ] && [ "$(cat "$2")" = "decoy lib.sh" ]; }    # decoyed <rc> <log>

# The cases, one per function: <script> <log>.
# The failure this guards: from here, in an interactive job started elsewhere.
case_here()           { run "$2" "$here" "$1" SLURM_SUBMIT_DIR="$elsewhere"; usage $? "$1" "$2"; }
case_by_path()        { run "$2" "$elsewhere" "$here/$1" SLURM_SUBMIT_DIR="$elsewhere"; usage $? "$1" "$2"; }
case_spool()          { run "$2" "$here" "$T/spool/$1/slurm_script" SLURM_SUBMIT_DIR="$here"; usage $? "$1" "$2"; }
case_spool_lost()     { run "$2" "$elsewhere" "$T/spool/$1/slurm_script" SLURM_SUBMIT_DIR="$elsewhere"; refused $? "$elsewhere" "$2"; }
case_export()         { run "$2" "$here" "$1" EX_EXAMPLE_DIR="$decoy" SLURM_SUBMIT_DIR="$here"; decoyed $? "$2"; }
case_export_spool()   { run "$2" "$here" "$T/spool/$1/slurm_script" EX_EXAMPLE_DIR="$decoy" SLURM_SUBMIT_DIR="$here"; decoyed $? "$2"; }
case_export_wrong()   { run "$2" "$here" "$1" EX_EXAMPLE_DIR="$elsewhere" SLURM_SUBMIT_DIR="$here"; refused $? "$elsewhere" "$2"; }

# every <case> <description>: the case for every entry script, one line.
every() {
    local s failed="" first=""
    for s in $scripts; do
        "case_$1" "$s" "$T/$1.$s.log" || { failed="$failed $s"; first="${first:-$s}"; }
    done
    if [ -z "$failed" ]; then ok "$2"; else bad "$2 - not:$failed" "$T/$1.$first.log"; fi
}

echo "every entry script finds this directory -> $T"
every here         "run from here with SLURM_SUBMIT_DIR elsewhere, each finds its own directory"
every by_path      "  and by path from anywhere"
every spool        "as sbatch's spooled copy submitted from here, each finds this directory"
every spool_lost   "  submitted from elsewhere, each stops with one clear error"
every export       "an exported EX_EXAMPLE_DIR wins over the script's own directory"
every export_spool "  and over the submit directory"
every export_wrong "  and a wrong one is refused by name, not replaced"

echo
echo "passed $pass, failed $fail   ($T)"
exit "$fail"
