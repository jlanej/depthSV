#!/bin/bash
# Quick portability probes, run under whatever bash invokes this file.
echo "bash: $BASH_VERSION"
S=/private/tmp/claude-501/-Users-Kitty-git-depthSV/1db72424-c4a1-45f4-b81e-5b466f3ec8a4/scratchpad/review_shell

echo "--- 1. BSD awk numeric compare as used in 01_prepare_inputs.sh:131 ---"
printf 'chr1\t0\t999999\nchr1\t999999\t3000000\nchr2\t0\t5\n' | awk -F'\t' '
    !($1 in max) { order[++n] = $1 }
    $3 > max[$1] { max[$1] = $3 }
    END { for (i = 1; i <= n; i++) print order[i] "\t" max[order[i]] }'

echo "--- 2. local arr=() and += on this bash ---"
f() { local a="$1" ff=(); ff+=("x y"); ff+=("$2"); printf '<%s>' "${ff[@]}"; echo; }
f one "two three"
g() { local sb=(); while [ $# -gt 0 ] && [ "$1" != "--" ]; do sb+=("$1"); shift; done; shift; echo "sb=${#sb[@]} rest=$*"; }
g --a --b -- --stage x

echo "--- 3. printf %q round-trip through source (paths.env) ---"
p="$S/dir with space/\$dollar/ünïcode'q\"uote"
printf 'EX_M_X=%q\n' "$p" > "$S/rt.env"
cat "$S/rt.env"
unset EX_M_X; source "$S/rt.env"
[ "$EX_M_X" = "$p" ] && echo "round-trip OK" || echo "round-trip BROKEN: <$EX_M_X>"
printf 'EX_M_EMPTY=%q\n' "" > "$S/rt2.env"; cat "$S/rt2.env"; source "$S/rt2.env"; echo "empty -> <${EX_M_EMPTY-unset}>"

echo "--- 4. EXIT trap on SIGTERM (mktemp cleanup in correct.sh/analyze.sh) ---"
rm -f "$S/trapflag"
bash -c 'trap "echo trap-ran > '"$S"'/trapflag" EXIT; sleep 30' &
pid=$!; sleep 0.5; kill -TERM $pid; wait $pid 2>/dev/null; sleep 0.3
[ -f "$S/trapflag" ] && echo "EXIT trap ran on SIGTERM" || echo "EXIT trap did NOT run on SIGTERM"

echo "--- 5. sort|uniq -d|head -3 under pipefail with many duplicates (join.sh:85-91) ---"
set -o pipefail
yes /x/SAMPLE$RANDOM.by1000.regions.bed.gz 2>/dev/null | head -20000 > "$S/dupman.txt"
awk '{print $0; print $0}' "$S/dupman.txt" | sed 's/SAMPLE/SAMPLE_/' > /dev/null
# 20000 identical lines -> 1 dup only; make many distinct dups:
awk 'BEGIN{for(i=0;i<20000;i++){print "/x/S" i ".by1000.regions.bed.gz"; print "/x/S" i ".by1000.regions.bed.gz"}}' > "$S/dupman.txt"
dsv_sample_name() { local b="${1##*/}"; b="${b%.gz}"; b="${b%.bed}"; b="${b%.regions}"; b="${b%.by1000}"; printf '%s\n' "$b"; }
dup="$( while IFS= read -r f; do [ -n "$f" ] || continue; dsv_sample_name "$f"; done < "$S/dupman.txt" | sort | uniq -d | head -3 )"; rc=$?
echo "rc=$rc dup=$(printf '%s' "$dup" | tr '\n' ' ')"
set +o pipefail

echo "--- 6. BSD sed usage filter as in dsv_usage ---"
printf '#!/x\n# a\n#\n# ---\n# b\n#\n\nbody\n' > "$S/us.sh"
awk 'NR>1 && /^#/ { sub(/^# ?/, ""); print; next } NR>1 { exit }' "$S/us.sh" | sed -e '/^-\{3,\}$/d' -e '/^$/{ $d; }' | cat -A 2>/dev/null || awk 'NR>1 && /^#/ { sub(/^# ?/, ""); print; next } NR>1 { exit }' "$S/us.sh" | sed -e '/^-\{3,\}$/d' -e '/^$/{ $d; }' | od -c | head

echo "--- 7. wc -l inside arithmetic on macOS ---"
printf 'a\nb\n' > "$S/wc.txt"; echo "rows=$(( $(wc -l < "$S/wc.txt") - 1 ))"

echo "--- 8. grep -c . on empty file, assignment under set -e ---"
: > "$S/empty.txt"
( set -e; n="$(grep -c . "$S/empty.txt")"; echo "survived n=$n" ) || echo "set -e killed the assignment (rc=$?)"
( set -e; n="$(grep -c . "$S/empty.txt" || true)"; echo "with ||true n=$n" )
