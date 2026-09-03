#!/usr/bin/env bash
# scripts/join.sh on 3,202 real-ID samples x 100k bins, with the example's
# resources (--jobs 8 --threads 4), timed with peak RSS; then the scratch
# footprint over time, sampled every 2 s.
cd "$(dirname "$0")"
ROOT=/Users/Kitty/git/depthSV/.claude/worktrees/agent-a265e6aa144253f7b
rm -rf g3202/join
mkdir -p g3202/join
( while true; do du -sk g3202/join 2>/dev/null | awk -v t="$(date +%s)" '{print t"\t"$1}'; sleep 2; done ) > g3202/join_du.log &
dupid=$!
/usr/bin/time -l -o g3202/join_time.txt bash "$ROOT/scripts/join.sh" --manifest g3202/manifest.txt --out g3202/join --jobs 8 --threads 4 > g3202/join.log 2>&1
rc=$?
kill $dupid 2>/dev/null
echo "join exit $rc"
cat g3202/join_time.txt | grep -E 'real|user|sys|maximum resident'
grep -E 'joining|batch 0000(01|10|20|30|40|50|57)/|row parity|wrote|manifest' g3202/join.log
ls -la g3202/join
echo "peak scratch KB: $(sort -k2 -n g3202/join_du.log | tail -1)"
