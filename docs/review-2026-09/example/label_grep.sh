#!/usr/bin/env bash
S=/private/tmp/claude-501/-Users-Kitty-git-depthSV/1db72424-c4a1-45f4-b81e-5b466f3ec8a4/scratchpad/review_example
W="$S/smoke"
echo "=== files mentioning synthetic / simulated / jitter anywhere under the smoke work dir:"
grep -rl -i -E 'synthetic|simulated|jitter' "$W" --include='*.md' --include='*.tsv' --include='*.env' --include='*.txt' 2>/dev/null | sed "s#$W/##"
echo "=== reports (eval/compare/profile) mentioning them:"
grep -l -i -E 'synthetic|simulated|jitter|smoke' "$W"/eval/*/summary.md "$W"/compare/summary.md "$W"/compare/*/summary.md "$W"/profile/profile_report.md 2>/dev/null | sed "s#$W/##"
echo "=== what those reports say about smoke:"
grep -h -i -E 'synthetic|simulated|jitter|smoke' "$W"/eval/*/summary.md "$W"/compare/summary.md "$W"/compare/*/summary.md "$W"/profile/profile_report.md 2>/dev/null | sort -u
