#!/usr/bin/env bash
S=/private/tmp/claude-501/-Users-Kitty-git-depthSV/1db72424-c4a1-45f4-b81e-5b466f3ec8a4/scratchpad/review_example
W="$S/smoke"
echo "=== run log (head 60)"; head -60 "$S/smoke.run.log"
echo; echo "=== paths.env"; for m in standard fast; do echo "-- $m"; cat "$W/inputs/$m/paths.env"; done
echo; echo "=== analyses.tsv (standard)"; cat "$W/inputs/standard/analyses.tsv"
echo; echo "=== prepare.summary (standard)"; cat "$W/inputs/standard/prepare.summary.txt"
echo; echo "=== prepare.summary (fast)"; cat "$W/inputs/fast/prepare.summary.txt"
echo; echo "=== cross_mode_samples"; cat "$W/inputs/cross_mode_samples.txt"
echo; echo "=== regions"; cat "$W/regions/standard.regions.txt"
echo; echo "=== assoc dir listing"; ls "$W/work/standard/" ; ls "$W/work/standard/assoc_ndim4" | head -80
echo; echo "=== eval standard"; cat "$W/eval/standard/summary.md"
echo; echo "=== eval fast (verdict + differing lines)"; grep -E 'verdict|FAIL|WARN' "$W/eval/fast/summary.md"
echo; echo "=== compare summary"; cat "$W/compare/summary.md"
echo; echo "=== standard_vs_fast summary"; cat "$W/compare/standard_vs_fast/summary.md"
echo; echo "=== concordance.tsv"; cat "$W/compare/standard_vs_fast/concordance.tsv"
echo; echo "=== profile report"; cat "$W/profile/profile_report.md"
echo; echo "=== smoke params"; cat "$W/smoke_mosdepth/fast/smoke.params.txt"
echo; echo "=== phenotypes head"; head -3 "$W/inputs/standard/phenotypes.tsv"
echo; echo "=== corrected log (standard chrM)"; cat "$W/work/standard/corrected/"corrected_ndim4.chrM*.log
echo; echo "=== top hits log2"; head -5 "$W/eval/standard/top_hits.log2_mtdna_cn.tsv"
