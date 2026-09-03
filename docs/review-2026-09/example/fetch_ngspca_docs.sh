#!/usr/bin/env bash
S=/private/tmp/claude-501/-Users-Kitty-git-depthSV/1db72424-c4a1-45f4-b81e-5b466f3ec8a4/scratchpad/review_example
D="$S/ngspca_docs"
mkdir -p "$D"
cd "$D" || exit 9
B=https://raw.githubusercontent.com/jlanej/NGS-PCA/master/example/1000G_highcov
curl -fsSL -o README.md "$B/README.md" && echo "README $(wc -l < README.md) lines"
curl -fsSL -o listing.json "https://api.github.com/repos/jlanej/NGS-PCA/contents/example/1000G_highcov" && echo "listing ok"
grep -o '"name": *"[^"]*"' listing.json | sed 's/"name": *//' | tr -d '"' | tr '\n' ' '; echo
curl -fsSL -o output_listing.json "https://api.github.com/repos/jlanej/NGS-PCA/contents/example/1000G_highcov/output" && echo "output listing ok"
grep -o '"name": *"[^"]*"' output_listing.json | sed 's/"name": *//' | tr -d '"' | tr '\n' ' '; echo
curl -fsSL -o ngspca_out_listing.json "https://api.github.com/repos/jlanej/NGS-PCA/contents/example/1000G_highcov/output/ngspca_output" && echo "ngspca_output listing ok"
grep -o '"name": *"[^"]*"' ngspca_out_listing.json | sed 's/"name": *//' | tr -d '"' | tr '\n' ' '; echo
curl -fsSL -o commit_7862645.json "https://api.github.com/repos/jlanej/NGS-PCA/commits/7862645" && echo "commit 7862645: $(grep -o '"date": *"[^"]*"' commit_7862645.json | head -1)" || echo "commit 7862645 not found"
curl -fsSL -o head_commit.json "https://api.github.com/repos/jlanej/NGS-PCA/commits/master" && echo "master head: $(grep -o '"sha": *"[^"]*"' head_commit.json | head -1) $(grep -o '"date": *"[^"]*"' head_commit.json | head -1)"
for f in 03a_mosdepth_coverage_summary.sh 03_collect_qc.sh 02_run_ngspca.sh 04_fast_mode_eval.sh 00_setup.sh; do
  curl -fsSL -o "$f" "$B/$f" && echo "fetched $f ($(wc -l < $f) lines)" || echo "MISSING $f"
done
