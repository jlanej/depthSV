#!/usr/bin/env bash
S=/private/tmp/claude-501/-Users-Kitty-git-depthSV/1db72424-c4a1-45f4-b81e-5b466f3ec8a4/scratchpad/review_example
D="$S/spectrum"
mkdir -p "$D"
cd "$D" || exit 9
B=https://raw.githubusercontent.com/jlanej/NGS-PCA/master/example/1000G_highcov
for f in svd.singularvalues.txt svd.samples.txt svd.bins.txt autosomal.median.txt svd.pcs.txt; do
  if curl -fsSL --retry 3 -o "$f" "$B/output/ngspca_output/$f"; then
    echo "std ok $f $(wc -l < "$f") lines"
  else
    echo "std MISSING $f"
  fi
done
for f in svd.singularvalues.txt svd.samples.txt svd.bins.txt autosomal.median.txt svd.pcs.txt; do
  if curl -fsSL --retry 3 -o "fast.$f" "$B/output_fast/ngspca_output/$f"; then
    echo "fast ok $f $(wc -l < "fast.$f") lines"
  else
    echo "fast MISSING $f"
  fi
done
curl -fsSL -o std.sample_qc.tsv "$B/output/qc_output/sample_qc.tsv" && echo "std qc $(wc -l < std.sample_qc.tsv)"
curl -fsSL -o fast.sample_qc.tsv "$B/output_fast/qc_output/sample_qc.tsv" && echo "fast qc $(wc -l < fast.sample_qc.tsv)" || echo "fast qc MISSING"
echo "--- singular values head"
head -3 svd.singularvalues.txt
echo "--- samples head"
head -2 svd.samples.txt
echo "--- bins head"
head -2 svd.bins.txt
echo "--- qc header"
head -1 std.sample_qc.tsv | tr '\t' '\n' | nl
echo "--- median head"
head -3 autosomal.median.txt
echo "--- pcs head (first 5 cols)"
head -2 svd.pcs.txt | cut -f1-5
