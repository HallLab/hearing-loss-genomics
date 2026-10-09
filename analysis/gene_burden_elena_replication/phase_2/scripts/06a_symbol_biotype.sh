#!/usr/bin/env bash
# Phase 2 / check 06a -- symbol -> biotype map from the release's own VEP, all 22
# chromosomes. Feeds the gene-level biotype measurement in check 06.
#
# Takes ~20 minutes and writes nothing until the end, because the whole scan is piped
# through one sort. Worth fixing if it is ever re-run often: emit per chromosome, sort
# afterwards, and progress becomes visible.
set -euo pipefail
R=/static/PMBB/PMBB-Release-2026-4.0/Exome/vep_annotations
OUT="$(cd "$(dirname "$0")/.." && pwd)/data/symbol_biotype.tsv"
: > $OUT
for c in $(seq 1 22); do
  awk -F'\t' -v OFS='\t' '
    NR==1{for(i=1;i<=NF;i++){if($i=="SYMBOL") s=i; if($i=="BIOTYPE") b=i}; next}
    $s!="" && $s!="-" && $b!="" {print $s, $b}
  ' $R/PMBB-Release-2026-4.0_genetic_exome.vep_annotations.chr$c.tsv
done | sort -u >> $OUT
echo "pares simbolo-biotipo distintos: $(wc -l < $OUT)"
