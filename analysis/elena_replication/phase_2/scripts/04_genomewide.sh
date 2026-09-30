#!/usr/bin/env bash
# Phase 2 / check 04 — is chr8 representative? Same measurement, all 22 chromosomes.
#
# Counts at the VARIANT level (not the annotation-row level): a variant is in the
# pLOF mask if any of its transcript rows has is_pLOF=True, and it has genuine
# loss-of-function support if any row carries a real LoF consequence.
#
# "Real" LoF, matched as an exact comma-separated term:
#   stop_gained frameshift_variant splice_acceptor_variant splice_donor_variant
#   start_lost stop_lost transcript_ablation
#
# The three terms in the pipeline's lof_terms list that are NOT loss-of-function
# (VEP IMPACT=LOW) are what this measures the cost of:
#   splice_polypyrimidine_tract_variant  splice_donor_region_variant
#   splice_donor_5th_base_variant
set -euo pipefail
CLS=/project/hall/analysis/hearing-loss-genomics/analysis/elena/rarevariantExWAS/variant_categories
OUT="$(cd "$(dirname "$0")/.." && pwd)/results/04_genomewide.tsv"
printf 'chrom\tplof_variants\twith_real_lof\tno_real_lof\tpct_no_real_lof\n' > "$OUT"
for c in $(seq 1 22); do
  awk -F'\t' -v C="$c" '
    BEGIN{ split("stop_gained frameshift_variant splice_acceptor_variant splice_donor_variant start_lost stop_lost transcript_ablation", a, " ")
           for (i in a) real[a[i]]=1 }
    NR==1{ for(i=1;i<=NF;i++) h[$i]=i; next }
    { v=$h["#Uploaded_variation"]
      if (tolower($h["is_pLOF"])=="true") plof[v]=1
      n=split($h["Consequence"], t, ",")
      for(i=1;i<=n;i++) if (t[i] in real) { rl[v]=1; break } }
    END{ p=0; r=0
         for (v in plof) { p++; if (v in rl) r++ }
         printf "%s\t%d\t%d\t%d\t%.1f\n", C, p, r, p-r, (p?100*(p-r)/p:0) }
  ' "$CLS/chr$c.classified.tsv" >> "$OUT"
done
awk -F'\t' 'NR>1{p+=$2; r+=$3} END{printf "TOTAL\t%d\t%d\t%d\t%.1f\n", p, r, p-r, 100*(p-r)/p}' "$OUT" >> "$OUT"
column -t "$OUT"
