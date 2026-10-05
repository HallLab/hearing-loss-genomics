#!/bin/bash
#
# Phase 4, step 1 of 3 -- split the rebuilt masks into per-chromosome group files.
#
# Why: SAIGE step 2 runs one chromosome at a time, but reads the whole group
# file every time.  Our ALL.txt is 317 MB and step 2 is 3 cohorts x 22 chr x
# 3 MAF cutoffs, so the unsplit file would be re-parsed 198 times.  No gene
# line in any of the four masks spans two chromosomes (checked), so splitting
# is purely an I/O saving -- it cannot change which genes get tested.
#
# Also converts tab to space.  Elena's mask files are space-delimited and ran
# through SAIGE 1.5.0 successfully; ours came out tab-delimited.  Both are
# legal, but we keep the format that is known to work here.
#
# Input : phase_2/results/masks_v2/{ALL,pLOF,pDM,pLOF_pDM}.txt   (rebuilt masks)
# Output: phase_4/data/masks_by_chrom/<mask>/chr<N>.txt
#
set -euo pipefail

REPL=/project/hall/analysis/hearing-loss-genomics/analysis/elena_replication
IN=${REPL}/phase_2/results/masks_v2
OUT=${REPL}/phase_4/data/masks_by_chrom

for MASK in pLOF pDM pLOF_pDM ALL; do
  mkdir -p ${OUT}/${MASK}
  echo "splitting ${MASK} ..."
  awk -v out="${OUT}/${MASK}" '
    {
      # chromosome comes from the first variant ID on the "var" line;
      # the "anno" line that follows inherits it
      if ($2 == "var") { split($3, a, ":"); chr = a[1] }
      if (chr == "") next
      line = $1
      for (i = 2; i <= NF; i++) line = line " " $i
      print line > (out "/chr" chr ".txt")
    }' "${IN}/${MASK}.txt"
done

echo
echo "=== genes per chromosome ==="
printf "%-5s %8s %8s %10s %10s\n" chr pLOF pDM pLOF_pDM ALL
for C in $(seq 1 22); do
  printf "%-5s" "chr${C}"
  for MASK in pLOF pDM pLOF_pDM ALL; do
    F=${OUT}/${MASK}/chr${C}.txt
    N=0; [[ -f ${F} ]] && N=$(( $(wc -l < ${F}) / 2 ))
    printf " %8d" "${N}"
  done
  echo
done
