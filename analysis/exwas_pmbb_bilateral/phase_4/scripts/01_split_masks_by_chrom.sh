#!/bin/bash
#
# Phase 4, step 0 -- split the masks into per-chromosome group files.
#
# SAIGE step 2 runs one chromosome at a time but reads the whole group file, and
# step 2 here is 3 cohorts x 3 MAF cutoffs, so an unsplit file is re-parsed 9
# times per mask. No gene line in any mask spans two chromosomes (checked
# upstream), so splitting is purely an I/O saving and cannot change which genes
# are tested.
#
# Also converts tab to space: the masks are tab-delimited, and the format known
# to work with this container is space-delimited.
#
# Input : phase_2/results/masks_v2/{pLOF,pDM,pLOF_pDM}.txt
# Output: phase_4/data/masks_by_chrom/<mask>/chr<N>.txt
#
set -euo pipefail

BASE=/project/hall/analysis/hearing-loss-genomics/analysis/exwas_pmbb_bilateral
IN=${BASE}/phase_2/results/masks_v2
OUT=${BASE}/phase_4/data/masks_by_chrom

for MASK in pLOF pDM pLOF_pDM; do
  mkdir -p ${OUT}/${MASK}
  echo "splitting ${MASK} ..."
  awk -v out="${OUT}/${MASK}" '
    {
      if ($2 == "var") { split($3, a, ":"); chr = a[1] }
      if (chr == "") next
      line = $1
      for (i = 2; i <= NF; i++) line = line " " $i
      print line > (out "/chr" chr ".txt")
    }' "${IN}/${MASK}.txt"
done

echo
echo "=== genes per chromosome ==="
printf "%-6s %9s %9s %10s\n" chr pLOF pDM pLOF_pDM
for C in $(seq 1 22); do
  printf "%-6s" "chr${C}"
  for MASK in pLOF pDM pLOF_pDM; do
    F=${OUT}/${MASK}/chr${C}.txt
    N=0; [[ -f ${F} ]] && N=$(( $(wc -l < ${F}) / 2 ))
    printf " %9d" "${N}"
  done
  echo
done
