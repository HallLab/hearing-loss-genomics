#!/bin/bash
#
# Phase 4, step 0 -- split the group file by chromosome.
#
# ONE FILE, THREE MASKS. SAIGE derives pLOF, pDM and pLOF:pDM from a single group file
# that carries a per-variant annotation, so only pLOF_pDM is split here.
#
# TWO GROUP FILES HERE, not one. The BRaVa arm also tests the SYNONYMOUS control set,
# which is a separate category (BRaVa 4) and so a separate file. It is split by the same
# pass, with the same lost-line refusal, because a control set that quietly lost half
# its genes would still produce a plausible-looking flat result -- which is exactly the
# failure it exists to rule out.
#
# A CONSEQUENCE OF THAT FORMAT, recorded because it is a definitional choice and not a
# technicality: 11,540 variants are both pLOF and pDM. The group file carries one label
# per variant and gives pLOF precedence, so those variants enter the pLOF test and are
# ABSENT from the pDM test. The pDM arm therefore means "probably damaging AND not
# truncating", not "probably damaging".
#
# Kept that way deliberately. A variant that truncates the protein is already being
# tested by the arm that exists for truncation; counting it again under pDM measures the
# same thing twice, which is the same objection as the duplicate-cell problem, one level
# down. Reversing it means emitting a second group file with pDM precedence and running
# the pDM arm separately -- three times the step 2 tasks, and overlapping masks that
# raise the correlation between the cells Cauchy combines.
#
# Splitting is an I/O saving only: no gene line spans two chromosomes (checked in phase
# 2's merge), so it cannot change which genes are tested. Tab becomes space because that
# is the delimiter this container is known to accept.
#
set -euo pipefail
B=/project/hall/analysis/hearing-loss-genomics/analysis/pmbb_hl_brava

split_one () {
  IN=$1; OUT=$2
  mkdir -p ${OUT}
  [ -s "${IN}" ] || { echo "MISSING: ${IN}" >&2; exit 1; }
  echo "--- $(basename ${IN})"

# Each gene is TWO lines: "<gene> var <variants...>" then "<gene> anno <labels...>".
# The chromosome is readable only from the var line -- on the anno line the third field
# is a label like pLOF -- so it is taken there and carried to the pair. Writing the two
# lines to different files, or dropping one, produces a group file SAIGE cannot use.
awk -F'\t' -v out="${OUT}" '
  $2 == "var" {
    split($3, a, ":"); chr = a[1]
    if (chr !~ /^[0-9]+$/ || chr+0 < 1 || chr+0 > 22) { bad++; pending=""; next }
    pending = chr
    gsub(/\t/, " "); print > (out "/chr" chr ".txt"); n[chr]++
    next
  }
  $2 == "anno" {
    if (pending == "") { bad++; next }
    gsub(/\t/, " "); print > (out "/chr" pending ".txt"); n[pending]++
    pending = ""
    next
  }
  { bad++ }
  END {
    if (bad) printf "  WARNING: %d lines with an unparseable chromosome\n", bad
    for (c in n) printf "  chr%-3s %7d lines\n", c, n[c] | "sort -V"
  }' "${IN}"

  TOT=$(cat ${OUT}/chr*.txt | wc -l)
  SRC=$(wc -l < "${IN}")
  echo "  split: ${TOT} of ${SRC} lines placed"
  [ "${TOT}" -eq "${SRC}" ] || { echo "REFUSING: $(( SRC - TOT )) lines lost from ${IN}" >&2; exit 1; }
}

split_one ${B}/phase_2/results/masks/pLOF_pDM.txt   ${B}/phase_4/data/masks_by_chrom
split_one ${B}/phase_2/results/masks/synonymous.txt ${B}/phase_4/data/syn_by_chrom
