#!/bin/bash
# Submit the ZNF175 VEP + ANNOVAR(dbNSFP4.5a) annotation to LSF.
set -euo pipefail
BASE=/project/hall/analysis/hearing-loss-genomics
SCRIPTS=$BASE/analysis/chapter_2/scripts/04
RESULTS=$BASE/analysis/chapter_2/results
LOGS=$SCRIPTS/logs
mkdir -p "$LOGS"

bsub -J znf175_vep -n 2 -M 16000 -W 60 \
     -o "$LOGS/znf175_vep.%J.out" -e "$LOGS/znf175_vep.%J.err" \
     bash "$SCRIPTS/run_vep_annovar_annotation.sh" \
          "$RESULTS/02/cmp_union_all_variants.csv" "$RESULTS/04"
