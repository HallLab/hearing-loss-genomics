#!/bin/bash
# Submit the ZNF175 Biofilter-4 annotation to LSF.
set -euo pipefail
BASE=/project/hall/analysis/hearing-loss-genomics
SCRIPTS=$BASE/analysis/chapter_2/scripts/03
RESULTS=$BASE/analysis/chapter_2/results
LOGS=$SCRIPTS/logs
mkdir -p "$LOGS"

bsub -J znf175_bf4 -n 2 -M 8000 -W 30 \
     -o "$LOGS/znf175_bf4.%J.out" -e "$LOGS/znf175_bf4.%J.err" \
     bash "$SCRIPTS/run_biofilter_annotation.sh" \
          "$RESULTS/02/cmp_union_all_variants.csv" "$RESULTS/03"
