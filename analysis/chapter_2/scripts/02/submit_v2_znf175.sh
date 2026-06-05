#!/bin/bash
# Submit the ZNF175 v2 (Hui / Release 2020 v2.0) extraction to LSF.
set -euo pipefail
BASE=/project/hall/analysis/hearing-loss-genomics
SCRIPTS=$BASE/analysis/chapter_2/scripts/02
RESULTS=$BASE/analysis/chapter_2/results/02
LOGS=$SCRIPTS/logs
mkdir -p "$LOGS"

V2=$BASE/data/pmbb_v2/Exome/pVCF/GL_by_chrom/PMBB-Release-2020-2.0_genetic_exome_chr19_GL.vcf.gz

bsub -J znf175_v2 -n 4 -M 8000 -W 60 \
     -o "$LOGS/znf175_v2.%J.out" -e "$LOGS/znf175_v2.%J.err" \
     bash "$SCRIPTS/extract_znf175_region.sh" \
          "$V2" "$RESULTS/v2_znf175_release" "19:51571283-51592510"
