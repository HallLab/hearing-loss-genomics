#!/bin/bash
# Submit the ZNF175 v1 (Park / Freeze One) extraction to LSF via the SAME pipeline as v2
# (bcftools region-extract + split + plink2 freq/geno-counts) for symmetric normalization.
set -euo pipefail
BASE=/project/hall/analysis/hearing-loss-genomics
SCRIPTS=$BASE/analysis/chapter_2/scripts/02
RESULTS=$BASE/analysis/chapter_2/results/02
LOGS=$SCRIPTS/logs
mkdir -p "$LOGS"

V1=/static/PMBB/PMBB_Freeze17/genotype/exome/all_variants/UPENN_Freeze_One_GRCh38.GL.pVCF.vcf.gz

bsub -J znf175_v1 -n 4 -M 8000 -W 60 \
     -o "$LOGS/znf175_v1.%J.out" -e "$LOGS/znf175_v1.%J.err" \
     bash "$SCRIPTS/extract_znf175_region.sh" \
          "$V1" "$RESULTS/v1_znf175_strict" "19:51571283-51592510"
