#!/bin/bash
# Annotate the ZNF175 union variant table with Biofilter 4 (report: annotation_master_variant).
# Returns, per variant: VEP consequence, LoF confidence (LOFTEE) = pLOF, AlphaMissense score/class,
# plus REVEL / CADD / SpliceAI / Pangolin / SIFT / PolyPhen and gnomAD frequencies.
#
# Usage: bash run_biofilter_annotation.sh <UNION_CSV> <OUTDIR>
#   UNION_CSV : results/02/cmp_union_all_variants.csv  (cols: key,CHROM,POS,REF,ALT,...)
#   OUTDIR    : where to write znf175_bf4_input.txt and znf175_bf4_annot.csv
set -euo pipefail

UNION_CSV="$1"
OUTDIR="$2"
mkdir -p "$OUTDIR"

SIF=/project/hall_shared/biofilter/images/bf4-hpc-4.1.4.sif
PG=/project/hall_shared/biofilter/databases/20260514/pgdata

# 1) build chr:pos:ref:alt input from the union table (skip '*' spanning-deletion alleles: not real alleles)
awk -F, 'NR>1 && $5!="*" {print $2":"$3":"$4":"$5}' "$UNION_CSV" > "$OUTDIR/znf175_bf4_input.txt"
echo "[$(date)] BF4 input variants: $(wc -l < "$OUTDIR/znf175_bf4_input.txt")"

# 2) run Biofilter; most_severe_only -> ~1 row per variant (most-severe transcript)
module load apptainer 2>/dev/null || true
T=$(mktemp -d); mkdir -p "$T/tmp" "$T/pg-run"
apptainer run --writable-tmpfs --pwd /tmp \
  --bind "$PG":/var/lib/postgresql/data \
  --bind "$T/tmp":/tmp \
  --bind "$T/pg-run":/var/run/postgresql \
  --bind "$OUTDIR":/workspace \
  "$SIF" biofilter report run \
    --report-name annotation_master_variant \
    --input-file /workspace/znf175_bf4_input.txt \
    --param most_severe_only=true \
    --output /workspace/znf175_bf4_annot.csv
rm -rf "$T"
echo "[$(date)] DONE -> $OUTDIR/znf175_bf4_annot.csv"
