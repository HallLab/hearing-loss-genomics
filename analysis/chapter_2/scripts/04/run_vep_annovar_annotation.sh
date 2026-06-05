#!/bin/bash
# Annotate the ZNF175 union variants by COORDINATE (covers ALL variants, incl. private ones BF4 missed):
#   VEP    -> consequence + IMPACT (HIGH = pLOF), per variant   (independent of any variant DB)
#   ANNOVAR dbNSFP4.5a (canonical) -> AlphaMissense score/pred (+ dbNSFP scores) for missense
# Does NOT use Biofilter. Window/variants come from the union table.
#
# Usage: bash run_vep_annovar_annotation.sh <UNION_CSV> <OUTDIR>
set -euo pipefail
UNION_CSV="$1"
OUTDIR="$2"
mkdir -p "$OUTDIR"
source /etc/profile.d/modules.sh 2>/dev/null || source /appl/Modules/5.3.1/init/bash

# 1) build a union VCF (ID = 19:POS:REF:ALT so VEP's Uploaded_variation = our key; skip '*' alleles)
{
  echo '##fileformat=VCFv4.2'
  echo '##contig=<ID=19>'
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
  awk -F, 'NR>1 && $5!="*" {print $2"\t"$3"\t"$2":"$3":"$4":"$5"\t"$4"\t"$5"\t.\t.\t."}' "$UNION_CSV" \
    | sort -t$'\t' -k2,2n
} > "$OUTDIR/znf175_union.vcf"
echo "[$(date)] union VCF variants: $(grep -vc '^#' "$OUTDIR/znf175_union.vcf")"

# 2) VEP — consequence + IMPACT (cache-only; this VEP can't index a bgzipped FASTA, so no --fasta).
#    --pick = one (most severe) consequence per variant.
module purge; module load ensemblvep/109.1
vep --offline --cache --dir_cache /project/ritchie/datasets/VEP --cache_version 109 \
    --assembly GRCh38 --force_overwrite --no_stats \
    --input_file "$OUTDIR/znf175_union.vcf" --format vcf --tab \
    --pick --symbol --canonical \
    --fields "Uploaded_variation,Location,Allele,Consequence,IMPACT,SYMBOL,Gene,BIOTYPE,CANONICAL" \
    --output_file "$OUTDIR/znf175_vep.tab"
echo "[$(date)] VEP done -> $OUTDIR/znf175_vep.tab"

# 3) ANNOVAR dbNSFP4.5a (canonical, uncompressed + .idx) -> AlphaMissense for missense
module purge; module load annovar/2025
table_annovar.pl "$OUTDIR/znf175_union.vcf" /project/ritchie/datasets/annovar_dbs/humandb \
    -buildver hg38 -out "$OUTDIR/znf175_annovar" \
    -protocol dbnsfp45a_canonical -operation f \
    -vcfinput -nastring . -polish
echo "[$(date)] ANNOVAR done -> $OUTDIR/znf175_annovar.hg38_multianno.txt"
