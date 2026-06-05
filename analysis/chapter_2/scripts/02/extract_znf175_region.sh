#!/bin/bash
# Extract the ZNF175 locus from a PMBB exome pVCF and compute allele freq + genotype counts.
# Independent of Daniel/Park: input = raw release pVCF; window = NCBI RefSeq coords; tools = bcftools + plink2.
#
# Usage:
#   bash extract_znf175_region.sh <VCF.gz> <OUT_PREFIX> [REGION]
# Args:
#   VCF.gz      : bgzipped, tabix-indexed pVCF (uses .tbi for fast random access)
#   OUT_PREFIX  : output path prefix (writes <prefix>_region.vcf.gz, <prefix>.afreq, <prefix>.gcount)
#   REGION      : default 19:51571283-51592510  (ZNF175, NCBI NC_000019.10, GRCh38, + strand)
set -euo pipefail

VCF="$1"
OUTPREFIX="$2"
REGION="${3:-19:51571283-51592510}"
# GRCh38 reference (Ensembl naming '19', matches the pVCF contig); needed for left-align + atomize
REF="${4:-/project/ritchie/datasets/VEP/vep115_cache/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz}"

source /etc/profile.d/modules.sh 2>/dev/null || source /appl/Modules/5.3.1/init/bash
module load htslib/1.21 bcftools/1.21
PLINK2=/appl/plink2-20240804/plink2

mkdir -p "$(dirname "$OUTPREFIX")"

echo "[$(date)] region extract: $REGION  from  $VCF"
# 1) region pull via .tbi random access; then CANONICAL normalization:
#    -a       : atomize MNVs/complex into atomic SNVs/indels (e.g. ATT:ACT -> T:C)
#    -m-any   : split multiallelics so each ALT is its own record
#    -f REF   : left-align indels + check REF against the GRCh38 reference
#    --set-id : unique chr:pos:ref:alt IDs (split records would otherwise share a compound ';' ID)
bcftools view -r "$REGION" "$VCF" -Ou \
  | bcftools norm -f "$REF" -a -m-any -Ou \
  | bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' -Oz -o "${OUTPREFIX}_region.vcf.gz"
tabix -f -p vcf "${OUTPREFIX}_region.vcf.gz"

NV=$(bcftools view -H "${OUTPREFIX}_region.vcf.gz" | wc -l)
echo "[$(date)] variants in region after split: $NV"

# 2) allele frequency + genotype counts (same recipe as the v1 extraction)
echo "[$(date)] plink2 freq + geno-counts"
$PLINK2 --vcf "${OUTPREFIX}_region.vcf.gz" \
        --vcf-half-call m \
        --freq --geno-counts \
        --memory 6000 --threads 4 \
        --out "$OUTPREFIX"

echo "[$(date)] DONE -> ${OUTPREFIX}.afreq / ${OUTPREFIX}.gcount"
