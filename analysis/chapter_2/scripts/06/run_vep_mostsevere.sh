#!/bin/bash
# VEP --most_severe on the union VCF -> authoritative most-severe consequence per variant
# (--pick can under-call pLOF by choosing a non-LoF transcript; --most_severe takes the worst across transcripts).
# Usage: bash run_vep_mostsevere.sh <UNION_VCF> <OUT_TAB>
set -euo pipefail
UNION_VCF="$1"; OUT="$2"
mkdir -p "$(dirname "$OUT")"
source /etc/profile.d/modules.sh 2>/dev/null || source /appl/Modules/5.3.1/init/bash
module purge; module load ensemblvep/109.1
vep --offline --cache --dir_cache /project/ritchie/datasets/VEP --cache_version 109 \
    --assembly GRCh38 --force_overwrite --no_stats \
    --input_file "$UNION_VCF" --format vcf --tab \
    --most_severe \
    --output_file "$OUT"
echo "DONE -> $OUT"
