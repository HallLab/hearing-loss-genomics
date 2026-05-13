#!/usr/bin/env bash
# submit_phase3.sh — submit Phase 3 (plink genotype extraction, light + chr21 pilot) to LSF.
#
# Phase 3 does:
#   - LIGHT MODE for all 22 chrs (decompress Daniel's per-chr files): ~1 min
#   - HEAVY PILOT on chr21 only (plink 1.9 + plink 2.0): ~10-30 min for both
#
# Total estimated wall: ~30-40 min. Request 60 min for safety.
# plink can use ~2-4 GB RAM when reading VCF; request 4 GB.

set -euo pipefail

PROJECT_ROOT="/project/hall/analysis/hearing-loss-genomics"
cd "$PROJECT_ROOT"

LOG_DIR="analysis/daniel/logs/phase3"
mkdir -p "$LOG_DIR"

TS=$(date +%Y%m%d_%H%M%S)
LSF_OUT="$LOG_DIR/lsf_${TS}.out"
LSF_ERR="$LOG_DIR/lsf_${TS}.err"

JOB_OUT=$(
    bsub \
        -J hl_phase3 \
        -o "$LSF_OUT" \
        -e "$LSF_ERR" \
        -W 60 \
        -M 4096 \
        -R "rusage[mem=4096]" \
        -n 1 \
        "bash $PROJECT_ROOT/analysis/daniel/scripts/run_phase3.sh"
)

echo "$JOB_OUT"
echo ""
JOB_ID=$(echo "$JOB_OUT" | grep -oP 'Job <\K[0-9]+')
echo "Monitor:    bjobs $JOB_ID"
echo "LSF stdout: $LSF_OUT"
echo "LSF stderr: $LSF_ERR"
echo "Run log:    $LOG_DIR/run_*.log  (created by the job itself)"
