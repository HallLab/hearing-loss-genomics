#!/usr/bin/env bash
# submit_phase2.sh — submit Phase 2 (SNP ID reconciliation, light mode) to LSF.
#
# Light mode reuses Daniel's pre-extracted pVCF IDs, so this is fast:
# ~30 s decompress + ~1-2 min for the Python join + filter + validation.
#
# Resources: 1 CPU, 1 GB RAM, 15 min wall (generous).

set -euo pipefail

PROJECT_ROOT="/project/hall/analysis/hearing-loss-genomics"
cd "$PROJECT_ROOT"

LOG_DIR="analysis/daniel/logs/phase2"
mkdir -p "$LOG_DIR"

TS=$(date +%Y%m%d_%H%M%S)
LSF_OUT="$LOG_DIR/lsf_${TS}.out"
LSF_ERR="$LOG_DIR/lsf_${TS}.err"

JOB_OUT=$(
    bsub \
        -J hl_phase2 \
        -o "$LSF_OUT" \
        -e "$LSF_ERR" \
        -W 15 \
        -M 1024 \
        -R "rusage[mem=1024]" \
        -n 1 \
        "bash $PROJECT_ROOT/analysis/daniel/scripts/run_phase2.sh"
)

echo "$JOB_OUT"
echo ""
JOB_ID=$(echo "$JOB_OUT" | grep -oP 'Job <\K[0-9]+')
echo "Monitor:    bjobs $JOB_ID"
echo "LSF stdout: $LSF_OUT"
echo "LSF stderr: $LSF_ERR"
echo "Run log:    $LOG_DIR/run_*.log  (created by the job itself)"
