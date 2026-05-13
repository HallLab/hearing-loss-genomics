#!/usr/bin/env bash
# submit_phase4.sh — submit Phase 4 (preparatory files) to LSF.
#
# Three Python scripts, all light: ~1-2 min total wall.
# Resources: 1 CPU, 256 MB RAM, 5 min wall.

set -euo pipefail

PROJECT_ROOT="/project/hall/analysis/hearing-loss-genomics"
cd "$PROJECT_ROOT"

LOG_DIR="analysis/daniel/logs/phase4"
mkdir -p "$LOG_DIR"

TS=$(date +%Y%m%d_%H%M%S)
LSF_OUT="$LOG_DIR/lsf_${TS}.out"
LSF_ERR="$LOG_DIR/lsf_${TS}.err"

JOB_OUT=$(
    bsub \
        -J hl_phase4 \
        -o "$LSF_OUT" \
        -e "$LSF_ERR" \
        -W 5 \
        -M 256 \
        -R "rusage[mem=256]" \
        -n 1 \
        "bash $PROJECT_ROOT/analysis/daniel/scripts/run_phase4.sh"
)

echo "$JOB_OUT"
echo ""
JOB_ID=$(echo "$JOB_OUT" | grep -oP 'Job <\K[0-9]+')
echo "Monitor:    bjobs $JOB_ID"
echo "LSF stdout: $LSF_OUT"
echo "LSF stderr: $LSF_ERR"
echo "Run log:    $LOG_DIR/run_*.log  (created by the job itself)"
