#!/usr/bin/env bash
# submit_phase7.sh — submit Phase 7 (ZNF175 carrier deep-dive) to LSF.
#
# Phase 7 = ZNF175 carrier extraction + second-hit hypothesis test.
# Light compute: 22 small plink jobs (each ~5s) + Python data manipulation.
# Resources: 1 CPU, 4 GB RAM, 30 min wall.

set -euo pipefail

PROJECT_ROOT="/project/hall/analysis/hearing-loss-genomics"
cd "$PROJECT_ROOT"

LOG_DIR="analysis/daniel/logs/phase7"
mkdir -p "$LOG_DIR"

TS=$(date +%Y%m%d_%H%M%S)
LSF_OUT="$LOG_DIR/lsf_${TS}.out"
LSF_ERR="$LOG_DIR/lsf_${TS}.err"

JOB_OUT=$(
    bsub \
        -J hl_phase7 \
        -o "$LSF_OUT" \
        -e "$LSF_ERR" \
        -W 30 \
        -M 4096 \
        -R "rusage[mem=4096]" \
        -n 1 \
        "bash $PROJECT_ROOT/analysis/daniel/scripts/run_phase7.sh"
)

echo "$JOB_OUT"
echo ""
JOB_ID=$(echo "$JOB_OUT" | grep -oP 'Job <\K[0-9]+')
echo "Monitor:    bjobs $JOB_ID"
echo "LSF stdout: $LSF_OUT"
echo "LSF stderr: $LSF_ERR"
echo "Run log:    $LOG_DIR/run_*.log  (created by the job itself)"
