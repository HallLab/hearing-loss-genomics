#!/usr/bin/env bash
# submit_phase1.sh — submit Phase 1 of Hui et al. 2023 replication to LSF (bsub).
#
# Wraps analysis/daniel/scripts/run_phase1.sh as an LSF job.
#
# Resource request:
#   - 1 CPU, 4 GB RAM, 30 min wall (very generous — Phase 1 takes ~2-4 min in practice)
#
# Outputs:
#   - LSF stdout/stderr per submission: analysis/daniel/logs/phase1/lsf_<timestamp>.{out,err}
#   - Pipeline log (inside the job): analysis/daniel/logs/phase1/run_<timestamp>.log
#   - Phase 1 outputs:                 analysis/daniel/outputs/phase1/

set -euo pipefail

PROJECT_ROOT="/project/hall/analysis/hearing-loss-genomics"
cd "$PROJECT_ROOT"

LOG_DIR="analysis/daniel/logs/phase1"
mkdir -p "$LOG_DIR"

TS=$(date +%Y%m%d_%H%M%S)
LSF_OUT="$LOG_DIR/lsf_${TS}.out"
LSF_ERR="$LOG_DIR/lsf_${TS}.err"

JOB_OUT=$(
    bsub \
        -J hl_phase1 \
        -o "$LSF_OUT" \
        -e "$LSF_ERR" \
        -W 30 \
        -M 4000 \
        -R "rusage[mem=4000]" \
        -n 1 \
        "bash $PROJECT_ROOT/analysis/daniel/scripts/run_phase1.sh"
)

echo "$JOB_OUT"
echo ""
echo "Monitor:    bjobs -l \$(echo \"$JOB_OUT\" | grep -oP 'Job <\\K[0-9]+')"
echo "LSF stdout: $LSF_OUT"
echo "LSF stderr: $LSF_ERR"
echo "Run log:    $LOG_DIR/run_*.log  (created by the job itself)"
