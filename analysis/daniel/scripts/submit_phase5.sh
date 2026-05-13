#!/usr/bin/env bash
# submit_phase5.sh — submit Phase 5 (heavy mode: IBD trick + plink filter + merge + biobin) to LSF.
#
# Resource estimates:
#   Step 5.0-5.1 (prep + IBD trick): ~10 s, <100 MB
#   Step 5.2 (22 plink filters): ~3-5 min sequential, <2 GB per chr
#   Step 5.3 (merge + VCF): ~1-3 min
#   Step 5.4 (biobin): ~5-30 min depending on cohort size
#   Total: ~15-45 min wall
#
# Request 32 GB / 4 cores / 90 min — generous; biobin is single-threaded by default
# (Daniel didn't use --threads), but extra cores cost nothing and 32 GB gives biobin
# headroom in case it builds large in-memory matrices for the regression.

set -euo pipefail

PROJECT_ROOT="/project/hall/analysis/hearing-loss-genomics"
cd "$PROJECT_ROOT"

LOG_DIR="analysis/daniel/logs/phase5"
mkdir -p "$LOG_DIR"

TS=$(date +%Y%m%d_%H%M%S)
LSF_OUT="$LOG_DIR/lsf_${TS}.out"
LSF_ERR="$LOG_DIR/lsf_${TS}.err"

JOB_OUT=$(
    bsub \
        -J hl_phase5 \
        -o "$LSF_OUT" \
        -e "$LSF_ERR" \
        -W 90 \
        -M 32768 \
        -R "rusage[mem=32768]" \
        -n 4 \
        "bash $PROJECT_ROOT/analysis/daniel/scripts/run_phase5.sh"
)

echo "$JOB_OUT"
echo ""
JOB_ID=$(echo "$JOB_OUT" | grep -oP 'Job <\K[0-9]+')
echo "Monitor:    bjobs $JOB_ID"
echo "LSF stdout: $LSF_OUT"
echo "LSF stderr: $LSF_ERR"
echo "Run log:    $LOG_DIR/run_*.log  (created by the job itself)"
