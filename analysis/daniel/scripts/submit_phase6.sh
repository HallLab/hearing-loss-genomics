#!/usr/bin/env bash
# submit_phase6.sh — submit Phase 6 (exome-wide all-genes burden test) to LSF.
#
# Architecture:
#   Step 1: LSF job array of 22 tasks (one per chr) — biobin per chromosome.
#   Step 2: dependent finalize job — concatenate + BH-correct + report top hits.
#
# Resources per chr biobin task:
#   - 8 GB RAM (Phase 5 used 177 MB for 967 bins; all-genes is ~5-10× more bins)
#   - 6 hours wall (Phase 5 was 46 min for 967 bins; all-genes chr could be 2-4h)
#   - 1 core (biobin is single-threaded)

set -euo pipefail

PROJECT_ROOT="/project/hall/analysis/hearing-loss-genomics"
cd "$PROJECT_ROOT"

LOG_DIR="analysis/daniel/logs/phase6"
mkdir -p "$LOG_DIR"

TS=$(date +%Y%m%d_%H%M%S)
ARRAY_OUT="$LOG_DIR/lsf_array_${TS}_chr%I.out"
ARRAY_ERR="$LOG_DIR/lsf_array_${TS}_chr%I.err"
FINAL_OUT="$LOG_DIR/lsf_finalize_${TS}.out"
FINAL_ERR="$LOG_DIR/lsf_finalize_${TS}.err"

# Step 1: Submit job array (22 tasks, one per chromosome)
echo "Submitting biobin job array (22 chromosomes)..."
ARRAY_OUT_LINE=$(
    bsub \
        -J "hl_phase6[1-22]" \
        -o "$ARRAY_OUT" \
        -e "$ARRAY_ERR" \
        -W 360 \
        -M 8192 \
        -R "rusage[mem=8192]" \
        -n 1 \
        "bash $PROJECT_ROOT/analysis/daniel/scripts/run_phase6_chr.sh"
)
echo "$ARRAY_OUT_LINE"
ARRAY_JOB_NAME="hl_phase6"

# Step 2: Submit dependent finalize (waits for ALL 22 array tasks to complete)
echo ""
echo "Submitting finalize job (depends on array completion)..."
FINAL_OUT_LINE=$(
    bsub \
        -J "hl_phase6_final" \
        -o "$FINAL_OUT" \
        -e "$FINAL_ERR" \
        -W 60 \
        -M 4096 \
        -R "rusage[mem=4096]" \
        -n 1 \
        -w "ended(${ARRAY_JOB_NAME})" \
        "bash $PROJECT_ROOT/analysis/daniel/scripts/run_phase6_finalize.sh"
)
echo "$FINAL_OUT_LINE"

echo ""
ARRAY_ID=$(echo "$ARRAY_OUT_LINE" | grep -oP 'Job <\K[0-9]+')
FINAL_ID=$(echo "$FINAL_OUT_LINE" | grep -oP 'Job <\K[0-9]+')
echo "Array job: $ARRAY_ID (22 tasks)"
echo "Final job: $FINAL_ID (waits for ended($ARRAY_JOB_NAME))"
echo ""
echo "Monitor: bjobs $ARRAY_ID $FINAL_ID"
echo "Per-chr LSF logs: $LOG_DIR/lsf_array_${TS}_chr*.out"
echo "Per-chr run logs: analysis/daniel/logs/phase6/chr/chr*_*.log"
echo "Final log:        $LOG_DIR/finalize_*.log (created by finalize job)"
