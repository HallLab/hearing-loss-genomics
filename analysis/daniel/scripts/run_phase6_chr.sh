#!/usr/bin/env bash
# run_phase6_chr.sh — single-chromosome biobin for Phase 6 (exome-wide all-genes burden).
# Called as one task of an LSF job array (chr from $LSB_JOBINDEX or $1).
#
# Daniel's per-chr filtered VCFs at data/PMBB_Exome/allGenes/ already have:
#   - keep-list applied (tokeep_moreHLcases, 41,748 individuals)
#   - --max-maf .01 filter applied (looser than Phase 5's .001)
#   - per-chr split done
# So we just decompress + run biobin (no --region-file → biobin uses all LOKI genes).

set -euo pipefail

CHR="${LSB_JOBINDEX:-${1:-}}"
[[ -n "$CHR" ]] || { echo "ERROR: CHR not provided (LSB_JOBINDEX or arg 1)"; exit 1; }

PROJECT_ROOT="/project/hall/analysis/hearing-loss-genomics"
cd "$PROJECT_ROOT"

if [[ -z "${VIRTUAL_ENV:-}" ]]; then
    # shellcheck disable=SC1091
    source venv/bin/activate
fi

# shellcheck disable=SC1091
source /etc/profile.d/modules.sh 2>/dev/null || true
module load rlsoftware/latest 2>/dev/null || true
export LD_LIBRARY_PATH="$PROJECT_ROOT/analysis/daniel/configs/lib-shims:${LD_LIBRARY_PATH:-}"

OUT_DIR="analysis/daniel/outputs/phase6"
VCF_DIR="$OUT_DIR/vcf"
BIOBIN_DIR="$OUT_DIR/biobin/HL_needAud"
LOG_DIR="analysis/daniel/logs/phase6/chr"

mkdir -p "$VCF_DIR" "$BIOBIN_DIR" "$LOG_DIR"
TS=$(date +%Y%m%d_%H%M%S)
LOG="$LOG_DIR/chr${CHR}_${TS}.log"

log() { printf '[phase6:chr%s] %s %s\n' "$CHR" "$(date +%H:%M:%S)" "$*"; }
fail() { log "ERROR: $*"; exit 1; }
exec > >(tee -a "$LOG") 2>&1

log "START Phase 6 chr${CHR}"

# Inputs
REF_VCF_GZ="data/PMBB_Exome/allGenes/allIndvs_burdenSNPs_allGenes_noRels_maf.01_chr${CHR}.vcf.gz"
CASES_CONTROL="analysis/daniel/outputs/phase4/cases_control.txt"
COVS="analysis/daniel/outputs/phase4/covs.txt"
LOKI_DB="/project/ritchie/datasets/loki/loki-20230816.db"

[[ -f "$REF_VCF_GZ" ]] || fail "Daniel's chr${CHR} VCF missing: $REF_VCF_GZ"
[[ -f "$CASES_CONTROL" ]] || fail "Phase 4 cases_control missing"
[[ -f "$COVS" ]] || fail "Phase 4 covs missing"
[[ -r "$LOKI_DB" ]] || fail "loki.db missing"

# Step 1: Decompress Daniel's per-chr VCF
LOCAL_VCF="$VCF_DIR/allIndvs_burdenSNPs_allGenes_noRels_maf.01_chr${CHR}.vcf"
log "Decompressing $REF_VCF_GZ → $LOCAL_VCF"
T0=$(date +%s)
zcat "$REF_VCF_GZ" > "$LOCAL_VCF"
T1=$(date +%s)
N_VARIANTS=$(grep -vc "^#" "$LOCAL_VCF" || true)
log "  decompressed: $(du -h $LOCAL_VCF | cut -f1), $N_VARIANTS variants, $((T1-T0))s"

# Step 2: Run biobin (NO --region-file → biobin uses LOKI genes for all the exome)
BIOBIN_PREFIX="$BIOBIN_DIR/HL_needAud_chr${CHR}"
log "Running biobin on chr${CHR}..."
log "  N variants: $N_VARIANTS"
log "  cases_control: $CASES_CONTROL"
log "  covariates: $COVS"
log "  loki: $LOKI_DB"
log "  no --region-file (LOKI gene boundaries)"

T0=$(date +%s)
biobin \
    -D "$LOKI_DB" \
    -V "$LOCAL_VCF" \
    -p "$CASES_CONTROL" \
    --covariates "$COVS" \
    --bin-regions Y \
    -G 38 \
    --test logistic \
    --report-prefix "$BIOBIN_PREFIX" \
    > "$BIOBIN_PREFIX.run_log.txt" 2>&1 \
    || fail "biobin chr${CHR} failed — see $BIOBIN_PREFIX.run_log.txt"
T1=$(date +%s)
ELAPSED=$((T1-T0))
log "biobin chr${CHR} done in ${ELAPSED}s"

BINS_CSV="$BIOBIN_PREFIX-bins.csv"
[[ -f "$BINS_CSV" ]] || fail "biobin did not produce bins.csv"
N_COLS=$(head -1 "$BINS_CSV" | awk -F',' '{print NF}')
N_ROWS=$(wc -l < "$BINS_CSV")
log "  bins.csv: $N_ROWS rows × $N_COLS cols ($(du -h $BINS_CSV | cut -f1))"

# Free disk: drop the decompressed VCF (Daniel's gz remains as backup)
rm -f "$LOCAL_VCF"
log "  removed local VCF (Daniel's gz remains in data/)"

log "DONE chr${CHR} — ${ELAPSED}s biobin wall"
