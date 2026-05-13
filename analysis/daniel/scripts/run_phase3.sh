#!/usr/bin/env bash
# run_phase3.sh — Phase 3 of Hui et al. 2023 replication: plink genotype extraction.
#
# Phase 3 takes the .extract from Phase 2 and uses plink to extract per-chr
# genotype matrices (bed/bim/fam) from the v2 pVCFs.
#
# This script does two things:
#  (a) LIGHT MODE for all 22 chrs: decompresses Daniel's pre-built per-chr
#      plink files (allIndvs_chr*.bed/bim/fam.gz) into our outputs. Validates
#      they're well-formed and consistent with our Phase 2 .extract.
#  (b) HEAVY PILOT on chr21 only: re-runs plink ourselves with BOTH plink 1.9
#      and plink 2.0 on the smallest autosomal pVCF (8.7 GB), then compares
#      outputs against Daniel's reference. This is the "can plink run here?"
#      smoke test AND the plink-version comparison.
#
# Output layout:
#   analysis/daniel/outputs/phase3/
#   ├── light/                     — Daniel's 22 per-chr files decompressed
#   ├── pilot_chr21/
#   │   ├── v19/                   — our plink 1.9 chr21 output
#   │   └── v20/                   — our plink 2.0 chr21 output
#
# Validation:
#   Light: bim row counts per chr, sum = our Phase 2 .extract count (9667)
#   Pilot: diff bim and fam files between v1.9 / v2.0 / Daniel reference;
#          md5sum the bed files.

set -euo pipefail

# ───────── Paths ─────────
PROJECT_ROOT="/project/hall/analysis/hearing-loss-genomics"
cd "$PROJECT_ROOT"

if [[ -z "${VIRTUAL_ENV:-}" ]]; then
    # shellcheck disable=SC1091
    source venv/bin/activate
fi

OUT_DIR="analysis/daniel/outputs/phase3"
LIGHT_DIR="$OUT_DIR/light"
PILOT_DIR="$OUT_DIR/pilot_chr21"
V19_DIR="$PILOT_DIR/v19"
V20_DIR="$PILOT_DIR/v20"
LOG_DIR="analysis/daniel/logs/phase3"

# Phase 2 extract list (input)
EXTRACT_LIST="analysis/daniel/outputs/phase2/matched_snp_IDs_annot_pVCF_noNA_noMultiallelic.extract"

# pVCFs (raw v2)
PVCF_DIR="data/pmbb_v2/Exome/pVCF/GL_by_chrom"
CHR21_PVCF="$PVCF_DIR/PMBB-Release-2020-2.0_genetic_exome_chr21_GL.vcf.gz"

# Daniel's reference (per-chr extracted bed/bim/fam — Step 3.1 output)
REF_GENO_DIR="data/PMBB_Exome/genotypes"

# plink binaries
PLINK_19="/appl/plink-1.90Beta6.18/plink"
PLINK_20="/appl/plink2-20240804/plink"

mkdir -p "$LIGHT_DIR" "$V19_DIR" "$V20_DIR" "$LOG_DIR"
TS=$(date +%Y%m%d_%H%M%S)
LOG="$LOG_DIR/run_${TS}.log"

# ───────── Logging ─────────
log() { printf '[phase3] %s %s\n' "$(date +%H:%M:%S)" "$*"; }
fail() { log "ERROR: $*"; exit 1; }
exec > >(tee -a "$LOG") 2>&1

log "START Phase 3 — plink genotype extraction (LIGHT + chr21 pilot HEAVY)"
log "Project root:   $PROJECT_ROOT"
log "Output dir:     $OUT_DIR"
log "Log file:       $LOG"

# ───────── Verify inputs ─────────
log ""
log "Verifying inputs..."
[[ -f "$EXTRACT_LIST" ]] || fail "Phase 2 .extract missing: $EXTRACT_LIST"
[[ -d "$PVCF_DIR" ]] || fail "v2 pVCF dir missing: $PVCF_DIR"
[[ -f "$CHR21_PVCF" ]] || fail "chr21 pVCF missing: $CHR21_PVCF"
[[ -x "$PLINK_19" ]] || fail "plink 1.9 binary missing: $PLINK_19"
[[ -x "$PLINK_20" ]] || fail "plink 2.0 binary missing: $PLINK_20"
[[ -d "$REF_GENO_DIR" ]] || fail "Daniel's genotypes dir missing: $REF_GENO_DIR"

N_EXTRACT=$(wc -l < "$EXTRACT_LIST")
log "  Phase 2 .extract:      $EXTRACT_LIST ($N_EXTRACT IDs)"
log "  v2 pVCF dir:           $PVCF_DIR"
log "  chr21 pVCF (pilot):    $CHR21_PVCF ($(du -h $CHR21_PVCF | cut -f1))"
log "  plink 1.9:             $($PLINK_19 --version 2>&1 | head -1)"
log "  plink 2.0:             $($PLINK_20 --version 2>&1 | head -1)"

# ═══════════════════════════════════════════════════════════════════════
# PART A — LIGHT MODE: decompress Daniel's per-chr files
# ═══════════════════════════════════════════════════════════════════════
log ""
log "════════════════════════════════════════════════════════════"
log "PART A — LIGHT MODE: decompressing Daniel's 22 per-chr files"
log "════════════════════════════════════════════════════════════"
T0=$(date +%s)
for i in {1..22}; do
    for ext in bed bim fam; do
        src="$REF_GENO_DIR/allIndvs_chr${i}.${ext}.gz"
        dst="$LIGHT_DIR/allIndvs_chr${i}.${ext}"
        [[ -f "$src" ]] || fail "Daniel's reference missing: $src"
        zcat "$src" > "$dst"
    done
done
T1=$(date +%s)
log "  decompressed 66 files (22 chrs × 3 exts) in $((T1-T0))s"

# ───────── Light validation ─────────
log ""
log "Light validation — per-chr counts + variant subset check"

# (a) Per-chr variant counts
log ""
log "  Per-chr variant counts (.bim row count):"
TOTAL_BIM_ROWS=0
for i in {1..22}; do
    n=$(wc -l < "$LIGHT_DIR/allIndvs_chr${i}.bim")
    TOTAL_BIM_ROWS=$((TOTAL_BIM_ROWS + n))
    printf '[phase3] %s     chr%-2d: %s variants\n' "$(date +%H:%M:%S)" "$i" "$n"
done
log "  TOTAL across 22 chrs: $TOTAL_BIM_ROWS variants"
log "  Phase 2 .extract had: $N_EXTRACT IDs"

if [[ "$TOTAL_BIM_ROWS" -eq "$N_EXTRACT" ]]; then
    log "  ✓ TOTAL bim rows == Phase 2 extract count"
else
    log "  ⚠ mismatch: $((TOTAL_BIM_ROWS - N_EXTRACT)) variants discrepancy"
fi

# (b) Sample count consistency across chrs (.fam row count)
log ""
log "  Sample counts (.fam row count) — should be identical across all 22 chrs:"
SAMPLE_COUNTS=$(for i in {1..22}; do wc -l < "$LIGHT_DIR/allIndvs_chr${i}.fam"; done | sort -u)
N_DISTINCT=$(echo "$SAMPLE_COUNTS" | wc -l)
if [[ "$N_DISTINCT" -eq 1 ]]; then
    log "  ✓ all 22 chrs have the same sample count: $SAMPLE_COUNTS"
else
    log "  ⚠ inconsistent sample counts across chrs:"
    echo "$SAMPLE_COUNTS" | head -10
fi

# (c) Variant subset check: every .bim variant ID should be in our Phase 2 .extract
log ""
log "  Variant subset check: do .bim variants come from Phase 2 .extract?"
ALL_BIM_IDS=$(cat $LIGHT_DIR/allIndvs_chr*.bim | cut -f2 | sort -u)
EXTRACT_IDS=$(sort -u "$EXTRACT_LIST")
N_BIM_NOT_IN_EXTRACT=$(comm -23 <(echo "$ALL_BIM_IDS") <(echo "$EXTRACT_IDS") | wc -l)
if [[ "$N_BIM_NOT_IN_EXTRACT" -eq 0 ]]; then
    log "  ✓ all bim variant IDs are in our Phase 2 .extract"
else
    log "  ⚠ $N_BIM_NOT_IN_EXTRACT bim IDs are NOT in our .extract — investigate"
fi

# ═══════════════════════════════════════════════════════════════════════
# PART B — HEAVY PILOT on chr21 (plink 1.9 vs plink 2.0)
# ═══════════════════════════════════════════════════════════════════════
log ""
log "════════════════════════════════════════════════════════════"
log "PART B — HEAVY PILOT on chr21 (plink 1.9 vs plink 2.0)"
log "════════════════════════════════════════════════════════════"
log "  Replicating Daniel's Step 3.1 command on chr21 with both plink versions."
log "  Daniel's command:"
log "    plink --vcf <chr21_pVCF> --vcf-half-call m --extract <our_extract> \\"
log "          --make-bed --out genotypes/allIndvs_chr21"

# ───────── B.1 — plink 1.9 ─────────
log ""
log "B.1 — plink 1.9 on chr21..."
T0=$(date +%s)
"$PLINK_19" \
    --vcf "$CHR21_PVCF" \
    --vcf-half-call m \
    --extract "$EXTRACT_LIST" \
    --make-bed \
    --out "$V19_DIR/allIndvs_chr21" \
    > "$V19_DIR/plink_stdout.txt" 2> "$V19_DIR/plink_stderr.txt" \
    || fail "plink 1.9 failed — see $V19_DIR/plink_stderr.txt"
T1=$(date +%s)
N19=$(wc -l < "$V19_DIR/allIndvs_chr21.bim")
S19=$(wc -l < "$V19_DIR/allIndvs_chr21.fam")
log "  plink 1.9 done in $((T1-T0))s: $N19 variants, $S19 samples"

# ───────── B.2 — plink 2.0 ─────────
log ""
log "B.2 — plink 2.0 on chr21..."
T0=$(date +%s)
"$PLINK_20" \
    --vcf "$CHR21_PVCF" \
    --vcf-half-call m \
    --extract "$EXTRACT_LIST" \
    --make-bed \
    --out "$V20_DIR/allIndvs_chr21" \
    > "$V20_DIR/plink_stdout.txt" 2> "$V20_DIR/plink_stderr.txt" \
    || fail "plink 2.0 failed — see $V20_DIR/plink_stderr.txt"
T1=$(date +%s)
N20=$(wc -l < "$V20_DIR/allIndvs_chr21.bim")
S20=$(wc -l < "$V20_DIR/allIndvs_chr21.fam")
log "  plink 2.0 done in $((T1-T0))s: $N20 variants, $S20 samples"

# ───────── B.3 — Pilot validation ─────────
log ""
log "B.3 — pilot validation: comparing v1.9, v2.0, and Daniel's reference"

# Disable set -e in validation block (diff returns non-zero on mismatch — expected here)
set +e

DANIEL_BIM=$(zcat "$REF_GENO_DIR/allIndvs_chr21.bim.gz" | sort)
V19_BIM=$(sort "$V19_DIR/allIndvs_chr21.bim")
V20_BIM=$(sort "$V20_DIR/allIndvs_chr21.bim")

# bim comparison
log ""
log "  .bim (variant list, sorted):"
log "    Daniel:    $(zcat $REF_GENO_DIR/allIndvs_chr21.bim.gz | wc -l) rows"
log "    v1.9:      $N19 rows"
log "    v2.0:      $N20 rows"

d_v19=$( { diff <(echo "$DANIEL_BIM") <(echo "$V19_BIM") || true; } | wc -l)
d_v20=$( { diff <(echo "$DANIEL_BIM") <(echo "$V20_BIM") || true; } | wc -l)
d_19_20=$( { diff <(echo "$V19_BIM") <(echo "$V20_BIM") || true; } | wc -l)
log "    diff lines Daniel vs v1.9: $d_v19"
log "    diff lines Daniel vs v2.0: $d_v20"
log "    diff lines v1.9 vs v2.0:   $d_19_20"

# fam comparison (just check sample sets are equal)
log ""
log "  .fam (sample list):"
DANIEL_FAM_IDS=$(zcat "$REF_GENO_DIR/allIndvs_chr21.fam.gz" | cut -f2 | sort)
V19_FAM_IDS=$(cut -f2 "$V19_DIR/allIndvs_chr21.fam" | sort)
V20_FAM_IDS=$(cut -f2 "$V20_DIR/allIndvs_chr21.fam" | sort)
log "    Daniel:    $(echo "$DANIEL_FAM_IDS" | wc -l) samples"
log "    v1.9:      $(echo "$V19_FAM_IDS" | wc -l) samples"
log "    v2.0:      $(echo "$V20_FAM_IDS" | wc -l) samples"
d_fam_v19=$( { diff <(echo "$DANIEL_FAM_IDS") <(echo "$V19_FAM_IDS") || true; } | wc -l)
d_fam_v20=$( { diff <(echo "$DANIEL_FAM_IDS") <(echo "$V20_FAM_IDS") || true; } | wc -l)
log "    diff lines Daniel vs v1.9 fam: $d_fam_v19"
log "    diff lines Daniel vs v2.0 fam: $d_fam_v20"

# bed md5sum
log ""
log "  .bed (binary genotype matrix, md5sum):"
DANIEL_BED_MD5=$(zcat "$REF_GENO_DIR/allIndvs_chr21.bed.gz" | md5sum | cut -d' ' -f1)
V19_BED_MD5=$(md5sum "$V19_DIR/allIndvs_chr21.bed" | cut -d' ' -f1)
V20_BED_MD5=$(md5sum "$V20_DIR/allIndvs_chr21.bed" | cut -d' ' -f1)
log "    Daniel:    $DANIEL_BED_MD5"
log "    v1.9:      $V19_BED_MD5"
log "    v2.0:      $V20_BED_MD5"

set -e

# ───────── Recommendation ─────────
log ""
log "════════════════════════════════════════════════════════════"
log "Verdict"
log "════════════════════════════════════════════════════════════"
log ""
RECOMMEND=""
if [[ "$DANIEL_BED_MD5" == "$V19_BED_MD5" && "$d_v19" -eq 0 && "$d_fam_v19" -eq 0 ]]; then
    log "  ✓ plink 1.9 reproduces Daniel's chr21 byte-for-byte"
    RECOMMEND="plink 1.9"
fi
if [[ "$DANIEL_BED_MD5" == "$V20_BED_MD5" && "$d_v20" -eq 0 && "$d_fam_v20" -eq 0 ]]; then
    log "  ✓ plink 2.0 reproduces Daniel's chr21 byte-for-byte"
    if [[ -z "$RECOMMEND" ]]; then RECOMMEND="plink 2.0"; fi
fi
if [[ "$V19_BED_MD5" == "$V20_BED_MD5" && "$d_19_20" -eq 0 ]]; then
    log "  ✓ plink 1.9 and plink 2.0 produce identical output (interchangeable)"
fi

if [[ -z "$RECOMMEND" ]]; then
    log "  ⚠ Neither plink version matches Daniel byte-for-byte."
    log "    .bim diffs (Daniel vs v1.9/v2.0): $d_v19 / $d_v20"
    log "    .bed md5sum matches: none"
    log "    Use semantic comparison (variant ID sets) to decide if differences are cosmetic."
    log "    See pilot outputs in $PILOT_DIR for forensics."
else
    log "  → RECOMMENDED for Phase 3+ heavy work: $RECOMMEND"
fi

log ""
log "DONE — Phase 3 (light mode + chr21 pilot) complete"
log "  Light outputs (all 22 chrs, ready for Phase 4): $LIGHT_DIR"
log "  Pilot outputs (chr21 plink runs):               $PILOT_DIR"
