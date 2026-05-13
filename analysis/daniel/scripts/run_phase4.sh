#!/usr/bin/env bash
# run_phase4.sh — Phase 4 of Hui et al. 2023 replication: preparatory files for biobin.
#
# Builds the 3 input files biobin will need (in Phase 5+):
#   - cases_control.txt          via case_control.py     (audbase + phecode 389)
#   - covs.txt                   via make_covs.py        (audbase + PMBB phenotype covariates)
#   - gene_list_regions.txt      via make_region_file.py (Phase 1 annot_genes_full_funcToInclude.txt)
#
# Step 4.4 (the actual biobin run) is deliberately NOT in this script — that's
# the territory of Phase 5 (combined with Phase 6's IBD trick) since Daniel
# overwrote his Phase 4 biobin outputs with the Phase 7 keepHLcases version.
#
# Outputs: analysis/daniel/outputs/phase4/
# Log:     analysis/daniel/logs/phase4/run_<timestamp>.log
#
# Validation: semantic set-equality of each of the 3 outputs against Daniel's
# pre-built references in data/PMBB_Exome/.

set -euo pipefail

# ───────── Paths ─────────
PROJECT_ROOT="/project/hall/analysis/hearing-loss-genomics"
cd "$PROJECT_ROOT"

if [[ -z "${VIRTUAL_ENV:-}" ]]; then
    # shellcheck disable=SC1091
    source venv/bin/activate
fi

OUT_DIR="analysis/daniel/outputs/phase4"
LOG_DIR="analysis/daniel/logs/phase4"
TMP_DIR="$OUT_DIR/.tmp"   # for decompressed inputs
SCRIPTS_DIR="analysis/daniel/scripts/pmbb_exome"

# Inputs from previous phases / raw v2
PHASE1_OUTPUT="analysis/daniel/outputs/phase1/annot_genes_full_funcToInclude.txt"
AUDBASE_GZ="data/PMBB_Exome/audbase_feb252021/RGC21_45k_aud_1.csv.gz"
PHECODE_GZ="data/PMBB_Exome/phecode_hl.txt.gz"
PHENOTYPE_COVS="data/pmbb_v2/Phenotype/2.0/PMBB-Release-2020-2.0_phenotype_covariates.txt"

# Daniel reference outputs (for validation)
REF_CASES_CONTROL="data/PMBB_Exome/cases_control.txt.gz"
REF_COVS="data/PMBB_Exome/covs.txt.gz"
REF_REGIONS="data/PMBB_Exome/gene_list_regions.txt.gz"

# Our outputs
CASES_CONTROL="$OUT_DIR/cases_control.txt"
COVS="$OUT_DIR/covs.txt"
REGIONS="$OUT_DIR/gene_list_regions.txt"

mkdir -p "$OUT_DIR" "$TMP_DIR" "$LOG_DIR"
TS=$(date +%Y%m%d_%H%M%S)
LOG="$LOG_DIR/run_${TS}.log"

# ───────── Logging ─────────
log() { printf '[phase4] %s %s\n' "$(date +%H:%M:%S)" "$*"; }
fail() { log "ERROR: $*"; exit 1; }
exec > >(tee -a "$LOG") 2>&1

log "START Phase 4 — preparatory files (case/control + covariates + region file)"
log "Project root:   $PROJECT_ROOT"
log "Output dir:     $OUT_DIR"
log "Log file:       $LOG"
log ""
log "Phase 4 scope: just the 3 preparatory files. Step 4.4 (biobin) deferred"
log "  to Phase 5 (which will fuse Phase 6's IBD trick + Phase 7's biobin run)."

# ───────── Verify inputs ─────────
log ""
log "Verifying inputs..."
[[ -f "$PHASE1_OUTPUT" ]] || fail "Phase 1 master list missing: $PHASE1_OUTPUT"
[[ -f "$AUDBASE_GZ" ]] || fail "Audbase missing: $AUDBASE_GZ"
[[ -f "$PHECODE_GZ" ]] || fail "phecode_hl.txt missing: $PHECODE_GZ"
[[ -f "$PHENOTYPE_COVS" ]] || fail "PMBB phenotype_covariates missing: $PHENOTYPE_COVS"
log "  Phase 1 master:       $PHASE1_OUTPUT ($(wc -l < $PHASE1_OUTPUT) rows)"
log "  audbase:              $AUDBASE_GZ ($(du -h $AUDBASE_GZ | cut -f1))"
log "  phecode_hl:           $PHECODE_GZ ($(du -h $PHECODE_GZ | cut -f1))"
log "  phenotype_covs (v2):  $PHENOTYPE_COVS ($(wc -l < $PHENOTYPE_COVS) rows)"

# ───────── Decompress inputs ─────────
log ""
log "Decompressing audbase + phecode_hl into $TMP_DIR..."
zcat "$AUDBASE_GZ"  > "$TMP_DIR/RGC21_45k_aud_1.csv"
zcat "$PHECODE_GZ"  > "$TMP_DIR/phecode_hl.txt"
AUDBASE_DECOMP="$TMP_DIR/RGC21_45k_aud_1.csv"
PHECODE_DECOMP="$TMP_DIR/phecode_hl.txt"
log "  audbase:    $(wc -l < $AUDBASE_DECOMP) rows"
log "  phecode_hl: $(wc -l < $PHECODE_DECOMP) rows"

# ═══════════════════════════════════════════════════════════════════════
# Step 4.1 — Case/control definition
# ═══════════════════════════════════════════════════════════════════════
log ""
log "Step 4.1 — case_control.py (audbase + phecode → cases_control.txt)"
log "  Logic: case if audiogram BL_SNHL=1; control if BL_SNHL=0 OR phecode=FALSE; NA otherwise."
T0=$(date +%s)
python "$SCRIPTS_DIR/case_control.py" "$AUDBASE_DECOMP" "$PHECODE_DECOMP" > "$CASES_CONTROL"
T1=$(date +%s)
N=$(wc -l < "$CASES_CONTROL")
log "  written: $CASES_CONTROL ($N rows, $((T1-T0))s)"
N_CASES=$(awk -F'\t' 'NR>1 && $2 == "1"' "$CASES_CONTROL" | wc -l)
N_CTRLS=$(awk -F'\t' 'NR>1 && $2 == "0"' "$CASES_CONTROL" | wc -l)
N_NA=$(awk -F'\t' 'NR>1 && $2 == "NA"' "$CASES_CONTROL" | wc -l)
log "  breakdown: cases=$N_CASES  controls=$N_CTRLS  NA=$N_NA  total_data_rows=$((N - 1))"

# ═══════════════════════════════════════════════════════════════════════
# Step 4.2 — Covariates
# ═══════════════════════════════════════════════════════════════════════
log ""
log "Step 4.2 — make_covs.py (audbase + phenotype_covariates → covs.txt)"
log "  Logic: age = audiogram_date - birth_year if audbase, else Age_at_Enrollment from PMBB."
log "         AgeSq = Age*Age. Sex coded 1=Male, 0=Female, NA otherwise. PC1-PC4."
T0=$(date +%s)
python "$SCRIPTS_DIR/make_covs.py" "$AUDBASE_DECOMP" "$PHENOTYPE_COVS" > "$COVS"
T1=$(date +%s)
N=$(wc -l < "$COVS")
log "  written: $COVS ($N rows, $((T1-T0))s)"

# ═══════════════════════════════════════════════════════════════════════
# Step 4.3 — Region file
# ═══════════════════════════════════════════════════════════════════════
log ""
log "Step 4.3 — make_region_file.py (Phase 1 master → gene_list_regions.txt)"
log "  Logic: for each gene, output chr + gene + min(start) + max(stop)."
log "  Daniel applies 'sort -gk1,1 -gk3,3' after the script (numeric sort by chr,start)."
T0=$(date +%s)
python "$SCRIPTS_DIR/make_region_file.py" "$PHASE1_OUTPUT" | sort -gk1,1 -gk3,3 > "$REGIONS"
T1=$(date +%s)
N=$(wc -l < "$REGIONS")
log "  written: $REGIONS ($N rows, $((T1-T0))s)"

# ═══════════════════════════════════════════════════════════════════════
# Validation — semantic set-equality
# ═══════════════════════════════════════════════════════════════════════
log ""
log "Validation — semantic equality against Daniel's references"

set +e   # tolerate non-zero diff exits

validate() {
    local ref_gz="$1" ours="$2" label="$3"
    local n_ref n_ours diff_count
    n_ref=$(zcat "$ref_gz" | wc -l)
    n_ours=$(wc -l < "$ours")
    diff_count=$({ diff <(zcat "$ref_gz" | sort) <(sort "$ours") || true; } | wc -l)

    if [[ "$diff_count" -eq 0 ]]; then
        log "  ✓ $label  (Daniel=$n_ref rows, ours=$n_ours rows, sorted-set ≡)"
        return 0
    else
        log "  ✗ $label  (Daniel=$n_ref, ours=$n_ours, $diff_count diff lines)"
        local n_ours_only n_ref_only
        n_ours_only=$(comm -13 <(zcat "$ref_gz" | sort) <(sort "$ours") | wc -l)
        n_ref_only=$(comm -23 <(zcat "$ref_gz" | sort) <(sort "$ours") | wc -l)
        log "      $n_ours_only rows only in ours, $n_ref_only rows only in Daniel"
        return 1
    fi
}

# Custom validator for covs.txt — rounds AgeSq (col 4) to 4 decimals to absorb
# Python 2 vs Python 3 float repr differences. Daniel's covs were generated under
# Python 2 (~12-char str(float), lossy); ours under Python 3 (shortest-roundtrip repr).
# All 43,722 AgeSq values are mathematically identical as floats (max abs diff = 0.0
# in earlier Python check); only textual repr differs. 4 decimals absorbs the lossy
# Python 2 string → float → round-half-even mismatches at the 6th decimal.
# AgeSq values are ~1000-7000, so 4 decimals = absolute precision 1e-4, relative ~2e-8.
validate_covs() {
    local ref_gz="$1" ours="$2" label="$3"
    local n_ref n_ours diff_count
    n_ref=$(zcat "$ref_gz" | wc -l)
    n_ours=$(wc -l < "$ours")

    # Round col 4 (AgeSq) to 4 decimals before comparing
    local awk_norm='NR==1 {print; next} {printf "%s\t%s\t%s\t%.4f\t%s\t%s\t%s\t%s\n", $1, $2, $3, $4, $5, $6, $7, $8}'
    diff_count=$({ diff <(zcat "$ref_gz" | awk -F'\t' "$awk_norm" | sort) \
                        <(awk -F'\t' "$awk_norm" "$ours" | sort) || true; } | wc -l)

    if [[ "$diff_count" -eq 0 ]]; then
        log "  ✓ $label  (Daniel=$n_ref rows, ours=$n_ours rows, ≡ after rounding AgeSq to 4 decimals)"
        log "      (note: raw byte-diff would show ~41k lines in AgeSq textual repr —"
        log "       Python 2 'str(float)' vs Python 3 shortest-roundtrip repr; values are identical)"
        return 0
    else
        log "  ✗ $label  (Daniel=$n_ref, ours=$n_ours, $diff_count diff lines after AgeSq rounding to 4 dec)"
        return 1
    fi
}

VAL_OK=0
validate      "$REF_CASES_CONTROL" "$CASES_CONTROL" "cases_control.txt    "  || VAL_OK=1
validate_covs "$REF_COVS"          "$COVS"          "covs.txt             "  || VAL_OK=1
validate      "$REF_REGIONS"       "$REGIONS"       "gene_list_regions.txt"  || VAL_OK=1

set -e

# Cleanup tmp
rm -rf "$TMP_DIR"

log ""
if [[ $VAL_OK -eq 0 ]]; then
    log "DONE — Phase 4 PASSED ✓"
    log "  All 3 preparatory files match Daniel's references (set-equality)."
    log "  Outputs in: $OUT_DIR"
    log "  Ready for Phase 5 (fused Phase 6 IBD trick + Phase 7 biobin run)."
    exit 0
else
    log "DONE — Phase 4 ran but at least one output DIFFERS from Daniel"
    log "Investigate the diff lines above before proceeding."
    exit 2
fi
