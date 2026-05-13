#!/usr/bin/env bash
# run_phase2.sh — Phase 2 of Hui et al. 2023 replication: SNP ID reconciliation.
#
# Reproduces steps 2.2-2.6 of Daniel's runbook (lines 47-62). Light mode:
# skips Step 2.1 (the zgrep on 781 GB of pVCFs) and reuses Daniel's
# pre-generated vcf_SNP_IDs_allchr file, which was a pure `cut -f 1-5`
# extract of the same pVCFs we still have access to today (no schema risk).
#
# Inputs:
#   - analysis/daniel/outputs/phase1/annot_genes_full_funcToInclude.txt (master list from Phase 1)
#   - data/PMBB_Exome/vcf_SNP_IDs/vcf_SNP_IDs_allchr.txt.gz (Daniel's pre-extracted IDs)
#
# Outputs: analysis/daniel/outputs/phase2/
# Log:     analysis/daniel/logs/phase2/run_<timestamp>.log
#
# Validation: semantic set-equality of VCF_IDs against Daniel's matched_*.txt
# and matched_*.extract reference files.

set -euo pipefail

# ───────── Paths ─────────
PROJECT_ROOT="/project/hall/analysis/hearing-loss-genomics"
cd "$PROJECT_ROOT"

# ───────── Env ─────────
if [[ -z "${VIRTUAL_ENV:-}" ]]; then
    # shellcheck disable=SC1091
    source venv/bin/activate
fi

OUT_DIR="analysis/daniel/outputs/phase2"
LOG_DIR="analysis/daniel/logs/phase2"
SCRIPTS_DIR="analysis/daniel/scripts/pmbb_exome"

# Phase 1 master list (input)
MASTER_LIST="analysis/daniel/outputs/phase1/annot_genes_full_funcToInclude.txt"

# Daniel's pre-extracted pVCF IDs (light-mode input — skips zgrep on 781 GB)
DANIEL_ALLCHR_GZ="data/PMBB_Exome/vcf_SNP_IDs/vcf_SNP_IDs_allchr.txt.gz"

# Daniel's reference outputs (for validation)
REF_MATCHED="data/PMBB_Exome/matched_snp_IDs_annot_pVCF.txt.gz"
REF_MATCHED_CLEAN="data/PMBB_Exome/matched_snp_IDs_annot_pVCF_noNA_noMultiallelic.txt.gz"
REF_EXTRACT="data/PMBB_Exome/matched_snp_IDs_annot_pVCF_noNA_noMultiallelic.extract.gz"

# Our outputs
ALLCHR_VCF_IDS="$OUT_DIR/vcf_SNP_IDs_allchr.txt"
MATCHED="$OUT_DIR/matched_snp_IDs_annot_pVCF.txt"
MATCHED_CLEAN="$OUT_DIR/matched_snp_IDs_annot_pVCF_noNA_noMultiallelic.txt"
EXTRACT="$OUT_DIR/matched_snp_IDs_annot_pVCF_noNA_noMultiallelic.extract"

mkdir -p "$OUT_DIR" "$LOG_DIR"
TS=$(date +%Y%m%d_%H%M%S)
LOG="$LOG_DIR/run_${TS}.log"

# ───────── Logging ─────────
log() { printf '[phase2] %s %s\n' "$(date +%H:%M:%S)" "$*"; }
fail() { log "ERROR: $*"; exit 1; }
exec > >(tee -a "$LOG") 2>&1

log "START Phase 2 — SNP ID reconciliation (LIGHT MODE)"
log "Project root:   $PROJECT_ROOT"
log "Output dir:     $OUT_DIR"
log "Log file:       $LOG"
log ""
log "Light mode: Step 2.1 (zgrep on 781 GB pVCFs) skipped — reusing Daniel's"
log "  pre-generated vcf_SNP_IDs_allchr.txt.gz from same v2 pVCFs (no schema risk)."

# ───────── Verify inputs ─────────
log ""
log "Verifying inputs..."
[[ -f "$MASTER_LIST" ]] || fail "Phase 1 master list missing — run Phase 1 first: $MASTER_LIST"
[[ -f "$DANIEL_ALLCHR_GZ" ]] || fail "Daniel's pre-extracted IDs missing: $DANIEL_ALLCHR_GZ"
[[ -f "$SCRIPTS_DIR/annot_IDs_vs_pVCF.py" ]] || fail "Script missing"
log "  master list:           $MASTER_LIST ($(wc -l < $MASTER_LIST) rows)"
log "  Daniel's allchr IDs:   $DANIEL_ALLCHR_GZ ($(du -h $DANIEL_ALLCHR_GZ | cut -f1))"

# ───────── Step 2.1 — bring in pre-extracted VCF IDs ─────────
log ""
log "Step 2.1 — bringing in Daniel's pre-extracted VCF IDs (allchr, gunzip)"
T0=$(date +%s)
zcat "$DANIEL_ALLCHR_GZ" > "$ALLCHR_VCF_IDS"
T1=$(date +%s)
N_ALLCHR=$(wc -l < "$ALLCHR_VCF_IDS")
log "  written:               $ALLCHR_VCF_IDS ($(du -h $ALLCHR_VCF_IDS | cut -f1), $N_ALLCHR rows, $((T1-T0))s)"
log "  this is all v2 pVCF variants across chr1-22 with chr/pos/ID/ref/alt"

# ───────── Step 2.3 — join master list with VCF IDs ─────────
log ""
log "Step 2.3 — joining master variant list with pVCF IDs (annot_IDs_vs_pVCF.py)"
log "  Reads $N_ALLCHR pVCF IDs and joins on (chr, pos) with the master list."
T0=$(date +%s)
python "$SCRIPTS_DIR/annot_IDs_vs_pVCF.py" "$MASTER_LIST" "$ALLCHR_VCF_IDS" \
    | sort -gk1,1 -gk2,2 > "$MATCHED"
T1=$(date +%s)
N_MATCHED=$(wc -l < "$MATCHED")
log "  written:               $MATCHED ($(du -h $MATCHED | cut -f1), $N_MATCHED rows, $((T1-T0))s)"

# ───────── Step 2.4 — drop NA + multi-allelic ─────────
log ""
log "Step 2.4 — filtering out NA rows (not in pVCF) and multi-allelic variants"
grep -v NA "$MATCHED" | grep -v ";" > "$MATCHED_CLEAN"
N_CLEAN=$(wc -l < "$MATCHED_CLEAN")
N_DROPPED=$((N_MATCHED - N_CLEAN))
log "  written:               $MATCHED_CLEAN ($(du -h $MATCHED_CLEAN | cut -f1), $N_CLEAN rows)"
log "  dropped $N_DROPPED rows: $(grep -c NA $MATCHED || true) NA + $(grep -c ';' $MATCHED || true) multi-allelic (overlap possible)"

# ───────── Step 2.5 — sanity check ref/alt alignment ─────────
log ""
log "Step 2.5 — sanity check: VCF.ref == annot.ref (col 4 vs col 7)"
N_MISMATCH=$(awk -F'\t' 'NR>1 && $4 != $7 {n++} END {print n+0}' "$MATCHED_CLEAN")
if [[ "$N_MISMATCH" -eq 0 ]]; then
    log "  ✓ no ref-allele mismatches between annotation and pVCF"
else
    log "  ⚠ $N_MISMATCH rows with mismatched ref allele — investigate"
fi

# ───────── Step 2.6 — extract VCF_IDs column ─────────
log ""
log "Step 2.6 — extracting VCF_ID column (col 3) into plink --extract format"
# Daniel's runbook says `cut -f 3 ...` but his actual .extract on disk has NO header
# (he must have stripped it elsewhere). plink --extract expects one ID per line with
# no header — including "VCF_ID" as a string would make plink hunt for that "SNP".
tail -n +2 "$MATCHED_CLEAN" | cut -f 3 > "$EXTRACT"
N_EXTRACT=$(wc -l < "$EXTRACT")
log "  written:               $EXTRACT ($N_EXTRACT lines, header skipped)"
log "  (this is the input plink will use in Phase 3's --extract flag)"

# ───────── Validation — semantic set-equality ─────────
log ""
log "Validation — semantic set-equality against Daniel's reference"
log "  (compares sorted unique VCF_ID sets, not byte-for-byte)"

# Disable set -e for validation block: `diff` returns 1 when files differ, which
# would abort the script under pipefail before we get to log the result.
set +e

VAL_OK=0

# (a) Compare matched_clean VCF_ID column
ref_clean_ids=$(zcat "$REF_MATCHED_CLEAN" | tail -n +2 | cut -f3 | sort -u)
our_clean_ids=$(tail -n +2 "$MATCHED_CLEAN" | cut -f3 | sort -u)
diff_count=$({ diff <(echo "$ref_clean_ids") <(echo "$our_clean_ids") || true; } | wc -l)
n_ref=$(zcat "$REF_MATCHED_CLEAN" | wc -l)
n_ours=$(wc -l < "$MATCHED_CLEAN")
if [[ "$diff_count" -eq 0 ]]; then
    log "  ✓ matched_clean VCF_ID set ≡ Daniel  (Daniel=$n_ref rows, ours=$n_ours rows)"
else
    log "  ✗ matched_clean: VCF_ID sets differ ($diff_count diff lines)"
    n_ours_only=$(comm -13 <(echo "$ref_clean_ids") <(echo "$our_clean_ids") | wc -l)
    n_ref_only=$(comm -23 <(echo "$ref_clean_ids") <(echo "$our_clean_ids") | wc -l)
    log "      $n_ours_only IDs only in ours, $n_ref_only IDs only in Daniel"
    VAL_OK=1
fi

# (b) Compare .extract files (line set)
ref_extract=$(zcat "$REF_EXTRACT" | sort -u)
our_extract=$(sort -u "$EXTRACT")
diff_count=$({ diff <(echo "$ref_extract") <(echo "$our_extract") || true; } | wc -l)
n_ref=$(zcat "$REF_EXTRACT" | wc -l)
n_ours=$(wc -l < "$EXTRACT")
if [[ "$diff_count" -eq 0 ]]; then
    log "  ✓ .extract line set ≡ Daniel  (Daniel=$n_ref lines, ours=$n_ours lines)"
else
    log "  ✗ .extract: sets differ ($diff_count diff lines)"
    VAL_OK=1
fi

set -e

log ""
if [[ $VAL_OK -eq 0 ]]; then
    log "DONE — Phase 2 PASSED ✓"
    log "  Pipeline produces same VCF_ID extract list as Daniel's reference."
    log "  Outputs in:  $OUT_DIR"
    log "  Key file for Phase 3 (plink --extract): $EXTRACT"
    exit 0
else
    log "DONE — Phase 2 ran but VCF_ID sets DIFFER from Daniel's reference"
    log "Possible causes:"
    log "  - master list (Phase 1 output) differs from Daniel's master"
    log "  - annot_IDs_vs_pVCF.py logic changed"
    log "  - sort order affected the dict iteration (unlikely with sort step)"
    exit 2
fi
