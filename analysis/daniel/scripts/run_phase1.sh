#!/usr/bin/env bash
# run_phase1.sh — Phase 1 of Hui et al. 2023 replication (gene list + annotation filter)
#
# Reproduces steps 1.2 and 1.3 of Daniel's runbook (data/PMBB_Exome/README.gz lines 25, 46)
# using PMBB v2 raw data via data/pmbb_v2/. Step 1.1 uses Daniel's pre-curated gene list
# directly (manual curation, not mechanically reproducible).
#
# Outputs: analysis/daniel/outputs/phase1/
# Log:     analysis/daniel/logs/phase1/run_<timestamp>.log
#
# Validation: diffs our outputs against Daniel's pre-built .gz intermediates.
# A clean diff confirms our pipeline + env reproduces the canonical Phase 1 outputs.

set -euo pipefail

# ───────── Paths ─────────
PROJECT_ROOT="/project/hall/analysis/hearing-loss-genomics"
cd "$PROJECT_ROOT"

# ───────── Env ─────────
# Activate venv if not already active (idempotent — safe if pre-activated)
if [[ -z "${VIRTUAL_ENV:-}" ]]; then
    # shellcheck disable=SC1091
    source venv/bin/activate
fi

OUT_DIR="analysis/daniel/outputs/phase1"
LOG_DIR="analysis/daniel/logs/phase1"
SCRIPTS_DIR="analysis/daniel/scripts/pmbb_exome"

GENE_LIST_SRC="data/PMBB_Exome/all_genes_including_ShadisList.txt.gz"
ANNOT_SRC="data/pmbb_v2/Exome/Variant_annotations/PMBB-Release-2020-2.0_genetic_exome_variant-annotation-counts.txt"

GENE_LIST_OUT="$OUT_DIR/all_genes_including_ShadisList.txt"
ANNOT_FULL_OUT="$OUT_DIR/annot_genes_full.txt"
ANNOT_FUNC_OUT="$OUT_DIR/annot_genes_full_funcToInclude.txt"

REF_ANNOT_FULL="data/PMBB_Exome/annot_genes_full.txt.gz"
REF_ANNOT_FUNC="data/PMBB_Exome/annot_genes_full_funcToInclude.txt.gz"

mkdir -p "$OUT_DIR" "$LOG_DIR"
TS=$(date +%Y%m%d_%H%M%S)
LOG="$LOG_DIR/run_${TS}.log"

# ───────── Logging helper ─────────
# All stdout/stderr is captured to $LOG via the exec redirect below;
# log() just adds a timestamp prefix to stdout (don't tee again — would duplicate).
log() { printf '[phase1] %s %s\n' "$(date +%H:%M:%S)" "$*"; }
fail() { log "ERROR: $*"; exit 1; }

exec > >(tee -a "$LOG") 2>&1

log "START Phase 1 replication"
log "Project root:   $PROJECT_ROOT"
log "Output dir:     $OUT_DIR"
log "Log file:       $LOG"

# ───────── Sanity: input files exist ─────────
log ""
log "Verifying inputs..."
[[ -f "$GENE_LIST_SRC" ]] || fail "Gene list missing: $GENE_LIST_SRC"
[[ -r "$ANNOT_SRC" ]]     || fail "v2 annotation file missing or unreadable: $ANNOT_SRC"
[[ -d "$SCRIPTS_DIR" ]]   || fail "Scripts dir missing: $SCRIPTS_DIR"

ANNOT_SIZE=$(du -h "$ANNOT_SRC" | cut -f1)
log "  gene list:        $GENE_LIST_SRC ($(du -h $GENE_LIST_SRC | cut -f1))"
log "  v2 annotation:    $ANNOT_SRC ($ANNOT_SIZE)"
log "  scripts:          $SCRIPTS_DIR"

# ───────── Step 1.1 — bring in Daniel's curated gene list ─────────
log ""
log "Step 1.1 — bringing in curated HL gene list (Daniel's manual curation, used as-is)"
zcat "$GENE_LIST_SRC" > "$GENE_LIST_OUT"
N_GENES=$(wc -l < "$GENE_LIST_OUT")
log "  written:          $GENE_LIST_OUT ($N_GENES genes)"

# Sanity: list size + spot-check known HL genes present + new-discovery gene absent
# Paper reports 173 genes; curated list has ~179 (extras likely X-linked, filtered later)
log "  sample genes:     $(head -5 $GENE_LIST_OUT | tr '\n' ' ')"

for expected_present in TCOF1 ESRRB OTOF; do
    if grep -qw "$expected_present" "$GENE_LIST_OUT"; then
        log "  ✓ $expected_present present (known HL gene — expected)"
    else
        fail "Known HL gene $expected_present NOT in list — list may be wrong"
    fi
done

# Negative control: ZNF175 was *discovered* in Phase 12 all-genes burden,
# it is NOT a member of the pre-burden known-HL-gene list.
if grep -qw "ZNF175" "$GENE_LIST_OUT"; then
    fail "ZNF175 unexpectedly present — known-HL list looks contaminated"
else
    log "  ✓ ZNF175 absent (correct — discovered gene, not in pre-burden list)"
fi

# ───────── Step 1.2 — filter PMBB annotation to HL genes ─────────
log ""
log "Step 1.2 — filtering v2 variant annotation to the HL gene set"
log "  (this reads the ${ANNOT_SIZE} annotation file; expect 1-3 min)"
T0=$(date +%s)
python "$SCRIPTS_DIR/only_HL_genes.py" "$GENE_LIST_OUT" "$ANNOT_SRC" > "$ANNOT_FULL_OUT"
T1=$(date +%s)
ELAPSED=$((T1 - T0))
N_ROWS=$(wc -l < "$ANNOT_FULL_OUT")
log "  written:          $ANNOT_FULL_OUT ($(du -h $ANNOT_FULL_OUT | cut -f1), $N_ROWS rows, ${ELAPSED}s)"

# ───────── Step 1.3 — filter to pLoF + missense REVEL>0.6 ─────────
log ""
log "Step 1.3 — filtering to functional classes (pLoF + missense REVEL>0.6)"
T0=$(date +%s)
python "$SCRIPTS_DIR/only_func_cats_to_include.py" "$ANNOT_FULL_OUT" > "$ANNOT_FUNC_OUT"
T1=$(date +%s)
ELAPSED=$((T1 - T0))
N_ROWS=$(wc -l < "$ANNOT_FUNC_OUT")
log "  written:          $ANNOT_FUNC_OUT ($(du -h $ANNOT_FUNC_OUT | cut -f1), $N_ROWS rows, ${ELAPSED}s)"

# ───────── Validation: semantic set-equality vs Daniel's outputs ─────────
#
# Why semantic, not byte-for-byte:
# Daniel's intermediates were generated from variant-annotations.txt (now removed
# from disk); we re-run against the replacement variant-annotation-counts.txt. The
# pipeline-relevant content (variant IDs, gene names, function categories, REVEL)
# matches, but metadata columns and header strings differ — schema migration that
# Daniel himself flagged at runbook line 42-43. See results/phase1/phase1_replication_report.md.
#
# A pass here means: same set of variants and same set of genes selected as Daniel.
log ""
log "Validation — semantic equality against Daniel's pre-built intermediates"
log "  (set of variant IDs + set of genes; not byte-for-byte — see report for rationale)"

validate_semantic() {
    local ref_gz="$1" ours="$2" label="$3"

    local n_ref n_ours
    n_ref=$(zcat "$ref_gz" | wc -l)
    n_ours=$(wc -l < "$ours")

    if [[ "$n_ref" -ne "$n_ours" ]]; then
        log "  ✗ $label: row count mismatch (Daniel=$n_ref, ours=$n_ours)"
        return 1
    fi

    # Variant ID set (col 1), skipping header
    local ids_diff_count
    ids_diff_count=$(diff <(zcat "$ref_gz" | tail -n +2 | cut -f1 | sort -u) \
                          <(tail -n +2 "$ours" | cut -f1 | sort -u) | wc -l)

    # Gene set (col 8 = Gene.refGene), skipping header
    local genes_diff_count
    genes_diff_count=$(diff <(zcat "$ref_gz" | tail -n +2 | cut -f8 | sort -u) \
                            <(tail -n +2 "$ours" | cut -f8 | sort -u) | wc -l)

    # Informational: byte-diff size + header comparison
    local ref_h1 ours_h1
    ref_h1=$(zcat "$ref_gz" | head -1 | cut -f1)
    ours_h1=$(head -1 "$ours" | cut -f1)

    if [[ "$ids_diff_count" -eq 0 && "$genes_diff_count" -eq 0 ]]; then
        log "  ✓ $label"
        log "      rows=$n_ours  variant-ID set ≡ Daniel  gene set ≡ Daniel"
        if [[ "$ref_h1" != "$ours_h1" ]]; then
            log "      (info: col 1 header differs — Daniel='$ref_h1', ours='$ours_h1' — schema migration, expected)"
        fi
        return 0
    else
        log "  ✗ $label: SETS DIFFER"
        [[ "$ids_diff_count" -ne 0 ]]   && log "      variant ID sets differ ($ids_diff_count diff lines)"
        [[ "$genes_diff_count" -ne 0 ]] && log "      gene sets differ ($genes_diff_count diff lines)"
        return 1
    fi
}

VAL_OK=0
validate_semantic "$REF_ANNOT_FULL" "$ANNOT_FULL_OUT" "annot_genes_full.txt"               || VAL_OK=1
validate_semantic "$REF_ANNOT_FUNC" "$ANNOT_FUNC_OUT" "annot_genes_full_funcToInclude.txt" || VAL_OK=1

log ""
if [[ $VAL_OK -eq 0 ]]; then
    log "DONE — Phase 1 PASSED ✓"
    log "  Pipeline produces same variant + gene sets as Daniel's reference."
    log "  Outputs in: $OUT_DIR"
    exit 0
else
    log "DONE — Phase 1 ran to completion but variant/gene sets DIFFER from Daniel's"
    log "Investigate possible causes:"
    log "  - gene list contamination or missing genes"
    log "  - changes to only_HL_genes.py or only_func_cats_to_include.py logic"
    log "  - PMBB v2 annotation file changed beyond schema migration (variants added/removed)"
    log "  - Python version drift (Daniel: ?, current: $(python --version 2>&1))"
    exit 2
fi
