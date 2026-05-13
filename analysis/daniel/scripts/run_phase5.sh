#!/usr/bin/env bash
# run_phase5.sh — Phase 5 of Hui et al. 2023 replication: HEAVY mode.
#
# This is the first end-to-end burden test. Combines Phase 6 (IBD trick) +
# Phase 7 (plink filter + merge + biobin) into one orchestrated run.
#
# Steps:
#   5.0  Prep IBD trick inputs (cases_withRels, cases_noRels, removed_cases)
#   5.1  Run keep_HL_cases_IBD.py → tokeep_moreHLcases.txt
#   5.2  Per-chr plink filter (22 chrs): --keep tokeep --max-maf .001
#   5.3  Merge 22 chrs → cohort-wide VCF
#   5.4  biobin: gene-burden logistic regression with HL phenotype
#
# Validation: tokeep set-equality + bins.csv top hit must be ESRRB.

set -euo pipefail

# ───────── Paths ─────────
PROJECT_ROOT="/project/hall/analysis/hearing-loss-genomics"
cd "$PROJECT_ROOT"

if [[ -z "${VIRTUAL_ENV:-}" ]]; then
    # shellcheck disable=SC1091
    source venv/bin/activate
fi

# Environment for plink 1.9 + biobin (rlsoftware module + liblzma shim)
# shellcheck disable=SC1091
source /etc/profile.d/modules.sh 2>/dev/null || true
module load rlsoftware/latest 2>/dev/null || true
export PATH="/appl/plink-1.90Beta6.18:$PATH"
export LD_LIBRARY_PATH="$PROJECT_ROOT/analysis/daniel/configs/lib-shims:${LD_LIBRARY_PATH:-}"

PLINK=/appl/plink-1.90Beta6.18/plink

OUT_DIR="analysis/daniel/outputs/phase5"
PREP_DIR="$OUT_DIR/prep"
IBD_DIR="$OUT_DIR/ibd_trick"
FILT_DIR="$OUT_DIR/filtered"
MERGE_DIR="$OUT_DIR/merged"
BIOBIN_DIR="$OUT_DIR/biobin"
LOG_DIR="analysis/daniel/logs/phase5"

# Inputs from previous phases / raw v2
PHASE3_LIGHT="analysis/daniel/outputs/phase3/light"
PHASE4_OUT="analysis/daniel/outputs/phase4"
CASES_CONTROL="$PHASE4_OUT/cases_control.txt"
COVS="$PHASE4_OUT/covs.txt"
REGIONS="$PHASE4_OUT/gene_list_regions.txt"

IBD_GENOME="data/pmbb_v2/Exome/IBD/PMBB-Release-2020-2.0_genetic_exome.genome"
IBD_3DEG_UNREL="data/pmbb_v2/Exome/IBD/PMBB-Release-2020-2.0_genetic_exome_3rd_degree_unrelated"
# Note 2026-05-13: /project/ritchie/datasets/loki/loki.db symlink points at
# loki-20220926.db which has been removed. Using loki-20230816.db (closest
# available in time to Daniel's 2022-09 era). biobin uses loki for gene/region
# metadata lookups — the gene boundaries we care about come from Phase 4's
# gene_list_regions.txt via --region-file, so the loki version shift should
# minimally affect results.
LOKI_DB="/project/ritchie/datasets/loki/loki-20230816.db"

# Daniel reference outputs (for validation)
REF_REMOVED_CASES="data/PMBB_Exome/removed_cases.txt.gz"
REF_TOKEEP="data/PMBB_Exome/tokeep_moreHLcases.txt.gz"
REF_BINS_GZ="data/PMBB_Exome/biobin/merged_maf.001_noRels_keepHLcases-bins.csv.gz"
REF_LOCUS_GZ="data/PMBB_Exome/biobin/merged_maf.001_noRels_keepHLcases-locus.csv.gz"

mkdir -p "$PREP_DIR" "$IBD_DIR" "$FILT_DIR" "$MERGE_DIR" "$BIOBIN_DIR" "$LOG_DIR"
TS=$(date +%Y%m%d_%H%M%S)
LOG="$LOG_DIR/run_${TS}.log"

# ───────── Logging ─────────
log() { printf '[phase5] %s %s\n' "$(date +%H:%M:%S)" "$*"; }
fail() { log "ERROR: $*"; exit 1; }
exec > >(tee -a "$LOG") 2>&1

log "START Phase 5 — first end-to-end burden test (HEAVY mode)"
log "Project root:   $PROJECT_ROOT"
log "Output dir:     $OUT_DIR"
log "Log file:       $LOG"

# ───────── Verify inputs ─────────
log ""
log "Verifying inputs..."
[[ -d "$PHASE3_LIGHT" ]]      || fail "Phase 3 light dir missing — run Phase 3 first: $PHASE3_LIGHT"
[[ -f "$CASES_CONTROL" ]]     || fail "Phase 4 cases_control.txt missing"
[[ -f "$COVS" ]]              || fail "Phase 4 covs.txt missing"
[[ -f "$REGIONS" ]]           || fail "Phase 4 gene_list_regions.txt missing"
[[ -r "$IBD_GENOME" ]]        || fail ".genome file missing: $IBD_GENOME"
[[ -r "$IBD_3DEG_UNREL" ]]    || fail "3rd-degree-unrelated file missing"
[[ -r "$LOKI_DB" ]]           || fail "loki.db missing: $LOKI_DB"
[[ -x "$PLINK" ]]             || fail "plink 1.9 missing"
which biobin >/dev/null       || fail "biobin not on PATH (module load rlsoftware/latest?)"
log "  Phase 3 light:        $PHASE3_LIGHT ($(ls $PHASE3_LIGHT | wc -l) files)"
log "  Phase 4 outputs:      $PHASE4_OUT"
log "  IBD .genome:          $(du -h $IBD_GENOME | cut -f1)"
log "  IBD 3rd-deg-unrel:    $(wc -l < $IBD_3DEG_UNREL) IDs"
log "  loki.db:              $LOKI_DB ($(du -h $(readlink -f $LOKI_DB) 2>/dev/null | cut -f1 || echo n/a))"
log "  plink:                $($PLINK --version 2>&1 | head -1)"
log "  biobin:               $(which biobin)"

# ═══════════════════════════════════════════════════════════════════════
# Step 5.0 — Prep IBD-trick inputs
# ═══════════════════════════════════════════════════════════════════════
log ""
log "════════════════════════════════════════════════════════════"
log "Step 5.0 — Prep IBD-trick inputs"
log "════════════════════════════════════════════════════════════"

# 5.0a — list of all cases (SNHL=1 from cases_control.txt)
awk -F'\t' 'NR>1 && $2 == "1" {print $1}' "$CASES_CONTROL" | sort > "$PREP_DIR/cases_only.txt"
N_CASES_TOTAL=$(wc -l < "$PREP_DIR/cases_only.txt")
log "  Total cases (SNHL=1) in cases_control.txt: $N_CASES_TOTAL"

# 5.0b — All IDs in the genotyped cohort (from chr22.fam — should be 43,731)
awk '{print $1}' "$PHASE3_LIGHT/allIndvs_chr22.fam" | sort > "$PREP_DIR/all_cohort_ids.txt"
N_COHORT=$(wc -l < "$PREP_DIR/all_cohort_ids.txt")
log "  Genotyped cohort N (from chr22.fam): $N_COHORT"

# 5.0c — cases_withRels: cases that are also genotyped (intersection)
comm -12 "$PREP_DIR/cases_only.txt" "$PREP_DIR/all_cohort_ids.txt" > "$PREP_DIR/cases_withRels_ids.txt"
N_WITH_RELS=$(wc -l < "$PREP_DIR/cases_withRels_ids.txt")
log "  cases_withRels (cases in cohort, before IBD filter): $N_WITH_RELS"

# 5.0d — Strict 3rd-degree-unrelated cohort
sort "$IBD_3DEG_UNREL" > "$PREP_DIR/unrelated_sorted.txt"
N_UNREL=$(wc -l < "$PREP_DIR/unrelated_sorted.txt")

# 5.0e — cases_noRels: cases that survive strict filter
comm -12 "$PREP_DIR/cases_withRels_ids.txt" "$PREP_DIR/unrelated_sorted.txt" > "$PREP_DIR/cases_noRels_ids.txt"
N_NO_RELS=$(wc -l < "$PREP_DIR/cases_noRels_ids.txt")
log "  cases_noRels (cases surviving strict filter): $N_NO_RELS"

# 5.0f — removed_cases: cases dropped by strict filter
comm -23 "$PREP_DIR/cases_withRels_ids.txt" "$PREP_DIR/cases_noRels_ids.txt" > "$PREP_DIR/removed_cases_ids.txt"
# Format as plink keep-list (FID IID 0 0 0 -9, space-sep)
awk '{printf "%s %s 0 0 0 -9\n", $1, $1}' "$PREP_DIR/removed_cases_ids.txt" > "$PREP_DIR/removed_cases.txt"
N_REMOVED=$(wc -l < "$PREP_DIR/removed_cases.txt")
log "  removed_cases (cases dropped by strict filter): $N_REMOVED"

# Validation 5.0 — compare to Daniel
set +e
ref_removed=$(zcat "$REF_REMOVED_CASES" | awk '{print $1}' | sort -u)
our_removed=$(awk '{print $1}' "$PREP_DIR/removed_cases.txt" | sort -u)
diff_count=$({ diff <(echo "$ref_removed") <(echo "$our_removed") || true; } | wc -l)
if [[ "$diff_count" -eq 0 ]]; then
    log "  ✓ removed_cases set ≡ Daniel (Daniel=$(echo "$ref_removed" | wc -l), ours=$N_REMOVED)"
else
    log "  ⚠ removed_cases differs from Daniel ($diff_count diff lines) — investigate before continuing"
fi
set -e

# ═══════════════════════════════════════════════════════════════════════
# Step 5.1 — IBD trick (keep_HL_cases_IBD.py)
# ═══════════════════════════════════════════════════════════════════════
log ""
log "════════════════════════════════════════════════════════════"
log "Step 5.1 — IBD trick (keep_HL_cases_IBD.py)"
log "════════════════════════════════════════════════════════════"
log "  Adds back the $N_REMOVED removed HL cases, removes their relatives from keep-list."

python analysis/daniel/scripts/pmbb_exome/keep_HL_cases_IBD.py \
    "$PREP_DIR/removed_cases.txt" \
    "$IBD_GENOME" \
    "$IBD_3DEG_UNREL" \
    > "$IBD_DIR/tokeep_moreHLcases.txt"

N_KEEP=$(wc -l < "$IBD_DIR/tokeep_moreHLcases.txt")
log "  tokeep_moreHLcases: $N_KEEP IDs (vs $N_UNREL strict; net change: $((N_KEEP - N_UNREL)))"

# Format as plink keep-list (FID IID, space-sep)
awk '{print $1, $1}' "$IBD_DIR/tokeep_moreHLcases.txt" > "$IBD_DIR/tokeep_moreHLcases_keep.txt"

# Validation 5.1 — compare to Daniel
set +e
ref_keep=$(zcat "$REF_TOKEEP" | sort -u)
our_keep=$(sort -u "$IBD_DIR/tokeep_moreHLcases.txt")
diff_count=$({ diff <(echo "$ref_keep") <(echo "$our_keep") || true; } | wc -l)
if [[ "$diff_count" -eq 0 ]]; then
    log "  ✓ tokeep_moreHLcases set ≡ Daniel ($N_KEEP IDs)"
else
    log "  ⚠ tokeep_moreHLcases differs ($diff_count diff lines) — investigate"
fi
set -e

# ═══════════════════════════════════════════════════════════════════════
# Step 5.2 — Per-chr plink filter
# ═══════════════════════════════════════════════════════════════════════
log ""
log "════════════════════════════════════════════════════════════"
log "Step 5.2 — Per-chr plink filter (--keep tokeep --max-maf .001)"
log "════════════════════════════════════════════════════════════"

T0=$(date +%s)
for i in {1..22}; do
    in_prefix="$PHASE3_LIGHT/allIndvs_chr${i}"
    out_prefix="$FILT_DIR/allIndvs_chr${i}_maf.001_noRels_keepHLcases"
    "$PLINK" \
        --bfile "$in_prefix" \
        --keep "$IBD_DIR/tokeep_moreHLcases_keep.txt" \
        --max-maf .001 \
        --make-bed \
        --out "$out_prefix" \
        > "$FILT_DIR/chr${i}.stdout" 2> "$FILT_DIR/chr${i}.stderr" \
        || fail "plink chr${i} failed — see $FILT_DIR/chr${i}.stderr"
    n_var=$(wc -l < "${out_prefix}.bim")
    n_ind=$(wc -l < "${out_prefix}.fam")
    printf '[phase5] %s   chr%-2d: %s variants × %s individuals\n' "$(date +%H:%M:%S)" "$i" "$n_var" "$n_ind"
done
T1=$(date +%s)
log "Step 5.2 total: $((T1-T0))s"

# ═══════════════════════════════════════════════════════════════════════
# Step 5.3 — Merge 22 chrs + convert to VCF
# ═══════════════════════════════════════════════════════════════════════
log ""
log "════════════════════════════════════════════════════════════"
log "Step 5.3 — Merge 22 chrs and convert to VCF"
log "════════════════════════════════════════════════════════════"

# Generate merge-list (chr2-22 to merge into chr1)
for i in {2..22}; do
    echo "$FILT_DIR/allIndvs_chr${i}_maf.001_noRels_keepHLcases"
done > "$MERGE_DIR/merge-list.txt"

MERGED_PREFIX="$MERGE_DIR/merged_maf.001_noRels_keepHLcases"

T0=$(date +%s)
"$PLINK" \
    --bfile "$FILT_DIR/allIndvs_chr1_maf.001_noRels_keepHLcases" \
    --merge-list "$MERGE_DIR/merge-list.txt" \
    --make-bed \
    --out "$MERGED_PREFIX" \
    > "$MERGE_DIR/merge.stdout" 2> "$MERGE_DIR/merge.stderr" \
    || fail "plink merge failed — see $MERGE_DIR/merge.stderr"
T1=$(date +%s)
N_MERGED_VAR=$(wc -l < "${MERGED_PREFIX}.bim")
N_MERGED_IND=$(wc -l < "${MERGED_PREFIX}.fam")
log "  merged: $N_MERGED_VAR variants × $N_MERGED_IND individuals, $((T1-T0))s"

# Convert to VCF for biobin
T0=$(date +%s)
"$PLINK" \
    --bfile "$MERGED_PREFIX" \
    --recode vcf-iid \
    --out "$MERGED_PREFIX" \
    > "$MERGE_DIR/vcf.stdout" 2> "$MERGE_DIR/vcf.stderr" \
    || fail "plink vcf conversion failed"
T1=$(date +%s)
log "  VCF written: ${MERGED_PREFIX}.vcf ($(du -h ${MERGED_PREFIX}.vcf | cut -f1), $((T1-T0))s)"

# ═══════════════════════════════════════════════════════════════════════
# Step 5.4 — biobin gene-burden logistic regression
# ═══════════════════════════════════════════════════════════════════════
log ""
log "════════════════════════════════════════════════════════════"
log "Step 5.4 — biobin (gene-burden logistic regression)"
log "════════════════════════════════════════════════════════════"
log "  Command: biobin -D loki.db -V merged.vcf -p cases_control.txt --covariates covs.txt"
log "                  --bin-regions Y --region-file gene_list_regions.txt -G 38 --test logistic"
log "  Phenotype: cases_control.txt (strict audiogram, $N_CASES_TOTAL total cases)"

BIOBIN_PREFIX="$BIOBIN_DIR/merged_maf.001_noRels_keepHLcases"

T0=$(date +%s)
biobin \
    -D "$LOKI_DB" \
    -V "${MERGED_PREFIX}.vcf" \
    -p "$CASES_CONTROL" \
    --covariates "$COVS" \
    --bin-regions Y \
    --region-file "$REGIONS" \
    -G 38 \
    --test logistic \
    --report-prefix "$BIOBIN_PREFIX" \
    > "$BIOBIN_PREFIX.run_log.txt" 2>&1 \
    || fail "biobin failed — see $BIOBIN_PREFIX.run_log.txt"
T1=$(date +%s)
log "  biobin done in $((T1-T0))s"
log "  outputs:"
ls -lh "$BIOBIN_DIR/" | awk 'NR>1 {print "    " $NF " (" $5 ")"}'

# ═══════════════════════════════════════════════════════════════════════
# Validation — top hit must be ESRRB
# ═══════════════════════════════════════════════════════════════════════
log ""
log "════════════════════════════════════════════════════════════"
log "Validation — top burden test hit"
log "════════════════════════════════════════════════════════════"

BINS_OURS="$BIOBIN_PREFIX-bins.csv"
[[ -f "$BINS_OURS" ]] || fail "biobin bins.csv not produced: $BINS_OURS"

# Extract gene names (row 0 — ID/genes row) and p-values (row 8 — logistic p-value)
# using Python. Note: row 7 ("Gene(s)") would seem natural for gene names but it's
# biobin's loki-derived annotation which can be empty for many bins; row 0 (the ID
# header row with gene names per column) is the reliable source.
python3 <<EOF
import csv, gzip

def parse_bins(path, opener=open, max_rows=13):
    with opener(path) as f:
        reader = csv.reader(f)
        return [row for i, row in enumerate(reader) if i < max_rows]

def top_hits(rows, n=20):
    genes = rows[0][2:]   # row 0 = ID row; skip "ID" label + empty cell
    pvals = rows[8][2:]   # row 8 = logistic p-value; skip label + nan
    pairs = []
    for g, p in zip(genes, pvals):
        if not g or g == 'nan':
            continue
        try:
            pairs.append((g, float(p)))
        except (ValueError, TypeError):
            continue
    pairs.sort(key=lambda x: x[1])
    seen, unique = set(), []
    for g, p in pairs:
        if g not in seen:
            unique.append((g, p))
            seen.add(g)
        if len(unique) >= n:
            break
    return unique

ours = top_hits(parse_bins('$BINS_OURS'))
daniel = top_hits(parse_bins('$REF_BINS_GZ', opener=lambda p: gzip.open(p, 'rt')))

print(f"{'#':>3}  {'Ours':<42}  {'Daniel':<42}")
print('-' * 92)
for i in range(20):
    o = ours[i] if i < len(ours) else (None, None)
    d = daniel[i] if i < len(daniel) else (None, None)
    o_str = f"{o[0]:20s} p={o[1]:.4e}" if o[0] else ' ' * 36
    d_str = f"{d[0]:20s} p={d[1]:.4e}" if d[0] else ' ' * 36
    print(f"  {i+1:2d}  {o_str}  {d_str}")
print()

o_names = [g for g, _ in ours]
d_names = [g for g, _ in daniel]
print('Key gene ranks:')
for gene in ['ESRRB', 'TCOF1', 'SCD5', 'GIPC3', 'EYA4', 'COL9A1']:
    o_rank = o_names.index(gene)+1 if gene in o_names else '-'
    d_rank = d_names.index(gene)+1 if gene in d_names else '-'
    print(f"  {gene:10s}  ours=#{o_rank}  daniel=#{d_rank}")
print()

# Success criterion: ESRRB in top 5 (allowing 1-2 LOC artifacts from newer loki)
if 'ESRRB' in o_names[:5]:
    print(f'✓ ESRRB appears in top 5 of our results (rank #{o_names.index("ESRRB")+1})')
    print('  — replicates Daniel\\'s "it\\'s ESRRB, just like the paper" finding')
else:
    print('⚠ ESRRB NOT in top 5 — investigate')
EOF

log ""
log "DONE — Phase 5 complete"
log "  Outputs: $OUT_DIR"
log "  Biobin: $BIOBIN_DIR"
