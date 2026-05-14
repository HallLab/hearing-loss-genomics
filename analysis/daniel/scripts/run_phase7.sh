#!/usr/bin/env bash
# run_phase7.sh — Phase 7: ZNF175 carrier deep-dive + second-hit hypothesis test.
#
# Project's main biological goal — extending beyond paper-published phases.
# Uses Daniel's preserved 142-carrier pLoF list as input. For each carrier,
# extracts variants in all 173 known HL genes; tabulates second-hit profile;
# compares HL-cases-among-carriers vs HL-controls-among-carriers.
#
# Outputs: analysis/daniel/outputs/phase7/
# Log:     analysis/daniel/logs/phase7/run_<timestamp>.log
#
# NOT a statistical "discovery" task — this is biological characterization of
# the carriers that Doug Epstein flagged via his mouse Zfp719 work.

set -euo pipefail

PROJECT_ROOT="/project/hall/analysis/hearing-loss-genomics"
cd "$PROJECT_ROOT"

if [[ -z "${VIRTUAL_ENV:-}" ]]; then
    # shellcheck disable=SC1091
    source venv/bin/activate
fi

# env (plink 1.9 explicit; lib shim if biobin needed later)
export PATH="/appl/plink-1.90Beta6.18:$PATH"

PLINK=/appl/plink-1.90Beta6.18/plink

OUT_DIR="analysis/daniel/outputs/phase7"
GENO_DIR="$OUT_DIR/genos"
RESULTS_DIR="$OUT_DIR/results"
LOG_DIR="analysis/daniel/logs/phase7"

# Inputs
DANIEL_CARRIERS="data/PMBB_Exome/ZNF175/ZNF175_carriers_pLOF_JoesList_change51587911_IDs.txt.gz"
CASES_CONTROL="analysis/daniel/outputs/phase4/cases_control.txt"
MASTER_ANNOT="analysis/daniel/outputs/phase1/annot_genes_full_funcToInclude.txt"
PHASE3_LIGHT="analysis/daniel/outputs/phase3/light"

mkdir -p "$OUT_DIR" "$GENO_DIR" "$RESULTS_DIR" "$LOG_DIR"
TS=$(date +%Y%m%d_%H%M%S)
LOG="$LOG_DIR/run_${TS}.log"

log() { printf '[phase7] %s %s\n' "$(date +%H:%M:%S)" "$*"; }
fail() { log "ERROR: $*"; exit 1; }
exec > >(tee -a "$LOG") 2>&1

log "START Phase 7 — ZNF175 carrier deep-dive + second-hit hypothesis"
log "Project root:   $PROJECT_ROOT"
log "Output dir:     $OUT_DIR"

# Verify inputs
[[ -f "$DANIEL_CARRIERS" ]]    || fail "Daniel's 142-carrier list missing"
[[ -f "$CASES_CONTROL" ]]      || fail "cases_control.txt (Phase 4) missing"
[[ -f "$MASTER_ANNOT" ]]       || fail "Phase 1 master annotation missing"
[[ -d "$PHASE3_LIGHT" ]]       || fail "Phase 3 light dir missing"
[[ -x "$PLINK" ]]              || fail "plink 1.9 not found"

# ───────── Step 7.0 — Prep carrier keep-list ─────────
log ""
log "════════════════════════════════════════════════════════════"
log "Step 7.0 — Prep carrier keep-list"
log "════════════════════════════════════════════════════════════"
zcat "$DANIEL_CARRIERS" | tail -n +2 > "$OUT_DIR/carriers_ids.txt"
awk '{print $1, $1}' "$OUT_DIR/carriers_ids.txt" > "$OUT_DIR/carriers_keep.txt"
N_CARRIERS=$(wc -l < "$OUT_DIR/carriers_ids.txt")
log "  $N_CARRIERS ZNF175 pLoF carriers (Daniel's final curated list)"

# ───────── Step 7.1 — Cross-ref with cases_control ─────────
log ""
log "════════════════════════════════════════════════════════════"
log "Step 7.1 — Cross-reference carriers with cases_control.txt (HL status)"
log "════════════════════════════════════════════════════════════"

python3 <<PYEOF
import csv

# Load cases_control
status = {}
with open('$CASES_CONTROL') as f:
    next(f)
    for line in f:
        parts = line.rstrip().split('\t')
        if len(parts) >= 2:
            status[parts[0]] = parts[1]

# Load carriers
with open('$OUT_DIR/carriers_ids.txt') as f:
    carriers = [line.rstrip() for line in f if line.strip()]

# Tabulate
counts = {'1': 0, '0': 0, 'NA': 0, 'not_in_cases_control': 0}
with open('$OUT_DIR/carriers_with_status.tsv', 'w') as f:
    f.write('PMBB_ID\tSNHL\n')
    for cid in carriers:
        s = status.get(cid, 'not_in_cases_control')
        f.write(f"{cid}\t{s}\n")
        counts[s] = counts.get(s, 0) + 1

print(f"Total carriers: {len(carriers)}")
for k, v in counts.items():
    print(f"  SNHL={k}: {v}")
PYEOF

log "  written: $OUT_DIR/carriers_with_status.tsv"

# ───────── Step 7.2 — Extract carrier genotypes per chr ─────────
log ""
log "════════════════════════════════════════════════════════════"
log "Step 7.2 — Extract carrier genotypes at all 9,667 HL-gene variants"
log "════════════════════════════════════════════════════════════"
log "  Running plink --keep carriers --recode A for each chr 1-22"

T0=$(date +%s)
for i in {1..22}; do
    in_prefix="$PHASE3_LIGHT/allIndvs_chr${i}"
    out_prefix="$GENO_DIR/carriers_chr${i}"
    "$PLINK" \
        --bfile "$in_prefix" \
        --keep "$OUT_DIR/carriers_keep.txt" \
        --recode A \
        --out "$out_prefix" \
        > "$GENO_DIR/chr${i}.stdout" 2> "$GENO_DIR/chr${i}.stderr" \
        || fail "plink chr$i failed"
done
T1=$(date +%s)
log "  done in $((T1-T0))s — 22 .raw files in $GENO_DIR/"

# ───────── Step 7.3-7.6 — Python analysis ─────────
log ""
log "════════════════════════════════════════════════════════════"
log "Step 7.3-7.6 — Build per-carrier × per-gene matrix + second-hit analysis"
log "════════════════════════════════════════════════════════════"

python3 <<PYEOF
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact, mannwhitneyu
import pathlib

OUT = pathlib.Path('$OUT_DIR')
GENO = pathlib.Path('$GENO_DIR')

# Load master annotation: variant_ID → gene
# NOTE: Phase 1 annotation uses ID format "1:6425014:C:A" (colon-separated);
# plink's --recode A produces "1_6425014_C_A" (underscore). Build lookup in both.
print("Loading Phase 1 master annotation...")
annot = pd.read_csv('$MASTER_ANNOT', sep='\t', usecols=['ID', 'Gene.refGene'])
annot.columns = ['variant_id', 'gene']
print(f"  {len(annot)} variants annotated to {annot['gene'].nunique()} genes")
# Build lookup keyed by plink-style ID (colons → underscores)
variant_to_gene = {}
for vid, gene in zip(annot['variant_id'], annot['gene']):
    variant_to_gene[vid] = gene                    # original "1:6425014:C:A"
    variant_to_gene[vid.replace(':', '_')] = gene  # plink-style "1_6425014_C_A"

# Load carriers + status
carriers_status = pd.read_csv(OUT / 'carriers_with_status.tsv', sep='\t')
print(f"  Loaded {len(carriers_status)} carriers from carriers_with_status.tsv")

# Load and combine all 22 per-chr genotype files
print("\nLoading per-chr .raw files (genotype counts per carrier per variant)...")
geno_dfs = []
for i in range(1, 23):
    raw = GENO / f'carriers_chr{i}.raw'
    if not raw.exists():
        print(f"  chr{i}: missing!")
        continue
    df = pd.read_csv(raw, sep=' ')
    # First 6 cols: FID IID PAT MAT SEX PHENOTYPE; rest = variant columns
    meta_cols = ['FID', 'IID', 'PAT', 'MAT', 'SEX', 'PHENOTYPE']
    variant_cols = [c for c in df.columns if c not in meta_cols]
    # Strip trailing '_<allele>' from variant IDs (plink suffix)
    rename = {}
    for c in variant_cols:
        if '_' in c:
            base = c.rsplit('_', 1)[0]
            rename[c] = base
    df = df.rename(columns=rename)
    df = df[['IID'] + list(rename.values())]
    geno_dfs.append(df)

# Combine on IID (PMBB_ID)
print("  Combining across 22 chrs...")
geno = geno_dfs[0]
for df in geno_dfs[1:]:
    geno = geno.merge(df, on='IID', how='outer')
print(f"  combined matrix: {len(geno)} carriers × {len(geno.columns)-1} variant columns")

# Drop variants that are missing for all carriers (NaN)
geno_data = geno.set_index('IID')
geno_data = geno_data.dropna(axis=1, how='all')
print(f"  after dropping all-NaN variants: {geno_data.shape}")

# For each (carrier, variant), determine if they carry the alt allele (>0)
carrier_mask = (geno_data > 0).fillna(False)
print(f"  total carrier-variant pairs (any alt allele): {carrier_mask.sum().sum()}")

# Map variants → genes
variants_in_data = [c for c in geno_data.columns]
genes_for_variants = [variant_to_gene.get(v, 'UNKNOWN') for v in variants_in_data]
print(f"  variants mapped to genes: {sum(1 for g in genes_for_variants if g != 'UNKNOWN')}/{len(variants_in_data)}")

# Per-carrier × per-gene aggregate (count of HL-gene-variants carried)
per_carrier_per_gene = {}
for variant, gene in zip(variants_in_data, genes_for_variants):
    if gene == 'UNKNOWN':
        continue
    per_carrier_per_gene.setdefault(gene, pd.Series(0, index=carrier_mask.index, dtype=int))
    per_carrier_per_gene[gene] += carrier_mask[variant].astype(int)

# Build long-format table: carrier × gene where carrier has ≥1 variant
rows = []
for gene, counts in per_carrier_per_gene.items():
    for carrier_id, count in counts.items():
        if count > 0:
            rows.append({'PMBB_ID': carrier_id, 'gene': gene, 'n_variants': int(count)})
hits_long = pd.DataFrame(rows)
print(f"\nTotal (carrier × HL-gene-with-hits) pairs: {len(hits_long)}")

hits_long.to_csv(OUT / 'per_carrier_HL_gene_hits.tsv', sep='\t', index=False)
print(f"  written: {OUT}/per_carrier_HL_gene_hits.tsv")

# Per-carrier summary: N_HL_genes_with_hits, N_variants_total
if len(hits_long) > 0:
    per_carrier_summary = hits_long.groupby('PMBB_ID').agg(
        n_HL_genes_with_hits=('gene', 'nunique'),
        n_HL_variants_total=('n_variants', 'sum'),
    ).reset_index()
else:
    per_carrier_summary = pd.DataFrame(columns=['PMBB_ID', 'n_HL_genes_with_hits', 'n_HL_variants_total'])

# Add carriers with 0 hits
all_carriers = carriers_status['PMBB_ID'].tolist()
zero_hit = pd.DataFrame({
    'PMBB_ID': [c for c in all_carriers if c not in per_carrier_summary['PMBB_ID'].values],
    'n_HL_genes_with_hits': 0,
    'n_HL_variants_total': 0,
})
per_carrier_summary = pd.concat([per_carrier_summary, zero_hit], ignore_index=True)

# Merge with case/control status
per_carrier_summary = per_carrier_summary.merge(carriers_status, on='PMBB_ID', how='left')
per_carrier_summary = per_carrier_summary.sort_values('n_HL_variants_total', ascending=False)
per_carrier_summary.to_csv(OUT / 'per_carrier_second_hit_summary.tsv', sep='\t', index=False)
print(f"  written: {OUT}/per_carrier_second_hit_summary.tsv")

# Statistical comparison: case-carriers vs ctrl-carriers
print("\n" + "=" * 70)
print("SECOND-HIT HYPOTHESIS TEST")
print("=" * 70)

cases = per_carrier_summary[per_carrier_summary['SNHL'] == '1']
ctrls = per_carrier_summary[per_carrier_summary['SNHL'] == '0']
na = per_carrier_summary[~per_carrier_summary['SNHL'].isin(['0', '1'])]

print(f"\nCarrier breakdown by SNHL status:")
print(f"  SNHL=1 (HL cases):    {len(cases)} carriers")
print(f"  SNHL=0 (HL controls): {len(ctrls)} carriers")
print(f"  SNHL=NA / not_in_cc:  {len(na)} carriers")

print(f"\nN_HL_genes_with_hits distribution:")
print(f"  Cases:    mean={cases['n_HL_genes_with_hits'].mean():.2f}, median={cases['n_HL_genes_with_hits'].median():.1f}, max={cases['n_HL_genes_with_hits'].max()}")
print(f"  Controls: mean={ctrls['n_HL_genes_with_hits'].mean():.2f}, median={ctrls['n_HL_genes_with_hits'].median():.1f}, max={ctrls['n_HL_genes_with_hits'].max()}")

# Fisher's exact: carriers with ≥1 second-hit by case/control
def has_hit(df):
    return (df['n_HL_genes_with_hits'] > 0).sum()
def no_hit(df):
    return (df['n_HL_genes_with_hits'] == 0).sum()

table = [
    [has_hit(cases), no_hit(cases)],
    [has_hit(ctrls), no_hit(ctrls)],
]
print(f"\nFisher's exact: ≥1 second-hit by case/control")
print(f"  cases:    {table[0][0]} with hit, {table[0][1]} without")
print(f"  controls: {table[1][0]} with hit, {table[1][1]} without")
odds_ratio, p_fisher = fisher_exact(table)
print(f"  OR={odds_ratio:.3f}, p={p_fisher:.4f}")

# Mann-Whitney: distribution of N_HL_genes_with_hits
if len(cases) > 0 and len(ctrls) > 0:
    u_stat, p_mw = mannwhitneyu(
        cases['n_HL_genes_with_hits'],
        ctrls['n_HL_genes_with_hits'],
        alternative='greater'  # cases > ctrls hypothesis
    )
    print(f"\nMann-Whitney U (one-sided, cases > controls in N_HL_gene hits):")
    print(f"  U={u_stat:.1f}, p={p_mw:.4f}")

# Top carriers by N_HL_genes_with_hits
print(f"\nTop 20 carriers by N_HL_genes_with_hits:")
print(per_carrier_summary.head(20).to_string(index=False))
PYEOF

log ""
log "DONE Phase 7 — ZNF175 carrier deep-dive complete"
log "  Key outputs in $OUT_DIR/:"
log "    carriers_with_status.tsv          (142 carriers + SNHL status)"
log "    per_carrier_HL_gene_hits.tsv      (long: carrier × HL gene × N variants)"
log "    per_carrier_second_hit_summary.tsv (per-carrier summary stats + statistical tests)"
