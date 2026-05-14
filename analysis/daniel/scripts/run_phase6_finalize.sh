#!/usr/bin/env bash
# run_phase6_finalize.sh — concatenate per-chr biobin outputs + BH correction + top hits.
# Runs AFTER the 22-task job array completes (via LSF -w dependency).

set -euo pipefail

PROJECT_ROOT="/project/hall/analysis/hearing-loss-genomics"
cd "$PROJECT_ROOT"

if [[ -z "${VIRTUAL_ENV:-}" ]]; then
    # shellcheck disable=SC1091
    source venv/bin/activate
fi

BIOBIN_DIR="analysis/daniel/outputs/phase6/biobin/HL_needAud"
OUT_DIR="analysis/daniel/outputs/phase6/results"
LOG_DIR="analysis/daniel/logs/phase6"

mkdir -p "$OUT_DIR" "$LOG_DIR"
TS=$(date +%Y%m%d_%H%M%S)
LOG="$LOG_DIR/finalize_${TS}.log"

log() { printf '[phase6:final] %s %s\n' "$(date +%H:%M:%S)" "$*"; }
fail() { log "ERROR: $*"; exit 1; }
exec > >(tee -a "$LOG") 2>&1

log "START Phase 6 finalize"

# Verify all 22 chr bins.csv exist
log "Verifying all 22 chr bins.csv files..."
missing=0
for i in {1..22}; do
    f="$BIOBIN_DIR/HL_needAud_chr${i}-bins.csv"
    [[ -f "$f" ]] || { log "  ✗ missing: $f"; missing=$((missing+1)); }
done
[[ "$missing" -eq 0 ]] || fail "$missing chr outputs missing — array tasks failed?"
log "  all 22 chr bins.csv present"

# Parse + concatenate + BH + report
python3 <<EOF
import csv
import pathlib

BIOBIN_DIR = pathlib.Path('$BIOBIN_DIR')
OUT_DIR = pathlib.Path('$OUT_DIR')

# Extract (gene, chr, n_variants, p_value) for every bin across all 22 chrs
all_bins = []
for chr_n in range(1, 23):
    f = BIOBIN_DIR / f'HL_needAud_chr{chr_n}-bins.csv'
    with open(f) as fh:
        reader = csv.reader(fh)
        # First 10 metadata rows
        rows = []
        for i, row in enumerate(reader):
            rows.append(row)
            if i >= 9:
                break

    # rows[0] = ID/gene names; rows[1] = Total Variants; rows[8] = logistic p-value
    genes = rows[0][2:]   # skip "ID" label + empty cell
    n_vars = rows[1][2:]  # skip label + nan
    pvals = rows[8][2:]   # skip label + nan

    for g, nv, pv in zip(genes, n_vars, pvals):
        if not g or g == 'nan':
            continue
        try:
            p = float(pv)
            n = int(float(nv))
        except (ValueError, TypeError):
            continue
        all_bins.append({'gene': g, 'chr': chr_n, 'n_variants': n, 'p': p})

print(f'Total bins across 22 chrs: {len(all_bins)}')

# Apply Benjamini-Hochberg correction
all_bins.sort(key=lambda b: b['p'])
m = len(all_bins)
for rank, b in enumerate(all_bins, start=1):
    b['rank'] = rank
    b['p_bh'] = b['p'] * m / rank
# Monotone-decreasing cumulative min from the top
running_min = 1.0
for b in reversed(all_bins):
    running_min = min(running_min, b['p_bh'])
    b['p_bh'] = running_min

# Save all results
out_all = OUT_DIR / 'all_chrom_meta_HL_needAud.txt'
with open(out_all, 'w') as fh:
    fh.write('gene\\tchr\\tn_variants\\tp_logistic\\tp_FDR\\trank\\n')
    for b in all_bins:
        fh.write(f"{b['gene']}\\t{b['chr']}\\t{b['n_variants']}\\t{b['p']:.6e}\\t{b['p_bh']:.6e}\\t{b['rank']}\\n")
print(f'Wrote {out_all} ({m} bins)')

# Top 30 hits — show LOC* and non-LOC separately
def is_loc(g):
    return g.startswith('LOC') or g.startswith('LINC')

# All top 30
print()
print('=' * 100)
print('TOP 30 — ALL (incl. LOC*/LINC artifacts from newer LOKI)')
print('=' * 100)
print(f"{'rank':>4}  {'gene':25s}  {'chr':>3}  {'n_var':>5}  {'p_logistic':>12}  {'p_FDR':>12}")
for b in all_bins[:30]:
    print(f"  {b['rank']:>3}  {b['gene']:25s}  {b['chr']:>3}  {b['n_variants']:>5}  {b['p']:>12.4e}  {b['p_bh']:>12.4e}")

print()
print('=' * 100)
print('TOP 30 — non-LOC/LINC (real characterized genes)')
print('=' * 100)
real = [b for b in all_bins if not is_loc(b['gene'])]
print(f"{'rank':>4}  {'gene':25s}  {'chr':>3}  {'n_var':>5}  {'p_logistic':>12}  {'p_FDR':>12}")
for b in real[:30]:
    # Find original rank
    orig_rank = b['rank']
    print(f"  {orig_rank:>3}  {b['gene']:25s}  {b['chr']:>3}  {b['n_variants']:>5}  {b['p']:>12.4e}  {b['p_bh']:>12.4e}")

# Look for paper's key genes (known HL + novel)
print()
print('=' * 100)
print('KEY GENES — paper hits look-up (top of ranking we expect)')
print('=' * 100)
print('Replicated known HL (paper Fig 2): ESRRB, TCOF1')
print('Novel discoveries (paper Fig 3):   ZNF175, COL5A2, HMMR, NNT, RAPGEF3')
print()
gene_idx = {b['gene']: b for b in all_bins}
for gene in ['ESRRB', 'TCOF1', 'ZNF175', 'COL5A2', 'HMMR', 'NNT', 'RAPGEF3', 'SCD5', 'GIPC3', 'OTOF']:
    if gene in gene_idx:
        b = gene_idx[gene]
        print(f"  {gene:12s}  chr={b['chr']:>2}  rank=#{b['rank']:<5}  p={b['p']:.4e}  p_FDR={b['p_bh']:.4e}  n_var={b['n_variants']}")
    else:
        print(f"  {gene:12s}  NOT FOUND")

# Genome-wide significance (Bonferroni for ~20k genes ≈ 2.5e-6 cutoff)
bonf_cut = 0.05 / m
n_bonf = sum(1 for b in all_bins if b['p'] < bonf_cut)
n_bh_05 = sum(1 for b in all_bins if b['p_bh'] < 0.05)
print()
print(f'Bonferroni threshold (0.05/{m}): {bonf_cut:.2e}')
print(f'  Bins below Bonferroni: {n_bonf}')
print(f'  Bins with BH-corrected p < 0.05: {n_bh_05}')
EOF

log ""
log "DONE Phase 6 finalize"
log "  results table: $OUT_DIR/all_chrom_meta_HL_needAud.txt"
