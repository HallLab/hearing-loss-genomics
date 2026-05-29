#!/usr/bin/env bash
# run_phase8.sh — Chapter 1 · Phase 7 (degree-HL burden, light-mode Etapa A)
#
# Goal: replicate Tables 3 & 4 of Hui et al. 2023 (PLOS Genetics) using Daniel's
# preserved supplementary table for the linear-regression-on-degreeHL burden
# analysis (20-PC covariates). Light-mode = no biobin re-run; just apply paper
# filters + BH-FDR to Daniel's preserved STable and compare to published values.
#
# Inputs:
#   - data/PMBB_Exome/allGenes/20PCs/degreeHL/results/allChr_STable_degHL.txt.gz
#         (Gene, Beta, SE, P, Carriers, Case_carriers — 18,547 genes)
#   - data/PMBB_Exome/all_genes_including_ShadisList.txt.gz (179 known HL genes)
#
# Outputs (in analysis/daniel/outputs/phase8/light_mode/):
#   - table3_known_hl_genes.tsv      — known HL genes, BH-FDR computed
#   - table4_novel_genes.tsv         — non-HL genes with >25 case carriers, BH-FDR
#   - paper_comparison.tsv           — side-by-side Daniel preserved vs paper Table 3+4
#   - run.log

set -euo pipefail

PROJECT_ROOT="/project/hall/analysis/hearing-loss-genomics"
cd "$PROJECT_ROOT"

source venv/bin/activate

OUT_DIR="$PROJECT_ROOT/analysis/daniel/outputs/phase8/light_mode"
LOG_DIR="$PROJECT_ROOT/analysis/daniel/logs/phase8"
TS=$(date +%Y%m%d_%H%M%S)
LOG="$LOG_DIR/run_${TS}.log"
mkdir -p "$OUT_DIR" "$LOG_DIR"

exec > >(tee -a "$LOG") 2>&1
echo "=== Phase 8 (Ch1 P7 light-mode Etapa A) start: $(date) ==="
echo "Project root: $PROJECT_ROOT"
echo "Python: $(which python3) — $(python3 --version)"
echo

python3 "$PROJECT_ROOT/analysis/daniel/scripts/pmbb_exome/degree_hl_burden_lightmode.py" \
    --stable     "$PROJECT_ROOT/data/PMBB_Exome/allGenes/20PCs/degreeHL/results/allChr_STable_degHL.txt.gz" \
    --hl-genes   "$PROJECT_ROOT/data/PMBB_Exome/all_genes_including_ShadisList.txt.gz" \
    --out-dir    "$OUT_DIR"

echo
echo "=== Phase 8 light-mode done: $(date) ==="
