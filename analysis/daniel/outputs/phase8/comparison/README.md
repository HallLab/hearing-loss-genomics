# Comprehensive p-values CSV

Two CSVs covering p-values from Phases 5, 6, 7 across all genes tested in Chapter 1.

## Files

| File | Rows | Use when |
|---|---:|---|
| **`all_real_genes_pvalues.csv`** ⭐ recommended | 18,555 | Normal analysis — only real protein-coding genes |
| `all_genes_pvalues.csv` | 45,558 | Full output — includes LOC/LINC/MIR artifacts from newer LOKI |

The two files are identical in columns; the only difference is whether LOC/LINC/MIR pseudogenes (added by the newer LOKI database to our Phase 6 run) are included.

## Columns

| Column | Description | Filled for |
|---|---|---|
| `Gene` | Gene symbol | All rows |
| `Is_real_gene` | "yes" for real genes, "no" for LOC/LINC/MIR artifacts | All |
| `In_HL_gene_set` | "yes" if in the 179 curated known HL genes | All |
| `In_paper_table_3_or_4` | "yes" if cited in Hui et al. 2023 Tables 3 or 4 | Only TCOF1, ESRRB, COL5A1, HMMR, RAPGEF3, NNT |
| `Phase5_HL_4PC_binary_p_ours` | Our Phase 5 p-value (biobin logistic, 4PC, HL genes only) | HL genes only (171) |
| `Phase5_HL_4PC_binary_p_daniel` | Daniel's Phase 5 p-value (same params) | HL genes only |
| `Phase6_WES_4PC_binary_p_ours` | Our Phase 6 p-value (biobin logistic, 4PC, exome-wide) | All genes tested in our Phase 6 (~44,879) |
| `Phase6_WES_4PC_binary_p_daniel` | Daniel's Phase 6 p-value (same params) | ~18,546 genes |
| `Phase7_WES_20PC_linear_p_daniel_stable` | Phase 7 p-value (R lm, 20PC, degree-HL) — light-mode, from Daniel's STable | 18,546 genes |
| `Phase7_WES_20PC_linear_beta_daniel_stable` | β coefficient | Same |
| `Phase7_carriers_daniel_stable` | Total carriers in Phase 7 STable | Same |
| `Phase7_case_carriers_daniel_stable` | HL-case carriers in Phase 7 STable | Same |
| `Paper_p_published` | p-value as published in Hui et al. 2023 | 6 paper genes only |
| `Paper_beta_published` | β as published | 6 paper genes only |
| `Paper_fdr_published` | FDR as published | 6 paper genes only |
| `Joe_PMBB_v1_p` | Joe Park's PMBB v1 p-value | **Empty — we don't have his per-gene results** |
| `Joe_PMBB_v1_notes` | Notes if/when Joe's data arrives | Empty |

## Quick examples (in Python)

```python
import pandas as pd
df = pd.read_csv('all_real_genes_pvalues.csv')

# All HL genes with Phase 5 data
hl = df[df['In_HL_gene_set'] == 'yes']

# Paper-cited genes only
paper = df[df['In_paper_table_3_or_4'] == 'yes']

# Top 20 by Phase 7 p-value (degree-HL linear, the paper's primary test)
top20 = df.dropna(subset=['Phase7_WES_20PC_linear_p_daniel_stable']) \
          .assign(p7=lambda d: d['Phase7_WES_20PC_linear_p_daniel_stable'].astype(float)) \
          .nsmallest(20, 'p7')

# Compare ours vs Daniel for Phase 6 on HL genes
hl_p6 = hl.dropna(subset=['Phase6_WES_4PC_binary_p_ours', 'Phase6_WES_4PC_binary_p_daniel']) \
          .assign(ours=lambda d: d['Phase6_WES_4PC_binary_p_ours'].astype(float),
                  daniel=lambda d: d['Phase6_WES_4PC_binary_p_daniel'].astype(float)) \
          .assign(ratio=lambda d: d['ours'] / d['daniel'])
```

## Quick examples (in Excel)

1. Open `all_real_genes_pvalues.csv` (use `all_genes_pvalues.csv` if you want pseudogenes too)
2. AutoFilter (Data → Filter)
3. Common filters:
   - `In_HL_gene_set = yes` → 173 HL genes
   - `In_paper_table_3_or_4 = yes` → the 6 paper-cited genes
   - Sort by `Phase7_WES_20PC_linear_p_daniel_stable` ascending → paper's primary test ranking

## What's NOT in this CSV

### Joe Park's PMBB v1 p-values
We don't have these. Joe ran the original ZNF175 burden analysis in PMBB v1 (~11k exomes) in 2018-2019. We only know the headline finding (ZNF175 pLOF burden with 8 driving cases) from Doug's kickoff narrative. **If Joe sends his per-gene results to Nikki**, we can populate the `Joe_PMBB_v1_p` column.

### Our independent Phase 7 run
Phase 7 was done in light mode — we adopted Daniel's preserved STable directly. The `Phase7_*_daniel_stable` columns are our Phase 7 values (= Daniel's STable values). If we ever do a heavy-mode Phase 7 (biobin → R lm() from scratch), we'd add `Phase7_*_ours` columns.

### Phase 6 with 20 PCs
Daniel ran a 20PC binary version as a sanity check (in `data/PMBB_Exome/allGenes/20PCs/binary/`). We did not include this column to keep the CSV focused. If needed, the source file is documented in [`chapter1_authoritative_pvalues.md`](../../../../results/chapter1_paper_replication/chapter1_authoritative_pvalues.md).

## How to regenerate

```bash
source venv/bin/activate
python3 analysis/daniel/scripts/pmbb_exome/build_all_genes_pvalues_csv.py \
    --project-root /project/hall/analysis/hearing-loss-genomics \
    --out          analysis/daniel/outputs/phase8/comparison/all_genes_pvalues.csv
```

Script: [`analysis/daniel/scripts/pmbb_exome/build_all_genes_pvalues_csv.py`](../../../scripts/pmbb_exome/build_all_genes_pvalues_csv.py).

## Authoritative narrative

For the narrative interpretation (why p-values differ across phases / between us and Daniel / between Daniel and Paper), see:

- [`results/chapter1_paper_replication/chapter1_authoritative_pvalues.md`](../../../../results/chapter1_paper_replication/chapter1_authoritative_pvalues.md) — single source of truth
- [`docs/andre/phase5_hl_genes_byte_equivalent_explained.md`](../../../../docs/andre/phase5_hl_genes_byte_equivalent_explained.md) — Phase 5 byte-equivalence
- [`docs/andre/phase6_hl_genes_drift_explained.md`](../../../../docs/andre/phase6_hl_genes_drift_explained.md) — Phase 6 LOKI drift
- [`docs/andre/phase7_wes_drift_explained.md`](../../../../docs/andre/phase7_wes_drift_explained.md) — Phase 7 iteration drift
