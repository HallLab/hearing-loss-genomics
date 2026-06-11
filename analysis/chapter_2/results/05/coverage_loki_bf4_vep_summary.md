# Coverage of the ZNF175 variant list — LOKI vs BF4 vs VEP

386 union variants (NB 02). Per-variant flags in `coverage_loki_bf4_vep.csv`. Date 2026-06-08.

## Coverage

| Source | Covered | How it matches | What it is |
|---|---|---|---|
| **VEP** | **384 / 386 (99%)** | by coordinate, **allele-specific** | computes consequence de novo for ANY variant (2 missing = `*` spanning-deletion alleles) |
| **LOKI** (`loki-20230816.db`, GRCh38) | **289 / 386 (75%)** | by **position only** (`snp_locus` has rs + chr + pos, **no ref/alt**) | knowledge base (dbSNP + genes + pathways) used by biobin/biofilter for gene/region boundaries |
| **BF4** (as loaded) | **97 / 386 (25%)** | allele-specific (gnomAD) | annotation DB; loaded **genome-only + AC ≥ 5** → drops rare (see `bf4_coverage_root_cause.md`) |

## Coverage by rarity (max cohort MAF)

| MAF bin | LOKI | BF4 | VEP |
|---|---|---|---|
| < 0.01% | 0.71 | 0.11 | 1.00 |
| 0.01–0.1% | 0.95 | 0.88 | 0.98 |
| 0.1–1% | 0.94 | 1.00 | 1.00 |
| ≥ 1% | 1.00 | 1.00 | 1.00 |

## Reading (⚠️ not apples-to-apples)

- **VEP is the most complete** for our task: by-coordinate, allele-specific → ~100% at every rarity, including private variants.
- **LOKI's 75% is position-based, so it overcounts** vs allele-specific matching: `snp_locus` records *positions* of known dbSNP variants (no ref/alt), so a **private allele sitting at a known dbSNP site counts as "covered"** even though that exact allele isn't in dbSNP. LOKI is a **knowledge base** (dbSNP/genes/pathways for biobin region collapsing), not a per-variant functional annotator — so this is a different question than BF4/VEP answer.
- **BF4's 25% is the genome-only + AC≥5 load artifact** — fixable by loading the joint exome+genome file (`bf4_coverage_root_cause.md`); not an inherent limit.

**Bottom line for annotating arbitrary/private variants:** VEP (by-coordinate, allele-specific) > LOKI (position-based dbSNP presence) > BF4-as-loaded (genome-only AC≥5). BF4 is fixable; LOKI answers a looser, position-level question by design.

*Outputs: `coverage_loki_bf4_vep.csv` (per-variant), this summary. LOKI query: `snp_locus WHERE chr=19 AND pos BETWEEN 51571283 AND 51592510` → 4,722 dbSNP positions in the locus.*
