# BF4 low-coverage — root cause (refined, 2026-06-08)

NB 03/05 found Biofilter 4 annotated only **97/384 (25%)** of our ZNF175 variants (75% `not_found`). We first attributed this to a **gnomAD AC ≥ 5** cutoff on the loaded DB.

## Root cause
The loaded BF4 database carried **genome-only** gnomAD frequencies, with an **AC ≥ 5** filter applied at load time. Because gnomAD **genomes** have far fewer samples than **exomes**, the genome-only allele count is low — so an AC ≥ 5 cutoff on *genome-only* counts drops variants whose evidence is mostly in the **exome** data.

## Example (gnomAD v4, GRCh38) — `19-51573333-C-A`
| | Exomes | Genomes | **Total (joint)** |
|---|---|---|---|
| Allele Count (AC) | 4 | 2 | **6** |
| Allele Number (AN) | 1,460,406 | 152,108 | 1,612,514 |

- **Genome-only AC = 2 < 5** → **excluded** by the BF4 load filter → comes back `not_found`.
- **Joint (exome+genome) Total AC = 6 ≥ 5** → would be **kept**.

## Fix
Load the **joint exome+genome** gnomAD file (use the **Total** AC), not the genome-only file. This recovers most of the `not_found` variants (the rare ones whose count lives in the exome data — exactly our burden candidates).

## Implication for the NB 05 verdict
BF4's low coverage here is a **DB-load configuration issue** (genome-only + AC ≥ 5), **fully fixable** — not an inherent BF4 limitation. This sharpens the NB 05 conclusion: VEP's coverage advantage over BF4 was **partly a config artifact** (genome-only load), though VEP's by-coordinate annotation and AlphaMissense remain genuine advantages for truly novel/private variants. A fair BF4-vs-VEP coverage comparison should re-run BF4 with the joint file loaded.

*Related: NB 03 (`03_annotate_union_biofilter.ipynb`), NB 05 (`05_reconcile_vep_vs_bf4.ipynb`), memory `reference-biofilter4-lpc`.*
