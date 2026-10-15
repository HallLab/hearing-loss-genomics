# Common-variant ExWAS — bilateral sensorineural hearing loss

**Status:** all phases complete · 3,164 cases · 0 significant variants
**Result:** the one promising signal (`BCAM`, chr19) was shown by its own
pre-registered power control to be selection, not subtype specificity
**Release:** PMBB-Release-2026-4.0 · **Test:** SAIGE single-variant
**Scope:** variants at **MAF >= 0.01** — *not* the rare-variant gene burden, which is
`../gene_burden_pmbb_bilateral/`

Phenotype: bilateral sensorineural only (reading C: `.2` and `.8` on the same date, twice).

Read [PREMISES.md](PREMISES.md) first. It is pre-registered: the significance threshold
and the two declared target genes were fixed before step 2 ran.

## Inherited, not rebuilt

| what | from |
|---|---|
| cohort | `../gene_burden_pmbb_bilateral/phase_1/` |
| covariates | `../gene_burden_pmbb_bilateral/phase_3/results/covariates/` |
| step 1 null models | `../gene_burden_pmbb_bilateral/phase_4/results/step1/` |
| genotypes | `analysis/elena/rarevariant_geneburden/plink_deduplicated/` |

This arm therefore contains only phases 2, 4 and 5. Phases 1 and 3 are deliberately absent
rather than duplicated — a copied cohort definition is a cohort definition that drifts.

## The 2x2

|  | broad phenotype | bilateral sensorineural |
|---|---|---|
| **gene burden** (MAF <= 0.01) | `gene_burden_elena_replication/` | `gene_burden_pmbb_bilateral/` |
| **ExWAS** (MAF >= 0.01) | `exwas_elena_replication/` | `exwas_pmbb_bilateral/` |
