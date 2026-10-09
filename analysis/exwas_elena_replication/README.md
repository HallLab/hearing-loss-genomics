# Common-variant ExWAS — broad hearing loss (replication)

**Status:** all phases complete · 6,752 cases · 0 significant variants
**Result:** replicates the delivered run at r = 0.98–0.99; the summer
poster's six significant results fail the poster's own call-rate filter
**Release:** PMBB-Release-2026-4.0 · **Test:** SAIGE single-variant
**Scope:** variants at **MAF >= 0.01** — *not* the rare-variant gene burden, which is
`../gene_burden_elena_replication/`

Phenotype: any hearing loss, Elena's phenotype as corrected in the replication.

Read [PREMISES.md](PREMISES.md) first. It is pre-registered: the significance threshold
and the two declared target genes were fixed before step 2 ran.

## Inherited, not rebuilt

| what | from |
|---|---|
| cohort | `../gene_burden_elena_replication/phase_1/` |
| covariates | `../gene_burden_elena_replication/phase_3/results/covariates/` |
| step 1 null models | `../gene_burden_elena_replication/phase_4/results/step1/` |
| genotypes | `analysis/elena/rarevariant_geneburden/plink_deduplicated/` |

This arm therefore contains only phases 2, 4 and 5. Phases 1 and 3 are deliberately absent
rather than duplicated — a copied cohort definition is a cohort definition that drifts.

## The 2x2

|  | broad phenotype | bilateral sensorineural |
|---|---|---|
| **gene burden** (MAF <= 0.01) | `gene_burden_elena_replication/` | `gene_burden_pmbb_bilateral/` |
| **ExWAS** (MAF >= 0.01) | `exwas_elena_replication/` | `exwas_pmbb_bilateral/` |
