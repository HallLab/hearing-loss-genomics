# Step Zero — Summary

## Cohorts
- v1 (Park / Freeze One): N=11,451, GRCh38, /static/PMBB/PMBB_Freeze17/genotype/exome/
- v2 (Hui / Release 2020 v2.0): N=43,731, GRCh38, data/pmbb_v2/Exome/pVCF/GL_by_chrom/

## ZNF175 locus (GRCh38)
- variant span: chr19:51,573,230-51,588,560  (padded tabix region in znf175_locus_grch38.txt)

## v2 variant sets already extracted by Daniel (counts differ by filter)
                        set  n_variants
        standard (allIndvs)          23
                      Joe's           9
 maf.001_noRels_keepHLcases          22
inclMultAllelic_noMAFfilter          26

## Per-sample concordance feasibility
- naive ID intersection (v1 vs v2): 0
- verdict: BLOCKED — IDs in different namespaces (UPENN_* vs PMBB*); need an ID crosswalk
- fallback available now: aggregate MAF comparison, stratified by ancestry

## Next steps
- NB 01: extract ZNF175 from v1 (Freeze One) via tabix region on the GRCh38 pVCF
- NB 02: harmonize variants (normalize, match by chr:pos:ref:alt) and compare MAF v1 vs v2 (by ancestry)
- Side-quest: locate a UPENN_* <-> PMBB ID crosswalk to unlock per-sample concordance
