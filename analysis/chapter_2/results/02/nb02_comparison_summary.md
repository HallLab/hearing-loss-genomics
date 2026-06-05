# NB 02 — ZNF175 MAF comparison: v1 (Park) vs v2 (Hui)

Window: chr19:51,571,283-51,592,510 (NCBI NC_000019.10, GRCh38, + strand).
Pipeline (symmetric, scientist-independent): bcftools view -r (.tbi random access)
-> norm -f <GRCh38> -a -m-any (left-align + atomize + split) -> plink2 --freq --geno-counts.

## Counts (whole-cohort)
- v1 (Freeze One, 11,451):       161 variants (142 rare <0.1%, 8 common >=1%)
- v2 (Release 2020 v2.0, 43,731): 380 variants (356 rare <0.1%)
- UNION: 386  ->  shared 155 | v1-only 6 | v2-only 225
- Pearson corr(MAF_v1, MAF_v2) on shared: 0.9449

## Scientific reading
The ZNF175 variant SUBSTRATE is stable between v1 and v2. The variants underpinning the original
signal are still present in v2 with correlated MAF (r=0.94). The apparent differences are explained:
- v2-only (225): EXPECTED -- 4x more samples reveal more rare variants.
- v1-only (6): NOT real disappearances -- 2 are '*' spanning-deletion markers
  (representation artifacts) + 4 ultra-rare private singletons (1-2 carriers).
- the 8 common (>=1%) variants are identical across freezes (stable backbone).

=> The 11K->45K signal loss is NOT explained by variant-substrate instability. It points instead to
   (a) winner's curse (a small-case signal regressing as N grows), or
   (b) association-pipeline differences flagged in the 2026-06-03 meeting (replications shifting to
       missense vs the pLOF discovery; possibly different phenotype-release versions).
CAVEAT: this step compares MAF; it does NOT re-test the association. It rules out ONE technical
        explanation (substrate). A definitive verdict needs the association re-run with matched pipeline.

## Outputs
- cmp_union_all_variants.csv  (UNION -- all 386 variants, blanks where absent in a cohort)
- cmp_shared_maf.csv, cmp_v1_only.csv, cmp_shared_maf_scatter.png, per_cohort_overview.csv

## Next steps
1. Annotate the union table with BF4 -> pLOF status + AlphaMissense scores, to focus on the variants
   that actually enter the burden (rare pLOF / deleterious missense).
2. Stratify MAF by ancestry (removes the remaining whole-cohort confounder: cohort size & ancestry mix).
3. (definitive) re-run the ZNF175-tinnitus association on v1 and v2 with a matched pipeline.
