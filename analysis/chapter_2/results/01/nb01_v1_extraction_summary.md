# NB 01 — v1 (Freeze One) ZNF175 extraction (independent, NCBI coords)

- Locus: chr19:51,571,283-51,592,510 (NCBI NC_000019.10, GRCh38, + strand)
- Source: Freeze One GL biallelic PLINK (PMBB raw); tool: plink2
- Variants in window (no MAF/pLOF filter): 151
- Observed position range: chr19:51,573,255-51,588,539
- MAF<0.1% (rare): 132
- with >=1 carrier: 151

Output table: results/v1_znf175_variants.csv

Next:
- annotate independently (VEP/ANNOVAR) to flag pLOF / missense (REVEL)
- ancestry-stratified MAF (needs v1 ancestry labels)
- NB 02: harmonize & compare with v2 by chr:pos:ref:alt
