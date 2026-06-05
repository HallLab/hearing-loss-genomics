# NB 06 — ZNF175 qualifying variants + carriers (Park filter, both cohorts)

Filter (Park discovery): pLOF (VEP --most_severe HIGH-impact LoF terms) + cohort MAF <= 0.1%.

## Qualifying variants
- v1 (Freeze One, 11,451): 10
- v2 (Release 2020 v2.0, 43,731): 24
- shared: 8

## Carriers (>=1 qualifying alt allele)
- v1: 34 / 11,451  (0.30%)
- v2: 72 / 43,731  (0.16%)

## Reading
This is the burden INPUT, identical pipeline both cohorts. Carrier counts + qualifying-set overlap
set up the regression (NB 07: tinnitus ~ carrier + age/age2/sex/PCs, EUR/AFR-stratified + meta).
Hui-definition counts (pLOF + missense REVEL>0.6) also computed for an optional comparison pass.

## Outputs
- znf175_qualified_variants.csv, qual_v{1,2}_park_ids.txt, carriers_v{1,2}.csv
