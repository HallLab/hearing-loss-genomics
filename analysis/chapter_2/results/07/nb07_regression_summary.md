# NB 07 — ZNF175 burden x tinnitus, Park pipeline on both cohorts

Model: tinnitus(rule-of-2, ICD 388.3x/H93.1x) ~ carrier + age + age2 + sex + PC1-10. Same pipeline v1 & v2.

## Cohort sizes (after rule-of-2 exclusion)
- v1 (Park, Freeze One): N=9161, carriers=26, tinnitus cases=131
- v2 (Hui, Release 2.0): N=42614, carriers=69, tinnitus cases=841

## Carrier-vs-tinnitus enrichment (2x2 / Fisher)
- v1: carriers 4/26 tinnitus vs non 127/9135 | OR=12.90 p=4.68e-04
- v2: carriers 4/69 tinnitus vs non 837/42545 | OR=3.07 p=4.76e-02

## Adjusted logistic (pooled, PC1-10)
- v1: OR=14.58 p=4.79e-06 (N=9161, cases=131)
- v2: OR=3.53 p=1.52e-02 (N=42605, cases=841)

## Stratified EUR/AFR + IVW meta
- v1: meta OR=33.71 p=2.05e-08
- v2: meta OR=4.59 p=3.74e-03

## Verdict (to interpret)
Same Park pipeline on both: if the carrier-tinnitus signal is present in v1 (11k) but absent in v2 (44k),
the loss is winner's curse / dilution (real null at scale), NOT a pipeline artifact. Small carrier/case counts
-> Fisher is the robust readout; logistic/meta are sensitivity checks.

## Outputs
- d1/d2 analysis frames saved as v1_analysis.csv / v2_analysis.csv
