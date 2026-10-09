# Phase 4 — merge report

66 of 66 cells merged. Checks in `02_merge_results.py`.

| cohort | tested | plink2 ceiling | shortfall | secondary bar | min p | hits at 5e-8 | hits at secondary |
|---|---:|---:|---:|---:|---:|---:|---:|
| combined | 70,350 | 70,374 | 24 | 7.107e-07 | 3.079e-05 | 0 | 0 |
| EUR | 74,566 | 74,585 | 19 | 6.705e-07 | 9.473e-06 | 0 | 0 |
| AFR | 118,049 | 118,083 | 34 | 4.236e-07 | 1.074e-06 | 0 | 0 |

Shortfall is the plink2 census minus what SAIGE tested. plink2 divides by observed
genotypes, SAIGE by the whole cohort, so the census is an upper bound. See
`phase_2/results/FINDINGS.md`.

