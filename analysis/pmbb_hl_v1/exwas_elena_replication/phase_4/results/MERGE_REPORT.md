# Phase 4 — merge report

66 of 66 cells merged. Checks in `02_merge_results.py`.

| cohort | tested | plink2 ceiling | shortfall | secondary bar | min p | hits at 5e-8 | hits at secondary |
|---|---:|---:|---:|---:|---:|---:|---:|
| combined | 70,366 | 70,383 | 17 | 7.106e-07 | 3.016e-06 | 0 | 0 |
| EUR | 74,536 | 74,553 | 17 | 6.708e-07 | 2.116e-05 | 0 | 0 |
| AFR | 117,953 | 117,996 | 43 | 4.239e-07 | 4.109e-06 | 0 | 0 |

Shortfall is the plink2 census minus what SAIGE tested. plink2 divides by observed
genotypes, SAIGE by the whole cohort, so the census is an upper bound. See
`phase_2/results/FINDINGS.md`.

