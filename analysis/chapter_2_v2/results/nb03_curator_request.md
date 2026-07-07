# PMBB curator request — ZNF175 crosswalk gap (v1 / Freeze One)

We are replicating Park et al. 2021 (ZNF175 pLOF → tinnitus) in Freeze One WES. Our qualifying-carrier × tinnitus
count is **4 carrier-cases**; a project-internal number is **8**. We have traced the gap to a **GENO_ID ↔ PT_ID
linkage gap** and need the master crosswalk to close it. Two lists (attached):

## List A — 7 WES carriers with NO phenotype linkage
GENO_IDs present in the WES `.fam` carrying a rare qualifying pLOF ZNF175 variant, but ABSENT from
`PMBB_Geno_Demographics_Deidentified_012020.csv` (no PT_ID → tinnitus status unknown).
File: `nb03_list_A_wes_carriers_no_phenotype.csv`.
Ask: the PT_ID for each, or confirmation they are unconsented / QC-dropped.

## List B — 77 tinnitus cases with NO genotype
PT_IDs meeting tinnitus rule-of-2 (ICD 388.3x / H93.1x, ≥2 dates) with NO GENO_ID in Demographics
(ZNF175 carrier status untestable — would need WES).
File: `nb03_list_B_tinnitus_no_genotype.csv`.
Ask: whether any of these have WES (a GENO_ID), i.e. the reverse crosswalk.

## Why it matters
A carrier who is also a tinnitus case would sit on List A (carrier) and, if unlinked, effectively on List B
(tinnitus, no linked genotype). The crosswalk lets us cross the two and settle 4 vs 8.
Note (expectation): the linked carrier tinnitus rate is 4/27 ≈ 15%, so recovering all 7 is expected to add ~1 case
(→5), not 4 — reaching 8 would require strong enrichment among the unlinked.
