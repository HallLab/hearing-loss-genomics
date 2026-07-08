# ZNF175 → tinnitus in PMBB Freeze One: why we count 4 carrier-cases (pre-conclusion)

**Author:** Andre Rico · **Date:** 2026-07-07 · **Status:** pre-conclusion (Chapter 2 v2)
**Scope:** reconcile our **4** ZNF175 tinnitus carrier-cases (v1/11K) with the project-internal **8**, and explain why the Park (2021) signal does not grow from 11K → 44K.

---

## Executive summary

Replicating Park et al. 2021 (ZNF175 pLOF → tinnitus) in Freeze One WES, we count **4 carrier-cases**. We have shown this **4 is robust** to every methodological choice we can vary — variant set, ICD code list, rule-of-2, event-vs-date counting, canonical phecode definition, **and annotation tool** (VEP vs ANNOVAR/REVEL/AlphaMissense). The same **4 individuals** drive the signal in both v1 (11K) and v2 (44K), proven by variant fingerprint. The gap to 8 is **not** a pipeline artifact; the most parsimonious reconstruction is **8 = 6 + 2** — a broader ear/rule-of-1 phenotype (6 on the linked carriers, which we reproduce) plus ~2 of the **7 unlinked** carriers, whose phenotype is blocked by a **GENO_ID ↔ PT_ID linkage gap** that cannot be closed with on-disk data. Separately, the failure to grow from 11K → 44K is **winner's curse**: the true effect is real but modest (carriers ~3× enriched), not the OR≈14.6 discovery estimate.

---

## The question

- Park 2021 (Nat Med): ZNF175 pLOF burden → tinnitus, **p = 3.24×10⁻¹⁰**, OR ≈ 14.6, **~35 carriers**, in ~11K WES (Freeze One).
- **Our replication (same Freeze One cohort): 34 carriers** (10 qualifying pLOF variants, MAF ≤ 0.1%) — matching Park's ~35 — of which **4 are tinnitus carrier-cases**; adjusted OR ≈ 14.6 (p ≈ 4.8×10⁻⁶), reproducing Park's effect a few orders weaker in p.
- The signal decays by ~44K (Hui-era) and by the exome-wide bar.
- Project-internal expectation: **8** carrier-cases in v1. We reproduce **4**. Why?

---

## What we established (chain of eliminations)

### 1. The variant set is faithful to Park
Our qualifying pLOF set (rare, MAF ≤ 0.1%) ≈ Park's preserved list: **10–11 pLOF variants / ~35 carriers** (Park) vs **10 / 34 carriers** (us). → the gap is **not** the variant set. Of the 34 carriers, **27 are phenotype-linked**, **7 are unlinked** (no PT_ID).

### 2. Canonical phecode phenotyping gives the same 4 (and fixed a real bug)
We rebuilt tinnitus with the canonical PheWAS pipeline (`createPhenotypes`, Phecode Map 1.2, rule-of-2, control exclusions). **Tinnitus = phecode 389.4.**

While validating, we found a format bug: feeding `createPhenotypes` a **pre-aggregated integer count** (instead of raw dates) **undercounts** cases. Corrected (raw dates as strings):

| phecode | pre-aggregated (bug) | corrected |
|---|---|---|
| 389.4 Tinnitus | 114 | **131** |
| 389 Hearing loss | 357 | **481** |
| 389.1 Sensorineural | 185 | **227** |

Cohort tinnitus rose 114→131, but the **carrier-cases stayed 4** — confirmed three ways (raw ICD, buggy phecode, corrected phecode). The 4 raw-ICD cases and the 4 phecode cases are the **same 4 people** (perfect overlap, not additive).

### 3. No counting rule reaches 8 (sensitivity grid, 27 linked carriers)

| definition | tinnitus 389.4 | hearing loss 389 |
|---|---|---|
| distinct-date ≥2 (rule-of-2, canonical) | **4** | 4 |
| distinct-date ≥1 (rule-of-1) | 5 | 6 |
| event-count ≥2 (not deduped by date) | **4** | 5 |
| event-count ≥1 (presence) | 5 | 6 |

Even counting raw events (not distinct dates) leaves tinnitus at **4** — no carrier has 2 tinnitus events on the same day. Broadening further to **all ear/hearing phecodes (388 + 389.x) at rule-of-1 saturates at 6** — it adds no one beyond tinnitus-OR-HL, because the ear-disease carriers are the same people (§6). Max reachable on linked carriers = **6**. **Nothing reaches 8 on the linked carriers.**

### 4. The 4 is robust to the annotation tool
There is a hard **ceiling**: only a tinnitus case carrying a *rare* ZNF175 variant could ever qualify — that is **11 people**. Our 4 are pLOF (frameshift); the other **7** carry variants no annotator upgrades:

| the 7 non-pLOF rare variants in tinnitus cases | REVEL | AlphaMissense |
|---|---|---|
| 4 × synonymous, 2 × intron, 1 × 5′UTR | — | — |
| 3 × missense (51586879, 51573428, 51587116) | 0.02–0.03 | Benign |

ANNOVAR gene-based LOF would not reclassify synonymous/missense as LOF; REVEL≥0.5 and AlphaMissense both leave the missense **benign**. → carrier-cases = **4**, robust to VEP vs ANNOVAR/REVEL/AlphaMissense.

### 5. Genotype ↔ phenotype bridge (one table, whole story)
Among each ear/HL phecode's cases, how many carry ZNF175 variants at each tier:

| phecode | cases | any ZNF175 | rare (MAF≤0.1%) | **qualifying pLOF** |
|---|---|---|---|---|
| 389.4 Tinnitus | 131 | 112 | 11 | **4** |
| 389 Hearing loss | 477 | 384 | 20 | 4 |
| 389.1 Sensorineural | 225 | 179 | 12 | 2 |

Reading tinnitus left→right is the burden filter live: 131 cases → 112 carry *some* ZNF175 variant (mostly common, uninformative) → 11 carry a *rare* one → **4** carry the *rare pLOF* signal.

### 6. The signal is tinnitus-specific
The 4 qualifying carrier-cases for hearing-loss (389) and tinnitus (389.4) are the **same 4 people**. Since 389.4 rolls up into 389 and 389 (=4) does not exceed 389.4 (=4), **no ZNF175 carrier has hearing loss without tinnitus**. The signal is phenotypically coherent with Park — specific to tinnitus, not diffuse ear disease.

### 7. Same 4 individuals in v1 and v2 — winner's curse made concrete
The v1 (11K) and v2 (44K) carrier-cases are the **same 4 people** — the ID scheme changed (UPENN → PMBB_ID) but the **variant fingerprint is identical**: `{51588428:C ×2, 51587727:C ×1, 51588214:C ×1}` in both freezes.

And "no new cases in v2" is **expected, not strange**:

| scenario | expected v2 carrier-cases |
|---|---|
| null (OR=1) | ~1.4 |
| our regressed OR ≈ 3.5 | ~4.5 |
| Park discovery OR ≈ 14.6 | ~16 |
| **observed** | **4** |

Carriers are still ~3× enriched for tinnitus (5.8% vs 2.0%) — a **real, modest** effect. Observed 4 matches OR≈3.5, sits above null, and is far below the ~16 that OR≈14.6 would predict. Quadrupling the cohort added ~0 cases because the true effect is modest and anchored on these 4 index individuals — textbook winner's curse.

---

## The residual: a crosswalk gap with two sides

The only way the count could exceed 4 is the **7 unlinked carriers** (real qualifying carriers — 3× frameshift 51581437, 2× stop_gained 51587814, 1× stop_gained 51588382, 1× frameshift 51588428 — whose phenotype is unknown). We confirmed this is a **data-access wall**, not analysis:

- v1 Freeze One: the 7 GENO_IDs are absent from the only crosswalk (`Demographics`).
- Published imputed V1.0 (~20K) and v2/v3/v4 all use `PMBB_ID`; **no UPENN↔PMBB_ID bridge exists on disk**.

The gap has a mirror side. Of **287** v1 tinnitus cases: **131** are WES-genotyped (our cohort, 4 carriers), **79** are array-genotyped only (blind to rare pLOF), and **77** have **no genotype at all**. So:

- **7 WES carriers** with no phenotype (carrier status known, tinnitus unknown) — *List A* (NB03).
- **77 tinnitus cases** with no genotype (tinnitus known, carrier status unknown) — *List B* (NB03).

A carrier who is also a tinnitus case would sit on both lists; only the master crosswalk lets us cross them.

---

## Leading reconstruction: **8 = 6 + 2**

The most parsimonious account of the project-internal **8** stacks two choices we deliberately do *not* make in our strict pipeline:

- **6** = a broader phenotype (all ear/hearing phecodes, **rule-of-1**) on the **27 linked** carriers — we reproduce this exactly (broadening beyond tinnitus saturates at 6; §3, §6).
- **+2** = two of the **7 unlinked** carriers being ear/tinnitus cases.

The arithmetic is internally consistent — no statistical luck required:

| base rate (linked carriers) | expected cases among the 7 unlinked |
|---|---|
| 4/27 = 15% (strict tinnitus, rule-of-2) | ~1.0 |
| **6/27 = 22% (broad ear, rule-of-1)** | **~1.5 → 2** |

Under the same broad definition that yields 6, the expected number of cases among the 7 unlinked is ~1.5 → **2**. So **6 + 2 = 8** falls out naturally. The alternative "4 + 4" would require 4/7 ≈ 57% of the unlinked to be cases — far above the observed rate, implausible.

**This reframes 4-vs-8 as not a discrepancy but two defensible choices:** strict tinnitus + linked-only = **4**; broad ear + rule-of-1 + full linkage = **8**.

**Two testable predictions:**
1. The 8 used **rule-of-1 + a broad ear/hearing phenotype** (not strict rule-of-2 tinnitus) — confirm by asking Park / the PI for the qualifying + phenotype rationale.
2. Recovering the 7 unlinked via the crosswalk yields **~2 ear/tinnitus cases** — confirm with the PMBB curators.

If both hold, 4-vs-8 is fully resolved.

---

## Paths forward (saídas)

1. **Get the master GENO_ID ↔ PT_ID crosswalk from the PMBB curators (Nikki).**
   Deliverables ready: `NB03` → *List A* (7 WES carriers, `nb03_list_A_...csv`) + *List B* (77 tinnitus cases, `nb03_list_B_...csv`) + a framed request (`nb03_curator_request.md`). Ask: (a) the PT_ID for each List-A carrier or confirmation they are unconsented/QC-drop; (b) whether any List-B tinnitus case has WES. *Expectation:* the linked carrier tinnitus rate is 4/27 ≈ 15%, so recovering all 7 likely adds ~1 case (→5); reaching 8 needs strong enrichment.

2. **Ask Park et al. (via the PI) for the exact rationale behind the 8.**
   Clarify the definition that produced the project-internal 8: (a) which qualifying-variant filter (MAF cutoff, pLOF term set, REVEL/missense inclusion, multiallelic handling); (b) which phenotype (tinnitus phecode vs broader HL, rule-of-2 vs other); (c) which cohort/freeze and sample set; (d) whether "8" counts carriers or carrier-cases, and whether it already includes the unlinked. This pins down whether 8 is a different definition we can reproduce, or requires the crosswalk.

3. **(Optional) Firth-penalized logistic replication** on the 4-anchored cell, to report the v1 effect with a discovery-faithful model and a confidence interval that makes the winner's-curse regression explicit.

---

## Key artifacts
- Notebooks: `scripts/01_extract_phecodes_11k.ipynb` (phenotype + §2b–2d bridge), `02_znf175_carrier_cases_phecode.ipynb` (carrier reconciliation, sensitivity, annotation robustness, v1↔v2), `03_lists_for_pmbb_curators.ipynb` (the two lists).
- Strategy: `strategy_8_vs_4_carriers.md`. Chapter-2 findings: `../chapter_2/findings_znf175_11k_vs_44k.md`.
- Reference papers: Park et al. 2021 (Nat Med), Hui et al. 2023 (PLOS Genet).
