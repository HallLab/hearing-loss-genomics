# Chapter 2 v2 — Strategy: why does the PI count 8 carriers/cases in v1 and we count 4?

**Author:** Andre Rico · **Date:** 2026-07-07 · **Status:** ✅ EXECUTED — see pre-conclusion `preconclusion_znf175_4_vs_8.md` (NB 01–03)

## Objective
Reconcile the **"8" (PI) vs "4" (us)** for ZNF175 in the v1 (Freeze One, ~11K) cohort. Find the exact methodological difference that produces 8, and either reproduce it or explain why ours is 4.

---

## ⭐ RESOLUTION (2026-07-07) — H1 refined: **8 = 6 + 2**

Our **4 is robust** to variant set, ICD list, rule-of-2, event-vs-date counting, canonical phecode, and annotation tool (VEP vs ANNOVAR/REVEL/AlphaMissense). The **leading reconstruction of the 8 is `6 + 2`**:
- **6** = broad ear/hearing phecodes + **rule-of-1** on the 27 linked carriers (broadening saturates at 6; tinnitus-specific).
- **+2** = ~2 of the 7 unlinked carriers being cases (expected 22%×7 ≈ 1.5→2; "4+4" would need 57%, implausible).

The 7 unlinked are bona-fide rare-pLOF carriers whose phenotype is **blocked by the missing UPENN↔PMBB_ID crosswalk** (not recoverable on disk — memory `project_upenn_pmbbid_no_bridge`). **Two testable predictions** remain: (1) the 8 used rule-of-1 + broad ear (ask Park/PI); (2) recovering the 7 yields ~2 cases (ask Nikki — NB 03 Lists A + B ready). Full detail + tables in the pre-conclusion.

---

## What we already established (Chapter 2 — do NOT re-derive)
- **Our variant set is faithful to Park.** Our qualifying pLOF set (10 variants, MAF≤0.1%) ≈ Joe Park's preserved list (`data/PMBB_Exome/ZNF175/Joe_analyses/znf175_variants.txt`, b37): **10–11 pLOF variants, 35 carriers** vs our **10 / 34 carriers**. → the gap is **NOT the variant set**.
- **Our v1 carrier breakdown:** 34 carriers → **27 phenotype-linked** (via `Demographics` GENO_ID→PT_ID) + **7 UNLINKED** (no PT_ID → tinnitus status unknown). Among the 27: **4 tinnitus cases** (rule-of-2, ICD 388.3x/H93.1x).
- **Ruled out:** broader tinnitus codes (H93.A, other ear codes) add **0** cases; relaxing rule-of-2 to ≥1 date adds only **1** (→5).

→ From our side, the missing cases are **not** in the variant set, the ICD code list, or the date rule. The leading suspects are **(a) the 7 unlinked carriers** and **(b) a different definition of "8" entirely** (Daniel's pipeline / phenotype).

---

## ⚠️ Step 0 — Disambiguate what "8" means (do this FIRST)
"8 carriers" is ambiguous and everything downstream depends on it. Pin down:
1. **8 carriers** (people carrying qualifying variants) **or 8 carrier-CASES** (carriers who are also phenotype-positive)?
   - If 8 **carriers**: our count is 34 — so the PI's set is a *narrower* variant/MAF filter. Which one?
   - If 8 **carrier-cases**: matches our "4 vs 8" framing → phenotype/linkage question.
2. **Which analysis** is the "8" from — **Park's tinnitus** (Nat Med) or **Daniel's unpublished ZNF175** (HL-focused, Hui pipeline)?
3. **Which phenotype** — tinnitus (phecode 389.x) or hearing loss (audiogram / phecode 389)?
4. **Which cohort/freeze** — 11K Freeze One, or a specific filtered subset?

Resolve via: (a) ask the PI directly, and (b) trace it in Daniel's preserved files/runbook (below) — the "8 signal-driving cases" is a project-internal number from **Daniel's runbook** (per memory `project-znf175-priority`).

---

## Hypotheses & how to test each

| # | Hypothesis | Test | Data |
|---|---|---|---|
| **H1 ✅ LEADING (refined: 8 = 6 + 2)** | 8 = broad-ear/rule-of-1 (**6** on linked) + **~2 of the 7 unlinked** carriers being cases | Reproduced the 6 exactly; recover phenotype for the 7 unlinked (needs master crosswalk — blocked on disk, ask Nikki). Predict ~2 cases (not 4). NB 03 Lists A+B ready | `Demographics`, master UPENN↔PMBB_ID crosswalk (ask Nikki/PMBB) |
| **H2** | "8" comes from **Daniel's exact ZNF175 pipeline** (his files/runbook), which differs from ours | Reproduce Daniel's ZNF175 carrier/case derivation from his preserved extracts + runbook; diff his variant set, sample set, phenotype vs ours | `data/PMBB_Exome/ZNF175/` (allIndvs, `_Joes`, `_maf.001_noRels_keepHLcases`, `Joe_analyses/`), `analysis/daniel/runbook_hui2023.txt` |
| **H3** | "8" = 8 **carriers** of a narrower set (fewer variants / stricter MAF / no multiallelic) | Sweep carrier counts across variant-set definitions; find which yields 8 | our NB06 + Daniel's `.bim` sets |
| **H4** | "8" uses a **different phenotype** (phecode via PheWAS map, or HL not tinnitus) | Re-derive cases with (a) tinnitus phecode via R `PheWAS`, (b) HL phenotype; count carrier-cases | R `PheWAS` map, ICD files |
| **H5** | "8" uses the **full 10,900** (we used 9,161 after covariate/phenotype merge) | Count carrier-cases on the full Park sample set (no covariate drop) | Freeze One + phenotype |

---

## Proposed sequence (NBs in `chapter_2_v2/scripts/`)
1. **NB 00 — trace the "8".** Read Daniel's runbook + preserved ZNF175 files; inventory his variant set, sample set, phenotype, and any explicit "8". Ask PI for the definition in parallel. *(H0/H2)*
2. **NB 01 — the 7 unlinked carriers.** Who are they? Try to recover their phenotype (fuller crosswalk / v2 linkage). *(H1)*
3. **NB 02 — sensitivity sweep.** Grid over {variant set × MAF cutoff × phenotype (ICD/phecode/HL) × sample set (9,161 vs full)} → table of carrier and carrier-case counts; find the cell(s) = 8. *(H3/H4/H5)*
4. **NB 03 — reconcile & report.** State the exact combination that produces 8, or the honest residual (e.g., "8 requires the 7 unlinked carriers' phenotype, which needs the full PMBB ID map").

---

## Success criteria
- **Best:** identify the exact definition/pipeline that yields 8, reproduce it, and document the delta vs our 4.
- **Acceptable:** show that 8 is only reachable with a specific choice we can't fully replicate (e.g., the full ID crosswalk for the 7 unlinked, or Daniel's exact phenotype) — and name precisely what's needed to close it.

## Key references
- Chapter 2 findings: `analysis/chapter_2/findings_znf175_11k_vs_44k.md`; carriers: `analysis/chapter_2/results/06/`, `07/`.
- Daniel's ZNF175: `data/PMBB_Exome/ZNF175/`; runbook: `analysis/daniel/runbook_hui2023.txt`.
- Memories: `project-znf175-11k-44k-resolution`, `project-znf175-priority`.
