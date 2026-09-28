# Cycle 2 — Rare Variation in Mendelian Hearing-Loss Genes and Adult-Onset Hearing Loss

**Status:** DRAFT charter — for review by Molly Hall, Douglas Epstein, Nikki Palmiero before any analysis runs
**Author:** Andre Rico · **Opened:** 2026-08-26
**Cohort:** PMBB Release 2026-4.0 (v4)

> Portuguese version (Andre's working copy): [`README.pt.md`](README.pt.md). This English file is canonical.

---

## 0. Why a new cycle

Everything outside this folder is **Cycle 1** (`analysis/`, `data/`, `docs/`, `results/`, `scripts/`). Cycle 1 was
organized around a *gene* — ZNF175 — and it did its job: it answered the question it was given, and the answer
was negative. Cycle 2 is organized around a *question* instead, so that the project stays meaningful whichever
way the result comes out.

Cycle 1 is not moved or archived — it is still live (Elena's v4 pipelines, two open data blockers). Cycle 2 sits
alongside it and reuses its infrastructure.

### What Cycle 1 established (settled — do not relitigate)

| Finding | Where |
|---|---|
| Hui et al. 2023 replicated on PMBB v2; all 6 paper genes recovered at FDR<0.05 | [`results/chapter1_paper_replication/`](../results/chapter1_paper_replication/) |
| The ZNF175→tinnitus signal at 11K reproduces (OR≈14.6) and decays by 44K | [`analysis/chapter_2/findings_znf175_11k_vs_44k.md`](../analysis/chapter_2/findings_znf175_11k_vs_44k.md) |
| The decay is **winner's curse**, not a pipeline artifact — the same 4 individuals drive it in both freezes | [`analysis/chapter_2_v2/preconclusion_znf175_4_vs_8.md`](../analysis/chapter_2_v2/preconclusion_znf175_4_vs_8.md) |
| The "8 vs 4" carrier discrepancy resolves as `8 = 6 + 2` (broader phenotype + unlinked carriers) | same |
| **ZNF175 is null in PMBB v4** across all masks × MAF × ancestry, with negative beta | [`analysis/elena/HL_only_rarevariant/`](../analysis/elena/HL_only_rarevariant/) |

**Consequence:** ZNF175 is no longer a subject. In Cycle 2 it is one row in a results table, nothing more.

---

## 1. Research question

### Primary — Q1

> **Does rare damaging variation in the known Mendelian hearing-loss genes contribute to adult-onset hearing
> loss in a hospital-based biobank — and if so, in which genes and at what effect size?**

### Secondary — Q2 (conditional on Q1 defining a carrier set)

> **Among carriers of those variants, what distinguishes carriers who develop hearing loss from carriers who
> do not?** (penetrance; age, sex, ancestry; a second hit elsewhere in the gene set; `CX3CR1` as the one
> *a priori* named modifier from the mouse work)

Q1 and Q2 are nested, not parallel: **Q1 builds the carrier set, Q2 asks what happens inside it.** Q1 is worth
doing regardless of its outcome; Q2 is only interpretable once Q1 exists.

### Why this question and not another

- **It has a real prior.** ~100 well-characterized non-syndromic congenital HL genes. A pre-specified set of
  ~100–200 genes replaces the exome-wide bar (~2.5×10⁻⁶) with FDR over ~10² tests — a large power gain for free.
- **The answer is genuinely open.** These genes are adjudicated for *congenital, mostly recessive, mostly severe*
  HL. Whether the heterozygous carrier state elevates risk of *adult-onset* HL is untested, not foregone. We do
  **not** expect a strong signal — if we did, the study would not be worth running. The expected outcome is null.
- **The negative is publishable.** "Carriers of rare damaging variants in Mendelian HL genes are not enriched
  for adult-onset HL in the EHR" directly contradicts the translational premise stated at kickoff (*linking
  adult HL to well-characterized congenital HL genes*). That is information, not failure.
- **It survives its own null.** Unlike a single-gene project, no outcome leaves us with nothing to report.

---

## 2. What we have (counted, not assumed)

### Phecode cohort — PMBB v4

| Phenotype | N with exome | cases | prevalence |
|---|---|---|---|
| Hearing impairment | 57,632 | **6,712** | 11.7% |
| Tinnitus | 53,096 | 2,732 | 5.1% |

*(counted from Elena's SAIGE covariate files, `analysis/elena/HL_only_rarevariant/new_SAIGE_covariates/`)*

Sample size is no longer the binding constraint. **Phenotype quality is.** The `hearing impairment` phecode
pools conductive, sudden, noise-induced and age-related loss into one label, and most adult HL is polygenic
and age-driven. This is why the audiogram matters — see below.

### Audiometric cohort — the underused asset

Brant's audiometry database (`audbase`) is already on disk from the Cycle 1 era, together with a
**PMBB linkage that was already built once, in February 2021**:

| File | Content |
|---|---|
| [`data/PMBB_Exome/brant audbase 1.7.21.TXT.gz`](../data/PMBB_Exome/) | raw audbase, 100,470 subjects / 216,542 records — **contains PHI** |
| [`data/PMBB_Exome/audbase_feb252021/RGC21_45k_aud_1.csv.gz`](../data/PMBB_Exome/audbase_feb252021/) | derived phenotypes: PTA air/bone per ear, `PTA`, `Bilateral_HL`, `Worse_ear`, `Degree_HL` (0–4), `BL_SNHL` |
| [`data/PMBB_Exome/audbase_feb252021/degree_HL_aud.txt.gz`](../data/PMBB_Exome/audbase_feb252021/) | 3,328 `PMBB_ID → Degree_HL_Aud`, de-identified |

**The key number:** that file holds **45,012 audiometry subjects by MRN, of which only 3,328 (7.4%) carry a
PMBB_ID.** The remaining ~41,700 are unlinked. The linkage was done in Feb 2021 by MRN + date of birth against
the PMBB of that era; PMBB has grown substantially since. Both sides of the join have grown.

So the item recorded as a blocker in the 2026-07-01 meeting — *"audiogram data not yet matched to PMBB IDs"* —
is more precisely: **a linkage exists, is five years stale, and covers 7.4% of the audiometry cohort.** We are
refreshing and extending a worked example, not building one from scratch.

> ⚠ **Data governance.** `audbase` raw and `RGC21_45k_aud_1.csv.gz` contain identifiers (name, address, phone,
> e-mail, DOB, MRN). Any re-linkage is an honest-broker / IRB action under the *Audiometric Phenotyping of PMBB
> Enrollees* project — not something to run on our own. Note the 2026-07-01 meeting recorded Elena as attached
> to the **wrong IRB** (Callback, not Audiometric Phenotyping); until Nikki's correction lands, these files are
> out of scope for her. Ratchet (Ritchie lab) is the assigned owner for the v4 refresh — coordinate, don't duplicate.

---

## 3. Design

Two tiers, each doing the job it is actually good at:

| Tier | Cohort | Phenotype | Role |
|---|---|---|---|
| **A — Discovery** | 57,632 exomes | phecode HL (binary), n=6,712 | powered; noisy label |
| **B — Validation** | ~3,328 exome + audiogram (growing with the refresh) | `Degree_HL` / PTA (quantitative), `BL_SNHL` | clean label; underpowered alone |

Tier B does **not** replace Tier A — 3.3K is far too small for rare-variant discovery. It validates and refines
it: a gene that survives Tier A and shows a dose-consistent shift in PTA in Tier B is a genuinely different
claim from one that only clears a phecode FDR. This is also the design that makes the linkage refresh pay off
scientifically rather than just administratively.

### Pre-specified analytic choices

Carried over from the 2026-07-01 decisions, plus the gaps that meeting left open:

| Choice | Value | Source |
|---|---|---|
| Phenotype | HL **and/or** tinnitus, combined | decided 2026-07-01 — *never built; Cycle 2 builds it* |
| Gene set (primary) | **ClinGen HL GCEP (40007), Definitive+Strong — n=100 genes**, snapshot 2026-08-26 | **provisional** — Andre 2026-08-26; tabled for ratification |
| Gene set (sensitivity) | same panel, +Moderate — n=120 genes | provisional, same decision |
| Masks | pLOF; pLOF+AlphaMissense; pLOF+REVEL; pLOF+AM+REVEL — **kept separate** | decided 2026-07-01; Cycle 1 collapsed them into one `pDM` |
| REVEL threshold | 0.5 primary, 0.6 sensitivity | decided 2026-07-01 (unresolved) |
| MAF | 0.01 / 0.001 / 0.0001 | decided 2026-07-01 |
| PCs | 5–6 (never 20) | decided 2026-07-01 |
| Multiple testing | **FDR within MAF group**, not across models | decided 2026-07-01 — *never applied to any result* |
| Method | SAIGE-GENE+ (burden + SKAT + Cauchy) | Elena's v4 pipeline, reused |
| Technical controls | non-HL gene–phenotype pairs with known large effects in PMBB (`BRCA1`/breast cancer, `TTN`/cardiomyopathy, `CFTR`/CF — the set Park 2021 used to validate), plus λ_GC and QQ on the exome-wide layer | new to Cycle 2 |

> **Not a control:** `GJB2`, `SLC26A4`, `MYO7A` and the rest of the primary set are the *hypothesis under test*,
> not a validity check on it. Treating them as positive controls would assume the answer. Pipeline validity is
> established on gene–phenotype pairs outside hearing.

### What counts as an answer

- **Q1 positive:** ≥1 gene from the primary set at FDR<0.05 within a MAF group, with technical controls behaving
  and λ_GC in range.
- **Q1 negative:** technical controls behave, calibration is clean, no gene survives → carrier-state variation in
  Mendelian HL genes does not detectably contribute to adult EHR-defined HL. **This is the expected outcome**, and
  it is reportable: it bounds the translational premise stated at kickoff.
- **Q1 uninformative:** technical controls fail or calibration is off → the finding is about the phenotype or the
  pipeline, and Tier B becomes the priority.

---

## 4. Scope guard — what Cycle 2 is *not*

- Not a ZNF175 project.
- Not a formal gene×gene interaction test (acknowledged underpowered at kickoff; Q2 stays descriptive).
- Not Menière's (too few cases; coordinate with Bogdan/Ian before entering).
- Not UK Biobank (access frozen since May 2026).
- Not a methods paper. The Cycle-1 winner's-curse work is discussion material, not the primary claim.

---

## 5. Open dependencies

| # | Item | Owner | Blocks |
|---|---|---|---|
| 1 | Which adjudicated HL gene list is authoritative | Doug / Andre | ~~blocks Q1~~ — **provisionally resolved 2026-08-26** (see §5.1); on the next meeting agenda |
| 2 | audbase ↔ PMBB v4 re-linkage | Ratchet (Ritchie lab); escalated by Doug | Tier B |
| 3 | Elena's IRB correction (Audiometric Phenotyping) | Nikki | Elena touching Tier B |
| 4 | GENO_ID ↔ PT_ID master crosswalk (7 Cycle-1 carriers) | Nikki / PMBB curators | Cycle 1 residual only |
| 5 | Final REVEL threshold (0.5 vs 0.6) | Molly / Nikki | mask construction |

### 5.1 Gene set — provisional decision (Andre, 2026-08-26)

**Adopted, pending ratification: the ClinGen Hearing Loss Gene Curation Expert Panel
(affiliate 40007), Definitive + Strong — 100 genes.** Sensitivity analysis on the same
panel + Moderate (120 genes). Snapshot `FILE CREATED: 2026-08-26`; sha256 pinned in
`cycle_2/data/clingen/manifest.json`. Regenerate with
[`cycle_2/analysis/clingen_fetcher/fetch_clingen_hl_genes.py`](analysis/clingen_fetcher/fetch_clingen_hl_genes.py).

Rationale — why this source and why not wait:

- **It is the only candidate that is itself adjudicated.** ClinVar is submitter-level and
  conflicted; DVD is a variant database, not a gene-validity authority; the Cycle-1
  173-gene set has no documented provenance. Only the GCEP publishes an evidence tier per
  gene-disease pair with a named panel and a date.
- **It matches the design.** §1 justifies the pre-specified set by trading the exome-wide
  bar for FDR over ~10² tests. 100 genes is exactly that regime — and it is a number the
  panel arrived at, not one we tuned.
- **Moderate is the sensitivity tier, not the primary.** Including it inflates the FDR
  denominator ~20% to admit genes the panel itself says may not hold up.
- **The decision is cheap to reverse.** The gene set enters as a file path. Swapping tiers
  is one flag on the fetch script and a re-run; no pipeline code changes. That is the
  actual argument for deciding now rather than waiting on the meeting — a provisional
  choice that is expensive to undo would not be worth making.

Still open, and *not* blocking:

- **Syndromic vs non-syndromic.** §1 motivates the set as "~100 non-syndromic congenital HL
  genes", but these 100 mix both. 16 genes carry pairs at different strengths (CIB2 =
  Refuted + Definitive, GJB6 = Definitive + Refuted, ADGRV1 = Definitive + Disputed); the
  script keeps each gene's *best* classification, so all 16 are in. Restricting to
  non-syndromic means filtering on `DISEASE LABEL`/MONDO, not on the gene — a different
  script, and a question for Doug.
- **The list drifts.** 164 pairs (2019) → 174 (Aug 2024) → 191 (Aug 2026). Report the
  snapshot date alongside any result. Cite Tshering KC et al., *Genet Med* 2025;27(5):101397
  (PMID 39987489), whose finding *is* the drift.

---

## 6. Next steps

1. Review this charter with Molly, Doug and Nikki — **before running anything** (the cadence set at kickoff).
2. ~~Resolve dependency #1 (the gene list).~~ Provisionally resolved — see §5.1. Raise at the next meeting; if the room disagrees, re-run the fetch script with a different tier and the downstream analysis is unaffected.
3. Build the combined HL-and/or-tinnitus phenotype for v4 — the one 2026-07-01 decision never executed.
4. Prepare the technical spec for the audbase↔PMBB refresh (join keys, expected yield, QC) so Doug and
   Ratchet receive a concrete ask rather than an open-ended one.
