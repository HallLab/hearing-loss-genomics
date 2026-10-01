# PMBB v4 rare-variant pipeline — replication status

**Andre Rico · 2026-10-01 · for the meeting with Nikki Palmiero**

An independent re-derivation of the pipeline, phase by phase from the institutional release, before
further work is built on its outputs. Phases 1–3 are done. Everything below is reproducible from
scripts in `analysis/elena_replication/`.

**This is not an audit of anyone's work.** It re-derives numbers before they are published, and it
runs in both directions — Phase 1 found the replication itself wrong and the pipeline right.

---

## What reproduces

| | |
|---|---|
| **Case/control definition** | exact, **70,925 of 70,925** person-for-person |
| Rule-of-2, ambiguous-middle exclusion, the v4 tinnitus relocation | all correct |
| Ancestry stratification | internally consistent |
| Combined-cohort PC count (5) | matches its own scree |
| VEP run, SpliceAI extraction, `pDM` and `ALL` masks | sound |

The v4 tinnitus fix deserves specific credit: 89% of tinnitus evidence sits in the OMOP `observation`
table rather than the phecode files, this is documented nowhere in the release, and missing it is what
made the replication's own first rebuild wrong by 558 people.

---

## What does not

| # | Finding | Measured | Direction |
|---|---|---|---|
| 1 | The analysis cohort is cut against the **imputed** genotype `.fam`, while every analytical step runs on **exome**. An exome LD-pruned set existed and was not used. | **431** analysable people absent, **40 of them cases** | toward the null |
| 2 | The phenotype file was regenerated after the covariates were built from it, and the hearing-loss branch was never rebuilt — though the tinnitus branch was, the next day. | **556** people the rules exclude ran as controls | toward the null |
| 3 | The exclusion in (1) is not a random 0.7%. | **15.6% of East Asian** participants against 0.73% cohort-wide | representativeness |
| 4 | The `pLOF` mask admits three VEP `IMPACT=LOW` splice terms unconditionally, and the documented SpliceAI ≥ 0.2 gate is never applied. | **63.8%** of the mask has no loss-of-function consequence — 694,802 of 1,089,876 variants, every chromosome between 60.9% and 65.7% | toward the null |
| 5 | PC counts exceed what the pipeline's own scree plots support. | EUR ran **9**, scree supports **4**; AFR ran **10**, supports **3** | over-correction |

Every defect dilutes rather than fabricates. **The risk is a lost finding, not a false one** — anything
that reached significance would reach it more easily once corrected. What cannot be said is that genes
absent from the results are genuinely without effect.

---

## Questions

**Nikki**
1. Where did the `lof_terms` list come from? Three of nine entries are `IMPACT=LOW`, named explicitly
   rather than caught by accident — a published definition or an inherited template is likelier than
   an oversight, and worth recording.
2. Is REVEL 0.5 vs 0.6 still open? Recorded as pending since 2026-07-01.

**Nikki and Elena**
3. Was the `.fam` filter deliberate? The comment reads *"Keep only samples with genotype data"*, which
   is the right intent against the wrong dataset.
4. The tinnitus covariates were rebuilt from the corrected phenotype on 1 Aug; the hearing-loss ones
   were not, and SAIGE ran on them on 3 Aug. Was that intentional, or did the branch get missed?

**Molly, Doug, Nikki**
5. Restore the 431? Correcting also returns **165 East Asian participants** — that stratum grows 17.6%
   while every other group moves by half a percent.
6. Rebuild the masks and re-run the gene burden? Corrected masks are already built.
7. Confirm **5 / 4 / 3** PCs for the corrected arm, or override.

---

## Already built, if the answers say go

| | |
|---|---|
| Two cohorts | `phase_1/results/cohorts/` — reproduction 57,632 · corrected 57,507 |
| Two mask sets | `phase_2/results/masks/` — corrected `pLOF` is 41% the size of the one that ran |
| PC decision | 5 / 4 / 3, from the eigenvalues |

The replication runs **two arms throughout**: a reproduction arm on what actually ran, and a corrected
arm. The difference between them is the result — it turns "there is a defect" into "the defect moves
this gene by this much", or into "it changes nothing", which is a different recommendation.

⚠️ One operational note: the covariate file the SAIGE run consumed is the **only surviving record of
what was analysed**, because the phenotype file it came from was overwritten. If it is deleted, the
reproduction arm cannot be rebuilt.

---

## Where to read more

| Page | Covers | Read if |
|---|---|---|
| [`06_plof_mask_defect.md`](06_plof_mask_defect.md) | Finding 4, with line references and a runnable check | you want to verify the mask defect |
| [`01_sample_frame_decision.md`](01_sample_frame_decision.md) | Findings 1–3, with the decision they need | you want the cohort story |
| [`03_step_1_explained.md`](03_step_1_explained.md) | Finding 1 in plain language | you want to hand it to someone without context |
| [`02_phase_1_reference.md`](02_phase_1_reference.md) | the full Phase 1 record, every number and its provenance | you want to check a figure |
| [`../phase_3/scripts/01_pc_selection.ipynb`](../phase_3/scripts/01_pc_selection.ipynb) | Finding 5, scree plots and PC scatters | you want to see the PC evidence |
| `phase_1/results/FINDINGS.md`, `phase_2/results/FINDINGS.md` | what is established vs inferred, per phase | you want the caveats |

---

## Not covered

Tinnitus as a phenotype, and the combined hearing-loss-and/or-tinnitus phenotype decided on
2026-07-01 — which was never built in this pipeline. Phases 4 and 5 (the association test and the
result tables) are not started: what these defects do to the published numbers is a Phase 4 question
and is deliberately not estimated here.
