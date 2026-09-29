# PMBB v4 hearing-loss cohort — 517 participants excluded for missing imputed PCs

**Author:** Andre Rico · **Date:** 2026-09-29 · **Status:** for review
**Artifacts affected:** files owned by both Nikki Palmiero and Elena (see *Where the rule lives*)
**Evidence:** `analysis/elena_replication/phase_1/` — script, outputs and per-person list

---

## Decision requested

1. **Was the exclusion deliberate, and which step introduced it?** 517 exome-sequenced
   participants are absent from the SAIGE sample list because they have no imputed-array principal
   components. The analysis is adjusted with exome PCs; the imputed ones are never used as
   covariates. The requirement appears in artifacts owned by both of you, and the files alone do not
   say which step introduced it or whether one inherited it from the other. → **Nikki and Elena**
2. **Should the 427 affected cases and controls be restored** for the hearing-loss analysis?
   → **Molly, Doug, Nikki**
3. **If restored, does anything already run need re-running?** → **Molly, Nikki, Elena**

---

## First, what checks out

The case/control definition was re-derived independently from
`/static/PMBB/PMBB-Release-2026-4.0/`, written from the release's own files without reference to
the existing code, and compared person by person:

| | existing pipeline | independent rebuild |
|---|---:|---:|
| cases | 6,752 | **6,752** |
| excluded, one date | 4,007 | **4,007** |
| excluded, related ear phenotype | 9,411 | **9,411** |
| controls | 50,755 | **50,755** |

**70,925 of 70,925 — 100%, person for person.** The phenotype logic is correct and now
independently confirmed.

Two details worth crediting specifically. The rule-of-2 and the related-ear exclusion behave
exactly as documented. And the pipeline correctly handles the PMBB v4 change that moved the
standard tinnitus ICD codes (388.3x, H93.1x) into the OMOP `observation` table — our first rebuild
missed it, disagreed on 558 people, and only reached 100% once that source was added. That
relocation is documented nowhere in the release and is easy to miss.

## What we found

The cohort shrinks between the phenotype file and the sample list given to SAIGE:

```
70,925   exome cohort in the release
57,507   case or control for hearing impairment
57,080   delivered to SAIGE
   427   dropped  --  40 cases, 387 controls
```

The 427 belong to a larger group of **517** lost between `HL_TIN_PMBBv4_keep.txt` (70,925) and
`HL_TIN_PMBBv4_SAIGE_samples.txt` (70,408). Checked against the release, all 517:

| | |
|---|---|
| have exome data | yes — all 517 are in the exome ancestry list |
| have complete `exome_PC1-6` | yes — zero missing |
| have any `imputed_PC1-20` | **no — all 20 missing, for all 517** |

They are exactly the set of people in the release covariates who have exome PCs but no imputed PCs:
517 of 70,925, or 0.7% of the cohort.

### The model is adjusted correctly — the filter is what went wrong

This is worth stating precisely, because the milder reading is the correct one.

The covariates the model consumes are **exome** PCs. Verified by value, not by column name: the
`PC1` column of `covariates_combined_5PCs_withBatch.txt` is numerically identical to the release's
`exome_PC1` (correlation 0.99965 over 70,404 people; `imputed_PC1` correlates 0.992 but the values
differ). The imputed PCs never enter the regression.

They act as a **gate**, not as a covariate. The sample list is exactly:

```
analysable  ∩  {has all 20 imputed PCs}
  57,507         (excludes 427)          =  57,080
```

— an exact set identity, not an approximation. The effect is what a completeness check spanning
every PC column in the file would produce, rather than only the six the model uses.

So this is not a case of an exome analysis being adjusted with ancestry components derived from
array data. That would be a methodological error; this is not it. The adjustment is right. What is
wrong is that a completeness requirement over unused columns silently reduced the cohort.

### Where the rule lives

| Artifact | Owner | State |
|---|---|---|
| `PMBBv4_phecodex/` phenotype, `all_statuses`, summary | Nikki | phenotype reproduces 100% |
| `PMBBv4_phecodex/hearing_impairment_PMBBv4_SAIGE.txt` | Nikki | 57,080 — the rule is exactly reproducible here |
| `HL_TIN_PMBBv4_keep.txt` | Elena | 70,925 — full cohort, rule not yet applied |
| `HL_TIN_PMBBv4_SAIGE_samples.txt` | Elena | 70,408 — rule applied, zero people without imputed PCs |

The requirement is present on both sides. Which step introduced it is question 1.

## Why it is worth a decision rather than a fix

**It is a reasonable default, not an obvious error.** Requiring complete covariates is standard, and
a covariate template carried over from an imputed-genotype analysis would produce exactly this. At
0.7% of the cohort it is well below the level that shows up in any summary statistic. We are not
assuming it was a mistake — that is the first question above.

**But the published summary and the analysed cohort disagree, and nothing says so.**

| | `PMBBv4_hearing_tinnitus_summary.csv` | file given to SAIGE |
|---|---:|---:|
| hearing-loss cases | 6,752 | **6,712** |
| controls | 50,755 | **50,368** |

Anyone citing the summary would report 6,752 cases. The analysis ran on 6,712. The summary is not
wrong about the phenotype; it describes a different sample frame than the one used, with no note
that the two differ. That is the part that will cause trouble later regardless of how the first
question is answered.

**The lost cases are not a random 0.7%.** Whether participants lacking array imputation differ
systematically from those who have it — by enrolment era, by site, by ancestry — we have not
checked. If they do, removing them is not neutral.

## Options

| | Action | Cost |
|---|---|---|
| **A** | Restore the 427 and re-run. Requires only that the covariate build stop requiring imputed PCs. | one SAIGE re-run |
| **B** | Keep the exclusion, and document it in the summary so the two numbers agree. | documentation only |
| **C** | Restore, and first check whether the 517 differ systematically from the rest of the cohort. | a short analysis, then A |

We are not recommending one. The choice depends on the answer to question 1, which we cannot
supply.

---

## How to check this yourself

```bash
cd /project/hall/analysis/hearing-loss-genomics/analysis/elena_replication/phase_1
python3 scripts/01_sample_frame.py           # the chain above, regenerated
python3 scripts/02_rebuild_cases_controls.py # the 100% phenotype agreement
```

Per-person list of the 427: `phase_1/results/01_dropped_participants.csv`.
Full write-up including what was not established: `phase_1/results/FINDINGS.md`.

---

## Context

This came out of an independent replication of the PMBB v4 rare-variant pipeline, which is being
re-derived phase by phase from the institutional release before further work is built on its
outputs. Phase 1 is the cohort definition. The replication is not an audit of any individual's
work — it is re-deriving numbers before publishing them, and it runs in both directions: this same
phase found the replication itself wrong about the tinnitus source, and the pipeline right.
