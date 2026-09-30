# PMBB v4 hearing-loss cohort — 517 participants excluded for missing imputed PCs

**Author:** Andre Rico · **Date:** 2026-09-29 · **Status:** for review
**Artifacts affected:** files owned by both Nikki Palmiero and Elena (see *Where the rule lives*)
**Evidence:** `analysis/elena_replication/phase_1/` — script, outputs and per-person list

---

## Decision requested

1. **Was the exclusion deliberate?** 517 exome-sequenced participants are absent from the SAIGE
   sample list because the analysis cohort was built by merging against the **imputed genotype
   `.fam`**. They have exome data and complete exome PCs; they have no array/imputed data. The
   analysis is adjusted with exome PCs throughout, so nothing it computes needs the imputed set.
   → **Nikki and Elena**
2. **Should the 427 affected cases and controls be restored** for the hearing-loss analysis?
   → **Molly, Doug, Nikki**
3. **If restored, does anything already run need re-running?** → **Molly, Nikki, Elena**
4. **Separately: the phenotype file was overwritten 36 minutes after the SAIGE covariates were built
   from it.** The run therefore consumed 57,632 people, not the 57,080 on disk, and 556 people the
   phenotype rules exclude were analysed as controls. Was this known? → **Nikki and Elena**

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

The cohort shrinks between the phenotype definition and the phenotype file handed downstream:

```
70,925   exome cohort in the release
57,507   analysable: case or control  -- what the phenotype supports
57,632   what the SAIGE run actually consumed
           -431  analysable people never entered it, 40 of them cases
           +556  people the phenotype rules exclude entered as controls
```

`57,507 - 431 + 556 = 57,632`. The run is **not a subset** of the analysable cohort: short in one
direction, long in the other.

A third number, **57,080**, appears in the phenotype file on disk. It was written after the run's
inputs had already been built from an earlier version and **was never consumed by anything** — it is
evidence for Finding 3, not a stage in the flow. Comparisons against it would measure nothing.

The 427 belong to a larger group of **517** lost between `HL_TIN_PMBBv4_keep.txt` (70,925) and
`HL_TIN_PMBBv4_SAIGE_samples.txt` (70,408). Checked against the release, all 517:

| | |
|---|---|
| have exome data | yes — all 517 are in the exome ancestry list |
| have complete `exome_PC1-6` | yes — zero missing |
| have any `imputed_PC1-20` | **no — all 20 missing, for all 517** |
| appear in the imputed genotype `.fam` | **no — none of the 517** |

They are exactly the set of people in the release covariates who have exome PCs but no imputed PCs:
517 of 70,925, or 0.7% of the cohort.

**The mechanism is a merge.** The ExWAS analysis blog logs it directly —
`FAM samples: 70493 | Phenotype samples: 70925 | Matched samples: 70408` — and the arithmetic closes:
the imputed `.fam` holds 70,493 people, 85 of whom have no exome, giving the 70,408 sample list. The
517 have no imputed PCs because they have no imputed genotypes at all; the missing PCs are a marker
of that, not a separate filter.

### The model is adjusted correctly — the filter is what went wrong

This is worth stating precisely, because the milder reading is the correct one.

The covariates the model consumes are **exome** PCs. Verified by value, not by column name: the
`PC1` column of `covariates_combined_5PCs_withBatch.txt` is numerically identical to the release's
`exome_PC1` (correlation 0.99965 over 70,404 people; `imputed_PC1` correlates 0.992 but the values
differ). The imputed PCs never enter the regression.

They act as a **gate**, not as a covariate. The phenotype file on disk is exactly:

```
analysable  ∩  {has imputed genotypes}
  57,507         (excludes 427)        =  57,080
```

— an exact set identity. That file is not the one the run consumed (see the section below); it is
shown here because it makes the gate visible.

— an exact set identity, not an approximation. The two conditions — absent from the imputed `.fam`,
and missing every imputed PC — pick out the same people, because both follow from having no imputed
genotypes.

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

| | `PMBBv4_hearing_tinnitus_summary.csv` | phenotype file on disk |
|---|---:|---:|
| hearing-loss cases | 6,752 | **6,712** |
| controls | 50,755 | **50,368** |

Anyone citing the summary would report 6,752 cases. The analysis ran on 6,712. The summary is not
wrong about the phenotype; it describes a different sample frame than the one used, with no note
that the two differ. That is the part that will cause trouble later regardless of how the first
question is answered.

**The exclusion is not a random 0.7% — this has now been checked, and it is concentrated.**

| Group | excluded | total | rate |
|---|---:|---:|---:|
| **EAS** | **208** | **1,333** | **15.60%** |
| AMR | 16 | 1,039 | 1.54% |
| EUR | 233 | 51,867 | 0.45% |
| AFR | 49 | 14,927 | 0.33% |

Cohort-wide the rate is 0.73%. **One in six East Asian participants was excluded — twenty-one times
that.** The same skew shows on sex (1.01% of women against 0.43% of men) and strongly on batch (1.09%
of batch 1 against 0.11-0.18% of batches 2 and 3), and the excluded are younger and enrolled earlier.

The batch and enrolment pattern suggests a technical coverage gap rather than anything about the
participants: exome-sequenced early, array genotyping never completed for a subset. That is the
shape of it, not a demonstrated cause.

It does not change the direction of bias — losing cases still biases toward the null. It does mean
the loss cannot be called neutral. East Asian participants are 1.6% of the analysed cohort and a
sixth of them are gone.

## A second, separate problem in the same area

Found while checking the first. The covariate files the SAIGE run consumed hold **57,632** people —
not the 57,080 in the phenotype file on disk. An inner join cannot increase N, and the build log
names what it read: `Phenotype samples: 57,636`.

```
2026-07-31 20:21   the covariate build reads the phenotype file, writes the SAIGE inputs
2026-07-31 20:57   the phenotype file is overwritten
```

**The file on disk is not the file that was used.** Re-running any later step against it today would
not reproduce what was run.

The 556 people present in the consumed file but not in today's all carry status
`excluded_related_ear_phenotype` — they were analysed as controls although the rule excludes them.
The 20:21 version had not applied that exclusion; the 20:57 version does, and was never used.

```
57,507 analysable  -  427  -  4  +  556  =  57,632 consumed
```

So the analysed cohort is short 427 legitimate people in one direction and carries 556 it should not
in the other. The 556 dilute the case-control contrast, which understates associations rather than
inventing them — 1.1% of controls, small but signed.

Detail and the reproducing script: `phase_1/results/FINDINGS.md` Finding 3.

---

## Options

| | Action | Cost |
|---|---|---|
| **A** | Restore the 427 and re-run. Requires only that the covariate build stop requiring imputed PCs. | one SAIGE re-run |
| **B** | Keep the exclusion, and document it in the summary so the two numbers agree. | documentation only |
| **C** | ~~Restore, and first check whether the 517 differ systematically.~~ The check is done — see above. C collapses into A, now with the knowledge that the loss is concentrated in the smallest ancestry group. | — |

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
