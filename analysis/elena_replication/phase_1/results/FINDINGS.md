# Phase 1 — who is in the study

**Status:** first check complete · **Date:** 2026-09-29
**Script:** [`../scripts/01_sample_frame.py`](../scripts/01_sample_frame.py)
**Outputs:** [`01_sample_frame.json`](01_sample_frame.json), [`01_dropped_participants.csv`](01_dropped_participants.csv)

---

## Finding 1 — 517 exome-sequenced participants were excluded for lacking PCs the analysis never uses

The cohort shrinks in a place nothing documents:

```
70,925   exome cohort in the release
70,925   phenotyped (all_statuses covers the whole cohort)
57,507   case or control for hearing impairment
57,080   the phenotype file on disk today
   427   dropped between the two -- 40 cases, 387 controls
```

The 427 are part of a larger group of **517** lost between `HL_TIN_PMBBv4_keep.txt` (70,925) and
`HL_TIN_PMBBv4_SAIGE_samples.txt` (70,408). Checked against the release, all 517:

- **have exome data** — every one is in the exome ancestry list;
- **have complete `exome_PC1-6`** — zero missing;
- **have no `imputed_PC1-20`** — all 20 missing, for all 517.

They are the entire set of people in the release with exome PCs but no imputed PCs: 517 of 70,925,
0.7% of the cohort.

**This is an exome study.** 517 participants — among them **40 hearing-loss cases** — were removed
from an exome analysis for missing array-imputation principal components that the analysis does not
consume.

### The model is adjusted correctly; the filter is what is wrong

The distinction matters, and the milder reading is the correct one.

The covariates the model consumes are **exome** PCs, verified by value rather than by column name.
The covariate files handed to SAIGE label their columns generically (`PC1 … PC5`), so the name
carries no information. Matched against the release: `PC1` is numerically identical to `exome_PC1`
(correlation 0.99965 across 70,404 people; `imputed_PC1` correlates 0.992, but the values differ —
e.g. 0.0003 vs 0.0004 for the same person). The imputed PCs never enter the regression.

They act as a **gate**. The delivered sample list is an exact set identity:

```
analysable  ∩  {has all 20 imputed PCs}   =   57,080
  57,507              (excludes 427)
```

Tested as set equality, not as a count match. The effect is what a completeness check spanning every
PC column in the covariates file would produce, rather than only the six the model uses.

So this is **not** an exome analysis adjusted with array-derived ancestry components — that would be
a methodological error, and it is not what happened. The adjustment is correct. A completeness
requirement over unused columns silently reduced the cohort.

Note that the release publishes `exome_PC1-20` and `imputed_PC1-20` — two separate PCA runs on two
separate assays, each numbered from 1. This is not a 6-versus-20 distinction within one family.

### Why it matters beyond the count

The published summary and the analysed cohort disagree, and nothing says so:

| | `PMBBv4_hearing_tinnitus_summary.csv` | phenotype file on disk |
|---|---:|---:|
| hearing-loss cases | 6,752 | **6,712** |
| controls | 50,755 | **50,368** |

Neither column is what the run consumed — that was a third set of numbers, 57,632 people with
50,920 controls. See Finding 3.

A reader citing the summary would state 6,752 cases. The analysis ran on 6,712. The summary is not
wrong about the phenotype — it is describing a different sample frame than the one used, with no
note that the two differ.

### What is not yet established

Whether the exclusion was deliberate. Requiring complete imputed PCs is a defensible default if a
covariate template was reused from an imputed-genotype analysis, and 0.7% is the kind of loss that
passes unnoticed. The check here establishes the mechanism and the cost, not the intent. Resolving
that means asking whoever built `HL_TIN_PMBBv4_SAIGE_samples.txt`.

### Consequence for the replication

The re-derived cohort should intersect with the **exome** sample list and require only the
covariates the model actually uses. On that frame the analysable set is 57,507, not 57,080. Any
later phase that compares against the pipeline has to hold this difference constant or it will
attribute a cohort difference to something else.

---

## Finding 2 — the case/control rule reproduces exactly

**Script:** [`../scripts/02_rebuild_cases_controls.py`](../scripts/02_rebuild_cases_controls.py)
**Output:** [`02_rebuild_cases_controls.json`](02_rebuild_cases_controls.json)

Cases and controls were re-derived from the release independently — written from the release's own
files without consulting the pipeline's code — and compared person-by-person.

| | pipeline | rebuilt |
|---|---:|---:|
| case | 6,752 | **6,752** |
| excluded_one_date | 4,007 | **4,007** |
| excluded_related_ear_phenotype | 9,411 | **9,411** |
| control | 50,755 | **50,755** |

**Agreement: 70,925 / 70,925 — 100%, person for person.** The phenotype definition is correct and
independently confirmed.

### Two things the check settled on the way

**The ear-family evidence needs the `observation` table, and this is where it bites.** A first pass
built ear-family evidence from `conditions_phecode_x` alone and disagreed with the pipeline on 558
people — all in the control / related-ear boundary, all in the direction of the rebuild seeing
*less* evidence than the pipeline. Every one of those 558 (100%) carries tinnitus in the OMOP
`observation` table. PMBB v4 relocated the standard tinnitus ICD codes (388.3x, H93.1x) there; the
phecode files retain only the ~3,016 pulsatile events. Adding that source took agreement from
99.21% to 100%.

This is a case where the pipeline was right and the replication was incomplete. That is worth
stating plainly: the check runs in both directions, and the v4 relocation is a genuine trap that
the pipeline had already avoided.

**Child-phecode rollup does not matter here.** PheWAS convention rolls a child phecode up into its
parent, and `SO_396` has five children with real volume (`.2` at 81,319 events, `.8` at 78,903,
plus `.1`, `.3`, `.9`). Both definitions were computed. They give **identical** results — no child
adds a case, which indicates the source records parent and child for the same encounter. The
question can be closed rather than argued.

### Independent corroboration of the extracts

Two counts produced here match the pipeline's own notebook exactly, from separate extractions:
ear-family PhecodeX rows **655,946**, and standard-tinnitus observation events **25,094**.

---

## Finding 3 — the phenotype file on disk is not the one the analysis used

**Script:** [`../scripts/03_what_saige_actually_used.py`](../scripts/03_what_saige_actually_used.py)
**Output:** [`03_what_saige_actually_used.json`](03_what_saige_actually_used.json)

### Scope of this check

Phase 1 asks who was in the study. That cannot be answered from the phenotype file alone, so this
check follows the cohort to the file the run consumed. It makes **no claim about Phase 4** — nothing
here examines the statistics computed from that file, only its membership. Phase 4 stays unexamined.

### What happened

| | N | cases | controls |
|---|---:|---:|---:|
| file the SAIGE run consumed | **57,632** | 6,712 | **50,920** |
| phenotype file on disk today | 57,080 | 6,712 | 50,368 |

An inner join cannot increase N, so these cannot be the same file. The build log settles it: the
merge that produced the consumed covariates records `Phenotype samples: 57,636`.

```
20:21   the covariate build reads the phenotype file (57,636 rows) and writes the SAIGE inputs
20:57   the phenotype file is overwritten (57,080 rows)
```

Thirty-six minutes. **The artifact on disk is not the artifact that was used**, and neither file says
so. This is a reproducibility problem independent of the counts: re-running any later step against
today's phenotype file would not reproduce what was run.

### Two defects, in opposite directions

**Finding 1 stands.** All 427 dropped participants, including all 40 cases, are absent from the
consumed file. The case sets of the two files are identical at 6,712. The loss is real.

**A second defect runs the other way.** The difference between the two files is not random:

```
only in the consumed file :  556  ->  all with status excluded_related_ear_phenotype
only in today's file      :    4  ->  legitimate controls
```

**556 people the pipeline's own rule excludes were analysed as controls.** They have ear-family
evidence that is not hearing impairment, so they are neither a clean case nor a clean control — the
rule exists precisely to set them aside. The 20:21 version of the phenotype file had not applied that
exclusion; the 20:57 version does, and was never used.

### What the analysis actually ran on

```
 6,712 cases      --  40 missing that the phenotype supports
50,920 controls   -- 556 present that the phenotype rules exclude, 4 legitimate ones absent
```

Arithmetic: 57,507 analysable - 427 - 4 + 556 = 57,632. The run is not simply a subset of the
analysable cohort; it is smaller in one direction and larger in the other.

Direction of bias from the 556: including people with ear disease among controls dilutes the
case-control contrast, biasing toward the null. At 1.1% of controls the effect is small, but it
understates associations rather than inventing them.

### Not established

Why the phenotype file was rewritten, and whether anyone knew the covariates had already been built
from the earlier version. A question for Nikki Palmiero and Elena, not something the files answer.

---

## Phase 1 verdict

The criterion in [README](../../README.md) §2 is that the re-derived case, control **and sample**
sets match person-for-person, or that every difference is explained by a documented rule.

- **Phenotype definition — passes.** Exact reproduction, 100% agreement (Finding 2).
- **Sample frame — does not pass.** 427 analysable participants, 40 of them cases, are dropped by
  an undocumented requirement for imputed PCs that the analysis never uses (Finding 1).
- **Provenance of what was run — does not pass.** The phenotype file was overwritten 36 minutes
  after the covariates were built from it, and 556 people the rules exclude were analysed as
  controls (Finding 3).

The split matters for what comes next. The rule for deciding who is a case is sound, so a
disagreement in a later phase should not be attributed to the phenotype. But the cohort that was
analysed is **not** simply a subset of what the phenotype supports — it is 431 people short in one
direction and 556 too many in the other, netting 57,632 against an analysable 57,507. Later phases
have to carry that as a known, signed difference rather than rediscover it.

---

## Note on ownership

`rarevariantExWAS/PMBBv4_phecodex/` is owned by `nikkipal`; `HL_TIN_PMBBv4_keep.txt`,
`HL_TIN_PMBBv4_SAIGE_samples.txt` and the covariate builds are owned by `elenas18`.

The PC filter of Finding 1 is present on both sides, and the files do not say which step introduced
it. Finding 3 spans both by construction: a file owned by one was overwritten after a build owned by
the other had read it. **Both findings are addressed to Nikki and Elena together**, not to either
alone — see [`../../docs/01_sample_frame_decision.md`](../../docs/01_sample_frame_decision.md).

---

## Note on ordering

This finding is the case for running the replication in pipeline order. Had the work started at
Phase 5 — rebuilding the summary table, the cheapest check — it would have confirmed that the table
faithfully reflects the per-stratum outputs and returned a pass. Those outputs were computed on a
cohort missing 40 cases, and nothing at Phase 5 can see that.
