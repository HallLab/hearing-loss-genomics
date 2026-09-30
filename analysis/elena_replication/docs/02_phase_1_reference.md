# Phase 1 — Cohort definition · reference record

**Author:** Andre Rico · **Date:** 2026-09-29 · **Status:** complete
**Type:** reference record. This page documents what was done and what it means. It asks for no
decision — the decision arising from it is [page 01](01_sample_frame_decision.md).
**Repo:** `analysis/elena_replication/phase_1/`

---

## 1. In one page

Phase 1 of the PMBB v4 replication checks the cohort: who is in the study, and who counts as having
hearing loss. It is the first phase because every later phase consumes its output.

Two things were checked, and they came out differently.

| | Question | Result |
|---|---|---|
| **The rule** | Who counts as a case? | **Passes.** Re-derived independently, 100% agreement, 70,925 of 70,925 people. |
| **The list** | Who made it into the study? | **Fails.** 427 analysable participants, 40 of them cases, were dropped by a filter that checks the wrong columns. |

The useful part is the split. The rule for deciding who is a case is now confirmed, so a
disagreement in a later phase cannot be blamed on the phenotype. The cohort that was actually
analysed is 427 people smaller than the phenotype supports, and that gap has to be carried forward
as a known quantity.

---

## 2. Why this phase exists at all

The question the whole pipeline asks is: *is there a gene where rare damaging variants are more
common in people with hearing loss than in people without?*

Before that can be asked, two things must be settled: **who has hearing loss** and **who is in the
study**. If either is wrong, every number downstream is wrong, and no amount of care in the
statistics fixes it. Phase 1 is that settling.

It is also why the replication runs in pipeline order rather than starting with the cheapest check.
Had this work begun at Phase 5 — rebuilding the final summary table, which takes minutes — it would
have reported a pass. The table does faithfully reflect the calculations. The calculations were run
on a cohort missing 40 cases, and nothing at Phase 5 can see that.

---

## 3. What Phase 1 read

Everything comes from the institutional release except the pipeline's own outputs, which are read
only as the reference being checked. Nothing in `analysis/elena/` was modified.

**From the release** — `/static/PMBB/PMBB-Release-2026-4.0/`

| Path | What it provides |
|---|---|
| `Phenotype/4.0/..._conditions_phecode_x.txt` | diagnosis codes with dates (12 GB) |
| `Phenotype/4.0/..._observation.txt` | the OMOP observation table (1.2 GB) — where v4 moved tinnitus |
| `Phenotype/4.0/..._covariates.txt` | age, sex, batch, and both PC families |
| `Exome/PCA/combined/..._samples_ancestries.tsv` | who has exome data (70,925 people) |

**From the pipeline** — `analysis/elena/rarevariantExWAS/`

| Path | Owner | What it is |
|---|---|---|
| `PMBBv4_phecodex/PMBBv4_hearing_tinnitus_all_statuses.csv.gz` | nikkipal | every person's case/control status |
| `PMBBv4_phecodex/PMBBv4_hearing_tinnitus_summary.csv` | nikkipal | the published counts |
| `PMBBv4_phecodex/hearing_impairment_PMBBv4_SAIGE.txt` | nikkipal | the file handed to SAIGE |
| `HL_TIN_PMBBv4_keep.txt` | elenas18 | full cohort, 70,925 |
| `HL_TIN_PMBBv4_SAIGE_samples.txt` | elenas18 | after the filter, 70,408 |
| `covariates/covariates_*_withBatch.txt` | elenas18 | the covariates SAIGE consumed |

---

## 4. Check 01 — the sample frame

**Script:** `phase_1/scripts/01_sample_frame.py`

Traces the cohort from the release to the file SAIGE actually received, and attributes every drop.

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

See §5b — this distinction was missed on the first pass and is a finding in its own right.

The 427 belong to a group of **517** lost between `HL_TIN_PMBBv4_keep.txt` (70,925) and
`HL_TIN_PMBBv4_SAIGE_samples.txt` (70,408). Checked against the release, all 517 have exome data,
have complete `exome_PC1-20`, and have **no** `imputed_PC1-20`.

### The mechanism, stated carefully

This is the part worth getting right, because a loose reading makes it sound worse than it is.

The release publishes **two separate families of principal components**, each numbered from 1:

```
exome_PC1   … exome_PC20     from exome sequencing
imputed_PC1 … imputed_PC20   from array genotyping + imputation
```

These are two different PCA runs on two different assays. It is **not** a "first 6 versus extended
20" distinction inside one family — both families go to 20.

**The model used the exome PCs.** Verified by value, not by column name: the covariate files handed
to SAIGE label their columns generically (`PC1 … PC5`), so the name carries no information. Matched
against the release, `PC1` is numerically identical to `exome_PC1` (r = 0.99965 over 70,404 people;
`imputed_PC1` correlates 0.992 but the values differ). **The imputed PCs never entered the
regression.**

They acted as a **gate**. The delivered sample list is an exact set identity:

```
analysable (57,507)  ∩  {has all 20 imputed PCs}  =  57,080
```

Tested as set equality, not as a matching count. The two conditions — absent from the imputed `.fam`, and missing every imputed PC — pick out the
same people, because both follow from having no imputed genotypes.

So: **the adjustment is correct; the filter is what is wrong.** An exome analysis adjusted with
array-derived ancestry components would be a methodological error. That is not what happened.

### 517 and 427 are different facts

All 517 have the same problem — no imputed PCs. Nobody in the release is missing exome PCs (checked:
zero). The 427 is about **phenotype status**, not about which PCs were absent:

```
517  lacking imputed PCs
 ├── 387  controls  ┐
 ├──  40  cases     ┘── 427  analysable       →  real loss
 ├──  64  already excluded, related ear phenotype  ┐
 └──  26  already excluded, one date               ┘── 90  already out anyway
```

The 90 had been excluded by the phenotype rules before any PC was considered, so removing them
changed nothing. **427 is the number that matters.**

### Why it shows up as two different published numbers

Both columns below are the pipeline's own numbers. Neither is this replication.

| | `summary.csv` | phenotype file on disk |
|---|---:|---:|
| cases | 6,752 | **6,712** |
| controls | 50,755 | **50,368** |

The summary was written at the phenotype step; the filter came afterwards and nobody rewrote it.
Anyone citing the summary reports 6,752 — a number the analysis did not use.

Neither column is what the run consumed either. That was a third set of numbers — §5b.

---

## 5. Check 02 — cases and controls

**Script:** `phase_1/scripts/02_rebuild_cases_controls.py`

Cases and controls were re-derived from the release, written from its own files without reference to
the existing code, then compared person by person.

### The rules, as the pipeline documents them

| Status | Condition |
|---|---|
| `case` | target evidence on **≥ 2 distinct dates** (the "rule of 2") |
| `excluded_one_date` | target evidence on exactly 1 date |
| `excluded_related_ear_phenotype` | no target evidence, but some other ear-family evidence |
| `control` | no ear-family evidence of any kind |

Target = hearing impairment (`SO_396`). Ear family = any `SO_39x` code, plus tinnitus from the
observation table (see below).

### Result

| | pipeline | rebuilt |
|---|---:|---:|
| case | 6,752 | **6,752** |
| excluded, one date | 4,007 | **4,007** |
| excluded, related ear | 9,411 | **9,411** |
| control | 50,755 | **50,755** |

**70,925 / 70,925 — 100%, person for person.** Not just matching totals: the same people.

### Sub-finding A — the pipeline was right and the rebuild was wrong

The first rebuild used only `conditions_phecode_x` and disagreed on **558** people, all in the
control / related-ear boundary, all in the direction of the rebuild seeing *less* evidence.

Cause: **PMBB v4 moved the standard tinnitus ICD codes (388.3x, H93.1x) into the OMOP `observation`
table.** The phecode files retain only the ~3,016 pulsatile-tinnitus events; the real bulk — 25,094
events — sits in `observation.txt`. All 558 (100%) carry tinnitus there.

Adding that source took agreement from 99.21% to 100%. The pipeline had already handled this. The
relocation is documented nowhere in the release, and it is the single easiest way to get a v4
phenotype wrong.

### Sub-finding B — child-phecode rollup makes no difference here

PheWAS convention rolls a child phecode up into its parent, and `SO_396` has five children with real
volume (`.2` = 81,319 events, `.8` = 78,903, plus `.1`, `.3`, `.9`). Both definitions were computed:

- `SO_396` exact
- `SO_396` plus all children

**Identical results.** No child adds a case, which means the source already records parent and child
for the same encounter. The question can be closed rather than argued in a meeting.

### Independent corroboration

Two extraction counts match the pipeline's own notebook exactly, from separately written code:

| | this replication | pipeline notebook |
|---|---:|---:|
| ear-family PhecodeX rows | 655,946 | 655,946 |
| tinnitus events in `observation` | 25,094 | 25,094 |

---

## 5b. Check 03 — what the run actually consumed

**Script:** `phase_1/scripts/03_what_saige_actually_used.py`

### Scope, and what kind of evidence this is

Phase 1 asks who was in the study. That cannot be answered from the phenotype file, because the
phenotype file was overwritten after being read. So this check follows the cohort to the file the run
consumed. It makes **no claim about Phase 4** — it inspects that file's membership, never the
statistics computed from it. Phase 4 remains unexamined.

**This check is an audit, not a replication, and it carries less weight than the other two.** Checks
01 and 02 re-derive from the institutional release: 02 rebuilds the cohort from scratch, and 01 uses
the release covariates to establish independently *why* the 517 were dropped. If the pipeline
vanished tomorrow, both conclusions would still stand. Check 03 re-derives nothing. Its substantive
inputs are three pipeline files compared against each other; the one release file it reads is used
only as a sanity check that everyone in the consumed file has exome data. So the evidence here is
entirely the pipeline's own — the question it answers is Phase 1's.

Separating what is proven from what is inferred:

| | Basis | Strength |
|---|---|---|
| The run consumed a different file than the one on disk | row counts (57,632 vs 57,080) and the build log recording `Phenotype samples: 57,636` | **solid** — no timestamps involved |
| 556 excluded people were analysed as controls | membership comparison plus their status in `all_statuses` | **solid** |
| The phenotype file was *overwritten*, 36 minutes later | file modification times | **inferred** — a copy, a `touch` or a restore would change mtime without changing content |

Phase 1's conclusion rests on the solid rows. The overwrite is the most likely explanation, not a
proven sequence, and the question of what actually happened goes to the people who ran it.

### What it found

| | N | cases | controls |
|---|---:|---:|---:|
| the file the run consumed | **57,632** | 6,712 | **50,920** |
| the phenotype file on disk | 57,080 | 6,712 | 50,368 |

An inner join cannot increase N, so these are not the same file. The build log resolves it —
`Phenotype samples: 57,636`:

```
20:21   the covariate build reads the phenotype file and writes the SAIGE inputs
20:57   the phenotype file is overwritten
```

**The artifact on disk is not the artifact that was used.** Re-running any later step against today's
phenotype file would not reproduce what was run, and nothing in either file says so.

### Two defects, opposite directions

Finding 1 holds: all 427, including all 40 cases, are absent from the consumed file, and the case
sets of both files are identical at 6,712.

But the 556 people present in the consumed file and not in today's all carry status
`excluded_related_ear_phenotype`. **They were analysed as controls although the pipeline's own rule
excludes them** — they have ear disease that is not hearing impairment, so they are neither a clean
case nor a clean control. The 20:21 phenotype file had not applied that exclusion; the 20:57 one does,
and was never used.

```
57,507 analysable  -  427  -  4  +  556  =  57,632 consumed
```

The run is not a subset of the analysable cohort. It is short in one direction and long in the other.
Bias from the 556 dilutes the case-control contrast, so it understates associations rather than
inventing them — 1.1% of controls, small but signed.

---

## 5c. Check 04 — ancestry stratification, and who the excluded are

**Script:** `phase_1/scripts/04_step5_and_dropped_profile.py`

Phase 1 is Step 1 **and Step 5** by the mapping in
[`pipeline_plan.md`](../pipeline_plan.md) §2. Checks 01–03 covered only Step 1. This closes the
phase.

### Step 5 verifies

The stratification is internally consistent. The EUR and AFR cohorts contain only participants the
release classifies as EUR and AFR; `combined` is exactly EUR + AFR + the smaller groups.

### Both arms, side by side

The strata were computed for each arm — the reproduction cohort that actually ran, and the corrected
cohort the phenotype supports. Nobody had computed the second.

| Stratum | reproduction (57,632) | corrected (57,507) | change |
|---|---|---|---|
| combined | 57,632 · 6,712 cases | 57,507 · **6,752** cases | −125 people, **+40 cases** |
| EUR | 43,016 · 5,160 | 42,786 · 5,183 | −230, +23 |
| AFR | 11,387 · 1,279 | 11,334 · 1,285 | −53, +6 |

**The corrected arm is smaller and has more cases.** It gains the 431 who were never let in (40 of
them cases) and loses the 556 controls the rules exclude. Case rate rises about 0.1 percentage
points in every stratum — a small, consistent shift, which is what one would expect from correcting
a defect that is real but not large.

Ancestry composition moves in one place and one place only:

| | reproduction | corrected | change |
|---|---:|---:|---:|
| EUR | 43,016 | 42,786 | −0.5% |
| AFR | 11,387 | 11,334 | −0.5% |
| **EAS** | **939** | **1,104** | **+17.6%** |
| SAS | 874 | 870 | −0.4% |
| AMR | 846 | 844 | −0.2% |

Correcting the filter restores **165 East Asian participants**, growing that stratum by nearly a
fifth while every other group changes by half a percent. That is Finding 4 seen from the other
direction: the exclusion was concentrated, so the correction is concentrated too.

### Finding 4 — the excluded are not a random 0.7%

Exclusion rate within each group, against 0.73% cohort-wide:

| Group | excluded | total | rate |
|---|---:|---:|---:|
| **EAS** | **208** | **1,333** | **15.60%** |
| AMR | 16 | 1,039 | 1.54% |
| SAS | 5 | 1,080 | 0.46% |
| EUR | 233 | 51,867 | 0.45% |
| AFR | 49 | 14,927 | 0.33% |

One in six East Asian participants, twenty-one times the cohort rate. The same skew appears on sex
(1.01% of women against 0.43% of men) and strongly on batch (1.09% of batch 1 against 0.11–0.18% of
batches 2 and 3); the excluded are younger and enrolled earlier.

Batch and enrolment point at a technical coverage gap — exome-sequenced early, array genotyping never
completed for a subset — rather than anything about the participants. That is the pattern; the cause
is not this replication's to establish.

### Declared out of scope

Step 1 produces **two** phenotypes. Checks 01–04 cover hearing impairment only. The tinnitus
phenotype, and the combined hearing-loss-and/or-tinnitus phenotype decided on 2026-07-01, are
unverified. Declared here rather than left implicit, because the combined phenotype will need the
same scrutiny when it is built.

---

## 6. Verdict

The criterion was: the re-derived case, control **and sample** sets match person-for-person, or every
difference is explained by a documented rule.

| Component | Verdict | Basis |
|---|---|---|
| Phenotype definition | **passes** | exact reproduction, 100% |
| Sample frame | **does not pass** | 427 dropped by an undocumented filter on unused columns |

Every other Phase 1 filter — rule of 2, related-ear exclusion — reproduced exactly. The PC
cohort was built by merging against the imputed `.fam`, and that is the only step that fails.

---

## 7. Every number in one place

These are Phase 1's numbers only. Each label says whose artifact the number comes from and whether
that artifact was actually used, because three different cohort counts are in circulation and the
difference between them is Finding 3.

**Cohort sizes, in the order they arose**

| Number | What it is | Used by the run? |
|---:|---|---|
| 70,925 | people with exome data in the PMBB v4 release | — the starting population |
| 57,507 | analysable for hearing impairment: case or control. **What the phenotype supports** — reproduced exactly by this replication | no |
| 57,632 | the covariate files the SAIGE run consumed. Built 2026-07-31 **20:21** (Elena) from the then-current phenotype file | **yes — this is what ran** |
| 57,080 | the phenotype file on disk today. Written **20:57** (Nikki), 36 minutes after the run's inputs had already been built from the earlier version | no — never reached SAIGE |

**Case and control counts**

| Number | What it is |
|---:|---|
| 6,752 | cases the phenotype supports — reproduced exactly |
| 6,712 | cases the run used; identical in both versions of the phenotype file |
| 50,920 | controls the run used |
| **431** | analysable people absent **from the run** — 40 cases + 391 controls. The operative number |
| 427 | analysable people absent from the **phenotype file on disk** — 40 cases + 387 controls (Finding 1). Not the same 431: it is 431 minus the 4 below |
| 4 | legitimate controls on the disk file but absent from the run |
| 556 | people the phenotype rules exclude, present in the run as controls (Finding 3) |
| 208 | East Asian participants excluded — 15.60% of the 1,333 in the cohort, against 0.73% cohort-wide (Finding 4) |
| 517 | people in the release with no imputed genotypes — the 427 plus 90 already excluded on phenotype grounds |
| 4,007 | excluded: target evidence on one date only |
| 9,411 | excluded: other ear-family evidence but not hearing impairment |

**Extraction and diagnostic counts**

| Number | What it is |
|---:|---|
| 558 | first-pass disagreements in check 02, all explained by the observation table |
| 655,946 | ear-family PhecodeX rows extracted from the release |
| 25,094 | tinnitus events in the OMOP observation table |
| 134,917 | `SO_396` hearing-impairment events |

The cohort sizes reconcile as `57,507 − 431 + 556 = 57,632`, or equivalently
`57,507 − 427 − 4 + 556`.

**The two arms, by stratum** (§5c)

| Stratum | reproduction 57,632 | corrected 57,507 | change |
|---|---|---|---|
| combined | 57,632 · 6,712 cases | 57,507 · 6,752 cases | −125 people, +40 cases |
| EUR | 43,016 · 5,160 | 42,786 · 5,183 | −230, +23 |
| AFR | 11,387 · 1,279 | 11,334 · 1,285 | −53, +6 |

The corrected arm is smaller and carries more cases. Case rate rises ~0.1 percentage points in every
stratum.

**Ancestry composition, the two arms**

| Number | What it is |
|---:|---|
| 939 → 1,104 | East Asian participants, reproduction arm → corrected arm |
| **165** | East Asian participants the correction restores — **+17.6%** of that stratum, against −0.5% or less for every other group |
| 1,333 | East Asian participants in the exome cohort in total |

### Which number goes to the next phase

What Phase 1 establishes is that **57,507** is the analysable cohort the phenotype supports, and that
what the run consumed was **57,632** — a different set, not a subset (§5b).

**Both numbers are carried forward, deliberately.** From Phase 3 onward the replication runs two
declared arms — reproduction on the 57,632 that actually ran, and corrected on the 57,507 the
phenotype supports — so that a later difference can be attributed to the defect rather than to a
cohort the replication changed. The design is in
[`pipeline_plan.md`](../pipeline_plan.md) §6; the decision about which cohort the *group's* analysis
should use is question 2 of [page 01](01_sample_frame_decision.md) and is not ours.

Which cohort each later phase inherits beyond that is that phase's question, not this one's, and is
deliberately not answered here. A provisional note peeking at variant-QC files suggested Phase 2 may be
cohort-independent; that was removed, because a claim about Phase 2 that has not been through Phase
2's own check would sit in this page unwatched and go stale. Phase 2 will establish it.

---

## 8. Scripts and outputs

All paths relative to `analysis/elena_replication/phase_1/`.

| File | Tracked in git | Contents |
|---|---|---|
| `scripts/00_extract_release_tables.sh` | yes | caches the two multi-GB release extracts |
| `scripts/01_sample_frame.py` | yes | traces the cohort chain, identifies the 517 |
| `scripts/02_rebuild_cases_controls.py` | yes | re-derives cases/controls, diffs against the pipeline |
| `scripts/03_what_saige_actually_used.py` | yes | compares the consumed file against the phenotype file on disk |
| `scripts/04_step5_and_dropped_profile.py` | yes | Step 5 for both arms; profiles the excluded |
| `results/FINDINGS.md` | yes | the full write-up, including what was not established |
| `results/01_sample_frame.json` | yes | every number in §4, machine-readable |
| `results/03_what_saige_actually_used.json` | yes | the §5b comparison, machine-readable |
| `results/03_improperly_included_controls.csv` | no — per-person IDs | the 556 |
| `results/04_step5_and_dropped_profile.json` | yes | the §5c tables, machine-readable |
| `results/02_rebuild_cases_controls.json` | yes | counts, agreement, confusion matrix |
| `results/01_dropped_participants.csv` | no — per-person IDs | the 427, with status |
| `results/02_disagreements_*.csv` | no — per-person IDs | now empty (header only): zero disagreements is the result |
| `data/_ear_family.tsv` | no — 24 MB extract | cached `SO_39x` rows, regenerable |
| `data/_tinnitus_obs.tsv` | no — 0.9 MB extract | cached observation tinnitus, regenerable |

Per-person files and cached extracts are deliberately untracked; the scripts regenerate them.

### How to re-run

```bash
cd /project/hall/analysis/hearing-loss-genomics/analysis/elena_replication/phase_1

bash scripts/00_extract_release_tables.sh          # only if data/ is empty; ~10 min
../../../venv/bin/python3 scripts/01_sample_frame.py
../../../venv/bin/python3 scripts/02_rebuild_cases_controls.py
../../../venv/bin/python3 scripts/03_what_saige_actually_used.py
../../../venv/bin/python3 scripts/04_step5_and_dropped_profile.py
```

Script 00 caches the two release extracts and prints the expected row counts (655,946 and 25,094) so
a silent change in the release shows up immediately. Scripts 01 and 02 take seconds once the cache
exists.

---

## 9. Concepts, for coming back to this later

| Term | Plain meaning |
|---|---|
| **case / control** | Someone the records say has the condition, versus someone they say does not. Getting the boundary right is most of Phase 1. |
| **rule of 2** | A diagnosis must appear on **two separate dates** to count. One mention could be a query, a rule-out, or a coding slip; two separate encounters is evidence. |
| **phecode** | A grouping of raw ICD billing codes into something clinically meaningful. `SO_396` = hearing impairment. `SO_39x` is the whole ear family. |
| **ear family** | Any ear-related code. Used for exclusion: someone with ear problems but *not* hearing impairment is neither a clean case nor a clean control, so they are set aside. |
| **principal components (PCs)** | Numbers summarising a person's genetic ancestry. They stop the analysis confusing *"this gene causes hearing loss"* with *"this gene is commoner in a population that happens to have more hearing loss."* Necessary — but they must come from the same assay as the data being tested. |
| **the two PC families** | `exome_PC*` from exome sequencing; `imputed_PC*` from array + imputation. Two separate PCA runs, both numbered 1–20. This study uses `exome_PC1-6`. |
| **covariate** | Something adjusted *for* in the model — age, sex, batch, PCs. |
| **gate** | Something used to decide who gets *in*. The bug here: imputed PCs were a gate, never a covariate. |
| **sample frame** | The set of people actually analysed, as opposed to the set the phenotype defines. When the two differ and nobody says so, published counts stop matching the analysis. |
| **OMOP `observation` table** | One of the standard tables in the PMBB data model. In v4, tinnitus ICD codes live here rather than with the other diagnoses — the trap in sub-finding A. |

---

## 10. Still open

| # | Question | Who |
|---|---|---|
| 1 | Was the PC filter deliberate, and which step introduced it? | Nikki, Elena — [page 01](01_sample_frame_decision.md) |
| 2 | Do the 427 come back? | Molly, Doug, Nikki |
| 3 | ~~Do the 517 differ systematically?~~ **Resolved, check 04 — they do.** 15.60% of East Asian participants were excluded against 0.73% cohort-wide, with the same skew on sex and strongly on batch. The loss is concentrated, not neutral. | — |
| 4 | Why was the phenotype file rewritten 36 minutes after the covariates were built from it, and did anyone know? | Nikki, Elena |

Cross-phase observations — things noticed about a later phase while working on this one — are not
recorded here. They go to the parent [`README.md`](../README.md) §6, so that no phase's page carries
claims its own checks did not establish.
