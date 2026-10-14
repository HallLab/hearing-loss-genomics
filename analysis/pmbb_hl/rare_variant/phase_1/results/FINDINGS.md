# Phase 1 — the two cohorts

**Run 2026-10-14** from `PMBB-Release-2026-4.0`. Both phenotypes built by one script,
`scripts/01_cohorts.py`.

## Result

| | broad | bilateral |
|---|---:|---:|
| cases | 6,752 | 3,164 |
| controls | 50,755 | 50,755 |
| **analysable** | **57,507** | **53,919** |
| excluded, has the phenotype on one date only | 4,007 | 2,695 |

By ancestry:

| | broad cases | bilateral cases | controls |
|---|---:|---:|---:|
| EUR | 5,183 | 2,547 | 37,603 |
| AFR | 1,285 | 490 | 10,049 |
| EAS | 102 | 43 | 1,002 |
| AMR | 75 | 37 | 769 |
| SAS | 70 | 32 | 800 |
| UNK | 37 | 15 | 532 |

Of 70,925 people with exome data, **20,170 carry some ear-family evidence** and are
therefore not controls.

## It reproduces the first pass

Written from scratch against the release, not ported:

| | rebuild | first pass |
|---|---:|---:|
| broad analysable | 57,507 | 57,507 |
| broad cases | 6,752 | 6,752 |
| bilateral analysable | 53,919 | 53,919 |
| bilateral cases | 3,164 | 3,164 |
| bilateral excluded on one date | 2,695 | 2,695 |

Two independent implementations of the same two rules agreeing to the person is the
evidence that the rules themselves are stated clearly enough to be reimplemented. That
was not guaranteed: the bilateral rule in particular is easy to get wrong, and the loose
reading of it gives 4,029 cases rather than 3,164.

## Nesting, asserted rather than assumed

All **3,164** bilateral cases are also broad cases. The script refuses to write anything
if that fails, because bilateral sensorineural hearing loss is a kind of hearing loss —
if a bilateral case did not satisfy the broad rule, one of the two rules would be wrong.

Every later comparison between the arms depends on this. In the first pass it was
checked afterwards, out of curiosity; here it is a precondition.

## The two counts, separated from the start

**57,507 is the analysable cohort.** It is not what SAIGE will test. Phase 3 drops
anyone missing a covariate — in the first pass that was 9 people with no `AGE` in the
broad combined cohort, all of them controls, leaving 57,498.

Both numbers are emitted at every phase from here on. In the first pass they drifted
into different documents and had to be reconciled after the fact.

## Definitions, in full

**Case, broad:** `SO_396` or any child (`.1` conductive, `.2` sensorineural, `.3` mixed,
`.5` sudden idiopathic, `.8` bilateral, `.9` unilateral) on **≥ 2 distinct dates**.

**Case, bilateral:** `SO_396.2` **and** `SO_396.8` recorded on the **same date**, on
**≥ 2 distinct dates**.

The same-date requirement is the whole rule. phecodeX splits a diagnosis across two
orthogonal axes — type and laterality — as separate rows, so a person coded
conductive-bilateral at one visit and sensorineural-unilateral at another carries both
`.2` and `.8` without ever having had bilateral sensorineural hearing loss.

**Control, both arms:** no ear-family evidence from any source. Not "no hearing loss" —
anyone with any `SO_39x` diagnosis, or tinnitus from the OMOP `observation` table, is
excluded rather than used as a control.

**Why `observation` is read.** PMBB v4 moved the standard tinnitus codes (388.3x,
H93.1x) out of the phecode files. Those keep 3,016 of 28,110 tinnitus events; the other
25,094 are only in `observation`. A control set built without that source silently
contains tinnitus patients.

## Files

| | |
|---|---|
| `results/cohort_broad.tsv` | IID, PHENO, ancestry |
| `results/cohort_bilateral.tsv` | the same, for the restricted phenotype |
| `results/01_cohorts.json` | counts, by ancestry, plus the nesting assertion |
| `data/_ear_family.tsv` | 655,946 rows, cached from the release |
| `data/_tinnitus_obs.tsv` | 25,094 rows, cached from `observation` |
