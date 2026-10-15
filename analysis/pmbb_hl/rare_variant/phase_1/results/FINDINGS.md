# Phase 1 — the three cohorts

**Run 2026-10-14; tinnitus arm added 2026-10-15** from `PMBB-Release-2026-4.0`. All
three phenotypes built by one script, `scripts/01_cohorts.py`, in a single pass — which
is why the control set is shared by construction rather than by coincidence.

## Result

| | broad | bilateral | tinnitus |
|---|---:|---:|---:|
| cases | 6,752 | 3,164 | 2,859 |
| controls | 50,755 | 50,755 | 50,755 |
| **analysable** | **57,507** | **53,919** | **53,614** |
| excluded, has the phenotype on one date only | 4,007 | 2,695 | 2,911 |
| excluded, other ear evidence | 9,411 | 14,311 | 14,400 |

By ancestry:

| | broad cases | bilateral cases | tinnitus cases | controls |
|---|---:|---:|---:|---:|
| EUR | 5,183 | 2,547 | 2,128 | 37,603 |
| AFR | 1,285 | 490 | 592 | 10,049 |
| EAS | 102 | 43 | 41 | 1,002 |
| AMR | 75 | 37 | 45 | 769 |
| SAS | 70 | 32 | 36 | 800 |
| UNK | 37 | 15 | 17 | 532 |

**AFR is the one place tinnitus is not simply the smallest arm.** It carries 592 cases
against bilateral's 490, so the AFR tinnitus analysis is better powered than the AFR
bilateral one despite tinnitus being smaller overall. Every other ancestry runs the
other way.

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

**There is no row for tinnitus, and that is the point.** No earlier pass ran it, so its
2,859 cases have nothing to be checked against. The rule was instead fixed in writing
before the arm ran — [`PREMISES_tinnitus.md`](../../../PREMISES_tinnitus.md) P1 states
the two-date rule, both sources, and the expected 2,859 / 53,614, so the count could be
compared against something stated in advance rather than against nothing.

## Nesting, asserted rather than assumed — and scoped to the pair it describes

All **3,164** bilateral cases are also broad cases. The script refuses to write anything
if that fails, because bilateral sensorineural hearing loss is a kind of hearing loss —
if a bilateral case did not satisfy the broad rule, one of the two rules would be wrong.

Every later comparison between the arms depends on this. In the first pass it was
checked afterwards, out of curiosity; here it is a precondition.

**Tinnitus breaks that relation, so the assertion must not be extended to it.** Tinnitus
is not a kind of hearing loss; a person can have one without the other. Running the same
refusal over the tinnitus arm would have stopped the pipeline on a correct cohort. The
check is therefore scoped to the one pair it describes, and the tinnitus relation is
**measured** instead:

| | | |
|---|---:|---|
| tinnitus cases | 2,859 | |
| also broad cases | **1,885** | **66%** — the arms overlap heavily |
| not broad cases | 974 | |
|   …with no hearing-loss evidence at all | 581 | tinnitus and nothing else |
|   …with hearing-loss evidence on one date only | 393 | below the two-date rule |

**The 974 and the 581 answer different questions and are easy to confuse.** 974 is "not
a broad *case*"; 581 is "no hearing-loss *evidence*". The 393 in between have hearing
loss recorded once, which is not enough to be a case but is enough that they are not
tinnitus-only people.

What follows from the 66%: **no comparison between the tinnitus arm and a hearing-loss
arm is independent.** Two thirds of the cases and every one of the 50,755 controls are
the same people. A result that looks similar across the two is the shared cases
speaking. Only a result that is notably *stronger* in tinnitus would say anything about
tinnitus, and even then a size-matched control would be needed first — which is
[`PREMISES_tinnitus.md`](../../../PREMISES_tinnitus.md) P6, written before the arm ran.

## The two counts, separated from the start

**57,507 is the analysable cohort.** It is not what SAIGE will test. Phase 3 drops
anyone missing a covariate — 9 people with no `AGE` in the broad combined cohort, all of
them controls, leaving 57,498.

The same gap exists in every arm, and it is always the same people:

| | analysable | tested | dropped |
|---|---:|---:|---:|
| broad combined | 57,507 | 57,498 | 9 |
| bilateral combined | 53,919 | 53,910 | 9 |
| tinnitus combined | 53,614 | 53,605 | 9 |

**It is the same 9 people every time**, not 27 — verified by set comparison across the
three covariate files, and all 9 are controls. That follows from the control set being
shared: the drop is a missing `AGE`, which does not depend on which phenotype is being
built. **No case is lost in any arm.**

Both numbers are emitted at every phase from here on. In the first pass they drifted
into different documents and had to be reconciled after the fact.

## Definitions, in full

**Case, broad:** `SO_396` or any child (`.1` conductive, `.2` sensorineural, `.3` mixed,
`.5` sudden idiopathic, `.8` bilateral, `.9` unilateral) on **≥ 2 distinct dates**.

**Case, bilateral:** `SO_396.2` **and** `SO_396.8` recorded on the **same date**, on
**≥ 2 distinct dates**.

**Case, tinnitus:** `SO_397` or any child, **or** the standard ICD codes (388.3x,
H93.1x) from the OMOP `observation` table, on **≥ 2 distinct dates**. The second source
is not optional — see below.

The same-date requirement is the whole rule. phecodeX splits a diagnosis across two
orthogonal axes — type and laterality — as separate rows, so a person coded
conductive-bilateral at one visit and sensorineural-unilateral at another carries both
`.2` and `.8` without ever having had bilateral sensorineural hearing loss.

**Control, all three arms:** no ear-family evidence from any source. Not "no hearing loss" —
anyone with any `SO_39x` diagnosis, or tinnitus from the OMOP `observation` table, is
excluded rather than used as a control.

**Why `observation` is read.** PMBB v4 moved the standard tinnitus codes (388.3x,
H93.1x) out of the phecode files. Those keep 3,016 of 28,110 tinnitus events; the other
25,094 are only in `observation`. A control set built without that source silently
contains tinnitus patients.

Once tinnitus became a case definition rather than only a control exclusion, the same
fact decides the arm itself: **an arm built on the phecode files alone would miss 89% of
its own phenotype.** The source was already being read for the controls, so this cost
nothing to add — but it would have been a quiet and total failure to get wrong.

## Files

| | |
|---|---|
| `results/cohort_broad.tsv` | IID, PHENO, ancestry |
| `results/cohort_bilateral.tsv` | the same, for the restricted phenotype |
| `results/cohort_tinnitus.tsv` | the same, for tinnitus |
| `results/01_cohorts.json` | counts, by ancestry, the nesting assertion, and the measured tinnitus overlap |
| `data/_ear_family.tsv` | 655,946 rows, cached from the release |
| `data/_tinnitus_obs.tsv` | 25,094 rows, cached from `observation` |
