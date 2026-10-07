# Phase 1 — who is in the study

**Status:** complete · **Scripts:** `00_extract_release_tables.sh` · `01_phenotype_bilateral.py`
**Output:** `cohort.tsv` · `01_phenotype.json` · **Premises:** P1, P2, P3, P4

This is the only phase that differs from the replication it inherits from. Everything downstream is
that pipeline, corrected, unchanged.

---

## The cohort

| | N |
|---|---:|
| exome samples in the release | 70,925 |
| **cases** — bilateral sensorineural, ≥ 2 dates | **3,164** |
| **controls** — no ear evidence of any kind | **50,755** |
| excluded — bilateral sensorineural, but on one date only | 2,695 |
| excluded — hearing loss, never bilateral sensorineural | 4,900 |
| excluded — other ear evidence only | 9,411 |

The first exclusion row is split out deliberately. Those 2,695 people **do** have the target
phenotype; they fail the rule of 2, which is a different reason from not having it. Lumping them
with the 4,900 — as a first version of this did — hides the largest single group the restriction
creates, and they are the obvious population to revisit if the rule of 2 is ever relaxed.

By ancestry, cases / controls: EUR 2,547 / 37,603 · AFR 490 / 10,049 · EAS 43 / 1,002 ·
AMR 37 / 769 · SAS 32 / 800 · unclassified 15 / 532.

---

## Finding 1 — the restriction costs 53% of cases

The replication's definition is `SO_396`, the phecodeX parent, which bundles conductive, mixed,
unilateral and unspecified hearing loss. That reproduces what the pipeline did. It is not the
phenotype the clinical lead considers the phenotype (premise **P1**).

Cases fall **6,752 → 3,164**.

This also settles a discrepancy raised in the 2026-10-02 meeting. Doug counted roughly 4,000 from
the PMBB phenotype browser using the bilateral sensorineural ICD code; the master code table shows
`H90.3 Sensorineural hearing loss, bilateral` on 3,964 people. He and the pipeline were counting
different things, and neither was wrong about its own question.

### Six ICD codes make a case

phecodeX puts two orthogonal axes in the children of `SO_396` — type (`.1` conductive, `.2`
sensorineural, `.3` mixed, `.5` sudden idiopathic) and laterality (`.8` bilateral, `.9` unilateral) —
and one ICD code carries one of each. Requiring `.2` and `.8` together resolves to:

| vocabulary | code | description | people with ≥1 event |
|---|---|---|---:|
| ICD10CM | `H90.3` | Sensorineural hearing loss, bilateral | 3,964 |
| ICD9CM | `389.18` | Sensorineural hearing loss, bilateral | 1,044 |
| ICD9CM | `389.16` | Sensorineural hearing loss, asymmetrical | 439 |
| ICD10CM | `H91.13` | Presbycusis, bilateral | 155 |
| ICD9CM | `389.12` | Neural hearing loss, bilateral | 85 |
| ICD9CM | `389.11` | Sensory hearing loss, bilateral | 26 |

`H91.13 Presbycusis` is worth noticing: age-related hearing loss is the phenotype this study is
ultimately about, and it qualifies.

What does **not** qualify, and would have under the broad definition: `H91.90 Unspecified hearing
loss, unspecified ear` (2,698 people), `H90.5 Unspecified sensorineural hearing loss` (1,011),
`H91.93 Unspecified hearing loss, bilateral` (1,464). Each is missing one of the two axes.

---

## Finding 2 — three readings of the rule, and they disagree

"Bilateral sensorineural on two distinct dates" can be built three ways:

| | reading | cases |
|---|---|---:|
| A | ≥ 2 dates carrying `.2` **or** `.8`, among people who have both | 4,029 |
| B | ≥ 2 dates of `.2` **and** ≥ 2 dates of `.8`, counted separately | 3,259 |
| **C** | **≥ 2 dates on which the same date carries `.2` and `.8`** | **3,164** |

**A is wrong.** A date carrying only `H90.0 Conductive hearing loss, bilateral` contributes `.8` and
counts toward the total, so conductive loss can carry a case over the line — under a definition
written to exclude it.

**B is closer but still wrong.** A sensorineural-unilateral visit and a conductive-bilateral visit
combine into a case who was never diagnosed with bilateral sensorineural loss at all.

**C is what the phrase means**: on at least two distinct dates, the diagnosis recorded was bilateral
sensorineural. C ⊂ B ⊂ A. The 865 people A admits and C does not are the measure of how much the
loose reading costs.

The first implementation here used A. It was caught by checking the three against each other rather
than by anything downstream — nothing further along would have noticed.

---

## Finding 3 — controls were deliberately not touched

A control is still someone with **no ear-family evidence at all**, from `conditions_phecode_x` or
from the OMOP `observation` table. 50,755, unchanged.

Narrowing the cases must not widen the controls. The 3,588 people who are no longer cases still have
hearing loss; moving them into the control group would blur the contrast worse than leaving them
among the cases. They are excluded, which is what the ambiguous middle is for.

The `observation` table matters here for the same reason it did in the replication (premise **P4**):
PMBB v4 moved the standard tinnitus codes (388.3x, H93.1x) out of the phecode files, and reading
only those files leaves 556 people who should be excluded sitting among the controls, invisible.

---

## Carried forward

**AFR has 490 cases.** The replication, with 1,285, already flagged 26 of its top 30 AFR hits as
fragile — each carried by a handful of alleles. At 490 that gets worse. The arm will be run, because
dropping it is its own kind of reporting bias, but it should be read as a null with very little
power rather than as a search.

**Audiograms would replace this phase entirely.** Diagnosis codes are a proxy for an audiogram, and
the clinical lead said plainly that audiograms are the safer bet. They are quantitative, they settle
laterality and type without inference, and they need a REDCap-to-PMBB ID bridge that does not exist
yet. When it does, only this phase changes.

**One axis is not used.** `SO_396.5` sudden idiopathic hearing loss sits on the type axis alongside
conductive, sensorineural and mixed. 203 of the broad cases carry it. Whether sudden idiopathic loss
belongs in a genetic study of age-related hearing loss is a clinical question, not one this phase
should answer on its own.
