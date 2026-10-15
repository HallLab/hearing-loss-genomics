# Phase 1 — Who is in the study

**Bilateral sensorineural hearing loss · PMBB Release 2026-4.0**
Part 1 of 5 · [Phase 2 — variants →](02_variants.md)

---

## In one line

**3,164 cases and 50,755 controls**, from 70,925 people with exome data — a phenotype restricted to
bilateral sensorineural hearing loss on the clinical lead's criteria.

---

## The cohort

| | N | rule |
|---|---:|---|
| **Cases** | **3,164** | bilateral sensorineural, on ≥ 2 distinct dates |
| **Controls** | **50,755** | no ear-family evidence of any kind |
| Excluded — right phenotype, one date only | 2,695 | fails the rule of 2 |
| Excluded — hearing loss, never bilateral sensorineural | 4,900 | wrong type or wrong laterality |
| Excluded — other ear evidence only | 9,411 | tinnitus, otitis, and so on |

By ancestry: EUR 2,547 / 37,603 · AFR 490 / 10,049 · EAS 43 / 1,002 · AMR 37 / 769 ·
SAS 32 / 800 · unclassified 15 / 532.

---

## Why this phenotype

From the meeting of 2026-10-02:

> "we typically exclude one-sided hearing loss, and we just focus on the bilateral sensorineural
> hearing loss […] unilateral hearing loss we exclude, because it's less likely genetic, more
> likely environment related" — **D. Epstein**

The broad phecode `SO_396` bundles conductive, mixed, unilateral and unspecified hearing loss
together. Restricting to bilateral sensorineural **halves the case count, 6,752 → 3,164**, and that
cost is the main thing to weigh when reading the results.

### Six ICD codes define a case

phecodeX puts two independent axes in the children of `SO_396` — **type** (`.1` conductive,
`.2` sensorineural, `.3` mixed) and **laterality** (`.8` bilateral, `.9` unilateral) — and one ICD
code carries one value from each. Requiring `.2` **and** `.8` together resolves to:

| code | description |
|---|---|
| `H90.3` | Sensorineural hearing loss, bilateral |
| `389.18` | Sensorineural hearing loss, bilateral *(ICD-9)* |
| `389.16` | Sensorineural hearing loss, asymmetrical |
| `H91.13` | **Presbycusis, bilateral** |
| `389.12` | Neural hearing loss, bilateral |
| `389.11` | Sensory hearing loss, bilateral |

`H91.13 Presbycusis` qualifying matters — age-related hearing loss is the phenotype this work is
ultimately about.

**What no longer qualifies**, and would have under the broad definition: `H91.90 Unspecified
hearing loss, unspecified ear` (2,698 people), `H90.5 Unspecified sensorineural` (1,011),
`H91.93 Unspecified hearing loss, bilateral` (1,464). Each is missing one of the two axes.

---

## Two implementation points worth knowing

**The two axes must meet on the same date.** `.2` and `.8` are separate axes, so counting dates that
carry *either* lets a conductive-bilateral visit and a sensorineural-unilateral visit combine into a
"case" never diagnosed with bilateral sensorineural loss. Requiring both on the same date, twice,
gives 3,164; the loose reading gives 4,029. **865 people** separate them.

**Tinnitus has to be read from the OMOP `observation` table.** PMBB v4 relocated the standard codes
(`388.3x`, `H93.1x`) out of the phecode files; reading only those leaves **556 people** carrying an
exclusion criterion sitting among the controls.

---

## Known limitations

- **AFR has 490 cases.** Run for completeness, but it should be read as a null with very little
  power rather than as a search.
- **Audiograms would replace this phase entirely.** Diagnosis codes are a proxy; `H91.90
  Unspecified hearing loss` alone covers 2,698 people who have hearing loss the record does not
  characterise. This is the measurable ceiling on phenotyping by code.
- **Sudden idiopathic loss (`SO_396.5`) is excluded.** 203 people in the broad definition carry it.
  Whether it belongs in a genetic study of age-related loss is a clinical question, not a
  statistical one.

---

*Detail: [`phase_1/results/FINDINGS.md`](../phase_1/results/FINDINGS.md) ·
code: [`phase_1/scripts/`](../phase_1/scripts/) ·
background on PMBB code systems: [`andre_notes/00_pmbb_codes.md`](../andre_notes/00_pmbb_codes.md)*
