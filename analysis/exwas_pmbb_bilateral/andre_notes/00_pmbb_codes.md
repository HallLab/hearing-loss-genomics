# PMBB diagnosis codes — ICD, OMOP, phecodes, and which table holds what

**Reference note.** Applies to **PMBB Release 2026-4.0**; the table layout changes between releases
and the traps at the end are mostly about exactly that.

Written because "is `H90.3` an ICD code or a phecode, and what is `SO_396`?" is a question that
costs an afternoon every time it is asked fresh.

---

## The short answer

A single diagnosis exists at **four levels at once**. They are not alternatives — they are the same
event described at four grains, and PMBB gives you a table for each.

| level | looks like | who made it | what it is for |
|---|---|---|---|
| **ICD** | `H90.3`, `389.18` | the clinician, at the visit | billing |
| **OMOP concept** | `SNOMED 60700002` | PMBB, on ingest | one vocabulary across hospitals |
| **Phecode 1.2** | `389.1` | research community, 2010s | the old standard for PheWAS |
| **PhecodeX** | `SO_396.2` | research community, 2023 | the current one |

**`SO_396` is a PhecodeX code, not an ICD code.** The `SO_` is a category prefix, not part of a
number.

---

## One diagnosis, traced through all four

A real person in the release, `PMBB9927890428710`, visit of 2015-12-23.

**1 — ICD, in `condition_occurrence`.** What was actually typed at the visit. Five codes that day:

```
H61.21    Impacted cerumen, right ear
H61.813   Other specified disorders of external ear, bilateral
H90.3     Sensorineural hearing loss, bilateral
H90.5     Unspecified sensorineural hearing loss
H90.8     Mixed conductive and sensorineural hearing loss, unspecified
```

**2 — OMOP concept, in the `condition_concept_id` column.** PMBB maps each ICD code to a standard
concept so that a code from one hospital's system means the same thing as a code from another's:

```
condition_concept_id  374366
  -> SNOMED 60700002   "Sensorineural hearing loss"
```

**3 — Phecode 1.2, in `conditions_phecode_12`.** The grouping used by most published PheWAS up to
about 2023. Looks like an ICD-9 code and is not one:

```
389.1     Hearing loss
```

**4 — PhecodeX, in `conditions_phecode_x`.** The current grouping, and what PMBB v4 leads with:

```
SO_396      Hearing loss
SO_396.2    ... sensorineural
SO_396.8    ... bilateral
```

Note that one ICD code (`H90.3`) produced **three** PhecodeX codes. That is deliberate and is
explained below.

---

## Why phecodes exist at all

ICD is built for billing, not research, and it shows. `H91.90 Unspecified hearing loss, unspecified
ear` and `H90.3 Sensorineural hearing loss, bilateral` are different billing events and, for most
genetic questions, the same disease. ICD-10 has roughly 70,000 codes; two clinicians seeing the same
patient routinely pick different ones.

A **phecode** collapses the ICD codes that mean one condition into one research-usable group. That is
all it is: a many-to-one map, maintained by researchers, from billing codes to disease concepts.

**Phecode 1.2 vs PhecodeX.** Phecode 1.2 was built on ICD-9 and shows it — the codes look like ICD-9
codes, and ICD-10-only concepts fit badly. PhecodeX (Shuey et al. 2023) was rebuilt for ICD-10, is
organised by body system, and splits conditions along clinically meaningful axes that 1.2 bundled
together. **Use PhecodeX for new work on v4**; use 1.2 only to compare against an older paper.

---

## Reading a PhecodeX code

```
SO_396.2
│  │   └── axis within the condition
│  └────── the condition
└───────── body system
```

### The 18 body-system prefixes in PMBB v4

With the number of distinct codes each carries in this release, largest first:

```
SO 311    CM 293    MS 279    GE 256    CA 244    NS 236
GI 206    CV 194    GU 177    EM 173    ID 150    RE 147
PP 129    DE 128    NB 104    BI 102    SS  78    MB  50
```

(Counted over the whole table. A sample of the first few million rows undercounts every one of
them by 5–50%, because the rare codes are not evenly spread through the file — worth knowing before
quoting a number you got in a hurry.)

**`SO` is Sense Organs** — eye and ear together. `SO_375.11` is a lacrimal gland disorder,
`SO_396` is hearing loss. The full category key is in the PhecodeX publication; the prefixes
above are what this release actually contains.

### The children are axes, not subtypes — the part that catches people

This is the single most useful thing in this note, and it is not obvious from the codes.

The children of a PhecodeX code are often **several independent axes at once**, and one ICD code
carries one value from each. For hearing loss:

| axis | codes |
|---|---|
| **type** | `.1` conductive · `.2` sensorineural · `.3` mixed · `.5` sudden idiopathic |
| **laterality** | `.8` bilateral · `.9` unilateral |

So `H90.3 Sensorineural hearing loss, bilateral` maps to `SO_396` **and** `SO_396.2` **and**
`SO_396.8` — the parent, the type, and the laterality.

**What this means in practice.** "Bilateral sensorineural hearing loss" is not one child code. It is
the intersection `.2 ∧ .8`, and the two have to co-occur **on the same date** to mean one diagnosis
rather than two visits that each had half of it. Treating the children as mutually exclusive
subtypes is the standard way to get this wrong, and it does not throw an error — it quietly builds
a different cohort.

### The parent is pre-rolled

A person with `SO_396.2` always carries `SO_396` as well, so you do not need to roll children up
into their parent yourself. Verified in v4: of 8,922 people with any `SO_396.x` child, **zero** lack
the parent.

---

## Which table to read

All under `/static/PMBB/PMBB-Release-2026-4.0/Phenotype/4.0/`, prefixed
`PMBB-Release-2026-4.0_phenotype_`.

| file | `condition_source_value` holds | use it for |
|---|---|---|
| `condition_occurrence` | raw ICD-9-CM **and** ICD-10-CM | the literal code at the visit |
| `conditions_icd10` | ICD-10-CM only | ICD-10-only work |
| `conditions_phecode_12` | Phecode 1.2 (`389.1`) | comparing against older papers |
| `conditions_phecode_x` | **PhecodeX** (`SO_396`) | **new phenotyping** |
| `conditions_primary_dx` | the visit's primary diagnosis | severity proxies |
| `observation` | ICD **and** CPT, in `observation_source_value` | see the trap below |

Every one of them has `person_id` and a date, which is what the rule of 2 needs.

---

## Traps

### 1 — codes move between tables across releases

**This one has already cost this group a cohort.** PMBB v4 moved the standard tinnitus codes
(`388.3x`, `H93.1x`) out of the phecode files into **`observation`**. The phecode files kept only
the ~3,016 pulsatile-tinnitus events; the other **25,094** are in `observation`.

Anyone building an ear-family phenotype from `conditions_phecode_x` alone silently loses **556
people** who then sit in the control group carrying an exclusion criterion.

**Rule: when a phenotype count looks low, check `observation` before checking your code.**

### 2 — `conditions_phecode_12` and `conditions_phecode_x` are different groupings

Not two spellings of the same thing. They disagree on what belongs together, and PhecodeX splits
conditions that 1.2 bundled. Never mix them in one phenotype definition, and say which you used.

### 3 — the parent bundles everything

`SO_396` is *all* hearing loss: conductive, mixed, unilateral, unspecified. A cohort built on the
parent is a different phenotype from one built on a specific child, and the difference is large —
in this project, 6,752 cases against 3,164.

### 4 — "unspecified" codes are common and they are not harmless

`H91.90 Unspecified hearing loss, unspecified ear` covers 2,698 people in v4. They have hearing
loss; the record does not say what kind. Any definition requiring a specific type or side drops all
of them, and some of them certainly have the thing you are looking for. This is the measurable
ceiling on phenotyping by diagnosis code, and the reason audiograms are worth the effort.

### 5 — a diagnosis on one date means little

The rule of 2 — a code on **two distinct dates** — is standard, and it is doing real work. A single
code can be a suspicion later ruled out, a typo, or a label attached to justify ordering a test.
It also costs: in this project 2,695 people have bilateral sensorineural hearing loss recorded
exactly once, and are excluded.

---

## Recipes

**All diagnoses in one PhecodeX family, with dates:**

```bash
awk -F'\t' '
  NR==1 { for (i=1;i<=NF;i++) h[$i]=i; next }
  $h["condition_source_value"] ~ /^SO_39[0-9]($|\.)/ {
    print $h["person_id"] "\t" $h["condition_start_date"] "\t" $h["condition_source_value"]
  }' PMBB-Release-2026-4.0_phenotype_conditions_phecode_x.txt
```

The `($|\.)` matters: without it, `SO_39` also matches `SO_390` through `SO_399` *and* anything
longer that happens to start the same way.

**Rule of 2, where two axes must meet on the same date:**

```python
pair = ear[ear.code.isin(["SO_396.2", "SO_396.8"])]
per_date = pair.groupby(["pid", "date"]).code.nunique()
qualifying = per_date[per_date == 2].reset_index().groupby("pid").date.nunique()
cases = set(qualifying[qualifying >= 2].index)
```

Counting dates that carry `.2` **or** `.8` instead gives a different and wrong answer. A real person
in v4, `PMBB2551063181373`, shows why:

```
2008-07-07   .1 .8      conductive, bilateral
2010-03-13   .1 .8      conductive, bilateral
2010-12-17   .2 .3 .9   sensorineural + mixed, UNILATERAL
2015-11-15   .1 .8      conductive, bilateral
2019-04-26   .1 .8      conductive, bilateral
2019-12-06   .1 .9      conductive, unilateral
2020-01-24   .3 .9      mixed, unilateral
2023-09-04   .8         bilateral, type not recorded
2024-11-02   .1 .8      conductive, bilateral
```

No visit carries `.2` and `.8` together, so this person was never diagnosed with bilateral
sensorineural loss. They have conductive bilateral loss six times and one unilateral sensorineural
visit. Counting dates that carry either axis gives **eight** and makes them a case; requiring both
on one date gives **zero**.

The axes only describe one diagnosis when they appear on the **same line of the record**. Taken
apart, they are two different diseases added together by accident. In this project the loose
reading admitted 865 such people — 27% more cases, every one constructed this way, and nothing
downstream would have objected.

**Which ICD codes sit behind a phecode:** no lookup table ships with the release. Join
`condition_occurrence` to `conditions_phecode_x` on `person_id` + date + `condition_occurrence_id`,
or ask whoever on the team already built the map for your phenotype.

---

## Sources

- **PhecodeX** — Shuey MM et al., *Bioinformatics* 2023. The category key and the axis design.
- **Phecode 1.2** — phewascatalog.org. Note the site 403s to scripted requests; the map also ships
  inside the PheWAS R package as `pheinfo`.
- **OMOP CDM** — ohdsi.org. What `condition_occurrence`, `observation` and `concept` are, and why
  they are shaped that way.
- **PMBB release notes** — the per-release documentation is the only authority on which table holds
  what *in that release*. Trap 1 exists because that changed.
