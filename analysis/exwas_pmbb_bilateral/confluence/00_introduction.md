# Rare-variant gene burden in bilateral sensorineural hearing loss

**PMBB Release 2026-4.0 · SAIGE-GENE+ · October 2026**
Hall Lab × Epstein Lab · Penn

---

## The question

Do people with bilateral sensorineural hearing loss carry more rare, damaging variants in any
particular gene than people without?

## The answer

**No gene reaches exome-wide significance**, in any of the three ancestry cohorts, under any
defensible correction. The test is well calibrated, so this is a real absence rather than an
inconclusive one, and the limiting factor is the size of the cohort rather than any analysis choice.

The rankings underneath that null do carry known deafness biology — which says the pipeline works
and the signal sits below the detection floor. **[Full reasoning in the conclusions →](06_conclusions.md)**

---

## Why this study exists

Adult-onset hearing loss is one of the most common sensory disorders, and its genetic contribution
is poorly characterised. Many genes are known to cause **congenital** deafness; whether rare variants
in those same genes also raise the risk of hearing loss acquired later in life is a different
question, and a biobank with both exome sequencing and decades of clinical records is one of the few
places it can be asked.

**Hui et al. 2023** (*PLOS Genetics*, [doi:10.1371/journal.pgen.1010584](https://doi.org/10.1371/journal.pgen.1010584))
did exactly that in PMBB and answered yes: rare deleterious variants in known congenital hearing-loss
genes do raise adult-onset risk, with four novel candidates emerging alongside. That study used
**40,627 people** and phenotyped on phecode 389 — hearing loss of any kind.

This work extends that line into **PMBB Release 2026-4.0**, which has **70,925 people with exome
data**, and sharpens the phenotype.

---

## Where this analysis sits

```
  Hui et al. 2023          PMBB ~40k · phecode 389 · burden in known HL genes
        │                  the precedent this line extends
        ▼
  ExWAS on PMBB v4         70,925 with exomes · single-variant + rare-variant burden
        │                  summer 2026 · results presented as a poster
        ▼
  Replication              the same pipeline re-run and audited
        │                  found 4 mask defects and 2 cohort defects
        ▼
  THIS ANALYSIS            the corrected pipeline + the clinical phenotype
        │                  bilateral sensorineural only
        ▼
  next                     audiograms · a larger cohort · external replication
```

Each step is a response to the one above it. The replication exists because results should be
checked before they are built on; **this analysis exists because the replication showed the pipeline
was sound once corrected, and the remaining question was whether the phenotype was right.**

`analysis/elena_replication/` holds the audit. The [conclusions](06_conclusions.md) set this
analysis against it directly.

---

## Why bilateral sensorineural

Hui et al. and the PMBB v4 ExWAS both used the broad definition — any hearing loss. Narrowing to
bilateral sensorineural **halves the case count, 6,752 → 3,164**, which is a real cost. Four reasons
it is still the right phenotype for a genetic study:

**1 · Unilateral loss is usually environmental.** Hearing lost in one ear points to something that
happened to that ear — trauma, infection, noise exposure on one side, a tumour. A genetic cause
affects both cochleae. Including unilateral cases adds people whose deafness has a known
non-genetic explanation.

> "unilateral hearing loss we exclude, because it's less likely genetic, more likely environment
> related" — **D. Epstein**, 2026-10-02

**2 · Conductive loss is a different organ.** Conductive loss is mechanical — the ossicles, the
middle ear, the eardrum — while sensorineural loss is the cochlea and the auditory nerve. They share
a symptom and almost nothing else. Of the **109 Definitive/Strong gene–disease pairs** in the ClinGen
hearing-loss panel, **none is labelled as conductive hearing loss**.

**3 · The target phenotype is bilateral sensorineural.** Presbycusis — age-related hearing loss, the
condition this line of work is ultimately about — is bilateral and sensorineural by definition. Its
ICD code `H91.13` qualifies under this restriction.

**4 · It is how the clinical side defines the condition.** A phenotype a clinician would not
recognise is hard to act on, whatever it does for statistical power.

### And the cost was measured, not assumed

Halving the cases could have thrown away the signal. Two controls were run to find out — the broad
phenotype drawn down to the same case count five times, and the discarded cases run as their own
case set. **The restriction's null is lost power, not lost biology.** Details in the
[conclusions](06_conclusions.md).

---

## The cohort in one table

| | N |
|---|---:|
| people with exome data | 70,925 |
| **cases** — bilateral sensorineural, ≥ 2 dates | **3,164** |
| **controls** — no ear evidence of any kind | **50,755** |

Analysed as three cohorts: **combined** (53,910 · 3,164 cases), **EUR** (40,143 · 2,547) and
**AFR** (10,539 · 490).

---

## The analysis pipeline

The five phases this publication documents, inside the box marked *THIS ANALYSIS* above.

```
   PMBB v4                PHASE 1              PHASE 2            PHASE 3
   release              who is in the      which variants      what is adjusted
                          study               count                away
┌──────────────┐      ┌─────────────┐     ┌────────────┐     ┌──────────────┐
│ phecode_x    │─────▶│   cases     │     │    VEP     │     │   ancestry   │
│ observation  │      │   3,164     │     │   masks    │     │     PCA      │
│ exome VEP    │──────────────────────────▶│  3 masks  │     │   5 / 4 / 3  │
│ exome PLINK  │      │  controls   │     └─────┬──────┘     └───────┬──────┘
└──────────────┘      │   50,755    │           │                    │
                      └──────┬──────┘           │            ┌───────▼──────┐
                             │                  │            │  covariates  │
                             └──────────────────┼───────────▶│  per cohort  │
                                                │            └───────┬──────┘
                                         PHASE 4│                    │
                                      ┌─────────▼────────────────────▼───┐
                                      │  SAIGE step 1   3 null models    │
                                      │  SAIGE step 2   66 jobs          │
                                      │    whole grid per call → Cauchy  │
                                      └────────────────┬─────────────────┘
                                              PHASE 5  │
                                      ┌────────────────▼─────────────────┐
                                      │  merge · Manhattan · QQ · tables │
                                      │  ClinGen enrichment · 2 controls │
                                      └──────────────────────────────────┘
```

| | what it decides | built or inherited |
|---|---|---|
| **[1 · Cohort](01_cohort.md)** | who is a case, a control, or neither | built here |
| **[2 · Variants](02_variants.md)** | which rare variants enter each gene's test | inherited — phenotype-independent |
| **[3 · Covariates](03_covariates.md)** | age, sex, batch, ancestry | PCA inherited, covariates built here |
| **[4 · The test](04_test.md)** | SAIGE-GENE+ burden, one omnibus p per gene | built here |
| **[5 · Results](05_results.md)** | what came out, and whether to believe it | built here |
| **[Conclusions](06_conclusions.md)** | the verdict, the three-arm comparison, next steps | — |

Phases 2 and 3 depend on the genome annotation and on who was sequenced, not on who was diagnosed,
so they were carried over from the preceding replication and verified byte-for-byte rather than
recomputed. Everything touching the phenotype was rebuilt.

---

## The decisions this analysis makes

Each is argued in full in [`PREMISES.md`](../PREMISES.md), with its source and the number it moves.

| | decision | why |
|---|---|---|
| **P1** | cases are **bilateral sensorineural only** | clinical: unilateral loss is less likely genetic — 6,752 → 3,164 cases |
| **P2** | controls unchanged | narrowing cases must not widen controls |
| **P3** | cohort framed on **exome** samples | the analysis runs on exome |
| **P4** | tinnitus read from OMOP `observation` | PMBB v4 moved the codes there |
| **P5** | **rebuilt** variant masks | SpliceAI gate, REVEL as a list, protein-coding only |
| **P6** | the **`ALL` mask is not run** | pools variants with no mechanism |
| **P7** | **5 / 4 / 3** principal components | what the eigenvalue scree supports |
| **P8** | `Batch` declared **categorical** | it is a label, not a quantity |
| **P9** | Bonferroni **and** FDR, both levels | the denominator has to be named, not implied |
| **P10** | **one SAIGE call per cohort and chromosome** | yields a per-gene omnibus from the tool itself |

---

## How to read the results

Two things make the numbers here different from a conventional burden table, and both are worth
knowing before looking at one.

**The headline p-value is a per-gene omnibus, not a best cell.** Each gene is tested nine times —
three masks × three frequency cutoffs — and taking the smallest of nine correlated tests against a
bar built for one test is undercharging. SAIGE's `Cauchy` row combines all nine with the search
already priced in, so the bar is simply `0.05 / genes`.

**The `search_penalty` column says how concentrated a gene's signal is.** It is the omnibus divided
by the best cell: 1× means the effect shows in all nine, 9× means it lives in exactly one. A gene
with the same p-value in all nine cells is more convincing than one that needed finding.

---

## Related work

| | |
|---|---|
| `analysis/elena_replication/` | the audit this analysis descends from — what the earlier pipeline did, and the four mask and two cohort defects found in it |
| `analysis/elena/HL_only_ExWAS/` | single-variant association on the same cohort, **common** variants (MAF ≥ 0.01) — complementary, and testing a disjoint variant set |

The three-arm comparison in the [conclusions](06_conclusions.md) sets this analysis against the
first of those.

---

*Repository: `analysis/exwas_pmbb_bilateral/` · every figure and number regenerable from the
scripts in each phase · what was copied rather than computed is listed in
[`PROVENANCE.md`](../PROVENANCE.md)*
