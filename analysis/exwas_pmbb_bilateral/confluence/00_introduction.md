# Rare-variant gene burden in bilateral sensorineural hearing loss

**PMBB Release 2026-4.0 · SAIGE-GENE+ · October 2026**

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

## The cohort in one table

| | N |
|---|---:|
| people with exome data | 70,925 |
| **cases** — bilateral sensorineural, ≥ 2 dates | **3,164** |
| **controls** — no ear evidence of any kind | **50,755** |

Analysed as three cohorts: **combined** (53,910 · 3,164 cases), **EUR** (40,143 · 2,547) and
**AFR** (10,539 · 490).

---

## The pipeline

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
