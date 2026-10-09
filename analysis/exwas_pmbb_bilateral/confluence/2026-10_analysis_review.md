# Hearing loss ExWAS — where we are

**Review for the group · October 2026 · Hall Lab × Epstein Lab**
Covers two rounds of work on the PMBB v4 exome analysis, set against the original run.

---

## In one slide

```
  original run (summer 2026)   ──▶   replication   ──▶   bilateral analysis
  PMBB v4, broad phenotype           same pipeline       corrected pipeline
                                     audited             + clinical phenotype
```

| | cases | significant genes |
|---|---:|---:|
| **Original** | 6,712 | 0 |
| **Replication** | 6,752 | 0 |
| **Bilateral** | 3,164 | 0 |

**No round finds a significant gene.** The conclusion is stable across every phenotype and method
combination we have tried, and the limiting factor is the cohort rather than any analysis choice.

**But the rankings underneath are not noise** — known deafness genes cluster at the top of two of
the three rounds more than chance allows. That is the one positive result, and it says the pipeline
works and the signal sits below the detection floor.

---

## Round 1 — the replication, and what it found

The original pipeline was re-run end to end against the institutional release. The case/control
rule reproduced **exactly**, 70,925 of 70,925. Six things did not.

### Cohort

| | effect |
|---|---|
| Cohort cut against the **imputed** `.fam` | 431 analysable people absent, **40 of them cases** — they have exome data but no imputed PCs, which the model never uses |
| Tinnitus read only from the phecode files | **556 people** the exclusion rule covers ran as controls; PMBB v4 moved those codes to the OMOP `observation` table |

### Variant masks

| | effect |
|---|---|
| SpliceAI gate calculated and never applied | **63.8%** of pLOF entries failed the pipeline's own documented criterion — 694,802 of 1,089,876 |
| Three IMPACT=LOW splice terms listed as pLOF | admitted variants near the splice site, not in it |
| REVEL parsed as one number, not a per-transcript list | **73.9%** of damaging missense variants lost — 5,391 kept where 20,618 qualify *(found by N. Palmiero)* |
| Non-coding genes in the masks | **1,101 of 19,038** — a loss-of-function test on a gene with no protein *(found by N. Palmiero)* |

Rebuilt, the masks move in **both** directions — three shrink because they held things that did not
belong, and `pDM` grows 23% because the REVEL bug had been hiding three of every four damaging
variants:

```
pLOF       1,002,120  →    462,144     −54%
pDM          720,983  →    890,332     +23%
pLOF_pDM   1,717,883  →  1,340,936     −22%
```

### Covariates and reporting

| | effect |
|---|---|
| PC counts 5 / 9 / 10 | the eigenvalue scree supports **5 / 4 / 3**; surplus components cost power and can introduce bias |
| 11 of 792 jobs produced no output | one task exceeded its memory limit and `errorStrategy = 'terminate'` cancelled ten healthy ones; the merge never noticed |
| Three result files lost their header | 5,391 rows landed in phantom columns, unusable, and one gene's result was consumed as the column names |

### The sharpest single example

**`TMC3-AS1` was the top gene of the original run**, in two of three cohorts. It is an `lncRNA` — a
gene with no protein — so a loss-of-function burden test on it cannot mean anything, and the whole
result rests on two variants. It is absent from both later rounds.

The name is part of why it matters: it is antisense to `TMC3`, and `TMC1` is an established deafness
gene. A reader scanning a top-gene table sees "TMC" and reads plausibility into it.

---

## Round 2 — the bilateral analysis

With the pipeline corrected, the open question became the **phenotype**. The original used phecode
`SO_396`, which bundles conductive, mixed, unilateral and unspecified hearing loss together.

> "we typically exclude one-sided hearing loss, and we just focus on the bilateral sensorineural
> hearing loss" — **D. Epstein**, 2026-10-02

Restricting to bilateral sensorineural takes cases **6,752 → 3,164**.

### Pipeline changes beyond the phenotype

| | |
|---|---|
| `ALL` mask dropped | pools variants with no mechanism — a position Doug reached independently |
| `Batch` declared categorical | it is a label; SAIGE was reading 1/2/3 as a quantity |
| One SAIGE call per cohort and chromosome | **66 jobs instead of 594**, and the tool returns a per-gene omnibus instead of nine p-values needing combination afterwards |
| Bonferroni **and** FDR, both aggregation levels | three denominators are defensible and differ ninefold, so the one in use has to be named |

---

## Results

Each cell names the best gene in that cohort. The p-value is the **per-gene omnibus** — one number
per gene, with the nine-way search across masks and frequency cutoffs already charged for.

| | combined | EUR | AFR |
|---|---|---|---|
| **Replication** | `PSG8` · 3.02 × 10⁻⁵ | `ATP13A4` · 6.91 × 10⁻⁵ | `CFH` · 2.17 × 10⁻⁵ |
| **Bilateral** | `PCBD1` · 2.05 × 10⁻⁵ | `FBXO25` · 1.37 × 10⁻⁵ | `IGSF9` · 9.14 × 10⁻⁶ |
| **significance threshold** | 2.79 × 10⁻⁶ | 2.79 × 10⁻⁶ | 2.81 × 10⁻⁶ |
| how far short | 7× to 11× | 5× to 25× | 3× to 8× |

**What the threshold is.** We test about 17,900 genes per cohort, so some will look impressive by
chance alone. Bonferroni divides the usual 0.05 by the number of tests —
`0.05 / 17,941 = 2.79 × 10⁻⁶` — which holds the probability of *any* false positive across the whole
scan at 5%. A gene has to beat that number, not 0.05, to count.

Nothing does, in either round, and Benjamini-Hochberg finds nothing either: the smallest q-values
are 0.15 to 0.22.

**None of the six genes above is a candidate.** They are the closest things to a signal, and each is
fragile in a way the p-value alone hides:

- `FBXO25` rests on **18 alleles**, 7 of them in cases, and its signal sits in exactly one of the
  nine mask-by-frequency cells — the maximum concentration possible.
- `CFH` tops AFR on **five alleles, all in cases.** Remove one carrier and it is gone.
- `PSG8` has a p-value of **0.40** on its burden component: the signal is all dispersion, with no
  consistent direction behind it.

λ runs 0.68 to 1.00 across both rounds — the test is calibrated and slightly conservative, so the
absence is informative rather than inconclusive.

### The one positive result

Tested against the **ClinGen Hearing Loss panel** — 100 Definitive/Strong genes, about 93 in each
round's tested set. If a ranking carries biology, those genes sit near its top more than chance
allows.

| | ClinGen genes in top 50 | expected | p |
|---|---:|---:|---|
| Original | 2 | 0.26 | 0.025 |
| **Replication** | **5** | 0.26 | **5.9 × 10⁻⁶** |
| Bilateral | 0 | 0.26 | 1 |

`SIX1`, `GJB3`, `COCH`, `MYO6`, `TMPRSS3` in the replication's top 50, where chance predicts a
quarter of one gene.

**What this does not license.** The jump from 2 to 5 and the fall to 0 both look like stories. They
are not supported: no pair of these counts is statistically distinguishable from another
(p = 0.436, 0.056, 0.678 …), and that test is already optimistic because the rounds share people and
genes. The honest claim is narrower — **two rounds carry known biology in their ranking, and no
round is shown to beat another.**

### Did the phenotype restriction cost us signal?

Tested, not argued. The broad phenotype was drawn down to 3,164 cases five times and rerun, and the
discarded cases were run as their own case set:

| | cases | ClinGen in top 50 |
|---|---:|---:|
| broad | 6,752 | 5 |
| broad at 3,164, five draws | 3,164 | 2 · 0 · 1 · 2 · 3 |
| bilateral sensorineural | 3,164 | **0** |
| the cases it discards | 3,588 | 1 |

The five draws are the yardstick: at 3,164 cases, 0 to 3 is what the broad definition finds. The
restricted arm's 0 sits inside that range, and the halves do not add up — 0 + 1 against 5.

> **Lost power, not lost signal.** The phenotype decision is not shown to discard biology — and also
> not shown to be better, only not worse for the reason suspected.

---

## Next actions

### Decisions we need from the group

| | question | who |
|---|---|---|
| **1** | Does **sudden idiopathic** hearing loss (`SO_396.5`, 203 people) belong in a genetic study of age-related loss? | clinical — Doug |
| **2** | Is a **null result** worth writing up as it stands, or do we hold for a larger cohort? | all |
| **3** | How should **Hardy–Weinberg** be handled across ancestries? It currently thins the combined GRM to 38,833 markers against 114,694 for AFR | methods — Molly / Nikki |

### Work, in order of what it would change

**1 · Audiograms.** The only item that changes the *question* rather than the answer. Diagnosis
codes are a proxy: `H91.90 Unspecified hearing loss, unspecified ear` alone covers **2,698 people**
whose loss the record does not characterise, and no phenotyping rule recovers what was never
written down. Blocked on a REDCap-to-PMBB ID bridge. **When it exists, only Phase 1 changes** —
everything downstream runs unaltered.

**2 · A larger cohort.** All of Us or UK Biobank. This is what a null at 3,164 cases asks for, and
it is the third sense of "replication" — reproducing a finding in independent data — that neither
round has addressed.

**3 · The X chromosome.** Currently uncovered, costing **6 of the 100** ClinGen Definitive/Strong
genes: `AIFM1`, `COL4A5`, `POU3F4`, `PRPS1`, `SMPX`, `TIMM8A`. Coverage is 94%. Needs X masks and
SAIGE handling hemizygosity in males.

**Deliberately not next:** sensitivity analyses on REVEL 0.5 vs 0.6 or SpliceAI 0.2 vs 0.5. They
would move the ranking but not the conclusion, and sensitivity analysis on a null result produces
more null results. They become worth doing if the cohort grows.

---

## Where everything lives

| | |
|---|---|
| `analysis/elena_replication/` | the audit — per-phase findings, both arms, all evidence |
| `analysis/exwas_pmbb_bilateral/` | the bilateral analysis, self-contained |
| `analysis/exwas_pmbb_bilateral/confluence/` | the seven-page write-up this review summarises |
| `…/phase_5/notebooks/01_three_way.ipynb` | the three-round comparison, reproducible |

Every figure and number in this document regenerates from the scripts in those folders. What was
copied rather than recomputed is listed in each folder's `PROVENANCE.md`.
