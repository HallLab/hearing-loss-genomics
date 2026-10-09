# Phase 4 — The association test

**Bilateral sensorineural hearing loss · PMBB Release 2026-4.0**
[← Phase 3 — covariates](03_covariates.md) · Part 4 of 5 · [Phase 5 — results →](05_results.md)

---

## In one line

SAIGE-GENE+ gene burden, **66 jobs** — one per cohort and chromosome, each covering the whole grid
of 3 masks × 3 max-MAF cutoffs and returning **one omnibus p-value per gene**.

---

## The question, finally

> Gene by gene: do people with bilateral sensorineural hearing loss carry more rare damaging
> variants in that gene than people without?

A specific rare variant may appear in 3 people out of 53,910, which supports no conclusion. So
variants are **summed per gene** — the *burden* — turning "is this variant associated?" into "is
*any* damaging variant in this gene associated?".

---

## Two steps, and why

**Step 1 — the null model.** Builds the expectation of who has hearing loss without looking at any
gene. Runs once per cohort.

```
logit P(case) = β₀ + β₁·AGE + β₂·AGE² + β₃·SEX + β₄·Batch + β₅₋ₙ·PCs + b
                                                                       ↑
                                        random effect, b ~ N(0, τ·GRM) — relatedness
```

Relatedness is the expensive part — 1.45 billion pairs for the combined cohort — and it does not
depend on the gene, which is why it is computed once and reused.

What crosses to step 2 is the **residual**, `yᵢ − μ̂ᵢ`: how much more or less likely each person was
to be a case than the model predicted.

**Step 2 — the test.** One question per gene: *are this gene's rare variants concentrated in the
surprising people?*

```
burden   S = Σᵢ Gᵢ · (yᵢ − μ̂ᵢ)       the gene's load × the person's surprise
SKAT     Q = Σⱼ wⱼ² · Sⱼ²            squared, so direction stops mattering
```

Separating the steps is not only economy. The baseline must be built **without ever seeing the
gene**; fitted together, covariates and gene compete for the same variance and the model can absorb
real signal into the covariates. Separated, every gene is measured against the same yardstick.

### The test is a score test, with a saddlepoint correction

Not a likelihood ratio test. A score test needs only the null model, which is what makes the
two-step structure possible — and it stays valid where the others break: a gene with 5 alleles all
in cases is perfect separation, where LRT and Wald would estimate an effect going to infinity.

The **SPA** (saddlepoint approximation) in `step2_SPAtests.R` fixes the tail. The standard score
test assumes a normal curve, which fails exactly when cases and controls are unbalanced and the
variant is rare — 3,164 against 50,746 here. Without it the p-values would be optimistic, and
precisely for the rarest genes.

---

## One call over the whole grid

Each gene is tested nine times — 3 masks × 3 MAF cutoffs — so the honest p-value has to charge for
nine looks. SAIGE's **Cauchy combination** does exactly that, but only when all nine sit in the same
execution. Both parameters take lists, and the MAF one defaults to our three cutoffs:

```
--annotation_in_groupTest="pLOF,pDM,pLOF:pDM"
--maxMAF_in_groupTest=0.0001,0.001,0.01
```

So **66 jobs** rather than 594, each emitting one `Cauchy` row per gene. What comes out follows a
simple rule:

```
omnibus ≈ best cell × ( 9 ÷ how many of the 9 cells carry the signal )

   signal in 1 cell   →  9.0×   exactly the Bonferroni price for nine looks
   signal in all 9    →  1.0×   no penalty
```

The design was tested before adoption: on chr21, the shared cells agree **100%** between one call
and nine, and SAIGE's native Cauchy matches a hand-computed ACAT for 165 of 208 genes.

### What SAIGE corrects, and what it does not

A common assumption, worth stating plainly:

```
within a gene, sweeping SKAT-O's ρ      SAIGE corrects
within a gene, the nine cells (Cauchy)  SAIGE corrects
across the ~18,000 genes                WE correct
```

Checked rather than assumed: if SKAT-O were a raw minimum its p-value would always equal the smaller
of burden and SKAT, and it is greater or equal in 471,567 of 471,615 tests. And the omnibus p-values
are uniform under the null — 0.048 below 0.05, 0.454 below 0.5 — which they could not be if any
gene-level adjustment had been applied.

The SPA is **not** a multiple-testing correction, despite the similar-sounding name.

---

## What ran

| cohort | samples | variance ratio | τ₂ |
|---|---:|---|---:|
| combined | 53,910 | 0.9847 | 0.155 |
| EUR | 40,143 | 0.9942 | 0.059 |
| AFR | 10,539 | 1.0000 | **0** |

**τ₂ = 0 in AFR** means the polygenic variance component converged to zero: at 490 cases there is no
detectable polygenic effect, and SAIGE reduces to logistic regression on the covariates. Not a
defect — the sample size showing up.

All 66 step-2 jobs completed with no failures. **Results in [Phase 5](05_results.md).**

---

*Detail: [`PREMISES.md`](../PREMISES.md) P8, P10 ·
code: [`phase_4/scripts/`](../phase_4/scripts/) ·
design test: [`phase_4/scripts/tests/README.md`](../phase_4/scripts/tests/README.md)*
