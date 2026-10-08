# Phase 3 — What is adjusted away

**Bilateral sensorineural hearing loss · PMBB Release 2026-4.0**
[← Phase 2 — variants](02_variants.md) · Part 3 of 5 · [Phase 4 — the test →](04_test.md)

---

## In one line

Covariate files for three cohorts, carrying age, sex, sequencing batch and **5 / 4 / 3** ancestry
principal components — the counts the eigenvalue scree supports.

---

## What has to be adjusted away, and why

| | if not adjusted |
|---|---|
| **Age** | cases skew older, so any variant commoner in the elderly looks associated |
| **Ancestry** | groups differ in millions of variants for historical reasons; if one group has more diagnoses, every variant commoner in it looks associated |
| **Sequencing batch** | three batches, run at different times; technical bias becomes signal if cases and controls are unevenly distributed |
| **Relatedness** | siblings share genes *and* tend to share diagnoses, for reasons that are not the gene being tested |

The first three enter as covariates here. Relatedness is handled inside SAIGE in the next phase.

---

## Principal components, and how many

Ancestry has no column to read — what exists is each person's genotypes. PCA compresses that into a
few numbers that capture the largest axes of variation, which in human genomes is ancestry. Nobody
tells it to look for ancestry; it finds it.

**Too few** components leave residual structure, which becomes false signal. **Too many** cost
degrees of freedom and can correlate with the phenotype by chance, introducing bias instead of
removing it.

The criterion is the eigenvalue scree — the elbow where the curve flattens:

```
EUR   eigenvalues  194.2  42.4  29.1  13.4  12.4  11.9
      drop            78%   31%   54%    8%    4%
                                          └── elbow: keep 4

AFR   eigenvalues   94.0   7.3   6.0   5.9   5.5   5.3
      drop            92%   19%    1%
                                   └── elbow: keep 3
```

**Adopted: 5 for combined, 4 for EUR, 3 for AFR.** The elbow is the first component whose drop falls
below 10%; moving that threshold to 8% or 12% changes no answer.

The **combined** cohort uses the release's global exome PCs, because there the structure to adjust
for is *between* ancestries. Within EUR or AFR a global PC1 is busy separating continents and says
little about internal structure, so those use within-ancestry PCA.

---

## Provenance

**The PCA was not rerun for this analysis.** It was computed over everyone with exome data per
ancestry — 51,867 EUR, 14,927 AFR — not over any analysis cohort, so changing the phenotype cannot
invalidate it. Carried over and verified byte-for-byte.

The same LD-pruned genotype sets double as the **GRM** that SAIGE uses to measure relatedness in
Phase 4: 38,833 markers for combined, 53,920 for EUR, 114,694 for AFR.

---

## The covariate files

Columns follow the pipeline convention: `IID PHENO AGE AGE2 SEX Batch PC1..PCn`

| cohort | written | cases | controls | PCs |
|---|---:|---:|---:|---:|
| combined | 53,910 | 3,164 | 50,746 | 5 |
| EUR | 40,143 | 2,547 | 37,596 | 4 |
| AFR | 10,539 | 490 | 10,049 | 3 |

**`AGE2`** is age squared: hearing loss does not progress linearly with age — little happens between
20 and 30, a great deal between 70 and 80 — and the squared term lets the model curve.

**`SEX` and `Batch` are declared categorical.** This matters for batch: the values 1, 2 and 3 are
labels, not quantities, and without declaring it SAIGE assumes batch 2 sits exactly halfway between
1 and 3. In AFR that would interpolate across groups of 7,853, 638 and 2,048 people.

**Nine people drop from combined, seven from EUR**, all for unrecorded age. This exclusion is
correct: age is a covariate the model consumes, so a row without it cannot be fitted.

---

## Known limitations

- **The combined GRM is thin** — 38,833 markers against 114,694 for AFR, because `--hwe 1e-6`
  applied across ancestries removes 176k variants where the same filter within EUR removes 19k.
  Hardy–Weinberg departure in a structured sample is expected and is not evidence of genotyping
  error, so the filter does the wrong job there. **Measured and left alone:** λ for this analysis
  runs 0.944–0.999, slightly conservative, not enough to matter. Rebuilding would mean deciding how
  to handle HWE across ancestries, which is a question for the group.
- **Small ancestry groups have no arm of their own.** EAS, SAS, AMR and the unclassified exist only
  in the combined cohort; the combined PCs are what adjust for the mixture.

---

*Detail: [`phase_3/results/FINDINGS.md`](../phase_3/results/FINDINGS.md) ·
code: [`phase_3/scripts/03_build_covariates.py`](../phase_3/scripts/03_build_covariates.py)*
