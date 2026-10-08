> **INHERITED DOCUMENT — this is `analysis/elena_replication/phase_3/results/FINDINGS.md`,**
> **copied here unchanged.** It describes the replication's audit, which has two arms and is
> written against the pipeline it was auditing. It is **not** a findings file for this analysis.
>
> It is kept because phase 3's outputs were produced by that work and carried over rather than
> recomputed, so this is the evidence behind them — see [`../../PROVENANCE.md`](../../PROVENANCE.md).
> Where it says "corrected arm", that is the arm whose outputs this folder inherited.
>
> This analysis's own account of phase 3 is in
> [`../../andre_notes/phase_03_analise_pt.md`](../../andre_notes/phase_03_analise_pt.md).

---

# Phase 3 — what is adjusted away

**Status:** complete · **Date:** 2026-10-05
**Scripts:** `01_pc_selection.ipynb` · `02_pca_by_ancestry.bsub` · `03_build_covariates.py`

Phase 3 covers the pipeline's Step 6 (covariate files) and Step 7 (LD pruning).

---

## Finding 1 — PC counts exceed what the pipeline's own scree supports

| cohort | ran with | scree supports | PC source |
|---|---:|---:|---|
| combined | 5 | **5** ✓ | the release's exome PCA |
| EUR | 9 | **4** | within-ancestry PCA |
| AFR | 10 | **3** | within-ancestry PCA |

**The combined cohort was right.** EUR and AFR used more than double what the eigenvalues support.

Elbow taken as the first PC whose drop from the previous falls below 10%. No case is borderline;
moving the threshold to 8% or 12% changes no answer.

### Confirmed twice, independently

The decision was first taken from the pipeline's own eigenvalue files. It has since been re-derived
from a PCA we ran ourselves:

```
EUR   ours   78.2%  31.2%  53.9%  ->  8.0%     elbow at PC5, keep 4
      theirs 54.6%  24.5%  12.0%  ->  7.1%     elbow at PC5, keep 4

AFR   ours   92.2%  18.6%  ->  1.4%            elbow at PC4, keep 3
      theirs 43.7%  19.2%  ->  3.8%            elbow at PC4, keep 3
```

Two routes, same answer. The conclusion no longer depends on a file whose provenance could not be
established — see [`PCA_NOTES.md`](PCA_NOTES.md) for the AFR profile discrepancy that remains
unexplained, and which does **not** appear in EUR: there PC1 is 194.24 against their 197.37.

### Why PCs of surplus are not free

The intuition that more components means more protection is wrong twice over. Each covariate spends
degrees of freedom, so components explaining nothing cost power to detect what is being looked for.
And a component from the plateau can correlate with the phenotype by chance — introducing bias rather
than removing it. This is the reasoning the 2026-07-01 meeting applied in rejecting 20 PCs.

---

## Finding 2 — the pipeline's PCA cannot cover the corrected cohort

Its within-ancestry PCA ran **after** the cohort was cut against the imputed `.fam`, so the
participants Phase 1 restores have no components: **203 in EUR, 34 in AFR**.

Three ways out, two of them wrong:

- **Use the release's global PCs instead.** Inside EUR a global PC1 is busy separating continents and
  carries almost nothing about structure within Europe. A methodological regression — and the pipeline
  got this right.
- **Drop the 237.** Excluding people for missing data, which is the defect being corrected, with 203
  of them from the group already excluded once.
- **Re-run the PCA without the pre-filter.** ← done.

Scope was everyone with exome data per ancestry (EUR 51,867, AFR 14,927) rather than the corrected
cohort, so one calculation serves both arms and a later difference stays attributable to cohort alone.
Parameters are the pipeline's own, unchanged.

**AFR reproduces its LD-pruned set to within 2 variants of 114,696**, which confirms the parameters
and leaves the extra samples as the only difference. Coverage of the corrected cohort is now complete:
0 without components in either stratum.

---

## Phase 3 output — covariates for both arms

`results/covariates/`, built by `03_build_covariates.py`.

| stratum | reproduction | corrected |
|---|---|---|
| combined | 57,632 · 5 PCs | 57,498 · 5 PCs |
| EUR | 43,016 · 9 PCs | 42,779 · 4 PCs |
| AFR | 11,387 · 10 PCs | 11,334 · 3 PCs |

Columns follow the pipeline exactly: `IID PHENO AGE AGE2 SEX Batch PC1..PCn`. The reproduction arm is
copied verbatim rather than rebuilt, for the same reason as Phase 1: the phenotype file it derives
from was overwritten.

### Nine people lost to missing age — and why that is different

The corrected arm loses **9** participants for having no `sample_age`: 9 in combined, 7 of them EUR,
0 AFR. They are the whole set of people in the release without an age, and all nine are among the 431
Phase 1 restores — none is in the pipeline's cohort, because they had already been excluded by the
`.fam` filter.

**This is what a correct exclusion looks like, and the contrast is the point.** Age is a covariate the
model actually fits. Dropping someone who lacks it is necessary. The 517 of Phase 1 were dropped for
lacking imputed PCs, which the model never consumes. Same mechanism, opposite justification.

The count is reported rather than silently dropped, which is itself the practice Phase 1 found missing.

---

## What is not established

- **The AFR scree profile discrepancy.** Ours gives PC1/PC2 = 12.8 against their 1.8. Outliers and
  normalisation were ruled out; their raw PCA output did not survive, so what their eigenvalues were
  computed on cannot be determined. It does not appear in EUR, which makes it specific rather than a
  general effect of the sample-set difference.
- **Whether the PCA itself was well constructed.** This phase checks whether PCs were used correctly,
  not whether the underlying PCA was sound. A PCA over a badly filtered variant set would produce
  eigenvalues that look reasonable and would not be caught here.
- **Step 7 beyond its PCA role.** LD pruning was re-run because the PCA needs it. Whether the
  pipeline's own LD-pruned sets are otherwise correct was not audited.
