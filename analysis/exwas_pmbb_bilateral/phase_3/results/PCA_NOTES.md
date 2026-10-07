# Within-ancestry PCA — what was re-run, and what could not be verified

**Script:** [`../scripts/02_pca_by_ancestry.bsub`](../scripts/02_pca_by_ancestry.bsub) · AFR job 49211936, 2026-10-02

## Why it was re-run

The pipeline's within-ancestry PCA ran *after* the cohort was cut against the imputed `.fam`, so the
participants Phase 1 restores have no components: 203 in EUR, 34 in AFR. Using global PCs instead
would be a methodological regression, and dropping the 237 would repeat the defect being corrected.

Scope is everyone with exome data per ancestry (EUR 51,867, AFR 14,927) rather than the corrected
cohort, so one calculation serves both arms and a later difference stays attributable to cohort
alone. Parameters are the pipeline's own, unchanged.

## AFR result

| | ours | the pipeline's |
|---|---:|---:|
| samples | 14,927 | 14,878 |
| variants after MAF/geno/HWE | 256,407 | — |
| variants after LD pruning | **114,694** | **114,696** |
| corrected-cohort members without a PC | **0** | 34 |

Reproducing the LD-pruned set to within **2 variants of 114,696** confirms the parameters are right
and that the only difference is the 49 extra samples, which shift allele frequencies and pairwise
correlation slightly.

**The elbow agrees: 3 PCs**, by both sets of eigenvalues.

## What does not agree, and cannot be resolved

The two scree *profiles* differ structurally:

```
ours            PC1 93.99   PC2  7.33     PC1 / PC2 = 12.8
the pipeline's  PC1 29.98   PC2 16.87     PC1 / PC2 =  1.8
```

Scaling would change magnitude, not ratio. Two things were checked and ruled out:

- **Not outliers.** Zero samples exceed |z| > 5 on PC1, PC2 or PC3. A dominant PC1 is expected here:
  the PMBB AFR group is largely African-American, where a continuous admixture gradient is real
  structure, not an artifact.
- **Not the normalisation.** Their `Variance` column is eigenvalue over the sum of the first 20, which
  reaches exactly 1.0 at PC20 — so "PC1 explains 12.1%" means 12.1% of what the top 20 capture, not of
  total variance. That changes neither the ratio nor the elbow.

**What could not be determined is what their eigenvalues were computed on.** Their `eigenvec` and
`eigenval` did not survive; only the summary table remains. Their `step7_ldprune.bsub` prunes the
exome, but nothing ties the surviving table to that output.

### What this changes

Earlier write-ups argued the PC count from *their* eigenvalues — that their own scree does not support
the number they used. That argument still holds for the elbow, and it is rhetorically strong because
the data is theirs. But it now rests on a file whose provenance cannot be established.

The conclusion no longer depends on it. **Our own PCA, on the exome, documented and reproducible,
gives 3 for AFR.** Where the two agree the decision is doubled; where they differ we use ours and say
why.
