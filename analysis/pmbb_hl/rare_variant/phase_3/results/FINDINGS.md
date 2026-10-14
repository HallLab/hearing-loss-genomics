# Phase 3 — covariates and the GRM

**Run 2026-10-14.** Covariates from the release; GRM job `49676669`, 2/2.

## Covariates

PCs come from the **release's own exome PCA**, at the count reaching 80% of the
variance in each cohort's eigenvalues.

| | PCs | variance | | broad analysable | tested | bilateral analysable | tested |
|---|---:|---:|---|---:|---:|---:|---:|
| combined | **11** | 80.3% | | 57,507 | 57,498 | 53,919 | 53,910 |
| EUR | **15** | 83.3% | | 42,786 | 42,779 | 40,150 | 40,143 |
| AFR | **16** | 84.1% | | 11,334 | 11,334 | 10,539 | 10,539 |

Cases: broad 6,752 / 5,183 / 1,285 · bilateral 3,164 / 2,547 / 490.

**Nobody is dropped for missing a PC.** That was the premise of using the release's
PCA, and the script refuses to write anything if it fails. The only people dropped are
9 (combined) and 7 (EUR) with no `AGE` value, all of them controls, so no case is lost
in any cell.

### This is the first run that meets the 80% rule

| | PCs | variance reached |
|---|---|---:|
| original run | 5 / 9 / 10 | 64.2% / 60.8% / 58.6% |
| first pass | 5 / 4 / 3 | 64.2% / 66.7% / 56.5% |
| **this rebuild** | **11 / 15 / 16** | **80.3% / 83.3% / 84.1%** |

The original's figures are computed on the same release eigenvalues it drew its PCs
from; the first pass's on its own PCA. Neither met the convention.

**A caveat carried rather than hidden:** "80% of the variance" means 80% of what the top
twenty eigenvalues capture, because plink publishes twenty. It is not 80% of total
genotypic variance, which twenty components come nowhere near. The rule is still the
rule — it is a proportion of a stated denominator rather than of everything.

### Why the release's PCA and not one of our own

The first pass re-ran a within-ancestry PCA, on the stated grounds that the original
pipeline's PCA had a coverage gap. That was wrong, checkably:

- the original computed **no PCA at all**. Its delivered PCs are identical to the
  release's `exome_PC*` columns at `max|difference| = 0.00e+00` across all 24 components
  in the three cohorts
- its `step7_ldprune.bsub` builds the **GRM**; there is no `--pca` in any of its scripts
  and no eigenvector file in its output
- the release's by-ancestry exome PCA covers **100%** of the analysable cohort

A GRM-building step was read as a PCA step, and a PCA was recomputed to close a gap that
did not exist there. That recomputation was the only real methodological divergence from
the original run, and it is removed here.

## The GRM

Rebuilt for EUR and AFR, because that gap **is** real:

| | rebuild | original | first pass |
|---|---:|---:|---:|
| EUR markers | 53,920 | 53,913 | 53,920 |
| EUR samples | **51,867** | 51,632 | 51,867 |
| AFR markers | 114,694 | 114,696 | 114,694 |
| AFR samples | **14,927** | 14,878 | 14,927 |

The original's pruned genotypes were built over a keep list that had already been cut.
SAIGE intersects its phenotype file with the GRM, so the 196 EUR and 34 AFR cohort
members absent from them would be dropped from the analysis with no message.

Reproducing the marker counts to within 7 and 2 confirms the parameters are the
original's; the difference is the 235 and 49 extra samples, which shift allele
frequencies enough to move a handful of variants across the filters.

Coverage checked: **0 of the tested cohort is outside the new GRM**, in either ancestry,
for either phenotype.

**Combined is not rebuilt.** The original's `ALL.common_ldpruned` covers all 70,925 and
has no gap. Verified: 0 of the 57,498 broad and 0 of the 53,910 bilateral combined
cohort are missing from it.

### One problem carried forward, not fixed

The combined GRM holds **38,833 markers for 70,925 people**, against EUR's 53,920 for
51,867. It is thinner than a subset of itself.

The cause is `--hwe 1e-6` applied **across ancestries**. Hardy-Weinberg departure in a
structured sample is a consequence of mixing populations, not evidence of genotyping
error, so the filter is removing population structure rather than bad genotypes.

Not fixed here, deliberately. Rebuilding it without the filter would be a third
difference from the original run and would change the null model of every combined
analysis. A thin GRM makes the relatedness correction somewhat worse, which is the
conservative direction, and nothing in any combined analysis has come close to a bar.

**If a combined result ever approaches significance, this is the first thing to redo.**

## Files

| | |
|---|---|
| `results/covariates/{broad,bilateral}_{combined,EUR,AFR}.txt` | IID PHENO AGE AGE2 SEX Batch PC1..PCn |
| `results/01_covariates.json` | counts, PC source and variance per cell |
| `data/grm/{EUR,AFR}.common_ldpruned.*` | rebuilt GRM genotypes |
| `data/{EUR,AFR}.keep` | the sample lists they were built over |
| combined GRM | `analysis/elena/rarevariantExWAS/exome_ldpruned/ALL.common_ldpruned`, reused |
