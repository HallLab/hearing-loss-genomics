# Phase 2 — MAF census

**Run 2026-10-09** · job `49564495`, 66/66 tasks · every histogram sums to its `.bim` count.

Nothing was built in this phase. The genotypes are the files both burden arms already
read and the files Elena's own ExWAS read, so there was nothing to rebuild and no
opportunity for this arm to diverge from the others on genotypes (premise P4).

What this phase produces is **one number per cohort**: how many variants SAIGE will
actually test. That number is the denominator of the secondary significance bar, and
it is fixed here, before step 2 exists.

## The census

13,223,743 variants on file per chromosome set. Of those, at MAF ≥ 0.01:

| cohort | samples | tested | share of file | secondary bar `0.05/n` |
|---|---:|---:|---:|---:|
| combined | 57,498 | 70,383 | 0.53% | 7.10e-07 |
| EUR | 42,779 | 74,553 | 0.56% | 6.71e-07 |
| AFR | 11,334 | 117,996 | 0.89% | 4.24e-07 |

**The primary bar of 5e-8 is roughly 10× stricter than Bonferroni on the actual test
space, in all three cohorts.** This is not a problem to be fixed — it is why both bars
are reported. A variant at p = 2e-7 would clear correction on the tests we really ran
and fail the threshold the poster used. Pre-registering both (premise P6) means that
case gets described rather than argued about.

## ALT1 frequency spectrum

| ALT1 frequency | combined | EUR | AFR |
|---|---:|---:|---:|
| [0, 0.01%) | 12,080,815 | 12,448,698 | 11,869,076 |
| [0.01%, 0.1%) | 826,073 | 572,134 | 944,450 |
| [0.1%, 1%) | 246,472 | 128,358 | 292,221 |
| **[1%, 5%)** | **35,251** | **37,825** | **86,371** |
| **[5%, 95%)** | **35,132** | **36,728** | **31,625** |
| ≥ 95% | 0 | 0 | 0 |

No variant anywhere has ALT1 frequency ≥ 0.95. REF is the major allele throughout,
so MAF ≥ 0.01 is exactly the two bolded bands and the tested count is exact rather
than a near-count. The mirror bins were requested in order to measure this instead of
assuming it.

## Two things in that table worth carrying forward

**1. The combined cohort tests fewer variants than EUR alone** — 70,383 against 74,553,
despite having 14,719 more samples. Pooling across ancestries pushes variants that are
common in one group below 1% overall. The combined cohort is not a superset of the
ancestry-specific ones in terms of what gets tested, which means a variant can be tested
in EUR and absent from the combined results without anything being wrong.

**2. AFR's extra variants sit entirely next to the threshold.** AFR tests 43,443 more
variants than EUR, and the whole excess is in the [1%, 5%) band (+48,546); in the truly
common band AFR tests *fewer* than EUR (−5,103).

Two causes, and this census cannot separate them:

- African-ancestry samples genuinely carry more low-frequency variation, which is
  expected and is visible in the rarer bands too.
- MAF is estimated with more noise in 11,334 samples than in 42,779. At a true MAF of
  1%, the standard error is ~0.07 percentage points for AFR against ~0.03 for EUR, so
  variants truly just below the threshold scatter above it more often.

The consequence is an asymmetry that matters for interpretation: **AFR gets the loosest
secondary bar (4.24e-7) precisely because it tests the most variants, but half of those
extra tests are in the 1–5% band of its smallest cohort — where power is worst.** The bar
is loosest exactly where the individual tests are weakest. That is an argument for
reporting the primary 5e-8 as the headline and the secondary bar as context, which is
how P6 was already written.

## Files

| | |
|---|---|
| `data/{combined,EUR,AFR}.keep` | plink2 sample lists, verified against the `.fam` |
| `results/census/<cohort>.chr<N>.afreq.alt1.bins` | 66 histograms |
| `results/census_summary.json` | the aggregate, including the full spectrum |

Histogram boundaries were placed at 0.0001, 0.001, 0.01, 0.05 and their mirrors, so
counting at any of those thresholds later is a recount of these files rather than a
re-read of the genotypes.
