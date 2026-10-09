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

**Correction to the reasoning originally recorded here.** The script's header says a
per-variant frequency file was avoided because plink2 can only compress output with
zstd and this environment lacks both the CLI and python's `zstandard`. That is wrong:
`zs` is an optional modifier, and plain `--freq` writes uncompressed text. The
histogram is still the right output for this phase — it is three orders of magnitude
smaller and answers the question asked — but it was chosen for a reason that does not
hold. Per-variant counts were in fact readable all along, and were used to diagnose
the discrepancy below.

---

## Reconciliation against what SAIGE actually tested

Phase 4 ran, and the counts above are **not** what SAIGE tested:

| cohort | plink2 census | SAIGE tested | difference |
|---|---:|---:|---:|
| combined | 70,383 | 70,366 | −17 (−0.024%) |
| EUR | 74,553 | 74,536 | −17 (−0.023%) |
| AFR | 117,996 | 117,953 | −43 (−0.036%) |

41 of the 66 cells differ, by 1 to 6 variants. **Every single difference is negative** —
SAIGE never tests more than the census predicts.

### The mechanism, verified

The two tools compute MAF over different denominators:

| | denominator |
|---|---|
| plink2 | `ALT_CTS / OBS_CT` — observed genotypes only |
| SAIGE | `ALT_CTS / 2N` — the whole cohort |

For a variant with missing genotypes, plink2's MAF is therefore the larger of the two,
and a variant can sit above 1% for plink2 and below it for SAIGE. The reverse is
impossible, which is why all 41 discrepancies point the same way.

Checked variant by variant on AFR chr1, the worst cell (−6):

| variant | ALT | OBS_CT | 2N | MAF obs | MAF 2N |
|---|---:|---:|---:|---:|---:|
| 1:18885456:G:A | 226 | 22,582 | 22,668 | 0.0100080 | 0.0099700 |
| 1:21561304:C:A | 226 | 22,544 | 22,668 | 0.0100248 | 0.0099700 |
| 1:114720447:T:C | 226 | 22,588 | 22,668 | 0.0100053 | 0.0099700 |
| 1:156596312:T:A | 225 | 22,442 | 22,668 | 0.0100258 | 0.0099259 |
| 1:161675397:A:G | 224 | 22,388 | 22,668 | 0.0100054 | 0.0098818 |
| 1:248639253:G:A | 221 | 22,088 | 22,668 | 0.0100054 | 0.0097494 |

All six straddle the threshold in exactly the predicted way, and zero variants went the
other direction (tested by SAIGE but not counted by plink2).

### What this means for the filter

**SAIGE's effective MAF floor is stricter than a nominal 1% for any variant carrying
missingness.** At 5% missing it requires an observed MAF of about 1.053% to pass. The
tested set is therefore slightly smaller and slightly cleaner than the census suggested
— boundary variants with the most missing data are the ones dropped.

### What this means for the bar

The denominator of the secondary bar is what SAIGE tested, not what plink2 counted:

| cohort | pre-registered bar | final bar |
|---|---:|---:|
| combined | 7.104e-07 | **7.106e-07** |
| EUR | 6.707e-07 | **6.709e-07** |
| AFR | 4.237e-07 | **4.239e-07** |

The bars move in the fourth significant figure. **The pre-registration held** — this is
the outcome that makes the exercise worth doing: the threshold was fixed before any
p-value existed, and reconciling it afterwards against the true test count changed
nothing that any conclusion rests on. The primary bar of 5e-8 is unaffected, being a
convention rather than a count.
