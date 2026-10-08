# Complement run — do the discarded cases carry the signal?

The power control showed that halving the cases explains most of the restricted arm's null, but not
whether the restriction also removed biology. This asks that directly.

The restricted cases are a **strict subset** of the broad ones — 3,164 of 6,752, checked — so the
complement is exactly the **3,588 the restriction throws away**: hearing loss that is unilateral,
conductive, mixed, or unspecified, plus bilateral sensorineural recorded on a single date.

Both readings were written into `scripts/01_build_cohort.py` before the run.

## What came back

| | cases | ClinGen Definitive/Strong in top 50 | expected | p |
|---|---:|---:|---:|---|
| broad, all | 6,752 | 5 | 0.26 | 5.9 × 10⁻⁶ |
| broad subsampled, five draws | 3,164 | 2 · 0 · 1 · 2 · 3 | 0.26 | — |
| restricted | 3,164 | 0 | 0.26 | 1 |
| **complement** | 3,588 | **1** — `SIX1` | 0.26 | 0.23 |

**Neither half carries it, and neither half is unusual.** Against the random-half distribution
(mean 1.6): restricted P(≤0) = 0.20, complement P(≤1) = 0.53. Both within range. And they do not
add up — 0 + 1 against 5 for the whole set.

**The enrichment needs all 6,752 cases together.** The restriction's null is lost power, and it is
not shown to discard biology.

## On size

The complement is 3,588 against 3,164 for every comparison point — 13% more. It was run at its
natural size rather than cut, because discarding 424 real cases for a round number is worse than
carrying the caveat. The direction is favourable: more cases biases **toward** finding enrichment,
so this null is stronger than a size-matched null would have been.

## What a positive result would have meant

Worth recording, since it did not happen. Known deafness genes associating with hearing loss that is
*not* bilateral sensorineural would have been biologically odd, and the likely cause would not have
been biology but ICD coding: `H91.90 Unspecified hearing loss, unspecified ear` alone covers 2,698
people, and someone with genuine bilateral sensorineural loss is often coded exactly that way. That
would have meant the restriction discards real cases that are merely badly labelled — and it would
have been a measured argument for audiograms rather than an asserted one.

## Reproducing

```
scripts/01_build_cohort.py   the 3,588 and the same 50,746 controls
scripts/02_saige.bsub        null model
scripts/03_step2.bsub        22 chromosomes, whole grid per call
scripts/04_readout.py        every arm side by side
```

Jobs 49546444 and 49546445, run 2026-10-08, 23 tasks, no failures.
