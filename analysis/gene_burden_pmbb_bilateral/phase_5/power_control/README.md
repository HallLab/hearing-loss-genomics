# Power control — is the restricted arm's null a lost signal or lost cases?

Phase 5 found ClinGen deafness genes clustering at the top of the broad analysis and nowhere near
the top of the restricted one. Two explanations survive that, and no amount of staring separates
them:

```
signal   the restriction removed the cases carrying the biology
power    the restriction removed half the cases, and 3,164 is too few for anything
```

They separate by holding N fixed. This draws the **broad** phenotype down to the restricted arm's
exact case count and runs it again, end to end.

**Five replicates, not one.** A single draw can be lucky in either direction, and with one number I
could report whichever answer I happened to draw. The decision rule was written into
`scripts/04_enrichment_by_replicate.py` before the result existed — that mattered here, because I
had already published the claim this test could confirm.

**Combined cohort only.** That is where the enrichment is. EUR and AFR would add 150 jobs to answer
a question neither is powered to settle.

## What came back

| | cases | ClinGen Definitive/Strong in top 50 | expected |
|---|---:|---:|---:|
| broad, all cases | 6,752 | 5 | 0.26 |
| draw 1 | 3,164 | 2 — `COCH`, `SIX1` | 0.26 |
| draw 2 | 3,164 | 0 | 0.26 |
| draw 3 | 3,164 | 1 — `ACTG1` | 0.26 |
| draw 4 | 3,164 | 2 — `COCH`, `TMPRSS3` | 0.26 |
| draw 5 | 3,164 | 3 — `GJB3`, `MYO6`, `SIX1` | 0.26 |
| restricted | 3,164 | **0** | 0.26 |

Pooled over the draws: **8 ClinGen genes in 250 top-50 slots against 1.3 expected, p = 6.4 × 10⁻⁵.**
The enrichment is not an artifact of sample size.

But the draws average 1.6, P(0 | Poisson 1.6) = 0.20, and one draw also returned zero. So the fall
from 5 to 1.6 is the cost of halving the cases, and the fall from 1.6 to 0 is not distinguishable
from sampling variation with one observation.

## What it did not settle, and what did

This left open whether the restriction removes signal *on top of* removing cases. The complement run
in `../complement/` closed it: the discarded cases do not carry the enrichment either.

## Reproducing

```
scripts/01_draw_subsamples.py        five covariate files, seed 20261008
scripts/02_step1.bsub                5 null models
scripts/03_step2.bsub                110 tasks, 5 replicates x 22 chromosomes
scripts/04_enrichment_by_replicate.py
```

Jobs 49545638 and 49545646, run 2026-10-08, 115 tasks, no failures.
