# Phase 5 — calibration, concordance, and the poster audit

**Run 2026-10-09.** Scripts: `01_calibration_and_concordance.py`, `02_poster_audit.py`,
`03_figures.py`. Figures in `figures/`.

Premise P9 — annotating whatever clears a bar — was registered before step 2 and is not
exercised here, because nothing cleared a bar. Recorded so that its absence reads as a
null result rather than as a step quietly dropped.

## Calibration

Lambda is legitimate in this arm in a way it was not in the burden arm: the genomic
inflation factor assumes a null statistic that is chi-square on 1 df, which is what a
single-variant test produces. The exact `qchisq(0.5, 1) = 0.4549364` is used, not the
0.456 rounded in her `step10_1_QQplot.R`.

| cohort | λ ours | λ hers |
|---|---:|---:|
| combined | 0.9893 | 0.9919 |
| EUR | 0.9991 | 0.9907 |
| AFR | 1.0207 | 1.0120 |

Both runs are well calibrated. The poster's claim of "adequate statistical calibration
and minimal evidence of systematic inflation" is correct and is reproduced.

## Concordance — the headline of this arm

| cohort | ours | hers | shared | Pearson r | Spearman ρ |
|---|---:|---:|---:|---:|---:|
| combined | 70,366 | 70,429 | 70,333 (99.9%) | 0.989 | 0.978 |
| EUR | 74,536 | 74,528 | 74,483 (99.9%) | 0.979 | 0.962 |
| AFR | 117,953 | 117,943 | 117,723 (99.8%) | 0.990 | 0.981 |

ρ is reported alongside r because with nothing significant, r is dominated by the bulk
near zero where both runs are measuring noise; rank agreement is the stricter test and
it holds.

**Every disagreement larger than one log unit is in the MHC** — 14 of 14 in EUR, 16 of
16 in AFR. Not most: all.

## A limitation of our own, found here

Splitting λ by region exposes a cost of our PC deviation that the whole-genome figure
hides:

| cohort | λ outside MHC (ours / hers) | λ inside MHC (ours / hers) |
|---|---|---|
| combined | 0.9886 / 0.9903 | 1.0410 / **1.1448** |
| EUR | 0.9950 / 0.9918 | **1.3386** / 0.9651 |
| AFR | 1.0195 / 1.0123 | **1.0921** / 0.9975 |

Outside the MHC the two runs are indistinguishable. Inside it, in EUR and AFR, **our 4
and 3 PCs leave residual inflation that her 9 and 10 do not.** EUR is the clear case:
λ = 1.34 against her 0.97, on 1,032 variants.

Our PC counts were chosen on scree plots across the genome, and that choice is still
defensible for the genome as a whole — outside the MHC it costs nothing. But the MHC is
the region where population structure is strongest, and there the extra PCs are doing
real work. Her choice is better in that window.

This changes no conclusion: λ = 1.34 on 1,032 variants does not manufacture a 5e-8 hit,
and the minimum p in EUR is 2.12e-05. It is recorded because it is a limitation of our
pipeline, found by looking, and because any future arm that cares about the MHC should
not inherit 4 PCs without reopening the question.

In combined the position reverses — hers is the inflated one, at 1.14 against our 1.04.

## The poster audit

`02_poster_audit.py` reads the claims out of the `.pptx` rather than retyping them,
measures each variant on the PMBB v4 release with plink2, and traces it through the
pipeline stages. It refuses to write a report if any variant cannot be measured, so that
a failed measurement can never be printed as an absent variant.

Six claimed results, three unique variants:

| gene | variant | cohort | poster p |
|---|---|---|---|
| KNL1 | chr15:40611486:A:G | Combined | 3.51e-12 |
| KNL1 | chr15:40611486:A:G | EUR | 7.21e-12 |
| ZFPM1 | chr16:88533289:A:C | EUR | 2.46e-10 |
| ZFPM1 | chr16:88533290:G:C | EUR | 2.52e-10 |
| ZFPM1 | chr16:88533290:G:C | Combined | 4.54e-09 |
| ZFPM1 | chr16:88533289:A:C | Combined | 6.18e-09 |

Measured on the release — **not on our pipeline, and not on hers**:

| variant | call rate | MAF | passes the poster's stated QC? |
|---|---:|---:|---|
| chr15:40611486:A:G | **0.2013** | 0.284 | no |
| chr16:88533289:A:C | 0.9617 | 0.295 | no |
| chr16:88533290:G:C | 0.9473 | 0.283 | no |

The poster's methods read: *"variants with a genotype call rate < 0.99 and a Minor Allele
Frequency (MAF) < 0.01 were removed."* None of the three clears 0.99. The variant
carrying the 3.51e-12 is genotyped in 20% of the cohort.

For scale, the call rate of what this arm actually tested: median 0.9999, minimum 0.9895
(combined), 0.9875 (EUR), 0.9534 (AFR). A small tail sits below 0.99 — 0.04%, 0.18% and
0.91% — because our rate is computed within a cohort while the QC was applied across all
70,925 samples, so a variant at 0.995 overall can read 0.985 inside AFR. Nothing tested
anywhere is below 0.95. The poster's worst variant is not in that tail; it is in a
different regime.

### Where each stops appearing

| stage | chr15 | chr16 |
|---|---|---|
| PMBB v4 release | present (546,760) | present (780,404) |
| `variant_qc/chrN_callrate.bim` | **gone** (476,227) | present (670,730) |
| `plink_deduplicated/chrN_allele_ids.bim` | gone (469,424) | **gone** (662,670) |
| `chrN_deduplicated.bim` — what step 2 reads | gone (469,422) | gone (662,667) |

A true 0.99 call-rate filter on the release leaves 475,376 variants on chr16. Her
`chr16_callrate.bim` holds 670,730. The stage named for call rate is therefore not
filtering at 0.99, which is how variants at 0.96 and 0.95 passed through it. The chr15
variant at 0.20 was removed there; the chr16 pair survived and was dropped later.

### Where the genes are present, they carry no signal

| | her delivered results | poster |
|---|---|---|
| KNL1 | 7 variants, best p **0.33** | 3.51e-12 |
| ZFPM1 | 7 variants, best p **0.26** | 2.46e-10 |

## What the arm concludes

Nothing is significant: min p 3.02e-06 / 2.12e-05 / 4.11e-06 against bars of 5e-8 and
~7e-07. Her delivered results say the same, and her own `Significance` column reads
"Not significant" throughout.

**This is a successful replication of the delivered pipeline and a failed replication of
the poster, and those are two different objects.** The arm Doug deprioritised as "just
contributing noise" produced the most consequential finding of the project — by failing
to reproduce a result rather than by finding one.

## Figures

| | |
|---|---|
| `figures/qq.png` | calibration, both runs, per cohort |
| `figures/manhattan.png` | whole exome, both bars drawn |
| `figures/concordance.png` | per-variant agreement, MHC marked |
| `figures/poster_variants.png` | the poster's three variants against the call rate of what was tested |
