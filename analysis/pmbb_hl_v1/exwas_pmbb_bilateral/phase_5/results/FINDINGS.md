# Phase 5 — bilateral sensorineural ExWAS, results

**Run 2026-10-09.** Step 2 job `49576771` (66/66) plus `49576838` for one retry.
Premise P9 — annotating whatever clears a bar — was registered before step 2 and is
not exercised, because nothing cleared a bar.

## Nothing is significant

| cohort | tested | secondary bar | min p | hits at 5e-8 | hits at secondary |
|---|---:|---:|---:|---:|---:|
| combined | 70,350 | 7.107e-07 | 3.08e-05 | 0 | 0 |
| EUR | 74,566 | 6.705e-07 | 9.47e-06 | 0 | 0 |
| AFR | 118,049 | 4.236e-07 | 1.07e-06 | 0 | 0 |

λ = 1.0159 / 1.0268 / 1.0130. Well calibrated.

## The two arms are nested, which governs everything below

Verified, not assumed:

- all **3,164** bilateral cases are also cases in the broad arm
- the **50,746** controls are the same individuals in both
- the restriction removes **3,588** cases and nothing else

So the bilateral analysis is the broad analysis with cases deleted. Three consequences:

1. **Agreement between the arms measures nothing.** It is the same data twice.
   Spearman ρ is 0.388 / 0.398 / 0.302 and is reported only to show it is *not* being
   used as a replication statistic.
2. **A variant getting weaker says nothing.** Losing 53% of cases does that by itself.
3. **A variant getting stronger cannot come from shared data**, so that is the only
   direction worth reading.

## The strengthening variants, and why they are not reported as findings

Applying "more than one log unit stronger than the broad arm, and p < 1e-4":

| cohort | count | notable |
|---|---:|---|
| combined | 2 | `chr19:44819634`, `chr19:44819705` — both `BCAM` |
| EUR | 4 | the same two `BCAM` variants, OR 1.17 → 1.35 |
| AFR | 8 | `chr7:158133844`, OR 1.80 → 2.90, the cohort minimum at 1.07e-06 |

`BCAM` (basal cell adhesion molecule, chr19q13.32) appears in two cohorts, one variant
missense and one synonymous. On its face this is the arm's most interesting result.

**It is not a result.** The power control in `phase_5/power_control/` drew five random
size-matched case subsets per cohort and applied the identical criterion:

| cohort | bilateral restriction | five random draws | draw mean | Poisson p |
|---|---:|---|---:|---:|
| combined | **2** | 8, 5, 7, 3, 5 | 5.6 | 0.976 |
| EUR | **4** | 5, 6, 3, 7, 3 | 4.8 | 0.706 |
| AFR | **8** | 18, 13, 10, 15, 15 | 14.2 | 0.972 |

The pre-registered rule fires on its first clause in all three: the observed count is
at or below the maximum of the draws. It is also **below the mean of the draws in every
cohort**. Removing 3,588 cases at random produces *more* apparent strengthening than
removing the 3,588 who lack bilateral sensorineural hearing loss.

The decision rule, including what would have counted as inconclusive, was written into
`power_control/README.md` before any of it ran. Without the control, `BCAM` was about
to be written up as this arm's finding.

## Why the null, quantitatively

Detectable odds ratio at 80% power, α = 5e-8, against the broad arm:

| cohort | cases | MAF 1% | MAF 5% | MAF 10% | MAF 50% |
|---|---:|---:|---:|---:|---:|
| bilateral combined | 3,164 | 2.27 | 1.45 | 1.31 | 1.18 |
| broad combined | 6,752 | 1.78 | 1.30 | 1.21 | 1.12 |
| bilateral EUR | 2,547 | 2.50 | 1.52 | 1.35 | 1.20 |
| bilateral AFR | **490** | **7.92** | 2.57 | 1.99 | 1.51 |

Common-variant effects on complex age-related traits sit at OR 1.05–1.15. This arm is
powered below that range everywhere except the most frequent variants in the two larger
cohorts. **AFR, with 490 cases, is effectively uninformative for common variants** — and
it is also the cohort whose random draws produced the most spurious strengthening, which
is the same fact seen twice.

## QC self-check

The same test this project applied to the summer poster's hits, applied here:

| cohort | ρ (call rate vs −log₁₀p) | p | median call rate, top 100 | rest |
|---|---:|---:|---:|---:|
| combined | −0.0001 | 0.98 | 0.99984 | 0.99989 |
| EUR | −0.0034 | 0.35 | 0.99989 | 0.99990 |
| AFR | −0.0031 | 0.29 | 0.99981 | 0.99991 |

No relationship in any cohort. A worse call rate buys no significance here.

## Premise P7 — the declared targets

`KNL1` and `ZFPM1` were named before this arm ran. P7 was amended on 2026-10-09 once
the three variants behind the poster's claims were shown to be absent from the
genotypes and to fail the poster's own call-rate filter — see
`../exwas_elena_replication/phase_4/results/FINDINGS.md`. There is no target in this
data to hit or miss, and the second caveat in P7 stands regardless: those genes were
found on the broad phenotype, so this arm could not have replicated them even had they
existed.

## Figures

| | |
|---|---|
| `figures/qq.png` | calibration, bilateral against broad |
| `figures/manhattan.png` | whole exome, both bars |
| `figures/crossarm.png` | the restriction's effect, strengthening variants circled |
| `figures/power.png` | what halving the cases costs |
| `figures/power_control.png` | observed against the five random draws |
