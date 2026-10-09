# Power control — is "stronger on half the cases" real, or what case removal does anyway?

**Written before the control ran. The decision rule below was fixed first.**

## The question

Phase 5 found variants whose association gets **stronger** when the phenotype is
restricted from any hearing loss to bilateral sensorineural, despite the restriction
removing 3,588 of 6,752 cases:

| cohort | variants stronger by >1 log unit, reaching p < 1e-4 |
|---|---:|
| combined | 2 |
| EUR | 4 |
| AFR | 8 |

The appealing reading is subtype concentration: the association lives in bilateral
sensorineural, so sharpening the phenotype sharpens the signal even on fewer cases.

The boring reading is selection. Roughly 70,000 to 118,000 variants are compared per
cohort. Remove half the cases at random and every p-value moves; some move the right
way by chance, and picking out the ones that did is choosing winners after the race.

Nothing in the numbers themselves separates these.

## The test

Hold the case count fixed and remove the phenotype. Draw cases **at random** from the
broad arm, down to the bilateral arm's exact count, and run the whole thing again:

| cohort | broad cases | drawn to | controls (unchanged) |
|---|---:|---:|---:|
| combined | 6,752 | 3,164 | 50,746 |
| EUR | 5,183 | 2,547 | 37,596 |
| AFR | 1,285 | 490 | 10,049 |

A random draw has no subtype meaning. Whatever strengthening it produces is the
floor — the amount that case removal alone buys.

**Five replicates.** One draw can be lucky in either direction, and with a single
number there would be a choice about which answer to report.

**All three cohorts**, unlike the burden arm's control which ran combined only. Here
the observed effect appears in all three, and EUR and AFR carry the strongest cases
(chr19:44819634 at 9.47e-06 in EUR, chr7:158133844 at 1.07e-06 in AFR).

Step 1 is refitted for every draw. The null model depends on who the cases are, so it
cannot be inherited the way the main arm inherits it.

## The decision rule, fixed now

For each cohort, let **k_obs** be the count in the table above and **k_1..k_5** the
counts from the five random draws, using the identical criterion (>1 log unit stronger
than the broad arm, and p < 1e-4 in the restricted run).

- **k_obs ≤ max(k_1..k_5)** → the strengthening is within what random case removal
  produces. **Conclusion: selection, not subtype.** The variants are not reported as
  findings.
- **k_obs > max(k_1..k_5)**, and a Poisson test of k_obs against the pooled draw rate
  gives **p < 0.05** → evidence the restriction concentrates signal beyond chance.
  Reported as suggestive, never as significant: no variant clears 5e-8 in either arm.
- **Anything between** → reported as inconclusive, with both counts shown.

The Poisson comparison pools the five draws rather than testing against the maximum,
because five counts estimate a rate better than they estimate an extreme.

**This rule is not revisited after the counts come back.** The burden arm's control
cost a claim I had already written down; the point of writing the rule first is that
the same thing can happen again and be survivable.

## What came back

Job `49576788` (15 step 1 fits) and `49576789` (330 step 2 tasks), plus `49576838`
re-running one task that exited before SAIGE started.

| cohort | bilateral restriction | five random draws | draw mean | Poisson p | verdict |
|---|---:|---|---:|---:|---|
| combined | **2** | 8, 5, 7, 3, 5 | 5.6 | 0.976 | selection |
| EUR | **4** | 5, 6, 3, 7, 3 | 4.8 | 0.706 | selection |
| AFR | **8** | 18, 13, 10, 15, 15 | 14.2 | 0.972 | selection |

**The rule fires on its first clause in all three cohorts: `k_obs ≤ max(draws)`.** The
strengthening variants are not reported as findings.

### It is not marginal, and it runs the other way

In every cohort the observed count is **below the mean of the random draws**, and in
combined and AFR it is below or equal to every individual draw. Removing 3,588 cases
*at random* produces more variants that "get stronger" than removing the 3,588 who do
not have bilateral sensorineural hearing loss.

So the restriction is not merely failing to concentrate signal — it is slightly less
productive of apparent strengthening than chance. Whatever `BCAM` looked like, there
is nothing there.

### Why the draws produce so many

Roughly 70,000 to 118,000 variants are compared per cohort. Remove half the cases and
every p-value moves; a handful move more than one log unit in the lucky direction and
land under 1e-4 purely because the estimate is noisier on fewer cases. AFR shows this
most strongly (mean 14.2) because 490 cases is where the estimates are noisiest.

Reading the top of that list as biology is choosing winners after the race. That is
exactly what this control exists to prevent, and it is what would have happened here:
`chr19:44819634` and `chr19:44819705` in `BCAM` were about to be written up as the
arm's most interesting result.

### What this does not say

It does not say bilateral sensorineural is the wrong phenotype. It says that **at this
cohort size the restriction cannot demonstrate subtype specificity**, because the
noise introduced by halving the cases is larger than any sharpening it might produce.
With 3,164 cases the question is unanswerable either way — the same conclusion the
burden arm's power control reached, by a different route.
