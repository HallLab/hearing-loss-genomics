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

*(filled in by `04_readout.py` after the runs complete)*
