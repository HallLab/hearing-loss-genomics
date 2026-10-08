# Phase 5 — reading the result

**Status:** complete · **Scripts:** `03_merge_results.py`
**Inputs:** 66 step-2 files, all present, all structurally whole

---

## The result

| cohort | genes | **omnibus p** | bar 0.05/genes | significant | best single cell |
|---|---:|---|---|---:|---|
| combined | 17,941 | 2.05 × 10⁻⁵ | 2.79 × 10⁻⁶ | **0** | 4.85 × 10⁻⁶ |
| EUR | 17,921 | 1.37 × 10⁻⁵ | 2.79 × 10⁻⁶ | **0** | 1.52 × 10⁻⁶ |
| AFR | 17,778 | 9.14 × 10⁻⁶ | 2.81 × 10⁻⁶ | **0** | 3.40 × 10⁻⁶ |

**Nothing is significant in any cohort.** The headline is SAIGE's own `Cauchy` row — one p-value per
gene with the nine-cell search already priced in (premise **P10**), so the bar is unambiguous.

Top by omnibus: combined `PCBD1`, `TJAP1`, `LSMEM2` · EUR `FBXO25`, `PCBD1`, `TJAP1` ·
AFR `IGSF9`, `OR6Y1`, `RAET1E`.

---

## Finding 1 — the restricted phenotype does not rescue the result

This is the finding the phase exists to produce, and it is negative.

The restriction was adopted on clinical grounds: unilateral and conductive hearing loss are less
likely to be genetic, so excluding them should leave a cleaner phenotype. It did not produce a
cleaner answer.

| | broad `SO_396` | bilateral sensorineural |
|---|---|---|
| cases | 6,752 | 3,164 |
| best omnibus, combined | 3.02 × 10⁻⁵ | 2.05 × 10⁻⁵ |
| anything significant | no | no |

Both arms are null. The restriction bought a marginally smaller best p-value at the cost of half the
cases, and crossed no threshold.

---

## Finding 2 — the top of the list reshuffles almost completely

| cohort | shared in top 50 |
|---|---:|
| combined | **9 / 50** |
| EUR | 11 / 50 |
| AFR | 6 / 50 |

Four fifths of each list changes. Against the replication's own comparison — 24, 17 and 15 of 50
between mask definitions — **the phenotype moves the ranking more than the mask defects did.**

That is worth holding onto. Phases 2 to 4 of the replication documented four mask defects and
surplus PCs and showed they reshuffled the top of the list. Changing who counts as a case reshuffles
it harder.

---

## Finding 3 — the ranking carries real biology, and the restriction is exonerated

**Written three times. The sequence is kept because it is the useful part:** "four known genes moved
the wrong way, probably nothing" → "the restriction removes signal" → this, which two controls
settled and which withdraws the middle version.

### The observation

Every ClinGen Hearing Loss GCEP gene at Definitive or Strong — 100 genes, 93 in our tested set.
In the broad combined analysis five sit in the top 50 against 0.26 expected:

```
SIX1 #4  ·  GJB3 #19  ·  COCH #23  ·  MYO6 #34  ·  TMPRSS3 #36        p = 5.9 × 10⁻⁶
```

In the restricted arm, none of the 93 reaches the top 250. The tempting reading — that restricting
the phenotype threw out the cases carrying the biology — took two controls to test, and did not
survive either.

### Control 1 — hold the sample size fixed

The broad cohort drawn down to 3,164 cases, the restricted arm's exact count, five independent
times, re-run end to end. Decision rule written into the script before the result existed.

### Control 2 — run the cases the restriction discards

The restricted cases are a strict subset of the broad ones, so the complement is exactly the 3,588
the restriction throws away: unilateral, conductive, mixed, and mostly unspecified hearing loss.

| | cases | ClinGen in top 50 | expected | p |
|---|---:|---:|---:|---|
| **broad, all cases** | 6,752 | **5** | 0.26 | 5.9 × 10⁻⁶ |
| broad subsampled, draws 1–5 | 3,164 | 2 · 0 · 1 · 2 · 3 | 0.26 | — |
| **restricted** (bilateral SN) | 3,164 | **0** | 0.26 | 1 |
| **complement** (what it discards) | 3,588 | **1** — `SIX1` | 0.26 | 0.23 |

### What the two controls establish

**The enrichment is real and survives size-matching.** Pooled over the five draws, 8 ClinGen genes
in 250 top-50 slots against 1.3 expected, p = 6.4 × 10⁻⁵. The broad phenotype's ranking carries
known deafness biology even at 3,164 cases. Worth having established against a phase whose headline
is null: the ordering of these lists is not pure noise.

**Neither half carries it, and neither half is unusual.** Against the random-half distribution
(mean 1.6):

```
restricted   0   P(<= 0) = 0.20    within range
complement   1   P(<= 1) = 0.53    within range
```

Neither is distinguishable from a random half of the same size. And the halves do not add up:
0 + 1 = 1 against 5 for the whole. **The enrichment needs the full 6,752 together.**

### The conclusion, and the claim withdrawn

**The restriction's null is lost power, not lost signal.** The phenotype decision is not shown to
discard biology, and the second version of this finding — which said it did — is withdrawn.

The reason the halves do not carry it is ordinary and worth stating plainly: at 3,164 cases nothing
in this study is detectable, so which 3,164 you pick barely matters. The signal is thin enough that
it only clears the noise when every case is in.

`SIX1` deserves a note as the one gene that shows up almost everywhere it could: top 50 in the full
broad arm, in two of five random draws, and the only one in the complement's top 50. Still 33× short
of significance, still not a finding — but it is the most consistent thing in this analysis.

---

## Finding 4 — AFR behaved as predicted, and that prediction should have carried more weight

`tau₂` converged to 0 for AFR, so SAIGE reduced to logistic regression with covariates: at 490
cases there is no detectable polygenic component. The replication saw the same at 1,285 cases.

The AFR arm was run because dropping it is its own kind of reporting bias. It should be read as a
null with almost no power, not as a search. `IGSF9` tops it at 9.14 × 10⁻⁶ against a bar of
2.81 × 10⁻⁶, which is the closest any cohort in either analysis has come — and in the cohort least
able to support it.

---

## What this phase does not answer

**Whether the restriction is right.** Settled as far as this data can settle it: the restriction is
not shown to discard signal, and its null is power. It stays defensible on clinical grounds and now
has a measurement behind it rather than only an argument. What this does **not** say is that the
restriction is better — only that it is not worse for the reason that was suspected. At 3,164 cases
neither definition can show anything, which is a statement about the cohort rather than about either
phenotype.

**Whether a bigger cohort would find something.** Both arms are null at this size. 3,164 cases is
small for rare-variant burden, and 6,752 was not large either.

**Audiograms.** Still the direction that would change the question rather than the answer.
