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

## Finding 3 — the ranking carries real biology; whether the restriction removes it is unresolved

**This finding has been rewritten twice, and the history is kept because it is the useful part.**
First version: four known genes moved the wrong way, "probably means nothing". Second: the ClinGen
test showed enrichment in the broad arm and none in the restricted one, so "the restriction removes
signal". The size-matched control says the second was overstated.

### What the ClinGen test found

Every ClinGen Hearing Loss GCEP gene at Definitive or Strong — 100 genes, 93 in our tested set.
In the broad combined analysis, five sit in the top 50 where chance predicts 0.26:

```
SIX1  #4  ·  GJB3  #19  ·  COCH  #23  ·  MYO6  #34  ·  TMPRSS3  #36      p = 5.9 × 10⁻⁶
```

In the restricted arm, none of the 93 reaches the top 250.

### What the size-matched control found

The broad cohort was drawn down to 3,164 cases — the restricted arm's exact count — five times, and
re-run end to end. 115 jobs. The decision rule was written into the script before the result existed.

| | ClinGen in top 50 | expected |
|---|---:|---:|
| broad, 6,752 cases | **5** | 0.26 |
| broad at 3,164 cases, five draws | **2 · 0 · 1 · 2 · 3** | 0.26 each |
| restricted, 3,164 cases | **0** | 0.26 |

Two things follow, and they point different ways.

**The enrichment is not an artifact of sample size.** Pooled over the five draws, 8 ClinGen genes
fall in 250 top-50 slots against 1.3 expected — p = 6.4 × 10⁻⁵. The broad phenotype's ranking still
carries known deafness biology when it has only 3,164 cases to work with. That is worth having
established: it means the ordering of these lists is not pure noise, which the phase's own null
result might otherwise suggest.

**But the restricted arm's zero is inside the range chance produces.** The draws average 1.6, and
P(0 | Poisson 1.6) = 0.20. One of the five draws also returned zero. So:

```
5  ->  1.6    attributable to halving the cases
1.6 ->  0     not distinguishable from sampling variation, with one observation
```

### The honest statement

The restriction's null is **consistent with lost power alone**. It may also have removed signal; the
data here cannot separate the two, because the restricted arm is one fixed set rather than a draw,
and its result sits where a low draw would sit.

My second version of this finding said the restriction removes signal. That claim is withdrawn.
What survives is the weaker and better-supported one: the ranking carries biology, and the
restricted arm is too small to show it.

### The test that would resolve it

The restricted cases are a strict subset of the broad ones — 3,164 of 6,752, with the complement
being exactly the **3,588 the restriction excludes**. Running that complement as its own case set
asks the question directly:

- complement shows the enrichment → the signal lives in the cases the restriction throws away, and
  the phenotype decision is costly
- complement shows nothing → the signal lives in the cases the restriction keeps, and its null here
  is power, full stop

One SAIGE run, same size as the others. It is the cleanest remaining question in this analysis.

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

**Whether the restriction is right.** Still open, and now open with numbers rather than opinion. The
broad ranking carries known deafness biology that survives size-matching; the restricted arm shows
none, but its null is consistent with power alone. The restriction remains defensible on clinical
grounds and is not condemned by this evidence. The complement run in Finding 3 is what would
settle it.

**Whether a bigger cohort would find something.** Both arms are null at this size. 3,164 cases is
small for rare-variant burden, and 6,752 was not large either.

**Audiograms.** Still the direction that would change the question rather than the answer.
