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

## Finding 3 — the broad phenotype carries real signal, and the restriction removes it

**This finding replaces an earlier version of itself.** The first reading of this phase said four
known deafness genes moving the wrong way "probably means nothing", on the grounds that both arms
are null and ranks wander inside a null. The systematic test says otherwise, and the earlier
conclusion was wrong.

### What was tested

Every ClinGen Hearing Loss GCEP gene at Definitive or Strong — **100 genes, 93 of them in our tested
set** — against both arms, by two tests that fail in different ways. The whole list, not four genes
picked by eye.

### Known deafness genes cluster at the top of the broad arm

| arm · cohort | top 50 | top 100 | top 250 | top 500 |
|---|---|---|---|---|
| **broad · combined** | **5 / 0.26 · p = 5.9 × 10⁻⁶** | 5 / 0.52 · p = 1.7 × 10⁻⁴ | 6 / 1.3 · p = 1.9 × 10⁻³ | 7 / 2.6 · p = 0.015 |
| broad · EUR | 2 / 0.26 · p = 0.028 | 2 / 0.52 · p = 0.095 | 3 / 1.3 · p = 0.14 | 4 / 2.6 · p = 0.26 |
| broad · AFR | 0 / 0.26 | 2 / 0.52 · p = 0.096 | 2 / 1.3 · p = 0.38 | 4 / 2.6 · p = 0.27 |
| **restricted · combined** | **0** | **0** | **0** | 2 / 2.6 · p = 0.74 |
| restricted · EUR | 0 | 0 | 1 / 1.3 · p = 0.73 | 4 / 2.6 · p = 0.26 |
| restricted · AFR | 0 | 0 | 2 / 1.3 · p = 0.38 | 4 / 2.6 · p = 0.27 |

Five of the fifty best genes in the broad combined analysis are established deafness genes, where
chance predicts a quarter of one:

```
SIX1     #4   p 9.3e-05     DFNA23 / branchio-oto-renal
GJB3     #19  p 1.3e-03     DFNA2B, connexin 31
COCH     #23  p 1.4e-03     DFNA9
MYO6     #34  p 2.3e-03     DFNA22 / DFNB37
TMPRSS3  #36  p 2.6e-03     DFNB8/10
```

**Under the restricted phenotype, not one appears in the top 250.**

The excess survives Bonferroni over all 60 tests run here (2 arms × 3 cohorts × 2 tiers × 5 cutoffs;
bar 8.3 × 10⁻⁴). It decays smoothly as the cutoff widens, which is the shape a genuine
top-concentrated signal has and a fluke usually does not.

### What this does and does not establish

**It does not make any gene significant.** `SIX1` at 9.3 × 10⁻⁵ is still 33× short of the bar. The
claim is about the *list*, not about any member of it: the ranking carries biology.

**It is not fully independent of what prompted it.** Four of the five drivers are the genes noticed
by hand, which is what motivated running this. The test does use all 93 tested ClinGen genes rather
than those four, and it adds `MYO6`, which was not among them — but a test built after seeing the
pattern it then confirms deserves the caveat stated rather than omitted.

**Power is the live alternative and it is not settled here.** The restricted arm has half the cases,
so its null could be lost power rather than lost signal. Against that: under pure power loss the
ClinGen genes should still sit high, just less sharply. Instead their median rank is 7,104 of 8,971
— no better than a random gene — and none reaches the top 250. That is a stronger absence than power
alone predicts, but it is an argument, not a measurement.

**The measurement that would settle it:** subsample the broad cohort to 3,164 cases, drawn to match
the restricted arm's size, and re-run. If the enrichment survives at matched N, the restriction
removed signal. If it vanishes, the restriction only removed cases. That is one more SAIGE run.

### Mann-Whitney finds nothing, and the disagreement is the point

Across all arms and cohorts, p > 0.12. ClinGen genes are **not** shifted as a body — their median
rank sits at or below the overall median everywhere.

Both tests are right. There is no mass shift of 93 genes; there is a handful at the very top. That
is exactly what a rare-variant study with a few real genes and no power for the rest looks like, and
it is why the top-N test was run alongside the one that uses every gene.

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

**Whether the restriction is right.** Finding 3 moved this from opinion toward measurement, and the
measurement is unfavourable: the broad phenotype's ranking carries known deafness biology and the
restricted one's does not. The restriction remains defensible on clinical grounds — unilateral loss
really is less likely to be genetic — but it now has evidence against it that did not exist when it
was adopted. The size-matched subsample is what would close the argument.

**Whether a bigger cohort would find something.** Both arms are null at this size. 3,164 cases is
small for rare-variant burden, and 6,752 was not large either.

**Audiograms.** Still the direction that would change the question rather than the answer.
