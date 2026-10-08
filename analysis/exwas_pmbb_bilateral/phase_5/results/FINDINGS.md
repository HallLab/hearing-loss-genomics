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

## Finding 3 — four known deafness genes all moved the wrong way, and it probably means nothing

Stated prominently because it is the first thing a reader will check, and stated with its caveat
because the caveat is the honest part.

| gene | combined, broad | combined, restricted |
|---|---|---|
| `SIX1` (DFNA23 / BOR) | p 9.3 × 10⁻⁵, **rank 4** | p 0.42, rank 6,703 |
| `COCH` (DFNA9) | p 1.4 × 10⁻³, rank 23 | p 0.08, rank 1,310 |
| `GJB3` (DFNA2B) | p 1.3 × 10⁻³, rank 19 | p 0.017, rank 312 |
| `TMPRSS3` (DFNB8/10) | p 2.6 × 10⁻³, rank 36 | p 0.021, rank 372 |

Four for four, all worse. The obvious reading — that the restriction threw out the cases carrying
the real signal — is **not supported**, for a reason that cuts both ways:

**Neither arm has any signal to throw out.** The p-value distributions are null in both:

```
                 lambda   p<0.01   p<0.001
combined broad    0.679   0.0084   0.00072
         restrict 0.805   0.0107   0.00061
EUR      broad    0.757   0.0104   0.00106
         restrict 0.846   0.0103   0.00112
AFR      broad    0.834   0.0095   0.00090
         restrict 0.677   0.0100   0.00079
```

Essentially identical, and both consistent with nothing. Within a null result, a gene's rank is
noise — `SIX1` at rank 4 in the broad arm was most likely noise too, and its fall is noise moving.
Four genes moving the same way is p ≈ 0.06 on a coin flip, and they are not independent: same
cohort, same controls, overlapping variants.

**So this is a flag, not a finding.** What would settle it is a systematic test against the full
ClinGen hearing-loss gene list rather than four genes picked by hand. That list lives in `cycle_2`;
it is a public endpoint and can be fetched independently into this folder. **It is the single most
useful next step**, and it is the only way to tell "the restriction removed real signal" from "there
was never any signal and ranks wander".

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

**Whether the restriction is right.** It is defensible on clinical grounds and it produced no
finding. Those are both true and neither settles the other. The ClinGen test in Finding 3 is what
would move this from opinion to measurement.

**Whether a bigger cohort would find something.** Both arms are null at this size. 3,164 cases is
small for rare-variant burden, and 6,752 was not large either.

**Audiograms.** Still the direction that would change the question rather than the answer.
