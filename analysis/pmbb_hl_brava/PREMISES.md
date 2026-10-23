# Pre-registered premises — the BRaVa annotation arm

**Written 2026-10-23, before anything ran.** Nothing in this file was chosen after
seeing a result.

## Why this arm exists

Two reasons, and the second is the one that justifies the compute.

**1. Comparability.** `analysis/pmbb_hl/confluence/hl_pmbb_v4_results.md` §10 proposes
replicating in All of Us and UK Biobank. Those consortia annotate to the **BRaVa**
standard (`docs/papers/BRaVa annotation recommendations.pdf`). Measured on chr21, our
`pDM` and BRaVa's overlap at **Jaccard 0.20** — they are not the same definition with a
different threshold, they are different variant sets. A result of ours placed beside a
BRaVa-annotated result from another biobank would differ for reasons no one could later
separate from biology.

**2. A robustness test on the project's only positive result.** The ClinGen enrichment
(5 of the top 50 against 0.26 expected, p = 5.91e-06, broad combined) has survived
changing the PCA source, raising PC counts from 5/4/3 to 11/15/16, and adding
deduplication before the combination. It has **never been tested against a different
annotation**. If it survives a mask set that shares only a fifth of its variants, it is
very hard to attribute to annotation artefact. If it does not survive, that is the most
important thing this project would learn this month.

## P1 — the rules, stated exactly

From the BRaVa document, categories in precedence order. **Each variant receives one
annotation**; lower numbers win.

```
pLOF        LOFTEE HC                                        (BRaVa category 1)

pDM         not category 1, AND any of:                      (BRaVa category 2)
              a. consequence in {missense, start_lost, stop_lost, inframe_insertion,
                 inframe_deletion} AND (REVEL >= 0.773 OR CADD_PHRED >= 28.1)
              b. SpliceAI DS >= 0.2, DS = max(DS_AG, DS_AL, DS_DG, DS_DL)
              c. LOFTEE LC

synonymous  synonymous_variant AND SpliceAI DS < 0.2         (BRaVa category 4)

pLOF_pDM    the union of pLOF and pDM, pLOF taking precedence
```

**Transcripts:** MANE Select where available, otherwise `CANONICAL`. This differs from
`pmbb_hl`, which reads every transcript row and takes the maximum REVEL across them.

**Frequency:** drop any variant whose maximum gnomAD broad-ancestry AF exceeds 0.01,
across `gnomADe_{AFR,AMR,EAS,NFE,SAS}_AF`. `pmbb_hl` applies no gnomAD filter — it
filters on the **internal** MAF at test time, which is a different thing: a variant rare
in PMBB and common in gnomAD is kept by us and dropped by BRaVa.

**Not built:** BRaVa category 3 (other missense). It is the complement of category 2
within missense, and no question here needs it.

## P2 — what is inherited, unchanged

| | |
|---|---|
| cohorts | `pmbb_hl/rare_variant/phase_1` — identical people, identical controls |
| covariates and PCs | `pmbb_hl/rare_variant/phase_3` — 11 / 15 / 16 at 80% variance |
| GRM | `pmbb_hl/rare_variant/phase_3` — relatedness does not depend on the annotation |
| **SAIGE step 1** | **`pmbb_hl/rare_variant/phase_4` — the null model does not depend on the mask, so the nine fitted models are reused rather than refitted** |
| genotypes | the same `chr*_deduplicated` files |

**Only the mask changes.** That is the whole point: anything that differs downstream is
attributable to annotation and to nothing else.

**The comparison is against `pmbb_hl`, not `pmbb_hl_v1`.** The v1 folder is the retired
first pass — its own PCA, 5/4/3 PCs, no tinnitus, no deduplication. Comparing against it
would confound annotation with ancestry correction, which is the exact confusion that
cost this project a rebuild.

## P3 — scope, and why it is not all eighteen cells

**Two cells: broad combined and broad EUR.** Those are the only arms where the ClinGen
enrichment exists (5 and 2 in the top 50; every other arm is 0 against 0.26 expected).
A robustness test is only meaningful where there is something to be robust.

That is ~44 step-2 tasks rather than 198. If the result is interesting, expanding is a
separate decision made with the numbers in hand.

### Addendum, 2026-10-23 — scope expanded to AFR, after the fact

**Recorded because it happened after the declared test was answered**, which is exactly
the circumstance pre-registration exists to constrain.

The AFR cohort was deliberately out of scope above: it carries 0 ClinGen genes in the
top 50 under our own annotation, so there was nothing there for a robustness test to be
robust about. It is now being run so the annotation comparison figure covers all three
cohorts rather than two.

**This cannot change the declared test and is not allowed to.** P5 names broad
*combined*, Definitive/Strong, top 50. That test was computed, reported, and is closed.
Nothing AFR produces may be substituted for it, added to it, or used to re-describe it.

If the AFR arm produces something interesting it is a **hypothesis**, generated after
seeing other results, and must be labelled as one wherever it appears.

## P4 — the synonymous mask is a control, and its expected result is nothing

Synonymous variants with no splice signal should not associate with hearing loss. If the
synonymous mask produces enrichment at the top of the ranking, or a λ far from 1, that
indicates a problem in the pipeline that **no other check in this project would catch** —
not annotation, not the phenotype, but something structural in how cases and controls
differ.

**Expected: a flat, boring result.** Recorded now so that a boring result is reported as
the confirmation it is, rather than quietly dropped for being uninteresting.

## P5 — the expected outcome, and what counts as "survived"

Stated in advance, because "the enrichment survived" is easy to rationalise afterwards.

**The masks will differ substantially.** chr21 gives Jaccard 0.20 on `pDM` and ~15% of
our `pLOF` rows falling outside LOFTEE HC. Expect the exome-wide figures to be close to
these.

**`pLOF` will shrink and `pDM` will change composition.** LOFTEE HC is stricter than our
term matching; CADD ≥ 28.1 admits variants we never consider; REVEL 0.773 excludes a
band we currently keep.

**The pre-declared test:** ClinGen Definitive/Strong genes in the **top 50** of the broad
combined ranking, hypergeometric p, exactly as `pmbb_hl/rare_variant/phase_5` computes it.

| outcome | reading |
|---|---|
| p < 0.05 | survives. The enrichment is not an artefact of our annotation choices |
| p >= 0.05 | does not survive. **This must be reported as prominently as the original finding**, and the original finding re-described as annotation-dependent |

**No other statistic will be substituted after the fact.** The top-50 depth and the
Definitive/Strong tier are fixed here, before the run, because both are choices that
could otherwise be made to favour the answer we want.

## P6 — what this arm cannot say

**It cannot say which annotation is better.** BRaVa is a consensus standard, not ground
truth. Ours uses AlphaMissense, which postdates the BRaVa document; theirs uses CADD and
LOFTEE, which ours does not. Neither is validated against the other here, and nothing in
this arm licenses calling one correct.

**It cannot rescue a null.** Nothing was significant under our annotation, and nothing is
expected to be significant under BRaVa's. If something does cross a bar it is a finding
about a single gene in a single annotation and will be treated as a hypothesis, not a
result.

**A different result is not necessarily a worse result.** If the enrichment weakens, the
honest conclusion is "this finding depends on annotation choices", not "BRaVa is wrong"
or "we were wrong". The two mask sets share a fifth of their variants; some movement is
expected on that basis alone.
