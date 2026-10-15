# Pre-registered premises — the tinnitus arm

**Written 2026-10-15, before the arm ran.** Nothing in this file was chosen after
seeing a p-value.

## Why this arm exists

This project began with `ZNF175`. Park et al. 2021 (*Nature Medicine*) reported it
associated with **tinnitus** in an earlier PMBB release, and the mouse ortholog `Zfp719`
produces hearing loss when disrupted. That pairing is what motivated extending the work
to PMBB v4.

**And no arm of this project has tested tinnitus.** Broad hearing loss, bilateral
sensorineural, and the common-variant versions of each — all four test hearing loss.
The hypothesis of origin has never been examined with a corrected pipeline.

## P1 — the phenotype

**Case:** tinnitus evidence on **≥ 2 distinct dates**, from either source:
- phecodeX `SO_397` or any child, from `conditions_phecode_x`
- the standard ICD codes (388.3x, H93.1x) from the OMOP `observation` table

The second source is not optional. PMBB v4 moved those codes out of the phecode files,
which keep 3,016 of 28,110 tinnitus events; the other 25,094 exist only in
`observation`. An arm built without it would miss 89% of its own phenotype.

**Control:** the same 50,755 used by the other arms — no ear-family evidence from any
source. Shared by construction, not by coincidence.

**Expected: 2,859 cases, 53,614 analysable.** Measured before this file was written, as
the basis for the power statement below; it is not a result.

## P2 — the arms are NOT nested, and this one breaks the phase 1 assertion

Every bilateral case is a broad case, and phase 1 refuses to proceed if that fails.
**Tinnitus is not a subset of hearing loss**, so that assertion must not be extended to
it. The check is scoped to the pair it describes.

**66% of tinnitus cases are also broad hearing-loss cases**, and only 581 have tinnitus
with no hearing-loss evidence at all. So this is an overlapping slice, not an
independent cohort, and any comparison between arms must say so.

## P3 — inherited unchanged

| | |
|---|---|
| masks | phase 2's, unchanged — a mask depends on annotation and genotypes, not on who is a case |
| GRM | phase 3's, unchanged — relatedness does not depend on the phenotype |
| PC counts | 11 / 15 / 16 from the release at 80% variance |
| bars | burden `0.05/genes` + FDR q < 0.10 · ExWAS `5e-8` + `0.05/variants` |

New work: the MAF census (the sample set differs, so MAF does), three null models, and
step 2 for both analyses.

## P4 — `ZNF175` is the declared target

Its p-value will be reported **whatever it is**, in its own named section, in all six
cells (2 analyses × 3 cohorts). Not left to be found in a top-50 table.

Unlike `KNL1` and `ZFPM1` in the ExWAS arm — whose variants turned out not to exist in
the delivered data — **this target exists**: `ZNF175` is in the mask and was tested by
the original run. There is something here to hit or miss.

For reference, its values in the original run's delivered tables:

| | p |
|---|---|
| tinnitus AFR | 0.064 |
| tinnitus EUR | 0.139 |
| tinnitus combined | 0.331 |
| hearing loss EUR | 0.022 |
| hearing loss combined | 0.040 |

Those come from the pipeline carrying the four mask defects this project corrected, with
5/9/10 PCs, so they are a landmark rather than a baseline. Note that they run the wrong
way for the hypothesis — `ZNF175` does better on hearing loss than on tinnitus there.
**That observation is recorded now so it cannot later be presented as a discovery.**

## P5 — the expected outcome, stated in advance

At 2,859 cases this arm will be in the same position as the bilateral one: the median
gene carries about five rare alleles among the cases, and the power analysis will say
the study cannot detect it. **The expected result is another measured null.**

It is still worth running. "We tested the hypothesis of origin with a corrected pipeline
and it does not hold at this cohort size" is a statement this project cannot currently
make, and it should be able to.

Writing the expectation down now means a null will not be presented as a disappointment
and a hit will not be presented as a surprise.

## P6 — what would make this arm interesting

A result here is interesting only if it is **better** than the hearing-loss arms, since
66% of the cases are shared. If `ZNF175` or anything else performs similarly across all
three phenotypes, that is the shared cases speaking, not tinnitus.

If something is notably stronger in tinnitus, the same caution applies as to the
bilateral arm's strengthening variants: the comparison is between overlapping sets, and
a size-matched control would be needed before claiming specificity.
