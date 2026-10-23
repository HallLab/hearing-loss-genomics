# Phase 5 (BRaVa) — the pre-declared test, answered

**Run 2026-10-23.** Step 2 job `49775068`, 88/88. The test and its threshold were fixed
in [`PREMISES.md`](../../PREMISES.md) P5 **before anything ran**.

## The declared test: it survives

> ClinGen Definitive/Strong genes in the **top 50** of the broad combined ranking,
> hypergeometric p. `p < 0.05` survives; `p >= 0.05` does not, and must be reported as
> prominently as the original finding.

| | observed | expected | p |
|---|---:|---:|---:|
| `pmbb_hl` annotation | 5 of 93 | 0.26 | **5.91e-06** |
| **BRaVa annotation** | **5 of 92** | **0.27** | **7.83e-06** |

**Essentially unchanged**, on a mask set that shares a fifth of its variants. The EUR arm
behaves the same way: 2 in the top 50 under both, p = 2.78e-02 ours against 3.04e-02
BRaVa.

This is the fourth thing that has failed to move this result. It has now survived
changing the PCA source, raising the PC counts from 5/4/3 to 11/15/16, adding
deduplication before the combination, and — now — an annotation that disagrees with ours
about **four fifths** of the variants in `pDM`.

## And the genes swap, which is the more interesting half

| rank | `pmbb_hl` | | rank | BRaVa |
|---:|---|---|---:|---|
| 4 | `SIX1` | | 4 | `COCH` |
| 16 | `GJB3` | | 13 | `SIX1` |
| 23 | `COCH` | | 32 | `RIPOR2` |
| 30 | `TMPRSS3` | | 39 | `CEACAM16` |
| 35 | `MYO6` | | 43 | `GJB3` |

**Three genes in common** — `SIX1`, `GJB3`, `COCH`. `MYO6` and `TMPRSS3` drop out;
`RIPOR2` and `CEACAM16` come in.

This is not a problem with the result. It is the clearest evidence yet for what the
result has always claimed: **it is a statement about the ranking, not about any gene.**

If any individual gene were carrying a detectable effect, it would hold its position
under a different annotation. None does. What holds is the *class* — known deafness
genes rise to the top at five times the rate chance allows, and which five arrive there
depends on which variants the mask happens to include.

A reader who took `COCH` and `MYO6` from the first result would have taken two genes,
one of which does not reproduce. **The ranking reproduces; the names do not.**

## The synonymous control is flat

BRaVa category 4 — synonymous variants with `SpliceAI < 0.2`, 16,781 genes, 1,593,025
variants. `pmbb_hl` has no equivalent. [`PREMISES.md`](../../PREMISES.md) P4 recorded in
advance that this is expected to show nothing.

**It shows nothing — and it shows nothing the way a real null does.** Across twelve
depth × tier × cohort combinations, ten find exactly zero ClinGen genes. The other two
are the EUR arm at the looser Moderate tier, which finds **one** gene at depth 100
(0.65 expected, p = 0.48) and **one** at depth 250 (1.63 expected, p = 0.81).

That is the right answer, and a cleaner one than zero everywhere would have been. A
control that returned a perfect zero in all twelve would itself be slightly odd —
chance alone puts about 0.65 genes in the top 100. Finding roughly the expected count,
with no p below 0.48, is what an uninformative ranking looks like.

**At the pre-declared tier and depth — Definitive/Strong, top 50 — the control is zero in
both cohorts.**

That is worth more than it reads. A synonymous mask carries no predicted functional
consequence, so it should not rank deafness genes. If it had, the problem would have been
**structural** — something separating cases from controls that is not biology, such as
batch, coverage or relatedness leaking through the null model — and no other check in
this project would have caught it. The enrichment would have looked exactly the same and
meant nothing.

The control was run in the same job array as the burden grid, so it could not be
forgotten and quietly assumed to have passed.

## Nothing is significant, in any of the six arms

| arm | genes | min omnibus p | best gene | min q | hits |
|---|---:|---:|---|---:|---:|
| burden, broad combined | 16,742 | 2.33e-05 | `PSG8` | 0.390 | 0 |
| burden, broad EUR | 16,697 | 1.77e-04 | `PSG8` | 0.914 | 0 |
| burden, broad AFR | 16,414 | 2.16e-05 | `RNF149` | 0.202 | 0 |
| synonymous, broad combined | 16,776 | 1.32e-04 | `DBT` | 0.896 | 0 |
| synonymous, broad EUR | 16,762 | 9.24e-05 | `RFX5` | 0.854 | 0 |
| synonymous, broad AFR | 16,699 | 2.52e-05 | `BBOX1` | 0.324 | 0 |

Bars as in `pmbb_hl`: `0.05/genes` ≈ 2.99e-06, plus FDR q < 0.10. **Zero hits under
either, anywhere.**

The BRaVa burden arm lands further from the bar than ours did (2.33e-05 against our
2.88e-05 in combined is comparable; EUR is markedly worse at 1.77e-04 against 3.31e-05).
That is expected and is not a verdict: BRaVa's masks are smaller, so each gene carries
fewer alleles and the arm has less to test.

## AFR, added after the fact and labelled as such

The AFR cohort was **out of scope** in `PREMISES.md` P3 — it carries 0 ClinGen genes in
the top 50 under our own annotation, so there was nothing for a robustness test to be
robust about. It was run on 2026-10-23, **after the declared test was answered**, so the
annotation comparison figure could cover all three cohorts. The
[addendum](../../PREMISES.md) records that.

**It finds no enrichment either:** 1 ClinGen gene in the top 50 against
0.28 expected, p = 0.25. Consistent with our own annotation,
which finds 0. Nothing here is a result, and nothing here touches the declared test —
which is still p = 7.83e-06 on broad combined, unchanged.

Its synonymous control behaves like the others: 0 at the declared depth and tier, 1 at
depth 100 against 0.55 expected (p = 0.43).

## What this does and does not establish

**It does not say BRaVa is better, or that ours is.** Neither annotation is validated
against the other here. Ours uses AlphaMissense, which postdates the BRaVa document;
theirs uses CADD and LOFTEE, which ours does not. `PREMISES.md` P6 fixed this in advance.

**It does not rescue a null.** Nothing was significant under our annotation and nothing
is under BRaVa's.

**What it does establish** is the thing the replication plan needs: the project's one
positive finding is not an artefact of our annotation choices. A BRaVa-annotated
replication in All of Us or UK Biobank is testing the same claim we are, and if it
disagrees, the disagreement will be about biology or power rather than about variant
definitions.

## Files

| | |
|---|---|
| `results/01_enrichment.json` | the declared test, flagged, plus every secondary depth and tier |
| `../phase_4/results/omnibus_*.tsv` | one row per gene, four arms |
| `../phase_4/results/03_merge.json` | counts, duplicates dropped, bars |
| `../phase_2/results/03_mask_comparison.json` | how far the two annotations are apart |
