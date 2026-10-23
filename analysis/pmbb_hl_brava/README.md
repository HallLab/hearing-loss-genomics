# `pmbb_hl_brava` — the same genotypes, annotated to the BRaVa standard

**Hall Lab × Epstein Lab · PMBB-Release-2026-4.0 · October 2026**

A second annotation of the `pmbb_hl` cohort, built to
[`docs/papers/BRaVa annotation recommendations.pdf`](../../docs/papers/).
**Not a correction** of `pmbb_hl` — that arm's masks are defensible and reproduce the
first pass byte for byte. This exists for two reasons, set out in
[`PREMISES.md`](PREMISES.md):

1. **Comparability.** The replication plan targets All of Us and UK Biobank, which
   annotate to BRaVa. Results annotated differently cannot be placed side by side.
2. **A robustness test on the project's only positive result.** The ClinGen enrichment
   has survived changing the PCA source, the PC counts and the deduplication. It has
   never been tested against a different annotation.

## Status

| phase | state |
|---|---|
| 2 — masks | **done**, 22/22, job `49775051` |
| 2 — comparison against `pmbb_hl` | **done** |
| 4 — SAIGE step 2 | **done**, 88/88, job `49775068`. Step 1 reused from `pmbb_hl` — the null model does not depend on the mask |
| 5 — the pre-declared enrichment test | **done** — [`phase_5/results/FINDINGS.md`](phase_5/results/FINDINGS.md) |

## The answer

**The ClinGen enrichment survives.** 5 Definitive/Strong genes in the top 50 of the broad
combined ranking, 0.27 expected, **p = 7.83e-06** — against 5.91e-06 under our own
annotation, on a mask set that shares a fifth of its variants.

**But three of the five genes are the same, not five.** `SIX1`, `GJB3` and `COCH` hold;
`MYO6` and `TMPRSS3` drop out and `RIPOR2` and `CEACAM16` take their places. That is the
clearest evidence yet for what the finding has always claimed — **it is about the
ranking, not about any gene.** A gene carrying a detectable effect would hold its place
under a different annotation. None does; the class does.

**The synonymous control is flat**, which rules out a structural artefact no other check
in this project would have caught.

**Nothing is significant in any of the four arms**, as expected.

Everything above step 2 is inherited unchanged from `pmbb_hl`: cohorts, covariates, PCs,
GRM, and the nine fitted null models. **Only the mask differs**, which is what makes any
downstream difference attributable.

## What the two annotations share

Measured, not asserted — `phase_2/results/03_mask_comparison.json`.

| mask | genes ours | genes BRaVa | Jaccard, all variants | median per-gene Jaccard |
|---|---:|---:|---:|---:|
| `pLOF` | 17,841 | 15,726 | 0.517 | 0.556 |
| `pDM` | 17,623 | 16,603 | **0.257** | **0.229** |
| `pLOF_pDM` | 17,945 | 16,755 | 0.438 | 0.458 |

**The per-gene figure is the honest one.** An exome-wide Jaccard pools a few huge genes
with 17,000 small ones; the burden test is per gene, so the per-gene distribution is
what the analysis experiences. For `pDM`, **10.3% of shared genes have no variant in
common at all**.

The five genes the enrichment rests on are present in both and keep about half their
variants:

| gene | ours | BRaVa | shared | Jaccard |
|---|---:|---:|---:|---:|
| `SIX1` | 44 | 33 | 31 | 0.674 |
| `GJB3` | 104 | 60 | 60 | 0.577 |
| `COCH` | 105 | 71 | 65 | 0.586 |
| `TMPRSS3` | 109 | 106 | 76 | 0.547 |
| `MYO6` | 242 | 178 | 147 | 0.538 |

So the pre-declared test in `PREMISES.md` P5 is answerable: the genes survive, their
contents do not.

## The four differences, and which way each cuts

| | `pmbb_hl` | BRaVa | effect |
|---|---|---|---|
| `pLOF` | seven consequence terms, plus splice rescued at `SpliceAI ≥ 0.2` | **LOFTEE HC** | BRaVa is stricter — ours admits LoF in last exons, rescued splice sites, single-exon genes |
| `pDM` | AlphaMissense P/LP, or REVEL ≥ 0.5 | REVEL ≥ 0.773, or **CADD ≥ 28.1**, or LOFTEE LC | different sets, not different thresholds. AlphaMissense postdates the BRaVa document |
| transcripts | **all rows**, max REVEL across them | **MANE Select**, else `CANONICAL` | ours is more permissive: a variant benign on MANE can enter via a minor transcript |
| frequency | internal MAF at test time | **gnomAD broad-ancestry AF ≤ 0.01** | a variant rare in PMBB and common in gnomAD passes ours, fails theirs |

## The synonymous mask

16,781 genes, 1,593,025 variants — BRaVa category 4, synonymous with `SpliceAI < 0.2`.
`pmbb_hl` has no equivalent.

It is a **negative control**, and `PREMISES.md` P4 records in advance that it is expected
to show nothing. If synonymous variants enrich at the top of the ranking, or λ sits far
from 1, the problem is structural — something about how cases and controls differ — and
no other check in this project would catch it.

## Layout

```
PREMISES.md                        written before anything ran
phase_2/scripts/01_build_masks.py  one chromosome per task, BRaVa categories
        scripts/01_build_masks.bsub    job name `bmasks`, not `masks`
        scripts/02_merge_masks.py  refuses on a partial mask; checks categories disjoint
        scripts/03_compare_masks.py    the distance from pmbb_hl, measured
        results/masks/             pLOF · pDM · pLOF_pDM · synonymous
        results/02_masks_manifest.json
        results/03_mask_comparison.json
```
