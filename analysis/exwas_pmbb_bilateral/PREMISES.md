# Premises

Every choice this analysis makes, where it came from, and what it moves. Nothing here is a default —
each line was argued somewhere, and the source is named so it can be challenged.

Two sources: **M** = the meeting of 2026-10-02 with Douglas Epstein and Nicole Palmiero.
**R** = a finding from `analysis/elena_replication/`, with its phase.

---

## P1 · Cases are bilateral sensorineural only · **M**

> "we typically exclude one-sided hearing loss, and we just focus on the bilateral sensorineural
> hearing loss" — D. Epstein, 2026-10-02

phecodeX encodes two orthogonal axes in the children of `SO_396`, and one ICD code carries one of
each:

| axis | codes |
|---|---|
| type | `.1` conductive · `.2` sensorineural · `.3` mixed · `.5` sudden idiopathic |
| laterality | `.8` bilateral · `.9` unilateral |

A case needs `.2` **and** `.8` **on the same date**, on two distinct dates. The same-date
requirement matters: without it a sensorineural-unilateral visit and a conductive-bilateral visit
combine into a "case" never diagnosed with bilateral sensorineural loss.

**Moves: cases 6,752 → 3,164.** A 53% cut. Of the 3,588 lost, all still have hearing loss — they are
excluded, never moved to controls.

This also resolves a discrepancy raised in the meeting. Doug counted roughly 4,000 from the PMBB
phenotype browser using the bilateral sensorineural ICD code; the master code table shows
`H90.3 Sensorineural hearing loss, bilateral` on 3,964 people. He and the pipeline were counting
different things, and neither was wrong about its own question.

## P2 · Controls stay as they were · **R, Phase 1**

No ear-family evidence at all, from `conditions_phecode_x` **or** the OMOP `observation` table.

Narrowing the cases must not widen the controls. The people P1 excludes still have hearing loss;
moving them into the control group would be worse than leaving them among the cases.

**Moves: nothing. Controls 50,755.**

## P3 · The cohort is framed on exome samples · **R, Phase 1**

The pipeline cut its cohort against an **imputed** LD-pruned `.fam`, dropping 431 analysable people —
40 of them cases — who have exome data but no imputed PCs. The analysis runs on exome; the frame
should be exome.

## P4 · Tinnitus exclusion reads the `observation` table · **R, Phase 1**

PMBB v4 relocated the standard tinnitus codes (388.3x, H93.1x) out of the phecode files into OMOP
`observation`. Reading only the phecode files leaves 556 people who should be excluded sitting in
the control group, invisible.

## P5 · Masks are the rebuilt ones · **R, Phase 2**

Four defects, four causes, all repaired:

| | |
|---|---|
| three IMPACT=LOW splice terms listed as pLOF | removed |
| SpliceAI gate never applied | enforced — it had failed on 63.8% of pLOF entries |
| REVEL parsed as one number, not a per-transcript list | fixed — a 73.9% undercount (N. Palmiero) |
| non-coding genes carrying burden tests | removed — 1,101 genes (N. Palmiero) |

Copied here as built. Phenotype-independent, so not recomputed.

## P6 · The `ALL` mask is not run · **M and R**

A burden test over every variant in a gene pools synonymous and intronic variants with the ones that
have a mechanism. Dropped on 2026-10-05 as a scope decision; the clinical lead reached the same
position independently three days earlier:

> "that's the confusing part to me, why you would ever use the all category... it's also going to
> give you a lot of noise... I wouldn't necessarily focus on the ones that don't have a likely causal
> effect" — D. Epstein, 2026-10-02

**Moves: 792 step-2 tasks → 594,** and removes the mask that held 15,951,289 entries against
1,340,936 for `pLOF_pDM`.

## P7 · PC counts are 5 / 4 / 3 · **R, Phase 3**

The pipeline used 5 / 9 / 10. The eigenvalue scree supports 5 / 4 / 3, confirmed by two independent
routes. Surplus components spend degrees of freedom and can correlate with the phenotype by chance.

The PCA itself is phenotype-independent — it ran over everyone with exome data per ancestry, not over
the analysis cohort — so it is copied here unchanged.

## P8 · `Batch` is categorical · **R, Phase 4**

Passed in `covarColList` without `qCovarColList`, SAIGE reads batch 1/2/3 as a quantity and assumes
batch 2 sits halfway between 1 and 3. In AFR that is interpolating across groups of 8,492 / 686 /
2,156 people. Batch is a label.

Worth noting for whoever picks this up: her single-variant pipeline *declares* `qcovar_cols =
'SEX,Batch'` with a comment saying batch is categorical, but the Nextflow script never passes it.
The intent was there; the wiring was not.

## P9 · Bonferroni and FDR, at both levels of aggregation · **R, Phase 5**

Three Bonferroni denominators are defensible and they differ ninefold, so the one being used has to
be named rather than implied. Benjamini-Hochberg is reported alongside, because it is the lenient
correction and where a weak real effect would appear first.

---

## What P1 costs in power, stated plainly

| cohort | cases | controls |
|---|---:|---:|
| combined | 3,164 | 50,755 |
| EUR | 2,547 | 37,603 |
| AFR | **490** | 10,049 |

**AFR at 490 cases is the number to worry about.** The replication, with 1,285 AFR cases, already
flagged 26 of its top 30 AFR hits as fragile — carried by a handful of alleles each. At 490 that
gets worse, not better. The AFR arm will be run, because not running it is its own kind of
reporting bias, but it should be read as a null with very little power rather than as a search.
