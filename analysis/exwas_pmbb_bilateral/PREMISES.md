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

### Tested after the fact, and it holds

This premise was adopted on clinical grounds alone. Phase 5 then found that known deafness genes
cluster at the top of the **broad** analysis (5 ClinGen Definitive/Strong genes in the top 50
against 0.26 expected, p = 5.9 × 10⁻⁶) and nowhere near the top of the restricted one — which looks
like evidence that this premise throws away real biology.

Two controls say it does not:

| | cases | ClinGen in top 50 |
|---|---:|---:|
| broad, all | 6,752 | 5 |
| broad subsampled to 3,164, five draws | 3,164 | 2 · 0 · 1 · 2 · 3 |
| **restricted (this premise)** | 3,164 | **0** |
| the cases this premise discards | 3,588 | 1 |

Neither half differs from a random half of the same size — P(0) = 0.20 and P(≤1) = 0.53 against the
draws — and the halves do not add up, 0 + 1 against 5. The enrichment needs all 6,752 together.

**So the restricted arm's null is lost power, not lost biology, and this premise is not shown to
discard signal.** What that does *not* say is that it is better; at 3,164 cases neither definition
can show anything, which is a fact about the cohort rather than about either phenotype.

Full write-up and both controls: [`phase_5/results/FINDINGS.md`](phase_5/results/FINDINGS.md)
Finding 3.

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

**Moves:** removes the mask that held 15,951,289 entries against 1,340,936 for `pLOF_pDM`. In the
replication's job structure that was 792 tasks down to 594; **P10** then takes it to 66.

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

## P10 · SAIGE runs the whole grid in one call · **R, measured 2026-10-07**

Not inherited — this analysis is the first to do it, and it came out of Nikki Palmiero asking why
the replication had no `Cauchy` rows.

SAIGE is built to take the grid whole: `--annotation_in_groupTest` takes a list, and
`--maxMAF_in_groupTest` **defaults** to exactly our `0.0001,0.001,0.01`. One call per cohort and
chromosome tests all nine combinations and emits **one Cauchy row per gene** covering all of them —
the per-gene omnibus, from the tool, with nothing computed on the side.

The replication overrode that default and sliced the grid into 594 jobs; the pipeline it replicates
sliced it into 792. Neither produced a per-gene omnibus as a result.

**Moves: 594 step-2 tasks → 66**, and the headline p-value stops depending on a calculation of mine
that runs after the fact.

### What it costs, measured rather than assumed

chr21, combined cohort, one call against the nine:

| comparison | genes | identical |
|---|---:|---|
| `pLOF:pDM` — same group file either way | 623 | **623 (100%)** |
| `pLOF` — from pLOF_pDM vs its own file | 620 | **620 (100%)** |
| `pDM` — from pLOF_pDM vs its own file | 522 | 381 (73%) |

Native Cauchy: one row per gene, matching the hand-computed ACAT exactly for 165 of 208 genes,
median ratio 1.000.

Only `pDM` moves, and the cause is label precedence: of 11,517 variants that are both, 11,516 are
written as `pLOF` in the combined file. So `pDM` read from there means "damaging missense that is
not also loss-of-function". These are variants that truncate on one transcript and are missense on
another, so either reading is defensible — and this one is arguably cleaner, since a variant that
truncates the protein is not a missense story.

Recorded so that a `pDM` result from this analysis is not compared naively against a `pDM` result
from the replication. They are answering slightly different questions.

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
