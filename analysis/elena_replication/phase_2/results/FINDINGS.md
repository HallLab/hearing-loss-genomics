# Phase 2 — which variants count

**Status:** chr8 pilot complete; genome-wide pending · **Date:** 2026-09-30
**Scripts:** [`../scripts/01_pilot_chr8_masks.py`](../scripts/01_pilot_chr8_masks.py),
[`../scripts/02_plof_diagnosis.py`](../scripts/02_plof_diagnosis.py)

---

## Finding 0 — which masks the run used, established rather than assumed

`rarevariant_geneburden/` holds **five** mask directories — `masks`, `masks_underscore`,
`masks_underscore2`, `masks_noduplicates`, `masks_deduplicated` — with nothing documenting which is
current. Phase 1's lesson was that an artifact on disk need not be the one that ran, so the answer
was taken from the Nextflow staged symlinks rather than from directory names or dates:

```
work/*/*/pLOF.txt -> rarevariant_geneburden/masks_deduplicated/pLOF.txt
```

All four masks resolve there, and it is the only one of the five containing all four files; the
other four hold `ALL` and `pLOF_pDM` only. **`masks_deduplicated` is what ran.**

---

## Finding 1 — half the pLOF mask meets neither half of its own documented definition

### The definition, from the pipeline's own analysis plan

Found in `elena_publishes/Hall Lab_ Analysis Plan Draft - Rare Variant ExWAS.docx`:

> **pLOF (predicted loss-of-function)**
> Frameshift, stop-gained, start-lost, stop-lost, or splice-site variants (SpliceAI ≥ 0.2)
> — *"VEP Consequence + SpliceAI"*

This matters more than it looks. SpliceAI scores **every** variant, synonymous ones included, so a
synonymous variant with a predicted splice effect is pLOF under this rule and not under a
consequence-only rule. An earlier draft of this finding treated synonymous-in-pLOF as prima facie
wrong. With the documented definition in hand, that inference does not hold, and the definition is
tested below rather than assumed.

Comparison is against the group files the release publishes itself
(`Exome/group_file_annotations/`), which carry `pLoF`, `damaging_missense`, `other_missense` and
`synonymous` per variant per gene. chr8 pilot:

| | pipeline `pLOF` | release `pLoF` |
|---|---:|---:|
| variants | **36,014** | 9,807 |

Shared: 9,198. The pipeline carries 26,816 the release does not, and the release 609 the pipeline
does not.

### Three explanations, separated

The gap could be a universe difference, a transcript policy, or a mis-assignment. Each was
quantified rather than argued:

```
26,816  in pipeline pLOF, not in release pLoF
 3,228    absent from the release universe entirely        (12%)
23,588    present in the release, annotated as something else
             other_missense     19,155
             damaging_missense   4,057
             synonymous            376
```

**Universe difference: minor.** Only 12% are unknown to the release.

**Transcript policy: minor.** A variant can be loss-of-function on a non-canonical transcript while
the group file reports the canonical call, and an any-transcript rule would legitimately include it.
Checking every transcript in the release's own VEP output: of the 23,588, only **1,919** have a
loss-of-function consequence on any transcript. **21,669 have none, on any transcript.**

**The documented definition explains 13% of it.** Of the 21,669 with no loss-of-function
consequence on any transcript, **2,824** carry SpliceAI ≥ 0.2 and therefore qualify under the plan's
splice-site clause. At stricter thresholds it is smaller still: 1,229 at 0.5, 547 at 0.8.

```
18,845  meet NEITHER criterion -- no LoF consequence on any transcript, and SpliceAI < 0.2
            other_missense     18,419
            synonymous            317
            damaging_missense     109
```

**That is 52.3% of the pipeline's chr8 pLOF mask.** The claim is not that these variants fail some
external standard — it is that they fail the pipeline's own stated rule, evaluated against the
release's VEP and SpliceAI output.

### One mechanism ruled out

A mask file pairs a `var` row with an `anno` row per gene, and a zip over lists of unequal length
would scramble the assignment. Checked across all four masks and all genes: **zero misaligned**, and
each file's annotations are homogeneous (`pLOF.txt` carries 1,002,120 annotations, all `pLOF`). The
mask files are structurally sound. Whatever happened, happened before they were written.


### Credit where it is due

The documented definition came from the analysis-plan documents collected in
`elena_replication/elena_publishes/`. Without it, this finding would have been written as
"synonymous variants cannot be loss-of-function, therefore the mask is broken" — a conclusion that
ignores a deliberate, documented design choice and would have been wrong in its reasoning even where
it happened to be right in its conclusion. The definition narrowed the claim from *most of the mask
is not LoF* to *half the mask fails the pipeline's own rule*, which is both smaller and far more
defensible.

---

## Finding 1b — the cause, established from the pipeline's own annotation

**Script:** [`../scripts/03_plof_root_cause.py`](../scripts/03_plof_root_cause.py)
**Output:** [`03_plof_root_cause.json`](03_plof_root_cause.json)

Finding 1 left one alternative open: that the pipeline's own VEP disagrees with the release's, making
the whole comparison a difference between two annotation runs rather than a defect. **It is a defect**,
and this is settled without reference to the release at all.

The mask builder reads `rarevariantExWAS/variant_categories/chr{N}.classified.tsv`, which carries the
pipeline's own `Consequence`, its own `SpliceAI_max`, and its own `is_pLOF` decision in the same row.
The question is therefore internal: does `is_pLOF` follow the documented rule?

### It does not, and the signature is exact

Of 197,002 annotation rows marked `is_pLOF = True` on chr8:

| | rows |
|---|---:|
| carry a genuine loss-of-function consequence | 50,252 — 25.5% |
| carry **no** LoF consequence, but contain the string `splice` | 146,750 — 74.5% |
| carry no LoF consequence and no `splice` | **0** |

Zero exceptions. Nothing is marked pLOF except through a real LoF consequence or a consequence
containing the word *splice*.

### The SpliceAI gate was never applied

The documented rule admits splice-site variants only at **SpliceAI ≥ 0.2**. Among splice-annotated
rows without a real LoF consequence:

| | rows |
|---|---:|
| `is_pLOF = True` with SpliceAI **< 0.2** | 129,118 |
| `is_pLOF = True` with SpliceAI ≥ 0.2 | 17,632 |
| `is_pLOF = **False**` although SpliceAI ≥ 0.2 | 9,094 |

SpliceAI is **uncorrelated** with the decision — it admits variants below the threshold and rejects
variants above it. The score was computed and written to the file; it was never used to gate.

### Three low-impact terms trigger it unconditionally

| VEP term | `is_pLOF` true | false | VEP impact |
|---|---:|---:|---|
| `splice_polypyrimidine_tract_variant` | 122,036 | **0** | LOW |
| `splice_donor_region_variant` | 17,580 | **0** | LOW |
| `splice_donor_5th_base_variant` | 7,134 | **0** | LOW |
| `splice_region_variant` | 46,693 | 54,531 | LOW — true only by co-occurrence |

The first three always trigger pLOF and never fail to. **None of them is a loss-of-function
consequence**: all are VEP `IMPACT=LOW` annotations for variants *near* a splice site rather than at
the donor or acceptor itself. `splice_polypyrimidine_tract_variant` alone accounts for 122,036 rows —
more than twice the entire genuine-LoF set.

### What this means

The defect is a **consequence-matching error in the classification step**, not an annotation
disagreement and not a mask-assembly bug. The pipeline's own VEP output is fine; the rule applied to
it is wrong in two independent ways — three low-impact terms admitted that should not be, and the
SpliceAI threshold that was supposed to gate them never consulted.

This also closes the alternative that Finding 1 flagged. The release was used to *detect* the problem;
the pipeline's own files *confirm* it. No appeal to an external reference is needed.

### Still not established

Whether chr8 is representative. Everything here is one chromosome. The classification code
(`step2_3_3_classifyvariants.bsub`) has also not been read — the behaviour is established from its
output, not from its source, so the exact expression that produces it is inferred rather than seen.

---

## Finding 1c — chr8 is representative; the genome-wide figure is 63.8%

**Script:** [`../scripts/04_genomewide.sh`](../scripts/04_genomewide.sh)
**Output:** [`04_genomewide.tsv`](04_genomewide.tsv)

Counted at the **variant** level across all 22 chromosomes — a variant is in the pLOF mask if any of
its transcript rows sets `is_pLOF`, and has genuine support if any row carries a real
loss-of-function consequence:

```
1,089,876   variants in the pLOF mask
  395,074   with a genuine LoF consequence          36.2%
  694,802   with none                               63.8%
```

Per chromosome the range is **60.9% to 65.7%** — every one of the 22 falls inside a five-point band.
chr8 sits at 64.5%, almost exactly the genome-wide figure. The pilot was representative, and the
defect is uniform rather than concentrated anywhere.

### A correction to how Finding 1b was stated

Finding 1b reported 74.5% for chr8. That is the **annotation-row** figure — 146,750 of 197,002 rows.
The mask contains **variants**, not rows, and a variant carries one row per transcript. At the variant
level chr8 is 64.5%.

Both numbers are correct for what they measure, and the row figure is the right one for diagnosing
the rule. **The variant figure is the one that describes the mask**, and it is the one to quote:
63.8% genome-wide, 694,802 variants.


---

## Finding 2 — CORRECTED: the pDM divergence is a parsing bug, not a threshold choice

> **This finding previously said the opposite.** It claimed the pDM divergence was a
> predictor-threshold judgement of the kind the 2026-07-01 meeting left open, and that pDM was
> "built from flags this defect does not involve". Nikki Palmiero found otherwise on 2026-10-05.
> The original reading is kept below the correction, because how it went wrong is the useful part.

**Script:** [`../scripts/06_revel_and_biotype.py`](../scripts/06_revel_and_biotype.py)

VEP writes REVEL per transcript, as a comma-separated list:

```
0.131,.,0.131,0.131,0.131
```

The classification script reads it with `pd.to_numeric(errors="coerce")`, which turns any
multi-valued entry into `NaN`. Only variants whose REVEL happened to arrive as a single number
survive.

| chr8 | |
|---|---:|
| missense rows with a usable score, as parsed | 65,806 — **12.1%** |
| with the list read | 517,274 — **95.3%** |
| variants passing REVEL ≥ 0.5, as parsed | 5,391 |
| variants passing REVEL ≥ 0.5, correct | **20,618** |
| **undercount** | **15,227 — 73.9%** |

So `pDM` is not intact. It is missing roughly three quarters of its REVEL-qualifying variants.

### Why this replication missed it

Checks 01 and 02 compared `pDM` against the release's group files **from outside**, saw a divergence
of the size a different threshold would produce, and recorded it as a judgement call.

**From outside, a different threshold and broken parsing are indistinguishable.** Both produce
"fewer variants than the reference". Only reading the field separates them.

That is a lesson about the method rather than about this defect: an external comparison detects a
discrepancy, it does not diagnose a cause. Phase 2 had an internal check available — the pipeline's
own classification file carries the raw `REVEL_score` — and did not use it for `pDM`, having already
concluded the mask was fine.

### What survives from the original reading

The contrast with `pLOF` still holds, and still matters. Disagreeing about *how damaging* a missense
variant is remains a judgement call, and the 2026-07-01 REVEL 0.5-versus-0.6 question is still open.
What is no longer true is that `pDM`'s divergence **is** that disagreement. Most of it is a bug.

---

## Finding 3 — non-coding genes are in the masks

Also raised by Nikki. Measured at gene level across all 22 chromosomes:

| mask | genes | not protein-coding | entries affected |
|---|---:|---:|---:|
| `pLOF` | 19,038 | **1,101 — 5.8%** | 20,078 — 2.0% |
| `pLOF_pDM` | 19,057 | **1,101 — 5.8%** | 20,342 — 1.2% |
| `pDM` | 17,336 | 8 — 0.0% | 264 |

`TMC3-AS1` is in the mask and is not protein-coding; the rest follow the same pattern — `A1BG-AS1`,
`ABCA9-AS1`, `ACTA2-AS1`.

**Large in genes, small in variants, and both are true.** Each of the 1,101 carries a burden test
that cannot mean anything — a lncRNA has no protein to lose function of — and each consumes
multiple-testing correction. But only about 2% of mask entries sit on them, so the dilution is minor.
Quoting either number alone misleads.

`pDM` is barely touched because a missense call requires a protein by definition.

---

## What is established, and what is not

| Claim | Basis | Strength |
|---|---|---|
| `masks_deduplicated` is what ran | resolved Nextflow symlinks | **solid** |
| 52.3% of chr8 pLOF meets neither clause of the documented definition | release VEP + SpliceAI, all transcripts | **solid for chr8, against this reference** |
| Mask files are structurally well-formed | var/anno lengths, annotation homogeneity | **solid** |
| pDM divergence is a threshold difference | 16,244 of 16,304 are release `other_missense` | **solid** |
| The cause of Finding 1 | the pipeline's own classification output (Finding 1b) | **established** |
| pDM is unaffected | — | **withdrawn** — it undercounts REVEL by 73.9% (Finding 2) |
| Non-coding genes in the masks | gene-level biotype, all 22 chromosomes (Finding 3) | **established** |
| Whether chr8 is representative | all 22 chromosomes, 60.9-65.7% (Finding 1c) | **established — it is** |

**Not established, and deliberately not guessed:**

- ~~**The cause.**~~ **Established — see Finding 1b.** It is a consequence-matching error in the
  classification step: three VEP `IMPACT=LOW` splice terms admitted unconditionally, and the
  documented SpliceAI ≥ 0.2 gate never applied. The classification *source* has still not been read,
  so the exact expression is inferred from its output.
- ~~**Whether the two annotation runs agree.**~~ **No longer an open question.** Finding 1b settles
  the defect from the pipeline's own annotation, without appeal to the release. The release detected
  the problem; the pipeline's own files confirm it.
- ~~**Whether it is genome-wide.**~~ **Established — see Finding 1c.** All 22 chromosomes fall
  between 60.9% and 65.7%; genome-wide, 694,802 of 1,089,876 pLOF variants have no loss-of-function
  consequence.
- **The effect on results.** That is Phase 4's question. A mask carrying non-LoF variants dilutes a
  burden test toward the null rather than inventing signal, but quantifying that is not Phase 2's
  job and is not attempted here.

## Next

1. Scale the comparison genome-wide (check 03) — establish whether chr8 is representative.
2. Read the pipeline's own VEP output (`vep_annotation_filtered/`) and compare it against the
   release's for the 18,845 residual variants. This distinguishes "the two VEP runs disagree" from
   "the mask builder mis-assigned", and must come before any statement about cause.
3. Read `step2_3_3_classifyvariants.bsub` to see how the documented definition was implemented.
4. Only then, assess what it means for the results already produced.
