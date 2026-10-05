# The pLOF mask admits variants its own definition excludes

**Author:** Andre Rico · **Date:** 2026-10-05 · **Status:** for review · **revised**
**For:** Nikki Palmiero · **Artifact:** `rarevariantExWAS/step2_3_3_classifyvariants.bsub`
**Evidence:** `analysis/elena_replication/phase_2/` — scripts, outputs, and corrected masks

---

## Decision requested

1. **Where did the `lof_terms` list come from?** Three of its nine entries are VEP `IMPACT=LOW`
   consequences. If the list was taken from a published definition or another pipeline, I would
   rather record that than assume it was an oversight. → **Nikki, Elena**
2. **Should the masks be rebuilt and the gene-burden step re-run?** A corrected mask is already
   built and ready. → **Nikki, Molly, Doug**
3. **Is the REVEL threshold still open?** The 0.5-versus-0.6 question has been pending since
   2026-07-01. It matters more now that the parsing bug below is fixed, since far more variants will
   be near the cut. → **Nikki**

---

## First, what checks out

This is a defect in one expression, not in the annotation work around it.

- **The VEP run is fine.** Consequences match what the release's own VEP assigns.
- **The SpliceAI computation is fine.** `SpliceAI_max` is extracted correctly at line 125 of the
  classification script.
- **`ALL` is fine.** No damage filter applies to it.
- **The mask files are structurally sound.** `var` and `anno` rows pair correctly for every gene, and
  annotations are homogeneous. Nothing is scrambled.
- **The documented definition is correct.** The analysis plan says exactly the right thing. The
  implementation departs from it.

---

## The defect

### The documented rule

From `Hall Lab_ Analysis Plan Draft - Rare Variant ExWAS.docx`:

> **pLOF (predicted loss-of-function)**
> Frameshift, stop-gained, start-lost, stop-lost, or splice-site variants **(SpliceAI ≥ 0.2)**
> — *"VEP Consequence + SpliceAI"*

### What the code does

`rarevariantExWAS/step2_3_3_classifyvariants.bsub`, **lines 46–56**:

```python
lof_terms = [
    "frameshift_variant",
    "stop_gained",
    "start_lost",
    "stop_lost",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "splice_donor_5th_base_variant",         # line 53  -- VEP IMPACT=LOW
    "splice_donor_region_variant",           # line 54  -- VEP IMPACT=LOW
    "splice_polypyrimidine_tract_variant"    # line 55  -- VEP IMPACT=LOW
]
```

**Line 138:**

```python
df["is_pLOF"] = df["Consequence"].apply(
    lambda x: any(term in str(x) for term in lof_terms)
)
```

Two departures:

**1. Three `IMPACT=LOW` terms are in the list.** None is a loss-of-function consequence. All three
annotate variants *near* a splice site rather than at the donor or acceptor — which is precisely the
population SpliceAI exists to adjudicate.

**2. The SpliceAI gate is never applied.** `SpliceAI_max` is computed at **line 125** and
`is_pLOF` is computed at **line 138**. Thirteen lines apart, and line 138 does not reference it.

---

## What it costs — measured

### chr8, at the annotation-row level, from the pipeline's own classification file

| | rows |
|---|---:|
| `is_pLOF = True` | 197,002 |
| …carrying a genuine LoF consequence | 50,252 — **25.5%** |
| …carrying none, but containing the string `splice` | 146,750 — **74.5%** |
| …carrying neither | **0** |

Zero exceptions: nothing enters the mask by any other route.

### SpliceAI is uncorrelated with the decision

Among splice-annotated rows without a real LoF consequence:

| | rows |
|---|---:|
| `is_pLOF = True` with SpliceAI **< 0.2** | 129,118 |
| `is_pLOF = **False**` with SpliceAI **≥ 0.2** | 9,094 |

It admits below the threshold and rejects above it.

### The three terms trigger unconditionally

| VEP term | triggers pLOF | never does | impact |
|---|---:|---:|---|
| `splice_polypyrimidine_tract_variant` | 122,036 | **0** | LOW |
| `splice_donor_region_variant` | 17,580 | **0** | LOW |
| `splice_donor_5th_base_variant` | 7,134 | **0** | LOW |
| `splice_region_variant` | 46,693 | 54,531 | LOW — only by co-occurrence |

`splice_polypyrimidine_tract_variant` alone accounts for more than **twice** the entire genuine-LoF
set.

### Genome-wide, at the variant level

| | variants |
|---|---:|
| in the pLOF mask | 1,089,876 |
| with a genuine LoF consequence | 395,074 — 36.2% |
| with none | **694,802 — 63.8%** |

Per chromosome the range is **60.9% to 65.7%** — all 22 inside a five-point band. The defect is
uniform, not concentrated.

### Theirs against ours

| | pipeline | corrected |
|---|---:|---:|
| `pLOF` entries | 1,002,120 | **409,179** |
| `pLOF_pDM` entries | 1,717,883 | 1,126,909 |
| `pDM` entries | 720,453 | 721,264 *(+811, see below)* |
| `ALL` entries | 21,415,507 | unchanged |

**The corrected pLOF mask is 41% the size of the one that ran.** Of what it keeps, 395,074 have a
genuine LoF consequence and roughly 14,000 enter through the SpliceAI gate — so that gate does admit
a real population. It was simply never consulted.

---

## Two further defects, found by Nikki

Raised 2026-10-05, after an independent review that reached the same pLOF diagnosis. Both were
confirmed here; neither had been found by this replication.

### REVEL is parsed as a single number

VEP writes it per transcript as a comma-separated list — `0.131,.,0.131,0.131,0.131` — and
`pd.to_numeric(errors="coerce")` turns any such entry into `NaN`. Only single-valued scores survive:
**12.1% of missense rows instead of 95.3%**. On chr8 the `pDM` set is **5,391 where it should be
20,618 — a 73.9% undercount**.

This replication had written that `pDM` was unaffected. **That was wrong**, and it is worth saying
why it was missed: checks 01–02 compared `pDM` against the release group files from outside, and from
outside a different threshold and broken parsing look identical — both give "fewer variants than the
reference". Only reading the field separates them.

### Non-coding genes are in the masks

Measured at gene level across all 22 chromosomes: **1,101 of 19,038 genes in `pLOF` are not
protein-coding**, carrying 20,078 mask entries. `TMC3-AS1` is among them, as are `A1BG-AS1`,
`ABCA9-AS1`, `ACTA2-AS1`.

Large in genes, small in variants — about 2% of entries — and both are true. Each of the 1,101 is a
burden test that cannot mean anything, and each consumes multiple-testing correction. `pDM` is barely
touched here (8 genes), because a missense call requires a protein.

### So all three damage masks are affected, by three independent causes

| defect | affects | size |
|---|---|---|
| low-impact splice terms + SpliceAI gate unapplied | `pLOF`, `pLOF_pDM` | 63.8% of variants |
| non-coding genes | `pLOF`, `pLOF_pDM` | 1,101 genes · 2% of entries |
| REVEL parsing | **`pDM`**, `pLOF_pDM` | 73.9% undercount |

---

## Why this is a defect and not a defensible disagreement

Worth stating, because two analysts can annotate the same variants and disagree legitimately.

**The confirmation is internal.** The classification file carries the pipeline's own `Consequence`,
its own `SpliceAI_max`, and its own `is_pLOF` in the same row. No external reference is needed: the
pipeline contradicts its own documented rule.

**The rule reproduces exactly.** Reimplementing the expression from lines 46–56 and 138 and comparing
against the file's `is_pLOF` column: **6,093,288 of 6,093,288 rows agree, zero divergences.** The
logic above is the logic that ran.

**The release detected it, but does not carry the claim.** The comparison began against
`Exome/group_file_annotations/`, which PMBB publishes with the same categories. That is how the gap
surfaced. The claim rests on the pipeline's own files.

---

## What is *not* being claimed

- **Not that results are invalidated.** The direction of the effect is dilution. Filling a gene's
  mask with inert variants makes a real signal harder to detect, not easier — a gene with 5
  protein-breaking variants and 15 inert ones is tested over all 20. Anything that reached
  significance under this mask would reach it more easily under the corrected one.
- **Not that false genes entered the analysis.** The defect adds variants to genes that were already
  there; it does not add genes. The risk is a **lost** finding, not a fabricated one.
- **Not that this was careless.** The documented definition is right, the SpliceAI score is computed
  correctly, and the three terms are explicitly named rather than caught by accident — which is why
  question 1 above asks where the list came from rather than assuming.
- **Not that every divergence is a defect.** Disagreeing about how damaging a missense variant is
  depends on predictor and cutoff, and the meeting left that open. That genuine judgement call still
  exists — what is no longer true is that `pDM`'s divergence *is* it. Most of that gap is the parsing
  bug above.

- **Not that the numbers here are final.** The corrected masks predate Nikki's two findings and are
  being rebuilt with the REVEL fix and a `BIOTYPE == protein_coding` restriction. The pLOF figures
  stand; the mask sizes will move.

---

## The corrected definition

```
pLOF =   an explicit list: frameshift / stop_gained / start_lost / stop_lost
         / splice_acceptor / splice_donor / transcript_ablation
  OR     any other splice annotation, but only with SpliceAI >= 0.2
  AND    BIOTYPE == protein_coding
```

Matched as exact terms against the comma-separated consequence list — **not** by VEP's `IMPACT`
field, which is per row: a row reading `stop_gained,splice_polypyrimidine_tract_variant` is `HIGH` as
a whole, so filtering on `IMPACT` would pull the low-impact splice terms straight back in.

`transcript_ablation` is added although the original list omits it: it is `IMPACT=HIGH` and
unambiguously loss of function, so leaving it out would be a second departure from the documented
definition rather than a fix.

---

## How to check this yourself

```bash
cd /project/hall/analysis/hearing-loss-genomics/analysis/elena_replication/phase_2

# the defect, from the pipeline's own classification file
../../../venv/bin/python3 scripts/03_plof_root_cause.py

# all 22 chromosomes
bash scripts/04_genomewide.sh

# build the corrected masks
../../../venv/bin/python3 scripts/05_emit_masks.py
```

Corrected masks: `phase_2/results/masks/`. Full write-up including what remains unestablished:
`phase_2/results/FINDINGS.md`.

---

## Context

This came out of an independent replication of the PMBB v4 rare-variant pipeline, re-derived phase by
phase from the institutional release before further work is built on its outputs. Phase 1 was the
cohort; this is Phase 2, the variant masks.

The replication is not an audit of any individual's work — it re-derives numbers before publishing
them, and it runs in both directions. Phase 1 found the replication itself wrong about the v4
tinnitus source and the pipeline right.
