# Pipeline plan — what the v4 pipeline does, and what this replication re-does

**Status:** reference document · **Author:** Andre Rico · **Written:** 2026-09-29
**Companion to:** [`README.md`](README.md), which states the scope and the defects. This file
explains the pipeline those defects sit in.

---

## 0. Why this document

The README names three defects by file path and line count. That is enough to act on but not enough
to understand: it does not say what the pipeline was trying to do, which step produced which
artifact, or why a defect in one step matters more than a defect in another.

This file supplies that. It is written to be read by someone who has not seen the pipeline before.

---

## 1. The question the pipeline answers

> **Is there any gene in the genome where rare, damaging variants appear more often in people with
> hearing loss than in people without?**

Two choices are built into that sentence, and both drive the pipeline's shape.

**Why rare variants.** Common variants have been swept exhaustively by association studies. Rare
ones are hard precisely because each is carried by very few people — no single rare variant has the
statistical power to show anything. The way around this is to **group by gene**: instead of testing
one variant at a time, ask "who carries *any* damaging variant in this gene?" That is a *burden*
test, and SAIGE-GENE+ is the implementation used here.

**Why exome-wide rather than a candidate list.** Decided 2026-07-01: scan all ~20,000 genes rather
than restricting to known hearing genes, using known genes only as positive controls. The cost is a
much harsher significance threshold; the benefit is that a novel gene can be found at all.

---

## 2. The pipeline — five phases, eleven numbered steps

Steps live in `analysis/elena/rarevariantExWAS/`. Phase 4 is where the defects are.

```
PHASE 1   Who is in the study
          step1   cases and controls from diagnosis codes (rule of 2 dates),
                  related individuals removed
          step5   check each ancestry group has enough cases and controls

PHASE 2   Which variants count
          step2   annotate every variant with VEP — what it does to the protein
          step3   variant-level quality control
          step4   group variants per gene into "masks"

PHASE 3   What has to be adjusted away
          step6   build covariate files: age, sex, batch, ancestry PCs
          step7   LD pruning — pick independent common variants so relatedness
                  between participants can be estimated

PHASE 4   The statistical test                      <-- defects 1 and 2 live here
          step8   SAIGE step 1: learn how hearing loss is distributed across the
                  cohort WITHOUT looking at any gene, already accounting for
                  relatedness and covariates
          step9   SAIGE step 2: test each gene against that null model
                  ---> this is what produces the 792 output files

PHASE 5   Reading the result
          step9_1  merge chromosomes        step10  QQ plots (calibration check)
          step9_2  map genes                step11  final tables
          step9_3-6  Manhattan plots
```

### The two directories, and how they chain

`rarevariantExWAS/` is the backbone. `HL_only_rarevariant/` is a re-run of Phases 4 and 5 against a
hearing-loss-only phenotype, and it **consumes the backbone's output directly** — verified in
[`step2_newcovariate.bsub`](../elena/HL_only_rarevariant/step2_newcovariate.bsub):

```
rarevariantExWAS/PMBBv4_phecodex/                        <- phenotype (Phase 1)
rarevariantExWAS/covariates/covariates_combined_5PCs_withBatch.txt
rarevariantExWAS/covariates/covariates_EUR_9PCs_withBatch.txt   <- covariates (Phase 3)
rarevariantExWAS/covariates/covariates_AFR_10PCs_withBatch.txt
```

This is why the replication cannot start at `HL_only_rarevariant/`: its inputs are produced
upstream, so an error in Phase 1 or Phase 3 propagates into it invisibly.

---

## 3. The masks — where biology enters

Step 4 is the only place in the pipeline where a judgement about *biological consequence* is made.
It labels each variant with what it does to the protein:

| Label | Meaning |
|---|---|
| `pLOF` | predicted loss of function — the variant **breaks** the protein (premature stop, frameshift) |
| `pDM` | predicted damaging missense — an amino-acid substitution scored as harmful (AlphaMissense / REVEL) |

SAIGE then forms four test groups from those two labels, via
`--annotation_in_groupTest="ALL,pLOF,pDM,pLOF:pDM"`:

| Test group | Variants included |
|---|---|
| `pLOF` | only protein-breaking |
| `pDM` | only damaging missense |
| `pLOF_pDM` | both together |
| `ALL` | every rare variant, no damage filter |

Keeping them separate matters: a gene can act through one mechanism and not the other, and a hit
that lands on the biologically expected mask is stronger evidence than one that does not.

### The release ships its own version of this

The institutional release does not only publish genotypes. It publishes the output of Phase 2, in
SAIGE's own group-file format:

```
/static/PMBB/PMBB-Release-2026-4.0/Exome/
  vep_annotations/         VEP annotations, 24 chromosomes
  group_file_annotations/  SAIGE group files, 24 chromosomes
  freqcounts/              allele frequencies per ancestry per chromosome (672 files)
```

A group file carries, per gene, a variant list and a parallel list of functional calls. Counting
those calls on chr8:

| Category | Variants |
|---|---:|
| `other_missense` | 142,037 |
| `synonymous` | 64,698 |
| `damaging_missense` | 23,139 |
| `pLoF` | 9,810 |

`pLoF` and `damaging_missense` map directly onto the `pLOF` and `pDM` labels that step 4 builds by
hand. Phase 2's output can therefore be checked against an independent authority at the cost of a
file comparison — see §7.

---

## 4. Why there are 792 outputs, not one

Phase 4 does not produce "a result". It produces one file per slice of the data:

```
3 ancestry groups  x  4 masks  x  22 chromosomes  x  3 MAF thresholds  =  792
```

Files actually on disk: **781**. The difference is exactly defect 2 (§5). The arithmetic closing to
the file is itself the cleanest confirmation that eleven results were lost rather than never
intended.

---

## 5. Where the three defects fall

```
PHASE 3  step6   ->  defect 3: PC counts 5 / 9 / 10; the recorded decision said 5-6
PHASE 4  step9   ->  defect 1: 3 files written with no header row (AFR / pDM / chr8)
PHASE 4  step9   ->  defect 2: 11 files never written (killed for memory, no retry)
PHASE 5  merge   ->  the summary table inherits both, and reports neither
```

Severity follows the phase. A Phase 5 defect corrupts how results are *read*. A Phase 4 defect means
results were never *computed*. A Phase 3 defect means everything downstream was computed against the
wrong adjustment — the most expensive kind to discover late, and the reason the replication checks
Phase 3 before it spends anything on Phase 4.

---

## 6. What the replication does to each phase

The replication walks the pipeline in its own order, 1 through 5. It does not invent a separate
sequence.

```
PHASE 1  rebuild cases/controls and the sample list      hours
PHASE 2  diff masks against the release's group files    minutes, no annotation compute
PHASE 3  rebuild ancestry PCs and covariates             hours
PHASE 4  re-run SAIGE on a sample of strata              days, ~1 TB
PHASE 5  rebuild the summary table                       minutes
```

**Why dependency order rather than cheapest-first.** An earlier draft ordered this work by cost, so
the cheap checks ran first. That was wrong for what the replication is for. Each phase consumes the
previous phase's output, so a check that passes out of order can be meaningless: confirming that the
summary table faithfully reflects the per-stratum outputs says nothing if those outputs were
computed on the wrong cohort. Worse, a disagreement found out of order cannot be attributed — a
difference visible at Phase 4 could originate at Phase 1, and there is no way to tell which from
Phase 4 alone.

Running in dependency order costs more up front and answers a better question. Every phase is
checked against a foundation already verified, so a disagreement belongs to the phase where it
surfaces.

**What that costs.** Phase 4 is days and ~1 TB, and it comes before the minutes-long Phase 5. There
is no quick partial answer early. That is the price of an attributable result rather than a fast one.

**One task sits outside this.** The summary table on disk is missing ~5,400 rows (§5), so what the
existing results say cannot presently be read. Re-parsing the raw outputs to recover them takes
minutes. That is triage — it makes existing output legible — and it is not evidence about whether
the output is correct. It is deliberately not a phase.

---

## 7. Phase 2 — how it is checked, and why this changed

**Earlier decision (2026-09-29): Phase 2 was out of scope.** The stated reasoning was that
re-running VEP across the exome is expensive and no known defect pointed at it. That is recorded
here rather than deleted, because the reversal is the informative part.

**Revised the same day.** Two things changed it.

First, Phase 2 is the highest-leverage unverified phase. Each of the three recorded defects affects
one stratum, one mask or one adjustment. An annotation difference moves *every gene at once*,
because annotation defines the masks that all four test groups are built from. This project also
has history here: a carrier-count discrepancy investigated in an earlier cycle was attributed in
part to a VEP-versus-ANNOVAR difference.

Second, and decisively, checking Phase 2 turned out not to require running an annotator at all. The
release publishes its own group files (§3), produced independently of this project. Comparing step
4's masks against them is a file diff.

### What the Phase 2 check asks

- Do the step-4 masks and the release group files contain the same variants per gene?
- Where both make a functional call, do they agree?
- Where they disagree, is it systematic (a threshold, a transcript choice) or scattered?

### What the Phase 2 check cannot settle

Agreement is not independent validation. Both the step-4 masks and the release group files derive
from VEP, so they can share an upstream assumption and agree while both being wrong in the same
direction. Two runs of one method are not two methods.

### Where Biofilter 4 fits

BF4 is structurally different: it looks a variant up in a knowledge base rather than computing its
consequence de novo. A third opinion that fails differently is worth more than a fourth that fails
the same way.

BF4 4.3.0 (module `biofilter/4.3.0`, bundle `20260914`) is a live candidate. Its bundle is built
from the **joint** gnomAD source (`dtp_variant_gnomad_joint`) and ships an `variant_alphamissense`
table — between them these address the load configuration that made an earlier BF4 evaluation
unusable for rare variants.

It is **not scheduled**. It becomes worth running if the Phase 2 check finds disagreement, or
afterwards as a deliberate third opinion. Two conditions apply if it is run:

1. Its coverage benchmark must be a **rarity-stratified sample drawn from the release itself**
   (`freqcounts/` carries the per-ancestry frequencies) — not a variant list inherited from other
   work, which would break the isolation rule in [README §4](README.md) and would not generalise
   beyond whatever locus it came from.
2. A lookup-based annotator is structurally disadvantaged on private variants, which an exome study
   has in abundance. The comparison can establish *where both can see, do they agree* and *how much
   does BF4 miss* — it cannot establish that BF4 substitutes for VEP.

---

## 8. How the facts in this document were established

| Claim | Checked by |
|---|---|
| Step list and phase grouping | file listing of `rarevariantExWAS/`, notebook headers |
| `HL_only_rarevariant` consumes backbone phenotype + covariates | paths in `step2_newcovariate.bsub` |
| Mask labels and test groups | `step4_burdenmasking.bsub`; `--annotation_in_groupTest` in the SAIGE command line |
| 792 expected vs 781 on disk | `find` count under `saige_results_newmasks/`, and 3x4x22x3 |
| Release ships VEP annotations, group files and freqcounts | directory listing of `/static/PMBB/PMBB-Release-2026-4.0/Exome/` |
| Group-file format and chr8 category counts | `awk` over `group_file_annotations/...chr8.txt` |
| BF4 4.3.0 bundle uses joint gnomAD | `bundle_plan.json` in `/project/hall_shared/datasets/biofilter/20260914/`; module `biofilter/4.3.0` |
| Defects 1-3 | see README §3, which records paths and counts |
