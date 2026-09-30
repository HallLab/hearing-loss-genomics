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

## 2. The pipeline, as its author documented it

Steps live in `analysis/elena/rarevariantExWAS/`. The numbering and sub-step titles below are taken
from the analysis blogs in [`elena_publishes/`](../elena_publishes/), not inferred from file names.
An earlier version of this section was inferred, and got Phase 4 wrong — see *Where the numbering
disagrees with itself*.

Our five phases are a grouping over her eleven steps. They are a reading aid, not a second
vocabulary: every claim in this replication is anchored to a step.

```
PHASE 1   Who is in the study
  Step 1  Phenotyping and Sample QC
            cases and controls from diagnosis codes (rule of 2)
            builds HL_TIN_PMBBv4_keep.txt, then the SAIGE sample list
  Step 5  Ancestry Stratification
            5.1  case/control counts per ancestry
            5.2  interpret the ancestry classification

PHASE 2   Which variants count
  Step 2  Variant Classification
            2.1  rewrite the raw PMBB v4 VEP annotation files
            2.2  keep a subset of annotation columns
            2.3  classify variants; compare REVEL against AlphaMissense
            2.4  count ZNF175 variants in cases vs controls
  Step 3  Variant QC
            3.1  restrict to rare variants meeting the thresholds
            3.2  create summary file
  Step 4  Burden Masking
            4.1  make the 4 masks
            4.2  convert ":" to "_" for SAIGE step 2 compatibility

PHASE 3   What has to be adjusted away
  Step 6  SAIGE Covariate File
            6.1  decide how many PCs for the combined dataset
            6.2  build the combined covariate file
            6.3  covariate QC — check for missing values
            6.4  decide PCs per ancestry (scree plots); build by-ancestry files
            6.5  add batch to all three cohorts
  Step 7  LD Pruning
            7.1  generate PLINK keep files for EUR and AFR
            7.2  LD prune all three datasets

PHASE 4   The statistical test
  Step 8  SAIGE
            8.1  gene burden      ---> the 792 outputs
            8.2  single variant, 4 MAF thresholds, no masks

PHASE 5   Reading the result
  Step 9  Manhattan Plots
            9.1  merge chromosomes across the 36 categories
            9.2  map genes to positions
            9.3-9.6  the plots
  Step 10 QQ Plots           10.1-10.4, calibration
  Step 11 Interesting Tables 11.1-11.4
```

### Where the numbering disagrees with itself

`step9` is overloaded on disk. `step9_saige_step2.bsub` runs the association test, while
`step9_1_mergechr.bsub` through `step9_6_manhattanMAC.R` produce Manhattan plots. The blogs place
SAIGE entirely in Step 8 and Manhattan plots in Step 9.

The file names and the documentation disagree, and this replication follows the documentation: **the
test is Step 8, Phase 4; everything from the chromosome merge onward is Step 9+, Phase 5.** An
earlier version of this page followed the file names, split Step 9 across two phases, and described
`step9` as "SAIGE step 2". Anyone tracing a finding back to a file should expect the mismatch.

### Decisions the documents record, which the files do not

Three choices are visible only in the prose or in code comments. Each would read as an omission from
the files alone:

| Decision | Where | Note |
|---|---|---|
| Related individuals are **not** removed | comment in `step1_phenotyping_sample_qc.ipynb`: *"Don't need to get rid of unrelated because SAIGE can account for relationships"* | deliberate and defensible — SAIGE models relatedness through the GRM |
| `pLOF` = frameshift, stop-gained, start-lost, stop-lost, **or SpliceAI ≥ 0.2** | analysis plan, "VEP Consequence + SpliceAI" | broader than a consequence-only rule; see §7 |
| PC count to be set by variance-explained plot, *"most likely 4 or 5"* | analysis plan | the runs used 5 combined, 9 EUR, 10 AFR |

### The sample list, and what Step 1 actually does

Step 1 does more than assign cases and controls. It also builds the list of people who enter the
analysis, and that is where the cohort is cut:

```python
fam_file = ".../Imputed/common_snps_LD_pruned/..._genetic_imputed.commonsnps.ldpruned.ALL.fam"
# Keep only samples with genotype data
matched = pheno[pheno["PMBB_ID"].isin(fam["IID"])].copy()
```

The intent in that comment is right — restrict to people with genotype data. The file consulted is
the **imputed** LD-pruned `.fam` (70,493 samples), while the analysis tests **exome** data. An exome
LD-pruned set exists and was built by the same pipeline (`exome_ldpruned/ALL.common_ldpruned`,
70,925 samples); the sample list was not gated on it. That single line is the mechanism behind the
517 (Phase 1, Finding 1).

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
PHASE 1  Step 1   ->  the sample list is cut against the imputed .fam; 517 dropped,
                      40 of them cases                            (Phase 1, Finding 1)
PHASE 1  Step 1   ->  the phenotype file was overwritten after the covariates were
                      built from it; 556 excluded people ran as controls (Finding 3)
PHASE 2  Step 4   ->  half the pLOF mask meets neither clause of its documented
                      definition                                  (Phase 2, chr8 pilot)
PHASE 3  Step 6   ->  PC counts 5 / 9 / 10; the plan said "most likely 4 or 5"
PHASE 4  Step 8   ->  3 files written with no header row (AFR / pDM / chr8)
PHASE 4  Step 8   ->  11 files never written (killed for memory, no retry)
PHASE 5  Step 9.1 ->  the merge inherits both, and reports neither
```

The first three were found by the replication; the last three by the review that preceded it.

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
| Step list, sub-steps and phase grouping | the analysis blogs in `elena_publishes/`, cross-checked against the file listing of `rarevariantExWAS/` |
| `step9` is overloaded; documentation and file names disagree | `step9_saige_step2.bsub` alongside `step9_1_mergechr.bsub`, against the blogs placing SAIGE in Step 8 |
| Relatedness is deliberately not filtered | code comment in `step1_phenotyping_sample_qc.ipynb` |
| The sample list is cut against the imputed LD-pruned `.fam` | the merge in `step1_phenotyping_sample_qc.ipynb`; `.fam` verified at 70,493 samples |
| `HL_only_rarevariant` consumes backbone phenotype + covariates | paths in `step2_newcovariate.bsub` |
| Mask labels and test groups | `step4_burdenmasking.bsub`; `--annotation_in_groupTest` in the SAIGE command line |
| 792 expected vs 781 on disk | `find` count under `saige_results_newmasks/`, and 3x4x22x3 |
| Release ships VEP annotations, group files and freqcounts | directory listing of `/static/PMBB/PMBB-Release-2026-4.0/Exome/` |
| Group-file format and chr8 category counts | `awk` over `group_file_annotations/...chr8.txt` |
| BF4 4.3.0 bundle uses joint gnomAD | `bundle_plan.json` in `/project/hall_shared/datasets/biofilter/20260914/`; module `biofilter/4.3.0` |
| Defects 1-3 | see README §3, which records paths and counts |
