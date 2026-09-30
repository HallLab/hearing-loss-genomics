# Step 1, explained — who gets into the study, and how 427 people did not

**Author:** Andre Rico · **Date:** 2026-09-30 · **Type:** explainer
Written to be read with no background. The formal record is
[`phase_1/results/FINDINGS.md`](../phase_1/results/FINDINGS.md); the decision it asks for is
[page 01](01_sample_frame_decision.md).

---

## What Step 1 is for

The pipeline asks: *is there a gene where rare damaging variants are more common in people with
hearing loss than in people without?*

Before that can be asked, two things must be settled — **who has hearing loss**, and **who is in the
study at all**. Step 1 settles both. Everything downstream inherits its answer, which is why it runs
first and why an error here cannot be repaired later.

---

## What Step 1 does, in four moves

### 1. Pull every ear-related diagnosis from the health records

**Why:** the phenotype has to come from somewhere, and in a hospital biobank that somewhere is
billing codes entered during ordinary care.

Two traps sit here. Diagnosis codes are grouped into *phecodes* — `SO_396` is hearing impairment,
`SO_39x` is the whole ear family — because a raw billing code is too fine-grained to be a phenotype.
And in PMBB v4 the tinnitus codes were **moved into a different table** than every other diagnosis.
Missing that relocation silently loses most tinnitus evidence. This pipeline handled it correctly.

### 2. Decide who is a case — the rule of 2

A person counts as a case only if the diagnosis appears on **two separate dates**.

**Why:** one mention can be a suspicion that was ruled out, a code entered to justify a test, or a
typo. Two separate encounters is evidence that a clinician kept thinking so. It is a deliberately
conservative rule: it costs real cases (4,007 people here) to avoid admitting false ones.

### 3. Set aside the ambiguous middle

People with ear-family evidence that is *not* hearing impairment are excluded — neither case nor
control.

**Why:** they are the people who would blur the comparison. Someone with ear disease might have
undiagnosed hearing loss, so counting them as a healthy control makes the two groups look more alike
than they are and hides real signal. 9,411 people are set aside this way.

### 4. Build the list of people who enter the analysis

**Why:** a phenotype is useless without genetic data to pair it with. Someone in the records but
never sequenced cannot contribute, and leaving them in would break the software.

**This is the move that goes wrong.**

---

## Where it goes wrong

The code is three lines:

```python
fam_file = ".../Imputed/common_snps_LD_pruned/..._genetic_imputed.commonsnps.ldpruned.ALL.fam"
# Keep only samples with genotype data
matched = pheno[pheno["PMBB_ID"].isin(fam["IID"])].copy()
```

The comment states the right intent: keep only people who have genetic data.

The problem is **which** genetic data. PMBB holds two different kinds, from two different
laboratory processes:

| | what it is |
|---|---|
| **Exome sequencing** | reads the protein-coding genes directly, letter by letter. Finds rare variants, including ones never seen before. |
| **Array + imputation** | reads a fixed panel of common positions, then statistically infers the rest. Cheap and broad; blind to genuinely rare variation. |

Most participants have both. Some have only one.

**This analysis is built on exome data from end to end.** The null model is fitted on exome
variants (`exome_ldpruned/*`), and every gene is tested on exome variants
(`plink_deduplicated/*`). The imputed data is never used for anything.

But the line above checks the **imputed** list. So the question it actually asks is *"does this
person have array data?"* when the question it needed to ask was *"does this person have exome
data?"*

An exome list existed — built by this same pipeline, in the same directory tree. It was not the one
consulted.

---

## What that cost

```
70,925   people with exome data
57,507   analysable: case or control for hearing impairment
57,080   survived the filter
   427   removed -- 40 cases, 387 controls
```

The 427 belong to a group of **517** who have exome data, have complete exome ancestry components,
and have no array data at all. Every one of them could have been analysed. None of them was.

Worth being precise about two things people get wrong on first reading:

- **All 517 have the same problem** — no array data. The gap between 517 and 427 is not about which
  data they lacked; it is that 90 of them had already been excluded on phenotype grounds, so
  removing them again changed nothing.
- **Nobody was removed for lacking exome data.** Zero people in the release are missing exome
  ancestry components. Had the filter checked the exome list, it would have removed no one.

---

## Why this is an internal inconsistency, not a disagreement

It would be a weaker finding if it were a difference of opinion — one analyst preferring the imputed
cohort, another the exome cohort. It is not.

**Every analytical step in this pipeline runs on exome data.** The filter gates on a dataset that no
downstream step consumes. The pipeline disagrees with itself, and the check that would have caught
it — does the sample list match the data being tested? — is one nobody runs, because it sounds
obviously true.

---

## What is not being claimed

- **Not that the phenotype is wrong.** It is right. Cases and controls were re-derived independently
  from the release and agree person-for-person, 70,925 of 70,925. The rule of 2, the ambiguous-middle
  exclusion, and the tinnitus-table fix all reproduce exactly.
- **Not that this was careless.** The intent is stated in the comment and it is the correct intent.
  The two `.fam` files differ by one word in a long path, and 0.7% of a cohort is below what any
  summary statistic would reveal.
- **Not that the results are void.** 40 cases out of 6,752 is a small loss of power, and the effect
  is to make associations slightly harder to detect, not to create false ones.

**But the loss is not evenly spread, and that has now been checked.** Cohort-wide the exclusion rate
is 0.73%. Among East Asian participants it is **15.60%** — 208 of 1,333, one in six. No other group
comes close (EUR 0.45%, AFR 0.33%). The same skew appears on sex (1.01% of women against 0.43% of
men) and strongly on batch (1.09% of batch 1 against 0.11–0.18% of the others), and the excluded are
younger and enrolled earlier.

That pattern looks technical — exome-sequenced early, array genotyping never completed for a subset —
rather than anything about the people. It does not change the direction of bias. It does mean the
loss cannot be called neutral: East Asian participants are 1.6% of the analysed cohort, and a sixth
of them are gone.

---

## The one-sentence version

The study runs on exome data, but the list of participants was filtered by who had array data — so
517 people who could have been analysed, 40 of them with hearing loss, were left out for missing
something the analysis never uses.
