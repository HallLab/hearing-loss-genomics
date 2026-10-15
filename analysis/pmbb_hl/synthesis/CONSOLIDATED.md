# Hearing loss in PMBB v4 — the rebuilt analysis, consolidated

**Hall Lab × Epstein Lab · PMBB-Release-2026-4.0 · October 2026**
Every number here is read from the file the phase that produced it wrote;
`scripts/01_collect.py` assembles it.

---

## The headline

![How close](figures/how_close.png)

**Twenty-four cells across four arms — two methods × four arms × three cohorts. Nothing
reaches its threshold in any of them.** The closest is 1.9× short; the furthest is 685×.

The second thing the figure says is the useful one: **the rare-variant arms land an
order of magnitude closer than the common-variant ones, in every run** — burden 1.9–24×,
ExWAS 15–685×.

| arm | phenotype | what it is |
|---|---|---|
| **original run** | broad | the PMBB v4 analysis delivered in October |
| **rebuild, broad** | broad | the same phenotype, this pipeline |
| **rebuild, bilateral** | bilateral sensorineural | a phenotype the original never ran |
| **rebuild, tinnitus** | tinnitus | the hypothesis this project started from |

> **The original run's burden figure is its omnibus recomputed by our method.** Its
> delivered tables report **per-cell** p-values corrected within one MAF family —
> `0.05/53,583` for combined, verified from its own `Pvalue_Bonferroni` column — across
> the 160,672 cell tests it ran. Placing that beside our omnibus would compare a minimum
> over nine times more tests. The figure instead uses the same ACAT over the same nine
> cells for all four, so one quantity is judged by one rule.
>
> Its ExWAS figure needs no such adjustment: single-variant p against 5e-8 is the same
> quantity in all four arms.

### The original run is the closest of the four in combined, and among the furthest in AFR

Burden, combined: **1.9×** against our 10×, 7.4× and 6.9×. Burden, AFR: **13×** against
our 5.0× and 5.2× for the hearing-loss arms — though the tinnitus arm at 14.6× is
further still, so the original is no longer the furthest there.

Neither is a verdict on pipeline quality. The gene sets differ — its 19,042 genes
include 651 antisense, 88 LINC and 14 microRNA genes that a `pLOF` mask should not
contain, and more genes means both a looser bar and more chances at a small p. The
comparison is here to show that **no run reaches its threshold**, not to rank them.

---

## What was run

| | broad | bilateral sensorineural | tinnitus |
|---|---|---|---|
| **rare-variant gene burden** | 6,752 cases | 3,164 cases | 2,859 cases |
| **common-variant ExWAS** | the same cohorts, MAF ≥ 0.01 | | |

All three phenotypes are built by one script, so the control set is shared by
construction. **Every bilateral case is also a broad case** — asserted in code, not
assumed. **Tinnitus is not nested**: 1,885 of its 2,859 cases (66%) are also broad
cases, and that overlap is measured rather than asserted, so the nesting check is
scoped to the pair it describes.

| cohort | broad: analysable / tested / cases | bilateral | tinnitus |
|---|---|---|---|
| combined | 57,507 / 57,498 / 6,752 | 53,919 / 53,910 / 3,164 | 53,614 / 53,605 / 2,859 |
| EUR | 42,786 / 42,779 / 5,183 | 40,150 / 40,143 / 2,547 | 39,731 / 39,724 / 2,128 |
| AFR | 11,334 / 11,334 / 1,285 | 10,539 / 10,539 / 490 | 10,641 / 10,641 / 592 |

*Analysable* is what the phenotype supports; *tested* is what reached SAIGE after
covariates. The difference is 9 and 7 people with no `AGE`, **all of them controls**, so
no case is lost anywhere.

---

## The masks — what a burden test is actually testing

Built from the release's own VEP annotation (106 GB), one chromosome per task. **One
build serves all three phenotypes**: a mask depends on the annotation and the genotypes,
not on who is a case, so building it per arm would only create a chance for them to
differ.

| mask | definition | genes | entries |
|---|---|---:|---:|
| `pLOF` | seven exact consequence terms, **or** any other splice annotation with `SpliceAI ≥ 0.2` | 17,841 | 462,144 |
| `pDM` | AlphaMissense pathogenic / likely pathogenic, **or** max REVEL across transcripts ≥ 0.5 | 17,623 | 890,332 |
| `pLOF_pDM` | the union, `pLOF` taking precedence | 17,945 | 1,340,936 |

All three require `BIOTYPE == protein_coding`. 139,846,060 annotation rows fall inside
the cohort's genotypes.

**`pDM` here means "probably damaging and not truncating."** The SAIGE group file carries
one annotation per variant and gives `pLOF` precedence, so the 11,540 variants that are
both enter the `pLOF` test and are absent from the `pDM` one. That is a definitional
consequence of the file format, not a choice, and it matters when reading a `pDM` result.

### Three corrections to the original run, each measured

| | count | what it fixes |
|---|---:|---|
| rows dropped, not `protein_coding` | 71,443,642 | the original's gene set carries 651 antisense, 88 LINC and 14 microRNA genes. Its top burden gene in two of three cohorts was `TMC3-AS1`, an antisense lncRNA |
| splice annotations rejected by the `SpliceAI` gate | 2,780,725 | the original matched `"splice"` as a substring, admitting three `IMPACT=LOW` terms unconditionally. The gate its own analysis plan specifies was never applied |
| `pDM` assignments recovered from REVEL | 1,366,320 | the original coerced VEP's multi-valued REVEL entries to `NaN`, silently dropping about three quarters of `pDM`. Found by Nikki Palmiero on 2026-10-05 |

For scale, the gate rejects 2,780,725 splice candidates and admits 426,418 — it removes
about **87%** of what substring matching would have let through.

**None of this rescues a result.** Nothing is significant with the corrected masks or
without them. The fixes are what make the gene list defensible, which is a different
claim and the one worth making.

### The masks reproduce the first pass byte for byte

Not the same counts — the same files, md5-identical on all three masks. The first pass
read the 106 GB serially; this reads it as 22 parallel tasks and merges. Byte-identical
output from a restructured implementation is the evidence that the restructuring
introduced nothing.

---

## Ancestry and relatedness — the PCs and the GRM

Two different corrections, often confused, and the confusion cost this project a
rebuild. **PCs** are covariates: they put continuous ancestry into the model as fixed
effects. **The GRM** is the random effect: it tells SAIGE who is related to whom. A
script that builds one is not building the other.

### PCs: the release's own exome PCA, at 80% of the variance

| | PCs | variance reached |
|---|---:|---:|
| original run | 5 / 9 / 10 | 64.2% / 60.8% / 58.6% |
| first pass | 5 / 4 / 3 | 64.2% / 66.7% / 56.5% |
| **this rebuild** | **11 / 15 / 16** | **80.3% / 83.3% / 84.1%** |

combined / EUR / AFR. **This is the first run that meets the 80% rule**, asked for by
Nikki Palmiero on 2026-10-14.

**Nobody is dropped for missing a PC** — that was the premise of using the release's,
and the script refuses to write if it fails. The only people dropped are 9 (combined)
and 7 (EUR) with no `AGE`, all controls, in all three phenotypes.

> **A caveat carried rather than hidden.** "80% of the variance" is 80% of what the top
> twenty eigenvalues capture, because plink publishes twenty. It is not 80% of total
> genotypic variance, which twenty components come nowhere near. The rule holds — it is
> a proportion of a stated denominator, not of everything.

### Why the release's PCA and not one of our own

The first pass re-ran a within-ancestry PCA to close a coverage gap in the original
pipeline. **That gap did not exist**, and it was checkable three ways:

- the original computed **no PCA at all** — its delivered PCs are identical to the
  release's `exome_PC*` columns at `max|difference| = 0.00e+00`, across all 24
  components in all three cohorts
- there is no `--pca` in any of its scripts and no eigenvector file in its output
- the release's by-ancestry exome PCA covers **100%** of the analysable cohort

Its `step7_ldprune.bsub` builds the **GRM**. A GRM step was read as a PCA step, and a
PCA was recomputed to fix something that was not broken. That recomputation was the only
real methodological divergence from the original run, and removing it is what this
rebuild exists for.

### GRM: rebuilt for EUR and AFR, because that gap is real

| | rebuild | original |
|---|---:|---:|
| EUR markers | 53,920 | 53,913 |
| EUR samples | **51,867** | 51,632 |
| AFR markers | 114,694 | 114,696 |
| AFR samples | **14,927** | 14,878 |

The original's pruned genotypes were built over a keep list already cut. **SAIGE
intersects its phenotype file with the GRM**, so the 196 EUR and 34 AFR cohort members
absent from them would have been dropped silently — no warning, no row in any log.

Reproducing the marker counts to within 7 and 2 is what confirms the parameters are the
original's; the difference is the extra samples. **0 of the tested cohort is outside the
GRM** in any of the nine cells.

Combined is **not** rebuilt: the original's `ALL.common_ldpruned` covers all 70,925 and
has no gap.

### One problem carried forward on purpose

The combined GRM holds **38,833 markers for 70,925 people**, against EUR's 53,920 for
51,867. **It is thinner than a subset of itself.**

The cause is `--hwe 1e-6` applied **across ancestries**. Hardy-Weinberg departure in a
structured sample is a consequence of mixing populations, not evidence of genotyping
error — so the filter is removing population structure rather than bad genotypes.

Left alone deliberately. Rebuilding it would be a third difference from the original run
and would change the null model of every combined analysis. A thin GRM makes the
relatedness correction somewhat worse, which is the **conservative** direction, and no
combined analysis has come close to a bar.

**If a combined result ever approaches significance, this is the first thing to redo.**

---

## Results

### Rare-variant gene burden

Bar: `0.05 / genes`. FDR (Benjamini-Hochberg, q < 0.10) reported alongside, as asked on
2026-10-14.

| phenotype | cohort | genes | min omnibus p | best gene | short by | min q | hits |
|---|---|---:|---:|---|---:|---:|---:|
| broad | combined | 17,943 | 2.88e-05 | `BPTF` | 10.3× | 0.198 | 0 |
| broad | EUR | 17,928 | 3.31e-05 | `GYS2` | 11.9× | 0.477 | 0 |
| broad | AFR | 17,789 | 1.40e-05 | `CFH` | 5.0× | 0.248 | 0 |
| bilateral | combined | 17,941 | 2.06e-05 | `PCBD1` | 7.4× | 0.239 | 0 |
| **bilateral** | **EUR** | 17,921 | **6.33e-06** | **`FBXO25`** | **2.3×** | **0.113** | 0 |
| bilateral | AFR | 17,778 | 1.46e-05 | `IGSF9` | 5.2× | 0.259 | 0 |
| tinnitus | combined | 17,940 | 1.91e-05 | `CTNNA3` | 6.9× | 0.343 | 0 |
| tinnitus | EUR | 17,923 | 6.75e-05 | `FXYD3` | 24.2× | 0.712 | 0 |
| tinnitus | AFR | 17,777 | 4.11e-05 | `POLD1` | 14.6× | 0.496 | 0 |

### Common-variant ExWAS

Bars: `5e-8` (primary, the convention) and `0.05 / variants` (secondary, ~7e-07). No
FDR — that was explicit, and the split is principled: FDR is a threshold over a count of
tests we performed, while 5e-8 counts nothing we did.

| phenotype | cohort | variants | min p | best variant | short of 5e-8 |
|---|---|---:|---:|---|---:|
| broad | combined | 70,366 | 2.13e-06 | `20:33064486:G:A` | 43× |
| broad | EUR | 74,536 | 2.47e-05 | `3:133755303:T:G` | 495× |
| broad | AFR | 117,953 | 4.08e-06 | `7:31775671:C:G` | 82× |
| bilateral | combined | 70,350 | 2.88e-05 | `7:150692505:T:C` | 576× |
| bilateral | EUR | 74,566 | 1.66e-05 | `19:44819634:C:T` | 332× |
| **bilateral** | **AFR** | 118,049 | **7.67e-07** | `7:158133844:C:T` | 15× |
| tinnitus | combined | 70,287 | 3.43e-05 | `1:94205228:T:C` | 685× |
| tinnitus | EUR | 74,563 | 1.25e-05 | `17:39105864:C:T` | 250× |
| tinnitus | AFR | 118,067 | 2.13e-05 | `5:14758467:C:T` | 426× |

---

## The one positive result, and the test that nearly undid it

Known deafness genes cluster at the top of the burden rankings far more than chance
allows:

| arm | ClinGen Definitive/Strong in top 50 | expected | p |
|---|---:|---:|---:|
| **broad combined** | **5** | 0.26 | **5.91e-06** |
| broad EUR | 2 | 0.26 | 2.78e-02 |
| everything else | 0 | 0.26 | 1.00 |

The five are `SIX1`, `GJB3`, `COCH`, `TMPRSS3`, `MYO6`.

**Identical to the first pass, to the digit, with the same five genes** — after changing
the PCA source, raising the PC counts from 5/4/3 to 11/15/16, and adding deduplication
before the combination.

**The tinnitus arm, added on 2026-10-15, is the nearest thing to a negative control.**
ClinGen's panel lists **deafness** genes, not tinnitus genes, so it should not enrich
there — and it finds zero in the top 50 of all three cohorts. Since 66% of the tinnitus
cases are also hearing-loss cases, an enrichment driven by the pipeline rather than by
the phenotype would have been carried across by those shared cases. It was not. The
inference is weak on its own, because that arm would likely miss a real enrichment too,
but it points the right way.

### The threat, and how it was resolved

ClinGen genes turn out to be **substantially better covered** than the average gene:
median MAC 453 against 190 in the broad combined arm, Mann-Whitney p < 0.001 in all six.
Well-studied genes tend to be larger and more variable, which is part of why they were
found. So the enrichment might have been measuring coverage rather than biology.

Two checks. Calibration is suggestive: across MAC sextiles the fraction of p below 0.05
is flat at ~0.05, and the *median* p rises with coverage — the opposite of the bias in
question. **The decisive test holds coverage fixed:**

| arm | observed | 2,000 MAC-matched draws: mean | max | empirical p |
|---|---:|---:|---:|---:|
| broad combined | **5** | 0.20 | 3 | **< 0.0005** |
| broad EUR | **2** | 0.24 | 3 | **0.0220** |

**Not one of 2,000 coverage-matched draws reached 5.** The enrichment survives intact.

It remains a statement about the **ranking**. None of the five is significant on its
own, and none may be named as a finding.

---

## Why nothing is significant

### The burden arm has almost nothing to test

| percentile | distinct variants | MAC | MAC in cases |
|---:|---:|---:|---:|
| 50 | 17 | 49 | **5** |
| 90 | 45 | 373 | 44 |

**64% of genes carry fewer than ten rare alleles among the 6,752 cases.** A burden test
on five case carriers cannot detect anything whatever the true effect.

Cohort needed for 80% power at `0.05/genes`:

| | OR 1.5 | OR 3 | OR 5 |
|---|---:|---:|---:|
| median gene | 2.1M | 287k | 134k |
| 90th percentile gene | 277k | 38k | 18k |

This cohort is 57,498. **For the median gene it is not powered even at OR 5.**

### The ExWAS is limited by effect size, not coverage

Everything it tests is above 1% by construction, so there is no untestable tail. Cohort
needed at `5e-8`:

| | OR 1.05 | OR 1.1 | OR 1.2 |
|---|---:|---:|---:|
| median tested variant (MAF 0.050) | 1.7M | 444k | 121k |
| common variant (MAF 0.25) | 428k | 112k | 31k |

Common-variant effects on complex age-related traits sit at **OR 1.05–1.15**. At the
median frequency this cohort would need 444,000 people to see OR 1.1.

### And the sharpest statement the project can make

In the four well-covered ClinGen genes — `GJB3`, `COCH`, `MYO6`, `TMPRSS3` — this study
**would have detected OR ≥ 1.7 and did not.** So if their effect exists in this
population, it is smaller than that: real but modest, consistent with polygenic risk in
adult-onset hearing loss rather than the Mendelian effects that put them on the panel.

`SIX1` is the exception worth noting: MAC 87, detectable only from OR 6.31, and still
the **4th gene** in the whole ranking.

---

## What a replication elsewhere would buy

Cohort needed to detect each effect in the four well-covered ClinGen genes:

| | OR 1.5 | OR 1.3 | OR 1.2 |
|---|---:|---:|---:|
| `GJB3` | 100k | **239k** | 494k |
| `COCH` | 108k | **258k** | 534k |
| `MYO6` | 115k | 274k | 568k |
| `TMPRSS3` | 119k | 285k | 591k |

A cohort at **All of Us's scale (~250k) lands exactly where OR 1.3 becomes visible** in
all four. UK Biobank's (~470k) approaches OR 1.2.

We already know the effect is **below 1.7**. Whether it is above 1.2 or 1.3 is what a
larger cohort would settle, and 250k sits in the middle of the plausible range.

Three honest caveats:

| | |
|---|---|
| `SIX1` is not rescued | 87 alleles; it would need 1.2M for OR 1.5 |
| the phenotype must be comparable | our prevalence is 11.7% on phecodeX with a two-date rule. **This has not been checked for either biobank and is the first thing to confirm.** |
| the power model is optimistic | it assumes every masked variant is causal, so real requirements are larger |

*Biobank figures are approximate and drawn to show scale.*

---

## What changed in the rebuild, and what it cost

| | |
|---|---|
| **PCs from the release**, at 80% variance (11/15/16) | the first pass re-ran a PCA to close a gap that did not exist — the original computed no PCA at all, its PCs are identical to the release's at `max\|difference\| = 0.00e+00`, and its `step7_ldprune.bsub` builds the GRM. A GRM step was read as a PCA step. |
| **GRM still rebuilt** for EUR and AFR | that gap is real: 196 and 34 cohort members are absent from the original's pruned genotypes |
| **Deduplication before the combination** | raised by Nikki Palmiero. Removes 1,542–6,911 cells per arm and changes the omnibus in 424–875 genes |
| **FDR alongside Bonferroni**, burden only | asked for on 2026-10-14; nothing clears q < 0.1, best is 0.113 |

**The MHC question is answered.** The first pass found λ = 1.3386 inside the MHC in broad
EUR with 4 PCs of its own PCA. With 15 release PCs the same arm reads **0.873** — from
inflated to conservative.

And it opened a new one: the **bilateral** arms now show that inflation (1.39, 1.23,
1.16) with the same PCs, GRM and region — only the phenotype differs. Most likely noise
at 3,164 cases on ~1,020 variants, but the first time the MHC has behaved differently
between the phenotypes with everything else fixed.

---

## What the set says

1. **Nothing is significant anywhere**, under any phenotype, either method, any cohort,
   against either of that cell's bars.
2. **The pipeline works.** Known deafness genes cluster at p = 5.9e-06, and that survives
   a test designed to break it.
3. **The limiting factor is the cohort.** The median gene carries five rare alleles among
   the cases.
4. **Rare-variant work is the better use of a larger cohort**, by an order of magnitude
   in how close it gets.
5. **The experiment that would answer the question is now specified**: ~250k people with
   comparable phenotyping resolves OR 1.3 in four known deafness genes.
6. **The hypothesis this project started from does not hold at this cohort size.**
   `ZNF175` — Park et al. 2021, *Nature Medicine*, tinnitus — was pre-registered as the
   declared target and reported in all eighteen cells whatever it said. It reaches nothing,
   and in the combined arm it is well enough covered (MAC 587 against a median of 178)
   that the statement is sharper than "underpowered": an effect of OR ≥ 2.8 had an 80%
   chance of being seen. See [`ZNF175.md`](../ZNF175.md).

### Next

| | |
|---|---|
| confirm phenotype comparability in All of Us / UK Biobank | the one thing that could invalidate the plan above |
| more cases, not more methods | every cell is power-limited; another analysis of the same 6,752 will not change it |
| size any replication for the rarest phenotype | prevalence enters power directly — tinnitus at 5.3% cases needs **twice** the cohort of broad hearing loss at 11.7% for the same variant and effect |
| revisit the combined GRM if a combined result ever approaches a bar | `--hwe` across ancestries leaves it thinner than EUR alone |
