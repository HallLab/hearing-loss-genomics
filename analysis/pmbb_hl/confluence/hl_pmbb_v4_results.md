# Hearing loss in PMBB v4 — results

**Hall Lab × Epstein Lab · PMBB-Release-2026-4.0 · October 2026**
Andre Rico · two analyses × three phenotypes × three cohorts = **18 cells**

---

> ### How to use this page
>
> **Lines marked `SAY` are the talk.** Read only those, in order, and you have delivered
> the whole thing in about ten minutes. Nothing else needs to be spoken.
>
> **Sections marked `BACKUP` are for questions.** Do not present them. Each one names
> the question it answers, so you can jump straight there.
>
> Every number on this page is read from the file that produced it. Source paths are at
> the bottom.

---

## 1 · The answer, before the detail

> **SAY —** *"We ran both analyses, on three phenotypes, in three cohorts. Eighteen
> cells. Nothing is significant in any of them — and we can now say exactly why, and
> exactly what would change it."*

> **SAY —** *"Three things came out of it that are worth your time: the pipeline
> demonstrably works, the limit is the cohort and not the method, and the hypothesis
> this project started from does not survive."*

**Nothing is significant — here is how close each analysis got.** The two are reported
separately because they count different things: a gene against a Bonferroni bar over
genes, a variant against the genome-wide convention.

| | best result | its bar | short by |
|---|---|---:|---:|
| **Rare-variant burden** | `FBXO25`, bilateral EUR, p = 6.33e-06 | Bonferroni `0.05 / 17,921 genes` = 2.79e-06 | **2.3×** |
| **Common-variant ExWAS** | `7:158133844:C:T`, bilateral AFR, p = 7.67e-07 | genome-wide `5e-8` | **15×** |

And the three things worth taking away:

| | |
|---|---|
| **The pipeline finds real biology** | known deafness genes cluster at the top, p = 5.9e-06, and that survived a test built to break it |
| **The limit is the cohort** | the median gene carries **five** rare alleles among the cases |
| **The origin hypothesis fails** | `ZNF175` reaches nothing, and not for lack of coverage |

---

## 2 · What we ran

**Figure:** `synthesis/figures/pipeline.png`

> **SAY —** *"Three phenotypes, built by one script, so they share a control set by
> construction rather than by coincidence."*

| cohort | HL broad | bilateral sensorineural | tinnitus |
|---|---|---|---|
| | *tested / cases* | *tested / cases* | *tested / cases* |
| combined | 57,498 / 6,752 | 53,910 / 3,164 | 53,605 / 2,859 |
| EUR | 42,779 / 5,183 | 40,143 / 2,547 | 39,724 / 2,128 |
| AFR | 11,334 / 1,285 | 10,539 / 490 | 10,641 / 592 |

**Controls: 50,755** in every arm — no ear evidence from any source.

> **SAY —** *"Every bilateral case is also a broad case, and the script refuses to run if
> that ever stops being true. Tinnitus is different — it is not a kind of hearing loss,
> so it overlaps rather than nests. Sixty-six percent."*

### BACKUP · *"Why is tinnitus not nested like the other one?"*

A person can have tinnitus without hearing loss. Running the nesting assertion over it
would have stopped the pipeline on a correct cohort, so the check is scoped to the pair
it describes and the tinnitus relation is measured instead.

| | | |
|---|---:|---|
| tinnitus cases | 2,859 | |
| also broad cases | **1,885** | 66% |
| not broad cases | 974 | |
| …no hearing-loss evidence at all | 581 | tinnitus and nothing else |
| …hearing loss on one date only | 393 | below the two-date rule |

**Consequence:** no comparison between tinnitus and a hearing-loss arm is independent.
Two thirds of the cases and *all* the controls are the same people.

### BACKUP · *"Why read the OMOP `observation` table?"*

PMBB v4 moved the standard tinnitus ICD codes (388.3x, H93.1x) out of the phecode
files. Those keep 3,016 of 28,110 tinnitus events; the other **25,094 exist only in
`observation`**. An arm built on the phecode files alone would miss 89% of its own
phenotype — and a control set built without that source silently contains tinnitus
patients.

---

## 3 · What a burden test is actually testing

> **SAY —** *"Each gene gets up to nine tests — three annotation masks crossed with
> three frequency ceilings — and those nine are combined into one p-value per gene."*

| mask | definition | genes |
|---|---|---:|
| `pLOF` | seven exact consequence terms, **or** any other splice annotation with `SpliceAI ≥ 0.2` | 17,841 |
| `pDM` | AlphaMissense pathogenic / likely pathogenic, **or** max REVEL ≥ 0.5 | 17,623 |
| `pLOF_pDM` | the union, `pLOF` taking precedence | 17,945 |

> **SAY —** *"We corrected three defects in how those masks were built. None of them
> rescues a result — nothing is significant with the corrections or without them. They
> are what makes the gene list defensible, which is a different claim and the one worth
> making."*

### BACKUP · *"What were the three corrections?"*

| | count | what it fixes |
|---|---:|---|
| rows dropped, not `protein_coding` | 71,443,642 | the original gene set carries 651 antisense, 88 LINC and 14 microRNA genes. Its top burden gene in two of three cohorts was `TMC3-AS1`, an antisense lncRNA |
| splice annotations rejected by the `SpliceAI` gate | 2,780,725 | `"splice"` was matched as a substring, admitting three `IMPACT=LOW` terms unconditionally. The gate the original analysis plan specifies was never applied — it removes about **87%** of what substring matching lets through |
| `pDM` recovered from REVEL | 1,366,320 | VEP's multi-valued REVEL entries were coerced to `NaN`, dropping about three quarters of `pDM`. Found by Nikki Palmiero, 2026-10-05 |

The masks reproduce the first pass **md5-identical** on all three files, from a
restructured implementation — which is the evidence the restructuring introduced
nothing.

### BACKUP · *"What does `pDM` mean exactly?"*

**"Probably damaging and not truncating."** The SAIGE group file carries one annotation
per variant and gives `pLOF` precedence, so the 11,540 variants that are both enter the
`pLOF` test and are *absent* from the `pDM` one. A consequence of the file format, not a
choice — but it changes how a `pDM` result reads.

---

## 4 · Ancestry and relatedness

**Figure:** `rare_variant/phase_3/figures/variance.png`

> **SAY —** *"Principal components are covariates — they put ancestry in the model.
> The GRM is a random effect — it tells SAIGE who is related to whom. They are different
> corrections, and confusing them cost this project a rebuild."*

> **SAY —** *"This is the first run of this analysis that meets the eighty percent
> variance rule."*

| | PCs | variance reached |
|---|---|---:|
| original run | 5 / 9 / 10 | 64.2% / 60.8% / 58.6% |
| first pass | 5 / 4 / 3 | 64.2% / 66.7% / 56.5% |
| **this rebuild** | **11 / 15 / 16** | **80.3% / 83.3% / 84.1%** |

### BACKUP · *"Why the release's PCA and not your own?"*

The first pass re-ran a PCA to close a coverage gap in the original pipeline. **That gap
did not exist**, and it was checkable three ways:

- the original computed **no PCA at all** — its PCs are identical to the release's at
  `max|difference| = 0.00e+00`, across all 24 components in all three cohorts
- there is no `--pca` in any of its scripts and no eigenvector file in its output
- the release's by-ancestry exome PCA covers **100%** of the analysable cohort

Its `step7_ldprune.bsub` builds the **GRM**. A GRM step was read as a PCA step. That
recomputation was the only real methodological divergence from the original run, and
removing it is what this rebuild exists for.

### BACKUP · *"Does '80% of the variance' mean 80% of everything?"*

No, and this should be said rather than hidden. It is 80% of what the **top twenty**
eigenvalues capture, because plink publishes twenty. Twenty components come nowhere near
80% of total genotypic variance. The rule holds — it is a proportion of a stated
denominator.

### BACKUP · *"Was the GRM gap real?"*

Yes, and it is the one thing we did rebuild. The original's pruned genotypes were built
over a keep list already cut: **196 EUR and 34 AFR** cohort members were absent. SAIGE
intersects its phenotype file with the GRM, so those people would have been dropped
**silently** — no warning, no log line. Marker counts reproduce to within 7 and 2, which
is what confirms the parameters are the original's. **0 of the tested cohort is outside
the GRM in any of the nine cells.**

### BACKUP · *"Anything you knowingly left broken?"*

One, deliberately. The combined GRM holds **38,833 markers for 70,925 people** against
EUR's 53,920 for 51,867 — thinner than a subset of itself. The cause is `--hwe 1e-6`
applied **across ancestries**, which removes population structure rather than bad
genotypes. Left alone because fixing it would be a third difference from the original
run and would change the null model of every combined analysis; a thin GRM errs
**conservative**, and no combined analysis has come close to a bar. **If a combined
result ever approaches significance, this is the first thing to redo.**

---

## 5 · Results

**Figure:** `rare_variant/phase_5/figures/manhattan.png`

> **SAY —** *"Rare-variant burden, nine cells. Each point is a gene, not a variant — so
> there are no peaks to look for. One high point alone is all a real signal would
> produce here. Nothing reaches the line."*

| phenotype | cohort | genes | min omnibus p | best gene | short by | min q |
|---|---|---:|---:|---|---:|---:|
| broad | combined | 17,943 | 2.88e-05 | `BPTF` | 10.3× | 0.198 |
| broad | EUR | 17,928 | 3.31e-05 | `GYS2` | 11.9× | 0.477 |
| broad | AFR | 17,789 | 1.40e-05 | `CFH` | 5.0× | 0.248 |
| bilateral | combined | 17,941 | 2.06e-05 | `PCBD1` | 7.4× | 0.239 |
| **bilateral** | **EUR** | 17,921 | **6.33e-06** | **`FBXO25`** | **2.3×** | **0.113** |
| bilateral | AFR | 17,778 | 1.46e-05 | `IGSF9` | 5.2× | 0.259 |
| tinnitus | combined | 17,940 | 1.91e-05 | `CTNNA3` | 6.9× | 0.343 |
| tinnitus | EUR | 17,923 | 6.75e-05 | `FXYD3` | 24.2× | 0.712 |
| tinnitus | AFR | 17,777 | 4.11e-05 | `POLD1` | 14.6× | 0.496 |

Bar: `0.05 / genes` ≈ 2.79e-06, plus FDR q < 0.10. **Zero hits under either.**

**Figure:** `exwas/phase_5/figures/manhattan.png`

> **SAY —** *"Common-variant ExWAS, same nine cells. Here each point is a variant. Best
> is bilateral AFR at 7.7e-07 — fifteen-fold short of the five-times-ten-to-the-minus-
> eight convention."*

**Figure:** `synthesis/figures/how_close.png`

> **SAY —** *"This is the one figure to look at. Twenty-four cells — our eighteen plus
> the original run's six — and the x-axis is how far each one fell short of its own bar.
> Everything is to the right of the line."*

> **SAY —** *"The closest cell on the whole figure is the original run's burden in
> combined, at one-point-nine. Ours is two-point-three. Neither reaches, and the gap
> between them is not a verdict on either pipeline — their gene set is bigger, which
> buys both a looser bar and more chances at a small p."*

> **SAY —** *"And the second reading is the useful one: the rare-variant arms land an
> order of magnitude closer than the common-variant ones, in every single run. That is
> the statistical version of the preference this group already held on biological
> grounds."*

### BACKUP · *"`FBXO25` is close — what is it?"*

Nothing yet, and the honest answer is thin. Its entire signal is `pLOF` in the ultra-rare
band: **16 variants, 18 alleles, 7 in cases against 11 in controls** among 2,547 cases
and 37,596 controls. SAIGE's `beta` implies OR 41; a Fisher test on the same counts gives
9.4 — that gap is what an effect estimated on seven alleles looks like. `pDM` alone is
null (p = 0.83) and adding it only dilutes.

There is also an open methodological question about `FBXO25` specifically — see §9.

### BACKUP · *"Why two different bars?"*

They answer different questions and the split is principled. `5e-8` is the genome-wide
convention and counts nothing we did. FDR is a threshold over a count of tests we
actually performed, which is why it applies to the burden arm and not to the ExWAS.

**The ExWAS carries two bars and they read very differently**, so have both:

| | bar | our best is short by |
|---|---:|---:|
| genome-wide convention | `5e-8` | **15×** |
| secondary Bonferroni, `0.05 / 118,049 variants` | 4.24e-07 | 1.8× |

The convention is the stricter and is the one reported. Quote 15×, not 1.8×.

---

## 6 · The one positive result

**Figure:** `rare_variant/phase_5/figures/clingen.png`

> **SAY —** *"Known deafness genes cluster at the top of the broad combined ranking far
> more than chance allows. Five ClinGen Definitive-or-Strong genes in the top fifty,
> against point-two-six expected. p equals five-point-nine times ten to the minus six."*

| arm | ClinGen in top 50 | expected | p |
|---|---:|---:|---:|
| **broad combined** | **5** | 0.26 | **5.91e-06** |
| broad EUR | 2 | 0.26 | 2.78e-02 |
| broad AFR | 0 | 0.26 | 1.00 |
| bilateral, all three | 0 | 0.26 | 1.00 |
| tinnitus, all three | 0 | 0.26 | 1.00 |

The five are `SIX1`, `GJB3`, `COCH`, `TMPRSS3`, `MYO6`.

> **SAY —** *"The obvious objection is that these genes are simply better covered — they
> are, more than twice over. So we tested it directly, and it failed: two thousand draws
> of ninety-three genes matched on coverage never once reached five. The result is
> stronger after that test than before it."*

> **SAY —** *"This is a statement about the ranking, not about any gene. None of the
> five is significant on its own, and none can be named as a finding."*

### BACKUP · *"How exactly did you test the coverage objection?"*

ClinGen genes have median MAC **453 against 190** (Mann-Whitney p = 1.1e-11), so the
objection is real and had to be answered rather than argued with.

2,000 draws of 93 non-ClinGen genes, matched on log₁₀ MAC within a window of 0.10:

| | |
|---|---|
| observed | **5** ClinGen genes in the top 50 |
| matched draws | mean 0.20, **max 3**, never reached 5 |
| empirical p | **< 0.0005** |

A weaker check points the same way: across MAC sextiles the fraction of p below 0.05 is
flat at ~0.05, and the *median* p **rises** with coverage — the opposite of the bias in
question.

### BACKUP · *"Did the rebuild change this result?"*

No — **identical to the first pass, to the digit, with the same five genes**, after
changing the PCA source, raising PC counts from 5/4/3 to 11/15/16, and adding
deduplication before the combination. At the looser Definitive+Strong+Moderate tier it is
stronger still: 6 of the top 50 against 0.31 expected, p = 6.57e-07.

### BACKUP · *"Tinnitus found zero — is that a failure?"*

The opposite; it is the closest thing to a negative control we have. ClinGen's panel
lists **deafness** genes, not tinnitus genes, so it should not enrich there. Since 66% of
tinnitus cases are shared with the broad arm, an artefact of the *pipeline* would have
carried across. It did not. Weak on its own — that arm would likely miss a real
enrichment too — but it points the right way.

---

## 7 · Why nothing is significant

**Figure:** `rare_variant/phase_5/figures/power.png`

> **SAY —** *"There is almost nothing to test. The median gene carries forty-nine rare
> alleles across fifty-seven thousand people — five of them in cases. Sixty-four percent
> of genes carry fewer than ten. A burden test on five case carriers cannot detect
> anything, whatever the true effect is."*

Cohort needed for 80% power at `0.05/genes`:

| | OR 1.5 | OR 3 | OR 5 |
|---|---:|---:|---:|
| median gene | 2.1M | 287k | 134k |
| 90th percentile gene | 277k | 38k | 18k |

**This cohort is 57,498. For the median gene it is not powered even at OR 5.**

**Figure:** `exwas/phase_6/figures/power.png`

> **SAY —** *"The common-variant side has the opposite problem — coverage is fine,
> effect size is the limit. Common-variant effects on age-related traits sit around odds
> ratio one-point-one. To see that at the median frequency we would need four hundred
> and forty-four thousand people."*

> **SAY —** *"So the honest headline is not 'there is no association'. It is 'no effect
> larger than this cohort could see' — and this cohort can see very little."*

### BACKUP · *"What is the sharpest thing you can actually claim?"*

In the four well-covered ClinGen genes — `GJB3`, `COCH`, `MYO6`, `TMPRSS3` — this study
**would have detected OR ≥ 1.7 and did not**. So if their effect exists in this
population it is smaller than that: real but modest, consistent with polygenic risk in
adult-onset hearing loss rather than the Mendelian effects that put them on the panel.

`SIX1` is the exception worth naming: MAC 87, detectable only from OR 6.31, and still the
**4th gene** in the whole ranking.

### BACKUP · *"Would All of Us or UK Biobank solve it?"*

For the rare-variant arm, substantially: ~250k people resolves OR 1.3 in those four
genes. Three caveats, and say all three —

1. `SIX1` is **not** rescued at that scale; 87 alleles would need 1.2M.
2. The power model is optimistic: it assumes every masked variant is causal.
3. **Phenotype comparability in those biobanks has not been checked.** That is the one
   thing that could invalidate the plan.

And a fourth if tinnitus comes up: prevalence enters power directly, so a cohort sized
for broad hearing loss is **half** what tinnitus needs.

---

## 8 · The hypothesis we started from

> **SAY —** *"This project exists because of `ZNF175` — Park 2021 in Nature Medicine
> reported it associated with tinnitus in an earlier PMBB release. No arm of this project
> had ever tested tinnitus. We added one, and we pre-registered the gene as the declared
> target before running it, so the result would be reported whatever it said."*

> **SAY —** *"It reaches nothing. And this one is not a coverage failure — the gene
> carries five hundred and eighty-seven rare alleles against a median of a hundred and
> seventy-eight. So the statement is sharper than 'underpowered': an effect of odds ratio
> two-point-eight or larger had an eighty percent chance of being seen, and was not."*

| | best of all 18 cells | bar |
|---|---:|---:|
| burden omnibus | p = 0.324 (tinnitus EUR) | 2.79e-06 |
| ExWAS | p = 0.014 (**bilateral** AFR) | 5e-8 |

> ⚠ **Note to self, not a line to read.** Say *"does not hold at this cohort size"*,
> never *"is wrong"*. Smaller effects are completely untouched by this result, and Park
> 2021 was a different release and a different pipeline.

### BACKUP · *"Was the direction a surprise?"*

No, and it was recorded in advance so it could not later be presented as a discovery. The
original run's own delivered tables put `ZNF175` **better on hearing loss** (p = 0.040
combined, 0.022 EUR) than on tinnitus (0.331 / 0.139 / 0.064) — the wrong way round for
the hypothesis. Our corrected run does not reproduce that ordering either way; both are
far too distant from any bar for the comparison to carry weight.

Full detail, all 18 cells: `ZNF175.md`.

---

## 9 · One open methodological question

> **SAY —** *"One thing I want a decision on rather than a nod."*

Each gene's nine cells are combined with ACAT. Two kinds of cell can be identical, and we
handle only one:

| | status |
|---|---|
| **same MAF cutoff, different mask** — a gene with no `pDM` variants makes `pLOF` and `pLOF;pDM` the same test | **removed before combining** (raised by Nikki, 2026-10-14) |
| **same mask, looser cutoff** — a gene with no variants between two ceilings makes those cells the same test | **not removed** |

`FBXO25` is the clean case: 16 `pLOF` variants, all ultra-rare, so its three `pLOF` cells
are **one test counted three times** and its best p enters with weight 3/9 instead of 1/5.
**56–68% of cells have at least one twin.**

Neither treatment is invalid — ACAT admits arbitrary weights. But the implicit weighting
is *"how many cutoffs did this test survive"*, which was never a stated choice and favours
genes whose signal is ultra-rare.

> **SAY —** *"No arm crosses any bar either way, so no conclusion turns on it. One
> sentence does: `FBXO25` is thirteen percent short of the FDR bar as we report it, and
> eighty-eight percent short if each distinct test counts once. Since it is our nearest
> miss, that sentence gets quoted — so I would like the choice made rather than
> inherited."*

| arm | min q, as reported | min q, distinct tests |
|---|---:|---:|
| broad combined | 0.198 | 0.198 |
| broad EUR | 0.477 | 0.446 |
| broad AFR | 0.248 | 0.248 |
| bilateral combined | 0.239 | 0.247 |
| **bilateral EUR** | **0.113** | **0.188** |
| bilateral AFR | 0.259 | 0.344 |
| tinnitus combined | 0.343 | 0.411 |
| tinnitus EUR | 0.712 | 0.801 |
| tinnitus AFR | 0.496 | 0.608 |

Both are already carried in `rare_variant/phase_4/results/step2/all_cells.tsv`, per cell
and per gene, so the decision can be made on the numbers.

---

## 10 · What I need next

> **SAY —** *"Three things, in order of how much they matter."*

| | |
|---|---|
| **1. Confirm phenotype comparability in All of Us / UK Biobank** | the single thing that could invalidate the replication plan |
| **2. More cases, not more methods** | every cell is power-limited; another analysis of the same 6,752 people will not change anything |
| **3. A decision on §9** | and on whether the poster correction (backup B1) goes anywhere |

---

## BACKUP · B1 — the poster, only if raised

The summer poster reports six genome-wide significant results — `KNL1` and `ZFPM1`,
three unique variants across two cohorts.

| variant | in the v4 release? | measured call rate | poster's own filter |
|---|---|---:|---|
| `chr15:40611486:A:G` | yes | **0.2013** | < 0.99 → should have been removed |
| `chr16:88533289:A:C` | yes | **0.9617** | < 0.99 → should have been removed |
| `chr16:88533290:G:C` | yes | **0.9473** | < 0.99 → should have been removed |

> *"variants with a genotype call rate < 0.99 and a Minor Allele Frequency (MAF) < 0.01
> were removed"* — the poster's own methods

All three are **absent from the pipeline's own QC output**: they exist in the release and
do not survive `variant_qc`, the stage named for call rate. They were never in the data
step 2 reads.

**We reproduce the delivered ExWAS files at Pearson r = 0.979–0.990** on −log₁₀p across
70k–118k shared variants. The poster disagrees with its own pipeline's output, not with
us.

**Framing:** this is "three variants did not clear a filter the poster says was applied",
not "the poster is wrong about the biology". The finding belongs to the pipeline, not to
a person — Elena is an undergraduate on a summer programme. Have the numbers; do not
volunteer the slide.

---

## Things not to say

| don't | do |
|---|---|
| "we found nothing" | "nothing reaches the threshold, and we measured how far away it is" |
| "we found `COCH` / `MYO6`" | "known deafness genes cluster at the top of the ranking" |
| "the bilateral restriction destroyed the signal" | "at 3,164 cases the restriction cannot answer the question either way" |
| "our pipeline is better than the original" | "no run reaches its threshold; the comparison is not a ranking" |
| "there is no association with hearing loss" | "no effect larger than this cohort could see" |
| "`ZNF175` was wrong" / "Park was wrong" | "it does not hold at this cohort size; in the combined arm we would have seen OR ≥ 2.8 and did not" |
| "tinnitus found nothing either, same story" | "it is the rarest phenotype — at 5.3% cases it needs twice the cohort" |

---

## Figures to attach, in order of use

| § | file |
|---|---|
| 2 | `synthesis/figures/pipeline.png` |
| 4 | `rare_variant/phase_3/figures/variance.png` |
| 5 | `rare_variant/phase_5/figures/manhattan.png` |
| 5 | `exwas/phase_5/figures/manhattan.png` |
| 5 | `synthesis/figures/how_close.png` |
| 6 | `rare_variant/phase_5/figures/clingen.png` |
| 7 | `rare_variant/phase_5/figures/power.png` |
| 7 | `exwas/phase_6/figures/power.png` |

Held in reserve, not attached unless asked:
`rare_variant/phase_3/figures/what_pcs_see.png` (what a PC is measuring),
`rare_variant/phase_5/figures/vs_first_pass.png` (did the rebuild move anything —
ρ = 0.93–0.97), `exwas/phase_5/figures/lambda_mhc.png` (is the MHC under control).

## Where the numbers come from

All paths relative to `analysis/pmbb_hl/`.

| | |
|---|---|
| cohorts | `rare_variant/phase_1/results/01_cohorts.json` |
| masks | `rare_variant/phase_2/results/02_masks_manifest.json` |
| PCs, GRM | `rare_variant/phase_3/results/01_covariates.json` |
| burden, per gene | `rare_variant/phase_4/results/omnibus_<pheno>_<cohort>.tsv` |
| burden, every cell | `rare_variant/phase_4/results/step2/all_cells.tsv` — 1.4M rows |
| ExWAS, per variant | `exwas/phase_4/results/merged_<pheno>_<cohort>.tsv` |
| ClinGen enrichment | `rare_variant/phase_5/results/02_enrichment.json` |
| power, matched test | `rare_variant/phase_6/results/` |
| the full written version | `synthesis/CONSOLIDATED.md` |
| `ZNF175` | `ZNF175.md` · pre-registration in `PREMISES_tinnitus.md` |

Each phase also writes a `FINDINGS.md` beside its results.
