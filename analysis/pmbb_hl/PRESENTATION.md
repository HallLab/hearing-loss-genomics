# Hearing loss in PMBB v4 — presentation outline

**Hall Lab × Epstein Lab · October 2026 · Andre Rico**

13 slides. Each one below gives what goes **on** the slide and, in italics, what to
**say**. All paths are relative to `analysis/pmbb_hl/`.

---

## 1 · Title

> **Rare-variant and common-variant analysis of hearing loss in PMBB v4**
> Two analyses × two phenotypes × three cohorts
> Andre Rico · Hall Lab × Epstein Lab · October 2026

*Nothing to say. Move on.*

---

## 2 · The question

> **Do rare or common coding variants associate with adult-onset hearing loss in 70,925
> PMBB participants with exome data?**
>
> And if not — **how far are we from being able to tell?**

*The second line is the one that makes this talk different from "we found nothing".
Say it out loud.*

---

## 3 · The pipeline

> **figure:** `synthesis/figures/pipeline.png`

*Point at two things only. Everything above step 2 is shared between the two analyses —
same cohorts, same covariates, same null models — so a difference downstream is
attributable. And the ExWAS does not refit the null model, because it does not depend
on the MAF filter; that is also what the original run did.*

---

## 4 · Phenotypes

> |  | definition | cases | analysable |
> |---|---|---:|---:|
> | **broad** | `SO_396` or any child, ≥ 2 distinct dates | 6,752 | 57,507 |
> | **bilateral sensorineural** | `SO_396.2` **and** `SO_396.8` on the **same date**, ≥ 2 dates | 3,164 | 53,919 |
>
> **controls: 50,755** — no ear evidence from any source. Shared by both arms.
>
> Every bilateral case is also a broad case — asserted in code, not assumed.

*The same-date rule is the whole bilateral definition and it is worth thirty seconds.
phecodeX splits a diagnosis across two axes — type and laterality — as separate rows, so
someone coded conductive-bilateral at one visit and sensorineural-unilateral at another
carries both codes without ever having had bilateral sensorineural hearing loss. The
loose reading gives 4,029 cases; this one gives 3,164.*

*Also: tinnitus lives in the OMOP `observation` table in v4, not in the phecode files.
Those hold 3,016 of 28,110 events. A control set built without that source silently
contains tinnitus patients.*

---

## 5 · Genotypes and variants

> **Rare-variant masks** — from the release VEP, 106 GB, 22 parallel tasks
>
> | mask | genes | entries |
> |---|---:|---:|
> | pLOF | 17,841 | 462,144 |
> | pDM | 17,623 | 890,332 |
> | pLOF + pDM | 17,945 | 1,340,936 |
>
> **Common-variant census** — MAF ≥ 0.01 → **70,366 to 118,049** variants per cohort,
> about 0.5% of the 13.2M on file
>
> Masks exclude non-coding genes. `pLOF` presupposes a protein to lose.

*If asked why the gene counts differ from the original run's 19,042: its set carries 651
antisense, 88 LINC and 14 microRNA genes, and its top burden gene in two of three
cohorts was `TMC3-AS1`, an antisense lncRNA.*

---

## 6 · Principal components

> **figure:** `rare_variant/phase_3/figures/variance.png`
>
> PCs from the **release's own exome PCA**, at **80% of the variance**:
> **combined 11 (80.3%) · EUR 15 (83.3%) · AFR 16 (84.1%)**
>
> | run | PCs | variance reached |
> |---|---|---:|
> | original | 5 / 9 / 10 | 64% / 61% / 59% |
> | our first pass | 5 / 4 / 3 | 64% / 67% / 57% |
> | **this run** | **11 / 15 / 16** | **80% / 83% / 84%** |

*This is the first run that meets the convention. And it fixed something measurable:
our first pass left λ = 1.34 inside the MHC in EUR; with the release's 15 PCs the same
arm reads 0.873.*

*Caveat if anyone presses: "80% of the variance" is 80% of what the top twenty
eigenvalues capture, because plink publishes twenty. It is a proportion of a stated
denominator, not of total genotypic variance.*

> **optional backup slide:** `rare_variant/phase_3/figures/what_pcs_see.png` — what a PC
> is actually measuring, if anyone asks

---

## 7 · Results — all 18 cells

> **figure:** `synthesis/figures/how_close.png`
>
> **Nothing reaches its threshold anywhere.** Closest 1.9× short, furthest 576×.
>
> **And the rare-variant arms land an order of magnitude closer than the common-variant
> ones — in every run.**

*The second line is the actionable one. Burden 1.9–13× short, ExWAS 15–576×. That is the
statistical version of the preference the group already held on biological grounds.*

*The original run's burden figure is its omnibus recomputed by our method — its delivered
tables report per-cell p corrected within one MAF family, which is a different quantity.
Say this before anyone asks.*

---

## 8 · The numbers behind that

> **Rare-variant gene burden** — bar `0.05/genes`, plus FDR q < 0.10
>
> | phenotype | cohort | min p | best gene | min q |
> |---|---|---:|---|---:|
> | broad | combined | 2.88e-05 | `BPTF` | 0.198 |
> | broad | EUR | 3.31e-05 | `GYS2` | 0.477 |
> | broad | AFR | 1.40e-05 | `CFH` | 0.248 |
> | bilateral | combined | 2.06e-05 | `PCBD1` | 0.239 |
> | **bilateral** | **EUR** | **6.33e-06** | **`FBXO25`** | **0.113** |
> | bilateral | AFR | 1.46e-05 | `IGSF9` | 0.259 |
>
> **Common-variant ExWAS** — bars `5e-8` and `0.05/variants`
>
> best: `7:158133844:C:T`, bilateral AFR, **7.67e-07** — 15× short of 5e-8

*`FBXO25` is the only thing worth watching: q = 0.113 against a 0.10 bar, 13% short. It
was the same arm's top gene in our first pass at 1.37e-05 and moved to 6.33e-06 here —
it **strengthened** when ancestry correction improved, which is the direction a real
signal moves and the opposite of residual stratification. It is not a result.*

---

## 9 · The one positive result

> **figure:** `rare_variant/phase_5/figures/clingen.png`
>
> Known deafness genes cluster at the top of the broad combined ranking:
> **5 ClinGen Definitive/Strong genes in the top 50 against 0.26 expected, p = 5.9e-06**
>
> `SIX1` · `GJB3` · `COCH` · `TMPRSS3` · `MYO6`
>
> **Identical to our first pass — same five genes, same p** — after changing the PCA
> source, the PC counts, and adding deduplication.

*A gene is significant where its point sits above a line. None does. This is a statement
about the **ranking**, not about any gene — none of the five is significant on its own
and none can be named as a finding.*

---

## 10 · The test that tried to break it

> ClinGen genes are **better covered** than average — median MAC 453 vs 190, p < 0.001
> in all six arms. So: biology, or just where the power is?
>
> **Hold coverage fixed.** 2,000 draws of 93 non-ClinGen genes matched on MAC:
>
> | | observed | matched draws: mean | max | empirical p |
> |---|---:|---:|---:|---:|
> | broad combined | **5** | 0.20 | 3 | **< 0.0005** |
> | broad EUR | **2** | 0.24 | 3 | **0.0220** |
>
> **Not one of 2,000 matched draws reached 5.**

*This is the slide to spend time on. The obvious alternative explanation was tested
directly and failed, so the result is stronger after the test than before it.*

---

## 11 · Why nothing is significant

> **figure:** `rare_variant/phase_5/figures/power.png`
>
> The median gene carries **17 variants summing to 49 alleles** across 57,498 people —
> **five of them in cases.** 64% of genes have fewer than ten rare alleles among cases.
>
> Cohort needed for 80% power at `0.05/genes`:
>
> | | OR 1.5 | OR 3 | OR 5 |
> |---|---:|---:|---:|
> | median gene | 2.1M | 287k | 134k |
> | 90th percentile gene | 277k | 38k | 18k |
>
> **This cohort is 57,498. For the median gene it is not powered even at OR 5.**

*A burden test on five case carriers cannot detect anything, whatever the true effect
is. That is the answer to "why nothing", and it is measured rather than assumed.*

---

## 12 · What would change it

> In the four **well-covered** ClinGen genes, this study **would have detected OR ≥ 1.7
> and did not.** So the effect, if real, is below that.
>
> | gene | OR 1.5 | OR 1.3 | OR 1.2 |
> |---|---:|---:|---:|
> | `GJB3` | 100k | **239k** | 494k |
> | `COCH` | 108k | **258k** | 534k |
> | `MYO6` | 115k | 274k | 568k |
> | `TMPRSS3` | 119k | 285k | 591k |
>
> **A cohort at All of Us's scale (~250k) lands exactly where OR 1.3 becomes visible in
> all four.**

*Three caveats, say them: `SIX1` is not rescued at that scale — 87 alleles, it would
need 1.2M. The power model is optimistic, it assumes every masked variant is causal.
And **the phenotype comparability in All of Us and UK Biobank has not been checked** —
that is the first thing that could invalidate this slide.*

---

## 13 · Where everything is

> **Read first:** `synthesis/CONSOLIDATED.md` — all 18 cells in one place
> **The pipeline:** `PIPELINE.md` — every decision with its reason
>
> **Full results, every p-value:**
>
> | | |
> |---|---|
> | burden, one row per gene | `rare_variant/phase_4/results/omnibus_<phenotype>_<cohort>.tsv` — 17,943 rows |
> | burden, every cell before combination | `rare_variant/phase_4/results/cells_<phenotype>_<cohort>.tsv` — 159,948 rows |
> | ExWAS, one row per variant | `exwas/phase_4/results/merged_<phenotype>_<cohort>.tsv` — 70,366 rows |
>
> Six files each, `{broad,bilateral}_{combined,EUR,AFR}`.
>
> Each phase writes a `FINDINGS.md` beside its results.

---

## Backup slides, if asked

| question | slide |
|---|---|
| "what is a principal component measuring?" | `rare_variant/phase_3/figures/what_pcs_see.png` |
| "did changing the PCs move anything?" | `rare_variant/phase_5/figures/vs_first_pass.png` — ρ = 0.93–0.97 |
| "is the MHC under control?" | `exwas/phase_5/figures/lambda_mhc.png` |
| "what about the common-variant power?" | `exwas/phase_6/figures/power.png` |
| "where are the genome-wide plots?" | `rare_variant/phase_5/figures/manhattan.png`, `exwas/phase_5/figures/manhattan.png` |

---

## Things not to say

| don't | do |
|---|---|
| "we found nothing" | "nothing reaches the threshold, and we measured how far away it is" |
| "we found `COCH` / `MYO6`" | "known deafness genes cluster at the top of the ranking" |
| "the bilateral restriction destroyed the signal" | "at 3,164 cases the restriction cannot answer the question either way" |
| "our pipeline is better than the original" | "no run reaches its threshold; the comparison is not a ranking" |
| "there is no association with hearing loss" | "no effect larger than this cohort could see" |
