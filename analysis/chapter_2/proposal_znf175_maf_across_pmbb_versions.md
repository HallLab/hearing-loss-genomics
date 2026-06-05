# Proposal — Comparing ZNF175 rare-variant MAF between Park (PMBB v1) and Hui (PMBB v2)

**Author:** Andre Rico
**Date:** 2026-06-05
**Status:** NB 00–05 done (substrate stable + annotated). Next: NB 06–07 = ZNF175 burden regression, Park pipeline on both cohorts — data assembly confirmed on-disk (see §2c).
**Project context:** Chapter 2 — resolving the ZNF175–tinnitus signal-loss puzzle (11K → 45K).
**All analysis lives in:** `analysis/chapter_2/` (notebooks in `scripts/0X/`, outputs in `results/0X/`)

---

## 1. Motivation — the 11K → 45K puzzle

The **ZNF175 → tinnitus** association was significant at ~11K exomes (Park et al., *Nat Med* 2021, pLOF discovery) and **vanished at ~45K** (Daniel's replication — negative). The central question from the 2026-06-03 meeting (Doug/Molly):

> **Is the signal loss a technical artifact or a real null?**

Molly's directive is explicit: **rule out the technical explanation before invoking biology**. Andre has already reproduced the ~45K negative result. This proposal supplies the other half of the mirror.

This is the **cheapest concrete first step** of Action Item #1: before re-running any association test, check whether the **ZNF175 variant substrate is stable** between the two cohorts.

> **Principle:** if the set of rare ZNF175 variants (or their MAFs) is not stable between freezes, re-running the association is pointless. Substrate stability is a **prerequisite**, not a consequence.

---

## 2. The proposal in one sentence

Pull **all variants at the ZNF175 locus** and compare their **presence and MAF** between **PMBB v1 (Park / Freeze One)** and **PMBB v2 (Hui)**, in order to **isolate the technical component** (pipeline / calling / QC / normalization) from the non-technical changes (cohort growth and ancestry composition).

---

## 2b. Results so far (NB 00–05, executed 2026-06-05)

Notebooks in `scripts/0X/`, outputs in `results/0X/`. All extraction/annotation via LSF jobs (`scripts/0X/*.sh`).

| NB | What | Headline result |
|---|---|---|
| **00** | Inventory & cohort-ID assessment | v1=Freeze One (11,451, GRCh38) ⟂ v2=Release 2020 v2.0 (43,731, GRCh38). **Per-sample concordance BLOCKED** — IDs in different namespaces (`UPENN_*` vs `PMBB*`), need a crosswalk. |
| **01** | Independent v1 extraction (NCBI coords `chr19:51,571,283-51,592,510`, plink2) | 151 variants (132 rare <0.1%) — *not* Daniel's curated set. |
| **02** | v1 × v2 MAF comparison (symmetric `bcftools norm -f -a` + plink2) | **Substrate STABLE**: union 386 → 155 shared (MAF r=**0.94**), 225 v2-only (expected, 4× cohort), 6 v1-only (all explained: `*` alleles + private singletons). The "24% disappearing variant" was an **MNV-vs-SNV representation** artifact, fixed by atomize. |
| **03** | Annotation via **Biofilter 4** | 97/386 found; 6 pLOF; **0 AlphaMissense** for ZNF175. **75% not_found because the loaded BF4 DB was filtered to gnomAD AC≥5** (AC<5 excluded — they exist in gnomAD; reversible by reloading without the cutoff). |
| **04** | Annotation via **VEP + dbNSFP4.5a** (by coordinate) | **384/386** annotated; **24 pLOF**; **177 AlphaMissense**. Covers the AC<5 variants BF4's cutoff excluded. |
| **05** | **VEP × BF4** reconciliation | 94% consequence agreement where both annotate; VEP superior on coverage + AlphaMissense + pLOF. BF4's gap is the AC≥5 DB filter, not an inherent limit. |

> **Conclusion (substrate question):** the ZNF175 variant substrate is **stable** between v1 and v2 → the 11K→45K signal loss is **NOT** a variant-substrate artifact. It points to **winner's curse** or **association-pipeline differences** (pLOF→missense in replications; phenotype-release version). ⚠️ This compares MAF, not the association — the definitive verdict needs the **association re-run** with a matched pipeline.

**Annotation-tool notes (LPC):** see memories `reference-biofilter4-lpc` and `reference-vep-annovar-lpc`. Open refinements: VEP `--most_severe` (vs `--pick`) for a pLOF-focused pass; LOFTEE HC/LC needs `human_ancestor`/gerp data (absent).

---

## 2c. Next phase (NB 06–07) — reproduce the ZNF175 burden, Park pipeline on BOTH cohorts

**The actual project question:** why did Park find ZNF175→tinnitus at 11K but Hui lose it at ~44K? Strategy: apply **ONE fixed pipeline (Park's)** to both cohorts, holding filters + phenotype + model constant, so any difference in result is attributable to the **cohort** (sample size/composition) — isolating winner's curse / dilution from a pipeline artifact.

### ✅ Data-assembly unblock (2026-06-05) — both cohorts are fully assemblable on-disk

The `UPENN_* ↔ phenotype` ID gap is **resolved**: the Freeze17 **Demographics** file is the crosswalk.

| Piece | v1 (Park / Freeze One) | v2 (Hui / Release 2020 v2.0) |
|---|---|---|
| ZNF175 genotype | region VCF (NB 02), ID `UPENN_*` | region VCF (NB 02), ID `PMBB_*` |
| **ID crosswalk + sex + age + race** | `/static/PMBB/PMBB_Freeze17/phenotype/PMBB_Geno_Demographics_Deidentified_012020.csv` — cols `PT_ID, GENO_ID(=.fam ID), GENDER_CODE, BIRTH_YEAR, RACE_CODE` | `data/pmbb_v2/Phenotype/2.0/..._phenotype_covariates.txt` — `PMBB_ID, sex, Age_at_Enrollment` |
| **PCs + ancestry (GIA)** | `/static/PMBB/PMBB_Freeze17/genotype/imputed/Additional_Files/OMNI_GSA_V1_V2_evec_gia.txt` — `IID, PC1–PC20, GIA` (AFR/EUR/…) | covariates `PC1–PC10` (+ a GIA-equivalent ancestry label TBD) |
| **Tinnitus phenotype** | `..._phenotype/PMBB_Geno_Nonsensitive_Diagnosis_Deidentified_012020.csv` — `PT_ID, ICD` (has `H93.1x` / `388.3x`) | `data/pmbb_v2/Phenotype/2.0/..._icd-10-matrix.txt` / `..._icd.txt` |

v1 join chain: `.fam GENO_ID → Demographics → PT_ID → Diagnosis (tinnitus)`; PCs/GIA via `evec_gia` keyed by the `UPENN_*` IID.

### NB 06 — qualifying variant set (Park filter, both cohorts)
From the NB 04 annotated table: keep **pLOF** (VEP consequence: frameshift / stop_gained / stop_lost / splice_donor / splice_acceptor) with **cohort MAF ≤ 0.1%**, applied **identically** to v1 and v2. Output: qualifying ZNF175 variants per cohort + **per-sample carrier status** (from the region VCFs). Also flag missense **REVEL > 0.6** separately (for an optional Hui-definition pass). Reveals whether the burden *input* changed between 11K and 44K.

### NB 07 — gene-focused burden regression (both cohorts)
Assemble per sample: **carrier (burden)** + **tinnitus** (phecode 389.x, rule-of-2: ≥2 dates = case, never = control) + **age, age², sex, PCs** + **GIA**. Run **logistic regression `tinnitus ~ burden + covariates`, stratified EUR/AFR + IVW meta-analysis** (Park's model). Run on v1 and v2 with the identical pipeline → the verdict: signal present at 11K but gone at 44K under the SAME pipeline ⇒ winner's curse / dilution; only gone under Hui's different pipeline ⇒ pipeline artifact.

---

## 3. Cohort & data mapping (corrected)

Important clarification established 2026-06-04/05 — "PMBB version" means different things across data modalities:

| Cohort (our label) | What it really is | Modality | N | Build |
|---|---|---|---|---|
| **v1 (Park)** | **"Freeze One"** — a pre-release RGC exome freeze (predates the numbered public releases) | WES | **11,451** | **GRCh38** |
| **v2 (Hui)** | PMBB **Release 2020 v2.0 WES** (= Daniel's "45K") | WES | 44,086 | GRCh38 |

Key facts:
- **There is no WES in public PMBB Release v1.0** — that release is imputed-array only (21,263 samples). WES first appears at **Release v2.0 (44,086)**. So Park's "11K" was never a numbered public release; it is the **Freeze One** RGC freeze.
- **The Freeze One `.fam` has exactly 11,451 samples** — matching the Park paper verbatim ("subset of 11,451 individuals... who have undergone WES"). Strong confirmation this is Park's discovery data.
- ✅ **Both v1 and v2 are GRCh38** → **no liftover needed between the two cohorts we compare.** This removes the build-difference confound for this specific comparison.
- ⚠️ Caveat: the Park *paper* cites GRCh37 ("mapped to GRCh37"), but the Freeze One files on disk are labeled **GRCh38** (`UPENN_Freeze_One_GRCh38.*`). So exact reproduction of the paper's published coordinates may differ; for our v1-vs-v2 MAF comparison both are GRCh38, which is what matters here.

### Data location & layout (v1 / Freeze One)

`/static/PMBB/PMBB_Freeze17/genotype/exome/`
- `all_variants/` — `UPENN_Freeze_One_GRCh38.{NF,GL}.pVCF{,.biallelic}.{vcf.gz,bed/bim/fam}` (+ annotation/VQSR companion files, no genotypes). Biallelic-split PLINK: **4,276,878 variants**.
- `PASS_variants/` — same, restricted to **sites passing VQSR**.

File flavors (per the freeze `README`):
- **NF** = "No Filter" — unfiltered, all positions / all calls.
- **GL** = genotype-quality filtered — low-quality genotypes set to no-call (SNP: QD<3 or DP<7; INDEL: QD<5 or DP<10) plus site-level Allele-Balance filters.
- `PASS_variants/` adds the **VQSR PASS** site filter on top.

### Data location & layout (v2 / Hui) — confirmed via `results/chapter1_paper_replication/`

PMBB **Release 2020 v2.0** Exome, **GL** flavor, per-chromosome pVCF. ZNF175 is on **chr19**:
- **Raw pVCF (chr19):** `data/pmbb_v2/Exome/pVCF/GL_by_chrom/PMBB-Release-2020-2.0_genetic_exome_chr19_GL.vcf.gz` (**56 GB**, has `.tbi` → region extraction via `tabix` is cheap; no need to read the whole file). Full set = 781 GB across 22 chrs.
- **Variant annotations:** `/static/PMBB/PMBB-Release-2020-2.0/Exome/Variant_annotations/PMBB-Release-2020-2.0_genetic_exome_variant-annotation-counts.txt` (5.4 GB).
- **Phenotype covariates:** `data/pmbb_v2/Phenotype/2.0/PMBB-Release-2020-2.0_phenotype_covariates.txt` (43,732 rows).
- **Daniel's pre-extracted chr19 PLINK (fastest route):** `data/PMBB_Exome/genotypes/allIndvs_chr19.{bed,bim,fam}.gz` (full, pre-IBD, 43,731 samples) and `allIndvs_chr19_maf.001_noRels_keepHLcases.{bed,bim,fam}.gz` (his Phase-6 filtered set).
- **N:** 44,086 released → **43,731** in the pVCF after IBD filtering (matches Daniel's runbook).
- **Validated tool:** plink **1.9** (v1.90b6.18) reproduces Daniel's chr19/chr21 outputs byte-for-byte.

> **Apples-to-apples note:** v2 here is the **GL** flavor (genotype-quality filtered). To match, use the **GL** flavor of v1 (Freeze One) too — and decide deliberately on `all_variants` vs `PASS_variants` (record VQSR status as a column rather than pre-filtering).

---

## 4. What this diagnoses

A variant-level comparison answers questions the burden re-run does **not**:

1. **Is the substrate stable?** As PMBB grew and was re-processed (joint-calling, QC, VQSR, normalization), the set of ZNF175 pLOF variants can change — variants present in v1 may disappear in v2 (re-genotyped as ref, QC-dropped) or appear new.
2. **Did the "rare" variants stay rare?** A variant at MAF 0.0008 in v1 may become 0.0012 in v2 and **cross the 0.1% threshold** → exits the burden set. A concrete mechanism for the burden set changing.
3. **Are the ~8 signal-driving cases** still present and still rare? (The original signal was *small-case* → **winner's-curse** territory.)
4. **Multiallelic / stop-loss:** the repo already carries `data/PMBB_Exome/addBack_multiallelic_stoploss/`, signalling that multiallelic/stop-loss normalization of ZNF175 was a **known prior problem**. A cross-freeze MAF comparison surfaces exactly this kind of artifact. **High prior of finding something.**

---

## 5. ⚠️ Design: don't conflate three sources of MAF change

A change in a variant's MAF between v1 and v2 can have **three causes**; only one is the target:

| Source of MAF change | Target? |
|---|---|
| **(1) Pipeline / calling / QC / VQSR / normalization** | ✅ **YES** — the "technical artifact" |
| **(2) Cohort growth** (11K→44K; more people → estimate refines; rare regresses to its true frequency) | ❌ winner's-curse mechanics, not an artifact |
| **(3) Ancestry composition** (different EUR/AFR mix per freeze) | ❌ pure confounder (MAF is ancestry-dependent) |

*(Note: build difference is NOT a source here — both cohorts are GRCh38.)*

### How to isolate the technical component (1):

- **Compare MAF WITHIN ancestry strata** (MAF_EUR and MAF_AFR separately), **never** global → removes (3).
- **If the cohorts are nested** (v1 ⊂ v2 in people, which is likely): run **per-sample genotype concordance** on the shared individuals. Since the person is identical, any genotype-call difference between v1 and v2 is **100% technical** → isolates (1), removing (2) and (3). **This is the strongest test.**
- **Track variants by normalized representation (HGVS / left-aligned, split multiallelics)**, not by raw coordinate — this is where multiallelic/stop-loss issues hide.

---

## 6. Variant ascertainment — pull ALL, then annotate (anti-circularity)

There are two different "which variants" questions with **opposite** answers:

| Context | Which variants | Why |
|---|---|---|
| Burden test (the association) | filter **hard & early**: pLOF + MAF≤0.1% | you want only the clean analysis set |
| **This diagnostic** (cross-freeze forensics) | take **ALL**, harmonize, filter **last** | the filter *status changing* is what we are hunting |

The trap: the pLOF/MAF filter is **exactly the thing whose stability we are investigating**. Pre-filtering each freeze independently would **blind** us to the interesting case — a variant that was pLOF+rare in v1 but crossed 0.1% or got re-annotated in v2 would simply vanish from both filtered sets. Pre-filtering is **circular**. Ascertainment must be **filter-agnostic**.

**Procedure:**
1. **Define the ZNF175 locus** (all exons + splice sites + small flank) on GRCh38.
2. **Pull ALL variants** in that region from each cohort, **unfiltered** (every site, every sample genotype).
3. **Harmonize**: normalize (left-align, split multiallelics), match by HGVS/normalized representation (no liftover needed — both GRCh38).
4. **Annotate the UNION with ONE consistent pipeline** (e.g. ANNOVAR/RefSeq or VEP) so any pLOF-status difference reflects the *variant* changing, not the annotation tool.
5. **Then classify/flag as columns** (pLOF? MAF_EUR/AFR? N_carriers? VQSR-PASS?) — keeping all variants visible. The pLOF+rare subset is *highlighted*; the rest stays in the table to explain *why* variants enter/leave.

### ⚠️ ZNF175-specific gotchas
- **Transcript definition:** RefSeq vs GENCODE differ at exon boundaries → fix one, record it.
- **KRAB-ZNF cluster (19q13.43):** ZNF175 sits in a dense paralog cluster, notoriously hard to map. Paralog reads can mismap → artifactual "variants". A pipeline/aligner change between freezes could shift this → a **strong** candidate for instability. Inspect closely.
- **VQSR PASS status:** record as a column rather than pre-filtering (whether v1/v2 used all_variants vs PASS could itself explain a MAF difference).

---

## 7. Step zero (quick & decisive)

1. **Confirm the v1↔v2 relationship:** intersect the sample-ID lists (Freeze One `.fam`, 11,451 vs PMBB v2 `.fam`, 43,731) → nested? What N is shared? *(Decides: per-sample concordance vs aggregate-MAF only.)*
2. ~~Locate the v2 (Hui) WES data~~ ✅ **Done** — see §3 (`data/pmbb_v2/Exome/pVCF/GL_by_chrom/`, GRCh38).
3. **Locate existing assets:** `data/PMBB_Exome/ZNF175/` and `data/PMBB_Exome/addBack_multiallelic_stoploss/` — what is already extracted, and in which build.
4. **Get the ZNF175 GRCh38 locus coordinates** (chr19) from the annotation file / gene model — to drive `tabix` region extraction on both pVCFs.

---

## 8. Deliverable

Master table of ZNF175 variants:

| Variant (HGVS, GRCh38) | v1 (Park/Freeze One) | v2 (Hui) |
|---|---|---|
| present? (Y/N) | | |
| MAF_EUR | | |
| MAF_AFR | | |
| N_carriers | | |
| pLOF? (Y/N) | | |
| VQSR PASS? | | |

Plus, if cohorts are nested:
- **Per-sample genotype concordance matrix** (shared IDs, v1 vs v2) for the key variants.
- **Highlight of the ~8 cases** that drove the original signal.

---

## 9. Outcome tree (every result is informative)

- **Variants / MAF unstable between freezes** → we found (part of) the **technical explanation** Molly asked for. 🎯 (likely foci: multiallelic/stop-loss normalization, MAF drift across 0.1%, QC/VQSR drops, KRAB-ZNF mismapping.)
- **Rock-stable across freezes** → substrate **ruled out** as culprit → the loss leans toward **winner's curse / real null**, and we narrow the problem to the *association* pipeline (e.g. the meeting's note that cross-biobank replications shifted to **missense**, not pLOF).

---

## 10. Scope — what this is NOT

- **Not** a re-test of the ZNF175–tinnitus association (this is variant forensics / QC; it comes first).
- **Not** an audiometry / quantitative-trait analysis (that phase is blocked by audiogram↔ID matching, meeting Action #7; lives in v4).
- **Not** the new exome-wide HL study on v4 (Nikki/Elena's lane, SKAT methods).
- **Does not assume nesting** — verifies it.

---

## 11. Open questions / to confirm

- [~] Is v1 nested within v2? → **Could not test by ID** (namespaces differ: `UPENN_*` vs `PMBB*`); per-sample concordance **blocked pending an ID crosswalk**. Comparison done on aggregate MAF instead (NB 02).
- [x] PMBB v2 (Hui) WES path on LPC, build, file layout? → `data/pmbb_v2/Exome/pVCF/GL_by_chrom/`, GRCh38, GL pVCF (see §3).
- [x] File flavor for v1: used **GL** (matches v2 GL). Both extracted via the same `bcftools norm -f -a` pipeline (NB 02).
- [x] Variant-set definition: pulled **ALL** variants (anti-circularity), then flagged pLOF via VEP IMPACT=HIGH + dbNSFP (NB 04). Reference MAF computed **internally** (plink2) per cohort, whole-cohort.
- [x] **ID crosswalk for v1** — RESOLVED: Freeze17 `Demographics` file maps `.fam GENO_ID ↔ PT_ID`, so v1 genotype links to v1 phenotype/covariates. (Note: this links v1-geno↔v1-pheno; a separate `UPENN_* ↔ PMBB` map would still be needed for v1↔v2 per-sample concordance — not required for the per-cohort regression.)
- [x] **Ancestry labels** — available: v1 via `evec_gia` GIA (+PC1–20); v2 via covariates PCs (GIA-equiv label TBD).
- [ ] **MAF by ancestry** — NB 02 used whole-cohort MAF; ancestry-stratified MAF/regression now possible via the labels above (NB 07).
- [ ] **Phenotype-release versioning** (PMBB phenotype v2.0→v2.3): did Park and Daniel use different phenotype releases? Check when assembling the tinnitus phecode (NB 07).
- [ ] **Burden regression (NB 06–07)** — Park pipeline (pLOF + MAF≤0.1%) on v1 and v2, tinnitus phecode, EUR/AFR-stratified + IVW meta. The definitive test (see §2c).

---

## Differences between the Park and Hui papers (for reference)

- **Phenotype:** Park's signal = **ZNF175 → tinnitus** (HL just under threshold), phecode-based; Hui defined **hearing loss** with an audiogram-hybrid definition. ZNF175 is **not** in the Hui paper.
- **Variant set:** Park = **pLOF only** for discovery (missense REVEL≥0.5 reserved for robustness); Hui = **pLOF + missense REVEL>0.6 in the primary burden**. Different threshold *and* role.
- **Scope:** Park = **exome-wide** (~1,518 genes); Hui = **HL-focused** (~173 known HL genes).
- **Cohort:** Park ~11K WES (Freeze One); Hui ~40K (Release v2.0 WES).

---

## Project references
- `docs/papers/paper_summary_park2021.md` — Park 2021 summary (source of ZNF175→tinnitus).
- `docs/papers/andre_notes_park2021.md` — method study notes (pLOF, MAF, burden, statistical model).
- `docs/meetings/2026-06-03_research-vision-and-elena-onboarding.md` — meeting minutes (11K→45K puzzle, action items).
- `data/PMBB_Exome/ZNF175/`, `data/PMBB_Exome/addBack_multiallelic_stoploss/` — Daniel's assets.
- `analysis/daniel/` — Hui replication runbook (sibling pipeline).
- **v1 data:** `/static/PMBB/PMBB_Freeze17/genotype/exome/` (Freeze One, GRCh38, 11,451 samples).
- **v2 data:** `data/pmbb_v2/Exome/pVCF/GL_by_chrom/` (Release 2020 v2.0, GRCh38, 43,731 samples); chr19 = ZNF175.
- `results/chapter1_paper_replication/` — chapter 1 reports (source of the v2 paths, plink-1.9 validation).
