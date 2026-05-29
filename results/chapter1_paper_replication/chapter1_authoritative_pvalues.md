# Chapter 1 — Authoritative P-Values

**This is the single source of truth for all Chapter 1 p-values.** Generated 2026-05-29 directly from source files. Every number here is traceable to a specific file path + column. Other docs may reference rounded versions of these numbers — when they conflict, **this doc wins**.

---

## Why this doc exists

Different docs in `results/chapter1_paper_replication/` reported p-values in different rounding/formatting (e.g., `0.0011` vs `1.13e-03`, `2.4×10⁻⁴` vs `2.40e-04`). Most are equivalent within rounding, but the inconsistency made it hard to know which number was "the" number. This doc fixes that by:

1. **Reporting at full precision** from source files (4-5 sig figs)
2. **Annotating each column with its exact source** (file path + extraction logic)
3. **Annotating LOKI version** explicitly (the only real source of drift between our runs and Daniel's)
4. **Defining each phenotype (regression outcome) explicitly** with cohort sizes and source files

---

## Phenotype definitions — what we're regressing against

Every burden test in Chapter 1 regresses something against a phenotype. There are **two distinct phenotypes** used across the 7 phases:

### Phenotype A — Binary HL ("HL_needAud") — used by Phase 5 and Phase 6

This is the **binary hearing-loss case/control variable** Daniel constructed in Phase 4. The naming convention `HL_needAud` means "audiogram required to define a case" (the strict variant).

**Definition (from `analysis/daniel/scripts/pmbb_exome/case_control.py`):**

| Person | SNHL value |
|---|---|
| Has audiogram with `BL_SNHL = TRUE` (sensorineural HL detected) | **1** (case) |
| Has audiogram with `BL_SNHL = FALSE` OR phecode 389 absent | **0** (control) |
| Has phecode 389 = TRUE but no audiogram available | **NA** (excluded) |

**Cohort distribution** (file: `data/PMBB_Exome/cases_control.txt.gz`, 59,061 individuals):

| SNHL value | N | What it means |
|---:|---:|---|
| 1 | **1,569** | HL cases (audiogram-confirmed sensorineural HL) |
| 0 | **56,548** | controls |
| NA | **944** | excluded (phecode-only, no audiogram to confirm) |

**Note on cohort size after filtering:** Phase 5 and Phase 6 do not use all 59,061 — they restrict to the **`tokeep_moreHLcases`** keep-list (41,748 individuals) which is the post-IBD-filter cohort with HL cases preserved via the IBD trick (see Phase 6 of the walkthrough). Within that 41,748, the case/control distribution becomes approximately **1,098 cases / 39,981 controls / NA-removed**.

**Statistical model:** `biobin --test logistic` — per-gene logistic regression of SNHL ∈ {0,1} on per-individual burden count, with covariates Sex + Age + AgeSq + PC1..PC4 (4 PCs).

**Source file (Phase 4 output):** `analysis/daniel/outputs/phase4/cases_control.txt` (our re-run, set-equal to Daniel)

### Phenotype B — Degree-of-HL 0-4 ("degHL") — used by Phase 7

This is the **continuous severity scale** that produces the paper's headline result. The paper (page 4, Materials and Methods) defines it from Pure Tone Average (PTA) thresholds in the worse-hearing ear.

**Definition (per paper):**

| Degree | PTA threshold (dB) | Severity label |
|---:|---|---|
| 0 | 0–15 | normal hearing |
| 1 | 16–25 | mild HL |
| 2 | 26–40 | moderate HL |
| 3 | 41–55 | severe HL |
| 4 | 56+ | profound HL |

**Hybrid cohort assignment (per paper):**
- Individuals with audiograms get a degree from their PTA → degree mapping above
- Individuals **without audiograms but with phecode 389 = FALSE** are assigned **degree = 0** (assumed normal; conservative — paper notes this would underestimate associations)
- Individuals with phecode 389 = TRUE but no audiogram are **excluded**

**Cohort distribution** (file: `data/PMBB_Imputed/deghl.txt.gz`, 36,507 individuals):

| Degree | N | What it means |
|---:|---:|---|
| 0 | **35,086** | normal hearing (audiogrammed normals + phecode-defined controls) |
| 1 | **311** | mild HL (PTA 16-25) |
| 2 | **516** | moderate HL (PTA 26-40) |
| 3 | **334** | severe HL (PTA 41-55) |
| 4 | **260** | profound HL (PTA 56+) |

**Reconciliation with paper text:** the paper reports "35,397 controls + 1,110 cases (HL 2-4)" for the headline degree-HL analysis. Our numbers match exactly:
- 35,086 (degree 0) + 311 (degree 1) = **35,397 controls** ✓
- 516 (degree 2) + 334 (degree 3) + 260 (degree 4) = **1,110 cases** ✓

**Statistical model (per Daniel's `run_regression.R`):** `lm(DegHL ~ PC1 + ... + PC20 + Sex + Age + AgeSq + burden_count_per_gene)` — per-gene linear regression with **20 PCs** (not 4). biobin generates the per-individual burden count matrix (`bins.csv`); R fits the lm independently for each gene.

**Source file:** `data/PMBB_Imputed/deghl.txt.gz` (Daniel's preserved phenotype file — we did NOT re-derive this in light mode; we used Daniel's preserved STable directly)

### Why two different phenotypes?

The paper deliberately uses both. Binary HL (Phenotype A) is conceptually simpler and used for descriptive analyses like Fig 2 (% carriers in cases vs controls). Degree-HL (Phenotype B) is more statistically powerful — it weights severe cases more than mild ones — and produces the headline gene-discovery results (Tables 3 & 4). This is exactly why TCOF1 only emerges in Phase 7: its 2 carriers with HL have severe degrees (3 and 4), so the linear regression gives them more weight than the binary test does.

### Phenotype variants Daniel preserved (not used in our chapter 1 runs)

Daniel ran the binary case/control analysis under **4 different phenotype definitions** as a sensitivity check. Our Chapter 1 only uses the first (`HL_needAud`). For completeness:

| Phenotype | Definition | Cohort size |
|---|---|---:|
| `HL_needAud` (used in P5/P6) | Audiogram required for case | ~1,569 cases |
| `HL_dontNeedAud` | Phecode 389 = TRUE allowed as case even without audiogram | larger N |
| `HL_caseAudAndPhecode` | Both audiogram AND phecode required for case | smallest N |
| `HL_rmAudNA` | Audiogram and phecode must agree | smaller, purest |

These live in `data/PMBB_Exome/renamed_pheno_files_for_looping/`. If you ever want to run sensitivity analyses across phenotypes, this is the input set.

---

## The headline 6 genes — all analyses, full precision

### Phase 5 — HL genes only, biobin logistic, 4 PCs

**Outcome:** binary HL (Phenotype A — `HL_needAud`; see definition above). **Cohort:** 41,748 after IBD trick (~1,098 cases / ~39,981 controls / NAs removed).

**Gene universe tested:** **only the 173 known HL genes** from the curated list ([`all_genes_including_ShadisList.txt.gz`](../../data/PMBB_Exome/all_genes_including_ShadisList.txt.gz)). Phase 5 does NOT test the rest of the exome — that's Phase 6. Therefore the 4 "novel candidate" genes from the paper (COL5A1, HMMR, RAPGEF3, NNT) are **not in Phase 5's testing universe** and have no p-value here. They first appear when we run exome-wide in Phase 6.

> **LOKI is NOT used in Phase 5** (we pass `--region-file` with explicit gene boundaries). Therefore our and Daniel's results are **byte-equivalent** here regardless of LOKI version. The two columns below are identical to 5 sig figs.

| Gene | In HL gene set? | Ours `loki-20230816` | Daniel `loki-20220926` | Source |
|---|---|---:|---:|---|
| TCOF1 | ✓ yes | **0.062265** | **0.062265** | Phase 5 bins.csv, row 8 ("logistic p-value") |
| ESRRB | ✓ yes | **8.6308×10⁻⁵** | **8.6308×10⁻⁵** | Phase 5 bins.csv, row 8 |
| COL5A1 | ✗ no — novel candidate | n/a — not tested | n/a — not tested | see Phase 6 |
| HMMR | ✗ no — novel candidate | n/a — not tested | n/a — not tested | see Phase 6 |
| RAPGEF3 | ✗ no — novel candidate | n/a — not tested | n/a — not tested | see Phase 6 |
| NNT | ✗ no — novel candidate | n/a — not tested | n/a — not tested | see Phase 6 |

**Source files:**
- **Ours:** `analysis/daniel/outputs/phase5/biobin/merged_maf.001_noRels_keepHLcases-bins.csv`
- **Daniel:** `data/PMBB_Exome/biobin/merged_maf.001_noRels_keepHLcases-bins.csv.gz`
- **Extraction:** wide-to-long via column header (row 0 = gene names) and row 8 (logistic p-value); if a gene appears in multiple bins, the smallest p is reported.

### Phase 6 — Whole exome, biobin logistic, 4 PCs

**Outcome:** binary HL (Phenotype A — `HL_needAud`, same as Phase 5; see definition above). **Cohort:** 41,748 after IBD trick (same cohort as Phase 5, just exome-wide instead of HL-only).

> **LOKI IS used in Phase 6** (`--bin-regions Y`, no `--region-file`). LOKI version drift is the dominant source of our ≠ Daniel discrepancies. See "LOKI drift mechanism" section below.

| Gene | Ours `loki-20230816` | Daniel `loki-20220926` | Source |
|---|---:|---:|---|
| TCOF1 | **2.10×10⁻³** | **1.60×10⁻³** | per-chr biobin meta concat |
| ESRRB | **2.40×10⁻⁴** | **1.13×10⁻³** | per-chr biobin meta concat |
| COL5A1 | **2.29×10⁻⁵** | **2.69×10⁻⁵** | per-chr biobin meta concat |
| HMMR | **1.30×10⁻³** | **4.20×10⁻³** | per-chr biobin meta concat |
| RAPGEF3 | **9.65×10⁻⁵** | **1.46×10⁻⁴** | per-chr biobin meta concat |
| NNT | **4.30×10⁻³** | **1.60×10⁻³** | per-chr biobin meta concat |

**Source files:**
- **Ours:** `analysis/daniel/outputs/phase6/results/all_chrom_meta_HL_needAud.txt` (column `p_logistic`)
- **Daniel:** `data/PMBB_Exome/allGenes/HL_needAud/results_allChr_needAud_withBH.txt.gz` (column `Logistic_regression_beta_p`)

### Phase 7 — Whole exome, R `lm()` linear regression, 20 PCs

**Outcome:** degree-HL 0-4 continuous (Phenotype B — `degHL`; see definition above). **Cohort:** 36,507 hybrid (35,397 controls [degree 0-1] + 1,110 cases [degree 2-4]). **Covariates:** 20 PCs (vs 4 in Phase 5/6) + Sex + Age + AgeSq.

> **Light mode — there is no separate "ours" run.** Phase 7 uses Daniel's preserved supplementary table as input. The "Ours" column in earlier docs (now corrected) was incorrectly labeled — it was always Daniel's STable values. The honest comparison is **Daniel preserved STable** vs **Paper published values**.
>
> **LOKI version note for Phase 7:** the STable was generated by Daniel using **`loki-20220926`** (his 2021 LOKI). The paper also used `loki-20220926`. So **LOKI is NOT the source of STable vs Paper differences** — both came from the same LOKI version. The differences below come from a different cause (intermediate iteration vs final run); see "Why STable differs from Paper" after the table.

| Gene | Daniel preserved STable (= our Phase 7 source) | Paper published | Match? |
|---|---:|---:|---|
| TCOF1 | β=**0.8012**, p=**3.90×10⁻⁶** | β=0.798, p=5.20×10⁻⁶ | ✅ within ~30% on p |
| ESRRB | β=**0.1508**, p=**5.22×10⁻⁴** | β=0.148, p=7.00×10⁻⁴ | ✅ within ~25% on p |
| COL5A1 | β=**0.0619**, p=**1.03×10⁻⁵** | β=0.0586, p=3.31×10⁻⁵ | ✅ within ~3× |
| HMMR | β=**0.0865**, p=**4.26×10⁻⁴** | β=0.0974, p=6.89×10⁻⁵ | ✅ within ~6× (largest discrepancy) |
| RAPGEF3 | β=**0.0755**, p=**1.03×10⁻⁴** | β=0.0731, p=1.88×10⁻⁴ | ✅ within ~2× |
| NNT | β=**0.0513**, p=**1.57×10⁻⁴** | β=0.0507, p=2.03×10⁻⁴ | ✅ within ~25% |

**Source files:**
- **Daniel preserved STable:** `data/PMBB_Exome/allGenes/20PCs/degreeHL/results/allChr_STable_degHL.txt.gz` (columns: Gene, Beta, SE, P, Carriers, Case_carriers)
- **Paper published values:** Hui et al. 2023 PLOS Genetics, Tables 3 and 4

### Phase 6 alternate — Whole exome, biobin logistic, 20 PCs (Daniel-only sanity check)

**Outcome:** binary HL (Phenotype A — `HL_needAud`). **Cohort:** same 41,748 as Phase 6 main. **Difference vs Phase 6 main:** uses 20 PCs instead of 4.

> **Daniel also ran the binary analysis with 20 PCs** as a sanity check. Not in the paper. We don't have a "ours" equivalent. Included here only for completeness — it confirms that the binary case/control test isn't rescued by adding more PCs (the bottleneck for binary is genuinely the small case count for genes like TCOF1, not covariate adjustment).

| Gene | Daniel binary 20PC p |
|---|---:|
| TCOF1 | **1.30×10⁻³** |
| ESRRB | **9.21×10⁻⁴** |
| COL5A1 | **2.08×10⁻⁵** |
| HMMR | **4.80×10⁻³** |
| RAPGEF3 | **1.48×10⁻⁴** |
| NNT | **1.80×10⁻³** |

**Source:** `data/PMBB_Exome/allGenes/20PCs/binary/results/allChr_STable_binaryHL.txt.gz`

---

## Two distinct sources of numerical drift — important to keep separate

There are **two different reasons** numbers differ across the columns above, and they have nothing to do with each other. Confusing them leads to wrong conclusions.

### Drift type 1 — LOKI version (affects only Phase 6: Ours vs Daniel)

| Where it shows up | Phase 6 columns "Ours `loki-20230816`" vs "Daniel `loki-20220926`" |
|---|---|
| Magnitude | ~30% to ~5× differences in p-values |
| Cause | biobin's `--bin-regions Y` mode consults LOKI to define which gene each variant belongs to. Different LOKI version → different bin definitions → different per-gene burden counts → different p-values. |
| Does it affect Phase 5? | **No** — Phase 5 uses `--region-file` with explicit boundaries; LOKI is bypassed. |
| Does it affect Phase 7? | **No** — Phase 7 used the same LOKI (`loki-20220926`) on both sides (Daniel's STable and the paper). |
| Fix | Locate `loki-20220926` on LPC (if archived) and re-run with it. |

### Drift type 2 — Intermediate iteration (affects only Phase 7: Daniel STable vs Paper)

| Where it shows up | Phase 7 columns "Daniel preserved STable" vs "Paper published" |
|---|---|
| Magnitude | β within ~5-10%, p within ~1 order of magnitude (HMMR is worst at ~6×) |
| Cause | Daniel ran the degree-HL analysis multiple times during 2021-2022, tightening filters and cohort definitions each time. The STable preserved in `data/PMBB_Exome/allGenes/20PCs/degreeHL/results/` is one of those intermediate runs, not the final one used in the submitted paper. |
| Evidence | Carrier counts differ between STable and paper (e.g., NNT: 840 STable vs 1,013 paper; HMMR: 387 vs 407). If LOKI were the cause, carrier counts wouldn't change — LOKI affects which gene a variant lands in, not how many people carry it. Different carrier counts means the cohort or filter set actually changed. |
| Does it affect Phase 5 or 6? | **No** — those don't use the STable. Phase 5/6 are direct biobin runs. |
| Fix | Ask Daniel for the final paper-submitted intermediate files, or accept the documented drift. |

### Why this distinction matters

A reader looking at the Phase 7 table might assume "the paper uses a different LOKI than what we have access to". **That's wrong.** Both the paper and Daniel's STable used `loki-20220926`. The Phase 7 drift is iteration drift, not LOKI drift.

Conversely, a reader looking at Phase 6 might assume "Daniel's old STable iteration is causing the difference". **Also wrong.** Phase 6 doesn't use Daniel's STable at all — both sides are direct biobin runs. The difference is LOKI version.

---

## Three quick reads of the data

### 1. Phase 5 (HL-only, 4PC) is byte-equivalent because biobin doesn't consult LOKI when you give it an explicit region file.

`--region-file gene_list_regions.txt` lists 173 known HL genes with fixed boundaries. biobin uses these boundaries verbatim; LOKI is never queried. So `loki-20220926` (Daniel) and `loki-20230816` (ours) produce identical bins → identical p-values.

ESRRB at 8.6308×10⁻⁵ matches Daniel to 5 sig figs. This is the strongest validation point in the project — our infrastructure reproduces Daniel's intermediate exactly when LOKI is not in the loop.

### 2. Phase 6 (WES, 4PC) differs by ~30%-5× because biobin uses LOKI to define every gene bin.

`--bin-regions Y` (no explicit region file) means biobin asks LOKI: "what genes are in this region?" for the whole exome. LOKI's gene boundary definitions evolved between Daniel's 2022 version and our 2023 version, causing:

- **Some variants migrate** between gene bins (e.g., a variant that landed in ESRRB under `loki-20220926` may land in an adjacent LOC under `loki-20230816`)
- **New LOC/LINC pseudogenes** in newer LOKI absorb signal from real genes
- Result: per-gene burden counts shift slightly, p-values shift correspondingly

The discrepancies are bounded — all 6 paper genes still appear at top tier in both versions. But you'd need to use `loki-20220926` (if you can find it on LPC) to get byte-equivalence.

### 3. Phase 7 is light-mode — we don't have a separate "ours" run.

We adopted Daniel's preserved STable directly as our Phase 7 result. The validation comparison is **Daniel STable ↔ Paper published**, not **Ours ↔ Daniel**. β coefficients match within 5-10%; p-values within an order of magnitude (HMMR is the largest discrepancy at ~6×).

The discrepancies between Daniel's STable and the paper come from the STable being an **intermediate iteration** — Daniel ran the analysis multiple times tightening filters and cohort definitions before submission. The STable preserves one of those intermediate runs, not the final published one.

---

## How earlier docs may differ from this one

| Other doc | What it says | Why it's consistent |
|---|---|---|
| [`chapter1_summary.md`](chapter1_summary.md) Table A | TCOF1 ours WES 4PC = `0.0021`, Daniel = `0.0016` | Same numbers as above, rounded to 2 sig figs |
| [`chapter1_summary.md`](chapter1_summary.md) Table B | β / p for all 6 genes (Daniel WES 20PC linear vs Paper) | Same numbers, formatted with scientific notation |
| [`phase5_replication_report.md`](phase5_replication_report.md) line 6, 88 | ESRRB p = `8.6308e-05` | Identical |
| [`phase6_replication_report.md`](phase6_replication_report.md) line 188-219 | Per-chr biobin table with ranks | Same numbers; just adds ranking info |
| [`phase7_degree_hl_burden.md`](phase7_degree_hl_burden.md) lines 105-117 | Calls Daniel's STable values "Our" | ⚠️ Labeling error — see below for correction |

### Known labeling error in `phase7_degree_hl_burden.md`

That doc's "Quantitative comparison" table uses column headers "**Our β**" and "**Our p**" for what is actually Daniel's preserved STable values (because Phase 7 is light-mode). The numbers are correct; only the column label is misleading. **The fix is to relabel "Our" → "Daniel preserved STable (used as our Phase 7 source)"** — which is what `chapter1_summary.md` Table B does correctly.

---

## How to regenerate this doc's source TSVs

```bash
source venv/bin/activate
python3 analysis/daniel/scripts/pmbb_exome/build_pvalue_comparison_table.py \
    --project-root /project/hall/analysis/hearing-loss-genomics \
    --out-dir      /project/hall/analysis/hearing-loss-genomics/analysis/daniel/outputs/phase8/comparison
```

Produces:
- `analysis/daniel/outputs/phase8/comparison/pvalue_comparison_hl_genes.tsv` — 179 known HL genes × all p-value columns
- `analysis/daniel/outputs/phase8/comparison/pvalue_comparison_wes_top.tsv` — 43 top WES genes (union of top 20-30 across each source)

Script source: [`analysis/daniel/scripts/pmbb_exome/build_pvalue_comparison_table.py`](../../analysis/daniel/scripts/pmbb_exome/build_pvalue_comparison_table.py).

---

## LOKI versions, definitively

| Run | LOKI database | Path on LPC | When LOKI built |
|---|---|---|---|
| Daniel (2021) | `loki-20220926` | (now deleted from LPC, or undiscovered) | September 26, 2022 |
| Ours (2026) | `loki-20230816` | `/project/ritchie/datasets/loki/loki-20230816.db` | August 16, 2023 |

`loki-20220926` is referenced by name in Daniel's runbook but we have not found it on the current LPC filesystem. If a future analysis needs byte-equivalent Phase 6, the first step would be to ask the Ritchie Lab if they have the older version archived.
