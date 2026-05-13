# Pipeline Walkthrough — Hui et al. 2023 (PMBB v2)

> Curated map of Daniel Hui's 1,262-line runbook ([`analysis/daniel/runbook_hui2023.txt`](../analysis/daniel/runbook_hui2023.txt)) into ordered phases. Each phase lists: cookbook line range, what it does, scripts invoked, inputs, outputs, and gotchas. Use this to know *where* in the runbook to look — then go to the runbook itself for exact commands.
>
> Cookbook line numbers refer to [`analysis/daniel/runbook_hui2023.txt`](../analysis/daniel/runbook_hui2023.txt) (decompressed from `data/PMBB_Exome/README.gz`). All scripts under [`analysis/daniel/scripts/`](../analysis/daniel/scripts/).
>
> **Scope:** Phase 1 of the project replicates this pipeline on **PMBB v2** (the release Daniel used). Raw v2 data is at [`data/pmbb_v2/`](../data/pmbb_v2/) → `/static/PMBB/PMBB-Release-2020-2.0/`. Daniel's intermediates are in [`data/PMBB_Exome/`](../data/PMBB_Exome/) etc.
>
> **Path substitution:** Daniel's runbook references `/project/PMBB/PMBB-Release-2020-2.0/...`. Replace with `data/pmbb_v2/...` or `/static/PMBB/PMBB-Release-2020-2.0/...` — the old `/project/PMBB/` path no longer exists on LPC.
>
> **Future-phase notes** at the end of each phase (labeled "v3+ porting notes") flag what will need attention when we later port the validated pipeline to PMBB v3 or v4. Ignore those during Phase 1.
>
> **PMBB v2 annotation file schema migration (discovered during Phase 1 re-run, 2026-05-12).** Daniel's runbook references two PMBB v2 annotation files at different points: the older `variant-annotations.txt` (used in early steps; **no longer exists on disk** as of 2026) and the newer `variant-annotation-counts.txt` (used in later steps; the only file currently present). Daniel flagged this migration at runbook line 42-43. Header column 1 differs (`Constant_ID` in the old file → `ID` in the new). For our re-runs, **use only `variant-annotation-counts.txt`** — Daniel's existing intermediates may have been built from the older file, so byte-diff against them will fail even when the pipeline is correct. Validate by set equality of variant IDs and genes instead. See [`results/phase1/phase1_replication_report.md`](../results/phase1/phase1_replication_report.md) for the detailed Phase 1 finding.

## Phase map (TL;DR)

| Phase | Cookbook lines | Topic |
|---|---|---|
| 1 | 1–46 | Build the 173-gene HL gene set + filter variant annotations |
| 2 | 47–62 | Reconcile SNP IDs between annotation and pVCF |
| 3 | 63–75 | Per-chromosome plink extraction + 3rd-degree IBD filter + merge |
| 4 | 77–93 | First case/control + covariates + biobin run |
| 5 | 96–100 | Ancestry split (EUR/AFR) biobin runs |
| 6 | 103–138 | New audbase; "keep HL cases despite relatedness" trick |
| 7 | 139–172 | Re-run biobin with `keepHLcases` cohort |
| 8 | 174–260 | **ZNF175 deep-dive** — the priority gene |
| 9 | 261–305 | Add back multiallelic + stoploss; ESRRB deep-dive |
| 10 | 307–345 | Sibling phenotypes (tinnitus, Meniere's) + 4 case/control definitions |
| 11 | 346–402 | Joe's 6 ZNF175 named individuals; matched-control sampling |
| 12 | 404–500 | **Exome-wide all-genes burden** (paper Fig 3) |
| 13 | 543–605 | Meta-analysis EUR × AFR per phenotype |
| 14 | 605–670 | Per-gene deep-dives (TCOF1, ADGRV1, ESRRB); FDR/BH correction; carrier counts |
| 15 | 672–728 | **UKBB replication** (200K exomes) |
| 16 | 764–840 | ClinVar pathogenic/likely-pathogenic carrier analysis |
| 17 | 782–875 | 20-PC covariate update + degHL (continuous) re-runs |
| 18 | 867–924 | Alternative regression models: Poisson, quasipoisson, negative binomial, randomized-pheno controls |
| 19 | 928–end | UKBB hearing-aid user phenotype; PRS/imputed parts (see also `data/PMBB_Imputed/`) |

---

## Phase 1 — Gene list curation + variant annotation filter (lines 1–46)

**Goal:** produce the final master variant×gene table [`annot_genes_full_funcToInclude.txt.gz`](../data/PMBB_Exome/annot_genes_full_funcToInclude.txt.gz), restricted to (a) the 173 HL gene set and (b) the variant functional classes used in the burden test (pLoF + missense REVEL>0.6).

### Re-run modes for Phase 1

Daniel's intermediates from this phase **already exist** in [`data/PMBB_Exome/`](../data/PMBB_Exome/):

| File | Size | Origin |
|---|---|---|
| [`all_genes_including_ShadisList.txt.gz`](../data/PMBB_Exome/all_genes_including_ShadisList.txt.gz) | 512 B | Step 1.1 output — the 173-gene HL list |
| [`annot_genes_full.txt.gz`](../data/PMBB_Exome/annot_genes_full.txt.gz) | 12 MB | Step 1.2 output — annotation table filtered to HL genes |
| [`annot_genes_full_funcToInclude.txt.gz`](../data/PMBB_Exome/annot_genes_full_funcToInclude.txt.gz) | 1.8 MB | Step 1.3 output — **master variant×gene table** used by every downstream phase |

**Light mode** (recommended for initial validation): skip Phase 1 entirely and use `annot_genes_full_funcToInclude.txt.gz` directly. Confirms we can read Daniel's outputs and matches his exact filter decisions.

**Heavy mode** (for end-to-end replication): regenerate all three files from raw v2 data. Confirms our infrastructure reproduces the same filter outputs Daniel got. Useful as a stepping stone toward the future v3/v4 port (where this phase WILL need to be re-run).

The remainder of this Phase 1 section documents **heavy mode**.

### Step 1.1 — Build the HL gene list (lines 1–24)

Three input sources are merged:

| Source | File | Origin |
|---|---|---|
| `autosomal_dominant.txt` + `autosomal_recessive.txt` | manually parsed from [hereditaryhearingloss.org](https://hereditaryhearingloss.org/) | DFNA + DFNB curated gene lists |
| `Doug_list_manual.txt` | Doug Epstein's email | additional candidates from Doug |
| `Hearing_loss_genes_justGenes.txt` | Shadi's prior paper | DFNA gene list from `data/DFNA/` |

Daniel concatenates these and dedupes:
```bash
cut -f 2 all_genes.txt | sort | uniq > justGenes_notFormatted.txt
# manual cleanup of notes/whitespace
sort justGenes_notFormatted.txt | uniq > justGenes_formatted.txt
cat Hearing_loss_genes_justGenes.txt justGenes_formatted_rename.txt | sort | uniq > all_genes_including_ShadisList.txt
```

**Gotcha — gene name reconciliation.** Several genes need renaming because the PMBB variant-annotation table uses different aliases (cookbook line 11-22, 28-37):

| Use this | Not this |
|---|---|
| LRTOMT | COMT2 |
| DIABLO | SMAC |
| ADGRV1 | VLGR1, GPR98 |
| USH1G | SANS |
| CDC14A | DFNB32 |
| RIPOR2 | FAM65B |
| KCNQ1 | KNCQ1 (typo) |
| GSDME | DFNA5 |
| PJVK | DFNB59 |
| MIRN96 | MIR96 |

X-linked genes excluded from autosomal burden: COL4A5, NDP (and MIRN96, which lives at an unconventional position; see runbook line 21).

### Step 1.2 — Filter PMBB variant annotation to HL genes (line 25)

Daniel's original command (runbook line 25):
```bash
python scripts/only_HL_genes.py \
    all_genes_including_ShadisList.txt \
    /project/PMBB/PMBB-Release-2020-2.0/Exome/Variant_annotations/PMBB-Release-2020-2.0_genetic_exome_variant-annotations.txt \
    > annot_genes_full.txt
```

**To re-run today** (with current paths and the still-existing annotation file):
```bash
zcat data/PMBB_Exome/all_genes_including_ShadisList.txt.gz > /tmp/genes.txt
python analysis/daniel/scripts/pmbb_exome/only_HL_genes.py \
    /tmp/genes.txt \
    data/pmbb_v2/Exome/Variant_annotations/PMBB-Release-2020-2.0_genetic_exome_variant-annotation-counts.txt \
    > analysis/daniel/outputs/phase1/annot_genes_full.txt
```

Two notable substitutions vs Daniel's original:
- Path remap: `/project/PMBB/...` → `data/pmbb_v2/...`
- **Filename:** the original `variant-annotations.txt` no longer exists in v2 — only `variant-annotation-counts.txt` (5.4 GB). Daniel himself migrated to the `-counts` version mid-runbook (line 42-43: "were replaced by these"). The two files have the same column structure, so the script works unchanged.

- **Script:** [`only_HL_genes.py`](../analysis/daniel/scripts/pmbb_exome/only_HL_genes.py) — reads column 8 (`line[7]`, 0-indexed) = `Gene.refGene`. Verified column position in the current `-counts.txt` header is still column 8 ✓.
- **Output:** [`annot_genes_full.txt.gz`](../data/PMBB_Exome/annot_genes_full.txt.gz) (Daniel's, 12 MB compressed) — all annotated variants in HL genes

QA check: any HL gene names NOT in the annotation get written to `not_in_annot_full.txt` and reviewed manually.

### Step 1.3 — Filter to functional classes used in burden test (line 46)

Daniel's command:
```bash
python scripts/only_func_cats_to_include.py annot_genes_full.txt \
    > annot_genes_full_funcToInclude.txt
```

**To re-run today:**
```bash
python analysis/daniel/scripts/pmbb_exome/only_func_cats_to_include.py \
    analysis/daniel/outputs/phase1/annot_genes_full.txt \
    > analysis/daniel/outputs/phase1/annot_genes_full_funcToInclude.txt
```

- **Script:** [`only_func_cats_to_include.py`](../analysis/daniel/scripts/pmbb_exome/only_func_cats_to_include.py)
- **What it keeps:** pLoFs (frameshift_*, stopgain, splicing) + missense with **REVEL > 0.6** (paper's burden criteria)
- **Output:** [`annot_genes_full_funcToInclude.txt.gz`](../data/PMBB_Exome/annot_genes_full_funcToInclude.txt.gz) (Daniel's, 1.8 MB) — the **master variant×gene table** used by every downstream phase
- **Sanity check** for our re-run: `diff <(zcat data/PMBB_Exome/annot_genes_full_funcToInclude.txt.gz) analysis/daniel/outputs/phase1/annot_genes_full_funcToInclude.txt` should be empty. If not, schema drift between Daniel's annotation file and the current `-counts.txt` — investigate.
- **Sibling scripts** (used in later phases): [`only_func_cats_to_include_UKBB.py`](../analysis/daniel/scripts/pmbb_exome/only_func_cats_to_include_UKBB.py) for UKBB Phase 15; [`only_func_cats_to_include_stoploss_newannot.py`](../analysis/daniel/scripts/pmbb_exome/only_func_cats_to_include_stoploss_newannot.py) for the multiallelic+stoploss add-back (Phase 9)

### v3+ porting notes (Phase 1) — future phase, ignore in Phase 1

- The PMBB v3 release at `/static/PMBB/PMBB-Release-2024-3.0/Exome/Variant_annotations/` will have a different (possibly more recent) annotation schema — column positions may have shifted. Inspect the v3 annotation header before running [`only_HL_genes.py`](../analysis/daniel/scripts/pmbb_exome/only_HL_genes.py) — the script's hard-coded `cut -f 8` may need updating.
- Gene name reconciliation list (Step 1.1 gotcha) should be re-applied — gene symbol changes between 2021 and 2024 may have happened (HGNC updates).
- REVEL threshold (0.6) and pLoF class definitions are paper-level decisions, not data-dependent — keep these constant.

---

## Phase 2 — SNP ID reconciliation (lines 47–62)

**Goal:** the variant annotation file and the pVCF use different SNP ID schemes — produce a lookup table to map between them, then a plink-compatible `--extract` file with 9,667 IDs.

**Replicated 2026-05-12** — see [`results/phase2/phase2_replication_report.md`](../results/phase2/phase2_replication_report.md). Both validation checks PASSED; `.extract` is bit-identical to Daniel's (`md5sum 5e80ebc0faa5e68277cfeb948af8b1da` when sorted unique).

### Re-run modes for Phase 2

| Mode | What it does | Cost | What it validates |
|---|---|---|---|
| **Light** (recommended) | Reuses Daniel's `data/PMBB_Exome/vcf_SNP_IDs/vcf_SNP_IDs_allchr.txt.gz` (77 MB) and runs only the join + filter + extract | ~20 s wall | Phase 1 → Phase 2 chain; the Python join logic + filter steps |
| **Heavy** | Re-runs `zgrep \| cut` on the raw v2 pVCFs (**781 GB total** — chr1 alone is 79 GB) via a 22-task LSF job array, then concatenates | ~15-20 min wall, requires job-array submission | Adds: the mechanical extract from raw pVCFs (low schema risk — pure VCF-spec cols 1-5) |

Light is the right call for Phase 1 v2 replication. Heavy becomes necessary only when porting to PMBB v3 (different pVCFs).

### Steps

1. **(Heavy only) Extract SNP IDs from each chr pVCF** (line 49):
   ```bash
   for i in {1..22}; do
     bsub "zgrep -v '#' /project/PMBB/PMBB-Release-2020-2.0/Exome/pVCF/GL_by_chrom/PMBB-Release-2020-2.0_genetic_exome_chr${i}_GL.vcf.gz \
       | cut -f 1-5 > vcf_SNP_IDs/vcf_SNP_IDs_chr${i}.txt"
   done
   ```
   In light mode, skip — use [`data/PMBB_Exome/vcf_SNP_IDs/vcf_SNP_IDs_allchr.txt.gz`](../data/PMBB_Exome/vcf_SNP_IDs/) directly.

2. **Join annotation to pVCF IDs** (line 53):
   ```bash
   python analysis/daniel/scripts/pmbb_exome/annot_IDs_vs_pVCF.py \
       <master list from Phase 1> \
       <allchr VCF IDs> \
       | sort -gk1,1 -gk2,2 > matched_snp_IDs_annot_pVCF.txt
   ```
   Script: [`annot_IDs_vs_pVCF.py`](../analysis/daniel/scripts/pmbb_exome/annot_IDs_vs_pVCF.py). Keys by `(chr, pos)`; collapses multi-allelic-in-annotation entries (11,661 → 11,236 unique positions). Outputs NA for variants in annotation but absent from pVCF.

3. **Drop NA + multi-allelic** (lines 55-57). 11,236 → 9,668 rows. Sanity check: `awk '$4 != $7'` (VCF ref vs annot ref) — expect 0 mismatches.

4. **Produce `.extract` file** (line 60):
   ```bash
   tail -n +2 matched_snp_IDs_annot_pVCF_noNA_noMultiallelic.txt | cut -f 3 \
       > matched_snp_IDs_annot_pVCF_noNA_noMultiallelic.extract
   ```
   **Important:** Daniel's runbook says just `cut -f 3 ...`, but his actual `.extract.gz` on disk has 9,667 lines (no header) — he must have stripped it elsewhere undocumentedly. `plink --extract` requires one SNP ID per line **with no header**; keeping `VCF_ID` as the first line makes plink hunt for a SNP literally named "VCF_ID" and silently drop the search. **Always strip the header.**

5. **Final attrition through Phase 2:** Phase 1 master 11,661 → 11,236 (chr+pos dedup) → 9,668 (NA + multi-allelic filter) → 9,667 (`.extract`, header stripped). 17% loss total; 1,173 multi-allelic come back later via Phase 9.

### Script trap encountered (relevant to future phase scripts)

When validating outputs with `diff <(...) <(...) | wc -l` under `set -euo pipefail`: if `diff` finds differences (returns 1), pipefail propagates the non-zero and `set -e` aborts the script before the `if [[ "$diff_count" -eq 0 ]]` branch can run. **Pattern for future phase scripts:** wrap validation blocks in `set +e` / `set -e` toggles, OR use `{ diff || true; } | wc -l`. See Phase 2 report Issue 2 for the full story.

### v3+ porting notes (Phase 2) — future phase, ignore in Phase 1

- Heavy mode mandatory — v3 pVCFs are different files.
- The "add back multi-allelic + stoploss" later add-back (Phase 9) was needed because the initial filter was too aggressive. For v3, decide upfront whether to include multi-allelic from the start.

---

## Phase 3 — plink genotype extraction + IBD filter (lines 63–75)

**Goal:** for each chromosome, extract just the burden-eligible SNPs from the pVCF, filter out 3rd-degree-related individuals (except HL cases — see Phase 6 trick), restrict to MAF < 0.001, and merge into a single cohort-wide bed/bim/fam.

### Steps

1. **Per-chr plink extraction** (line 62):
   ```bash
   for i in {1..22}; do
     bsub "plink --vcf .../PMBB-Release-2020-2.0_genetic_exome_chr${i}_GL.vcf.gz \
         --vcf-half-call m --extract matched_snp_IDs_annot_pVCF_noNA_noMultiallelic.extract \
         --make-bed --out genotypes/allIndvs_chr${i}"
   done
   ```
   Output: [`data/PMBB_Exome/genotypes/allIndvs_chr{1..22}.bed/bim/fam`](../data/PMBB_Exome/genotypes/).

2. **IBD-based individual filter** (line 69):
   - Source N (chr22): 43,731 individuals in the raw pVCF
   - After 3rd-degree-related filter: 41,757 (file: `PMBB-Release-2020-2.0_genetic_exome_3rd_degree_unrelated`)
   - Format as plink keep-list: `awk '{print $0, $0}' ... > PMBB_no3rdRels.txt`

3. **Apply keep-list + MAF filter** (line 70-71):
   ```bash
   for i in {1..22}; do
     plink --bfile genotypes/allIndvs_chr${i} \
         --keep PMBB_no3rdRels.txt --max-maf .001 \
         --make-bed --out genotypes/allIndvs_chr${i}_maf.001_noRels
   done
   ```

4. **Merge all chromosomes** (line 73):
   ```bash
   plink --merge-list genotypes/merge-list.txt --make-bed \
       --out genotypes/allIndvs_maf.001_noRels_merged
   ```

### Key files produced

- [`data/PMBB_Exome/PMBB_no3rdRels.txt.gz`](../data/PMBB_Exome/PMBB_no3rdRels.txt.gz) — 41,757 sample keep-list
- `genotypes/allIndvs_chr{1..22}_maf.001_noRels.{bed,bim,fam}` — per-chr filtered
- `genotypes/allIndvs_maf.001_noRels_merged.{bed,bim,fam}` — cohort-wide

### v3+ porting notes (Phase 3) — future phase, ignore in Phase 1

- v3 IBD file location: `/static/PMBB/PMBB-Release-2024-3.0/Exome/IBD/` (check exact filename). Sample N will differ.
- MAF cutoff `.001` is a paper-level constant. Cohort-MAF cutoff `.01` (for later "addBack" runs in Phase 9) is the other constant.

---

## Phase 4 — First case/control + covariates + biobin (lines 77–93)

**Goal:** first end-to-end gene-burden test on the HL gene set, with audiogram-derived cases, on the merged cohort.

### Case/control definition (line 78-81)

```bash
# audiogram column from Brant's audbase
cut -d "," -f 5,16 rgc21_45k_aud_1.csv | sed -e 's/"//g' | sed -e "s/,/ /g" \
    | awk '$2 == 1 {print $1}' > BL_SNHL_audiogram.txt

# phecode 389 (HL) column from PMBB diagnosis file (column 723)
cut -f 1,723 PMBB_Diagnosis_Deidentified_072020_Phecodes.txt > phecode_hl.txt

# join into case/control file
python scripts/case_control.py \
    audbase_feb252021/RGC21_45k_aud_1.csv phecode_hl.txt > cases_control.txt
```

- **Script:** [`scripts/pmbb_exome/case_control.py`](../analysis/daniel/scripts/pmbb_exome/case_control.py)
- **Definition (from script comments):** case if audiogram BL_SNHL = TRUE; control if BL_SNHL = FALSE OR phecode = FALSE; NA if phecode 389 = TRUE but audiogram missing
- **Output:** [`data/PMBB_Exome/cases_control.txt.gz`](../data/PMBB_Exome/cases_control.txt.gz)

### Covariates (line 84-88)

```bash
python scripts/make_covs.py \
    audbase_feb252021/RGC21_45k_aud_1.csv \
    /project/PMBB/PMBB-Release-2020-2.0/Phenotype/PMBB-Release-2020-2.0_phenotype_covariates.txt \
    > covs.txt
```

- **Script:** [`scripts/pmbb_exome/make_covs.py`](../analysis/daniel/scripts/pmbb_exome/make_covs.py)
- **Covariates:** sex, age (audiogram date − birth year), `AgeSq`, PCs
- **Note (line 86):** Binglan's UKBB biobin code at `/project/ritchie07/personal/binglan/hearing_loss/src/ukbb_hl/...` was the template

### Region file (line 89)

```bash
python scripts/make_region_file.py annot_genes_full_funcToInclude.txt \
    | sort -gk1,1 -gk3,3 > gene_list_regions.txt
```

biobin's `--region-file` requires chr/start/end/gene rows; this script produces them from the annotation table.

### biobin run (line 92-93)

```bash
biobin -D ~/group/datasets/loki/loki.db \
    -V genotypes/allIndvs_maf.001_noRels_merged.vcf \
    -p cases_control_aud_and_phecode.txt \
    --covariates covs.txt \
    --bin-regions Y --region-file gene_list_regions.txt \
    -G 38 --test logistic \
    --report-prefix biobin/allIndvs_merged > biobin/allIndvs_merged_log.txt 2>&1
```

- **Tool:** `biobin` (Ritchie Lab — needs `loki.db` at `~/group/datasets/loki/loki.db`)
- **Test:** logistic regression of HL case status ~ aggregated burden per gene + covariates
- **Output:** [`data/PMBB_Exome/biobin/`](../data/PMBB_Exome/biobin/) — `*-bins.csv` (per-individual carrier matrix) + `*-locus.csv` (per-variant) + log

---

## Phase 5 — Ancestry split (lines 96–100)

Split into EUR-only and AFR-only biobin runs for ancestry-stratified analysis. Files: [`EUR.txt.gz`](../data/PMBB_Exome/EUR.txt.gz), [`AFR.txt.gz`](../data/PMBB_Exome/AFR.txt.gz). Re-runs biobin with same params on each subset.

---

## Phase 6 — New audbase + "keep HL cases despite relatedness" (lines 103–138)

**Critical phase — the "IBD trick" that lets HL cases stay in the cohort.**

After running with the strict 3rd-degree-unrelated filter, Daniel realized cases were being filtered out (only 1,087 cases remained from 2,247 audiogram individuals). He relaxes the IBD filter specifically for HL cases:

```bash
# Cases with relatives (from full cohort)
awk '$2 == 1' cases_control.txt | cut -f 1 \
    | grep -w -f - genotypes/allIndvs_chr22.fam > cases_withRels.txt

# Cases that survived the strict filter
awk '$2 == 1' cases_control.txt | cut -f 1 \
    | grep -w -f - genotypes/allIndvs_maf.001_noRels_merged.fam > cases_noRels.txt

# Cases removed by the strict filter
grep -w -v -f cases_noRels.txt cases_withRels.txt > removed_cases.txt

# Add removed cases back IF they're not >3rd-degree to any kept individual
python scripts/keep_HL_cases_IBD.py removed_cases.txt \
    /project/PMBB/PMBB-Release-2020-2.0/Exome/IBD/PMBB-Release-2020-2.0_genetic_exome.genome \
    /project/PMBB/PMBB-Release-2020-2.0/Exome/IBD/PMBB-Release-2020-2.0_genetic_exome_3rd_degree_unrelated \
    > tokeep_moreHLcases.txt
```

- **Script:** [`scripts/pmbb_exome/keep_HL_cases_IBD.py`](../analysis/daniel/scripts/pmbb_exome/keep_HL_cases_IBD.py)
- **Output:** [`tokeep_moreHLcases.txt.gz`](../data/PMBB_Exome/tokeep_moreHLcases.txt.gz) — the **final cohort** used in the paper

**Note (line 129-138):** 3 individuals (PMBB2821699971795, PMBB9129875068625, PMBB9394987694139) have HL family members and were further excluded. Final case N matches paper's reported numbers.

### v3+ porting (Phase 6) — future phase, ignore in Phase 1

The IBD "genome" file and 3rd-degree-unrelated file for v3 may use different formats. This phase is the most likely place for cohort N to drift between v2 and v3.

---

## Phase 7 — Re-run biobin with `keepHLcases` cohort (lines 139–172)

Re-runs Phase 4's biobin on the `tokeep_moreHLcases` cohort. **Top hit: ESRRB** (line 151 — "it's ESRRB, just like the paper"). Same pipeline structure, different keep-list.

Key output: [`data/PMBB_Exome/biobin/merged_maf.001_noRels_keepHLcases-bins.csv`](../data/PMBB_Exome/biobin/).

---

## Phase 8 — ZNF175 deep-dive (lines 174–260) **★ project priority**

**Goal:** identify ZNF175 carriers and the variants they carry.

### Step 8.1 — ZNF175-only annotation + region (lines 174-181)

Same as Phase 1-4 but restricted to `ZNF175/`:
- [`scripts/pmbb_exome/only_HL_genes.py`](../analysis/daniel/scripts/pmbb_exome/only_HL_genes.py) on a single-gene list
- [`scripts/pmbb_exome/only_func_cats_to_include.py`](../analysis/daniel/scripts/pmbb_exome/only_func_cats_to_include.py)
- [`scripts/pmbb_exome/make_region_file.py`](../analysis/daniel/scripts/pmbb_exome/make_region_file.py)
- plink extraction on chr19 only
- biobin with EUR-AFR covariates

### Step 8.2 — Carrier identification (lines 188-202)

**Important** — this is where the "8 carriers" number comes from:

```bash
# Extract ZNF175 region from the full pVCF, including multi-allelic, no MAF filter
plink --bfile .../PMBB-Release-2020-2.0_genetic_exome_GL \
    --extract ZNF175/ZNF175_annot_genes_full_funcToInclude.extract \
    --make-bed --out ZNF175/ZNF175_allIndvs_inclMultAllelic_noMAFfilter
plink --bfile ZNF175/ZNF175_allIndvs_inclMultAllelic_noMAFfilter \
    --recode vcf-iid --out ZNF175/ZNF175_allIndvs_inclMultAllelic_noMAFfilter

# Identify carriers at specific positions
python scripts/ZNF175_carrier.py ZNF175/ZNF175_allIndvs_inclMultAllelic_noMAFfilter.vcf \
    | grep 51587727 | cut -f 1 | grep -w -f - cases_control.txt
# Result: 8 carriers, 1 phecode-case (PMBB2106731298975)

python scripts/ZNF175_carrier.py ZNF175/ZNF175_allIndvs_inclMultAllelic_noMAFfilter.vcf \
    | grep 51581437 | cut -f 1 | grep -w -f - cases_control.txt
# Result: 90 carriers, 2 phecode-cases (PMBB7501686571326, PMBB1245988577461)
```

- **Script:** [`scripts/pmbb_exome/ZNF175_carrier.py`](../analysis/daniel/scripts/pmbb_exome/ZNF175_carrier.py)
- **Positions of interest:** chr19:51587727 (rare, 8 carriers), chr19:51581437 (more common, 90 carriers) — both filtered later because 51581437 was high frequency

### Step 8.3 — Joe's pLoF list (lines 212-258)

Joe Park provided a list of **10 lifted-over positions** (lines 213-233) that should be considered ZNF175 pLoFs. After QA against PMBB v2 annotation, several didn't survive (multi-allelic, not annotated, etc.). The surviving subset goes into [`ZNF175/ZNF175_allIndvs_chr19_Joes.{bed,bim,fam,vcf}`](../data/PMBB_Exome/ZNF175/) and a carrier file [`ZNF175/ZNF175_pLOF_Joes.txt`](../data/PMBB_Exome/ZNF175/) — this is the **canonical "Joe's ZNF175 pLoF list"** carrier set.

### Step 8.4 — Separate carriers into pLoF vs missense (line 257)

```bash
python scripts/missense_or_pLOF.py \
    ZNF175/ZNF175_annot_genes_full_funcToInclude.txt \
    ZNF175/ZNF175_matched_snp_IDs_annot_pVCF.txt \
    ZNF175/ZNF175_carriers_allIndvs.txt | cut -f 1,3 \
    > ZNF175/ZNF175_carriers_allIndvs_wCategory.txt
```

- **Script:** [`scripts/pmbb_exome/missense_or_pLOF.py`](../analysis/daniel/scripts/pmbb_exome/missense_or_pLOF.py)

### v3+ porting + project relevance — future phase, ignore in Phase 1

This is **the** phase for the Phase 1 priority. Two questions to resolve before re-running:
1. **Which "8" do we deep-dive?** Position 51587727 produces 8 carriers but only 1 is phecode-case. The kickoff meeting says "8 signal-driving cases" — most likely interpretation: deep-dive all 8 carriers regardless of phecode status. Confirm with Doug/Daniel.
2. **Which pLoF list to use?** Joe's 10 positions filtered down (line 235-246), giving the final "Joe's pLoFs" set. For v3, redo this filter against v3 annotation.

---

## Phase 9 — Add back multiallelic + stoploss (lines 261–305)

Initially the filter dropped multi-allelic and stoploss variants. Daniel adds them back and re-runs the burden test. New workspace: [`data/PMBB_Exome/addBack_multiallelic_stoploss/`](../data/PMBB_Exome/addBack_multiallelic_stoploss/).

Two MAF cutoffs explored: `.001` (paper rare) and `.01` (more inclusive). The `.01` version uses **logistic regression** for the binary HL phenotype.

**ESRRB pLoF-only carriers** identified in [`addBack_multiallelic_stoploss/ESRRB/`](../data/PMBB_Exome/addBack_multiallelic_stoploss/ESRRB/). Same workflow as ZNF175.

---

## Phase 10 — Sibling phenotypes + 4 case/control definitions (lines 307–345)

Daniel runs the burden test against **4 phenotype variants** as a specificity check:
- `caseAudAndPhecode` — case if BOTH phecode AND audiogram (strict, smaller N)
- `needAud` — case if phecode AND (audiogram=HL OR audiogram missing); base from Phase 4
- `dontNeedAud` — case if phecode positive, audiogram ignored (largest N)
- `rmAudNA` — case if phecode AND audiogram both available and agree (purest, smallest N)

Also runs negative-control phenotypes: tinnitus (phecode column 727) and Meniere's disease (column 717). These should NOT show ZNF175 / ESRRB enrichment.

Scripts:
- [`case_control_allowOnlyPhecode.py`](../analysis/daniel/scripts/pmbb_exome/case_control_allowOnlyPhecode.py)
- [`case_control_allowOnlyPhecode_rmAudNA.py`](../analysis/daniel/scripts/pmbb_exome/case_control_allowOnlyPhecode_rmAudNA.py)
- [`HL_case_aud_and_phecode.py`](../analysis/daniel/scripts/pmbb_exome/HL_case_aud_and_phecode.py)

Outputs: `renamed_pheno_files_for_looping/cases_control_{needAud,dontNeedAud,caseAudAndPhecode,rmAudNA}.txt`.

---

## Phase 11 — Joe's 6 ZNF175 named individuals (lines 346–402)

Joe Park emails a list of **6 specific ZNF175-carrier PMBB IDs** (line 348-353):
- PMBB8949342677388
- PMBB4664277909557
- PMBB2106731298975
- PMBB7083231520332
- PMBB9701760809542
- PMBB9508968070076

Daniel cross-references against:
- `cases_control.txt` (audiogram phenotype) — only PMBB2106731298975 is a phecode-case
- `cases_control_allowOnlyPhecode.txt` (looser) — 5 are cases, 1 is NA (PMBB9701760809542)
- `tinnitus_phecode_case_control.txt` — only PMBB4664277909557 and PMBB9508968070076 are tinnitus cases
- `menieres_phecode_case_control.txt` — all 0/NA, no Meniere's enrichment

This is **also a candidate set for the 8-case deep-dive** — it's 6 here, not 8. Question for Doug/Daniel: does the deep-dive cover these 6 + 2 more, or a different selection entirely?

Then matches **4 controls per carrier** by age/sex/PCs:
```bash
python scripts/get_ctrls_4max.py ZNF175/ZNF175_pLOF_Joes.txt \
    /project/PMBB/.../PMBB-Release-2020-2.0_genetic_exome_ancestries.txt \
    /project/PMBB/.../PMBB-Release-2020-2.0_phenotype_covariates.txt \
    | sort -k2 > ZNF175/ZNF175_matched_controls_4ctrls.txt
```

Script: [`get_ctrls_4max.py`](../analysis/daniel/scripts/pmbb_exome/get_ctrls_4max.py).

---

## Phase 12 — Exome-wide all-genes burden (lines 404–500) **★ Fig 3 of paper**

Repeats Phase 1-7 but **without the HL gene-list filter** — burden across all ~20K genes. Uses MAF<0.01 (less rare than HL-gene run at .001). Output: per-gene biobin results combined into one table per phenotype variant.

Per-chr biobin runs (parallel via `bsub`), then concatenated, then BH-corrected (Benjamini-Hochberg, R `p.adjust`).

Outputs in [`data/PMBB_Exome/allGenes/HL_{needAud,dontNeedAud,rmAudNA,caseAudAndPhecode}/`](../data/PMBB_Exome/allGenes/).

This is the run that **discovered the novel candidates** (HMMR, COL5A2, NNT, RAPGEF3) at FDR < 0.05.

---

## Phase 13 — Meta-analysis EUR × AFR (lines 543–605)

Per-phenotype, per-chromosome, per-ancestry biobin runs, then meta-analysed across EUR + AFR using inverse-variance weighting in [`scripts/pmbb_exome/meta.R`](../analysis/daniel/scripts/pmbb_exome/meta.R) and [`parse_results_to_meta.R`](../analysis/daniel/scripts/pmbb_exome/parse_results_to_meta.R).

Outputs: `allGenes/HL_meta_{phenotype}/...`.

---

## Phase 14 — Per-gene deep-dives + FDR (lines 605–670)

ADGRV1, TCOF1, ESRRB single-gene region extracts. BH/FDR correction on the meta-results. Count-of-carriers per gene.

---

## Phase 15 — UKBB replication (lines 672–728)

UK Biobank 200K exomes. Different annotation file (`UKB_WES_200K_Summary_All_Variants_121520.txt`), different sample-keep file. Same biobin pipeline. Phenotype = self-reported hearing-difficulty fields from UKBB tab.

Cohort: **148,970 individuals (39,272 cases / 109,698 controls)** — matches the paper.

---

## Phase 16 — ClinVar pathogenic carrier analysis (lines 764–840)

Filter the variant annotation to ClinVar `Pathogenic` / `Likely_pathogenic` / `Pathogenic/Likely_pathogenic`. Re-extract genotypes. Identify carriers and tabulate by inheritance pattern (AD vs AR), zygosity (heterozygous vs homozygous), and HL gene status. Paper numbers: 6.80% controls / 7.93% cases ClinVar carriers (p = 0.147, not significant alone).

Scripts: [`ClinVar_carrier.py`](../analysis/daniel/scripts/pmbb_exome/ClinVar_carrier.py), [`ClinVar_carrier_addGene_HLstatus.py`](../analysis/daniel/scripts/pmbb_exome/ClinVar_carrier_addGene_HLstatus.py), [`ClinVar_and_nonClinVar_counts.py`](../analysis/daniel/scripts/pmbb_exome/ClinVar_and_nonClinVar_counts.py).

---

## Phase 17 — 20-PC update + degHL continuous phenotype (lines 782–875)

Switch from 4 PCs to **20 PCs** in covariates. Rerun the all-genes burden with the updated covariates against the degree-of-HL (0–4 ordinal) phenotype.

Two transformations explored:
- Raw degree-HL (0-4)
- "1-baseline" version subtracting baseline
- Log of "1-baseline"

Uses linear regression (`--test linear` in biobin) instead of logistic.

Outputs in [`data/PMBB_Exome/allGenes/20PCs/{binary,degreeHL,degreeHL_1baseline,degreeHL_1baseline_log}/`](../data/PMBB_Exome/allGenes/20PCs/).

---

## Phase 18 — Alternative regression models (lines 867–924)

Daniel sanity-checks the burden signal with three count-model alternatives:
- Poisson ([`run_poisson.R`](../analysis/daniel/scripts/pmbb_exome/run_poisson.R))
- Quasi-poisson ([`run_quasipoisson.R`](../analysis/daniel/scripts/pmbb_exome/run_quasipoisson.R))
- Negative binomial ([`run_negbinom.R`](../analysis/daniel/scripts/pmbb_exome/run_negbinom.R))

And **randomized-phenotype controls** for each — shuffle case/control labels and re-run, to confirm signal doesn't appear under null. Outputs in `allGenes/20PCs/degreeHL/loglinear/` and `results_randPheno/`.

---

## Phase 19 — UKBB hearing-aid phenotype + imputed/PRS (lines 928–end)

Alternative UKBB phenotype: self-reported hearing aid use (UKBB field 1522). Same biobin pipeline.

The **imputed-data side** (GWAS, PRS-CS, PRSice, heritability) lives in [`data/PMBB_Imputed/`](../data/PMBB_Imputed/) with its own scripts at [`analysis/daniel/scripts/pmbb_imputed/`](../analysis/daniel/scripts/pmbb_imputed/). Not covered in this walkthrough yet — see [`paper_replication.py`](../analysis/daniel/scripts/pmbb_imputed/paper_replication.py), [`format_ss_PRScs.py`](../analysis/daniel/scripts/pmbb_imputed/format_ss_PRScs.py), and the GCTA / PRS-CS output dirs.

---

## What to do next (Phase 1 — v2 replication)

1. **Resolve the "8 cases" definition** with Doug/Daniel. The four candidate interpretations are:
   - 8 carriers at chr19:51587727 (regardless of phecode)
   - 8 individuals from a different position
   - The 6 named in Joe's email + 2 from another source
   - Aggregate ZNF175 pLoF carriers from the multiallelic+stoploss add-back run
2. **Verify tools are available on LPC:** `biobin` + `loki.db`, `plink` (v1.9), R packages, SLURM (vs. LSF).
3. **Pilot re-run** of Phase 4 (the first biobin run) on the existing intermediates in [`data/PMBB_Exome/biobin/`](../data/PMBB_Exome/biobin/) to confirm we reproduce a known number (e.g., top hit = ESRRB, line 151 of the runbook).
4. **End-to-end replication:** re-run the full pipeline on raw v2 data via [`data/pmbb_v2/`](../data/pmbb_v2/) (path-substituting `/project/PMBB/PMBB-Release-2020-2.0/` → `data/pmbb_v2/` in Daniel's commands). Confirm the published cohort numbers (1,110 cases / 35,397 controls; h² = 4.53%; ZNF175 + ESRRB + TCOF1 hits).

**After Phase 1 is validated** (future phase, not now): port to PMBB v3 (`PMBB-Release-2024-3.0`) and/or v4 (`PMBB-Release-2026-4.0`). Start with Phase 1 gene-list curation against the v3+ annotation schema — most likely place for v2→v3 schema drift.
