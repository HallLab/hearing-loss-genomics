# Data Inventory — `data/`

> First-pass inventory of `/home/alrico/hall/analysis/hearing-loss-genomics/data/`. Identifies what each subdirectory and key file is for, based on file names, READMEs, and the original pipeline notes. Updated 2026-05-12.
>
> Companion docs:
> - [`pipeline_walkthrough.md`](pipeline_walkthrough.md) — the curated phase-by-phase guide to Daniel's runbook
> - [`../analysis/daniel/README.md`](../analysis/daniel/README.md) — the replication workspace where decompressed scripts and runbook live

## Top-level

| Path | Origin | Purpose |
|---|---|---|
| [`PMBB_Exome/`](../data/PMBB_Exome/) | Daniel Hui (Hui et al. 2023) | **Main rare-variant burden workspace.** All pipeline scripts, phenotype/covariate files, gene lists, intermediate plink/VCF outputs, and per-gene deep-dives (ZNF175, TCOF1, ESRRB, ADGRV1). 5.6 GB, 6,258 files. |
| [`PMBB_Imputed/`](../data/PMBB_Imputed/) | Daniel Hui | **GWAS & PRS workspace.** Common-variant GWAS, PRS-CS, PRSice, heritability, paper-SNP replication. 385 GB (mostly imputed genotypes), 1,406 files. |
| [`DFNA/`](../data/DFNA/) | Shadi (separate, prior paper) | Adjacent dataset — Shadi's DFNA analysis. NOT the focus of our project but referenced (HL gene list is partially built from Shadi's list). 910 MB, 188 files. |
| [`pmbb_v3/`](../data/pmbb_v3/) | Andre (Apr 2026) | **Symlinks into the current PMBB v3 release** (`/static/PMBB/PMBB-Release-2024-3.0/`). Three files so far: covariates, PheCode12_long, PheCodeX_long. This is the *new* data we will run the replication against. |
| [`README.gz`](../data/README.gz) | Daniel | One-line legend for the three sibling dirs. |

## PMBB v3 — current release (target of replication)

`pmbb_v3/` contains only **three symlinks** so far — these are the only PMBB v3 files actively pulled into the project. The full release sits at `/static/PMBB/PMBB-Release-2024-3.0/` with subdirs `Exome/`, `Genotype/`, `Imaging/`, `Imputed/`, `Phenotypes/`, `README.txt`.

| Symlink | Target |
|---|---|
| [`pmbb_v3/covariates.txt`](../data/pmbb_v3/covariates.txt) | `/static/PMBB/PMBB-Release-2024-3.0/Phenotypes/3.0/PMBB-Release-2024-3.0_covariates.txt` |
| [`pmbb_v3/phecode12_long.txt`](../data/pmbb_v3/phecode12_long.txt) | `..._phenotype_conditions_PheCode12_long.txt` (ICD9-mapped phecodes, classic version — what Hui et al. used) |
| [`pmbb_v3/phecodeX_long.txt`](../data/pmbb_v3/phecodeX_long.txt) | `..._phenotype_conditions_PheCodeX_long.txt` (PheCodeX extended version) |

**To-do for inventory:** when replication starts, add symlinks for the v3 exome pVCFs, ancestries file, IBD/relatedness file, and audiogram linkage (audiogram data lives outside `/static/PMBB/Phenotypes/`, per kickoff notes — need to confirm where).

## `PMBB_Exome/` — rare-variant burden workspace

### Top-level files (gene lists, phenotype, covariates, annotations)

**Gene lists**
- [`all_genes.txt.gz`](../data/PMBB_Exome/all_genes.txt.gz), [`Hearing_loss_genes.txt.gz`](../data/PMBB_Exome/Hearing_loss_genes.txt.gz), [`Hearing_loss_genes_justGenes.txt.gz`](../data/PMBB_Exome/Hearing_loss_genes_justGenes.txt.gz) — manually curated HL gene set (DFNA + DFNB from hereditaryhearingloss.org + Doug's manual additions)
- [`DFNA_genes.txt.gz`](../data/PMBB_Exome/DFNA_genes.txt.gz), [`DFNB_genes.txt.gz`](../data/PMBB_Exome/DFNB_genes.txt.gz) — split by inheritance
- [`autosomal_dominant.txt.gz`](../data/PMBB_Exome/autosomal_dominant.txt.gz), [`autosomal_recessive.txt.gz`](../data/PMBB_Exome/autosomal_recessive.txt.gz), [`Doug_list_manual.txt.gz`](../data/PMBB_Exome/Doug_list_manual.txt.gz) — sub-lists
- [`all_genes_including_ShadisList.txt.gz`](../data/PMBB_Exome/all_genes_including_ShadisList.txt.gz) — merged with Shadi's DFNA gene list; this is the **173-gene set used in Fig 2 of the paper**

**Phenotypes / case-control definitions**
- [`PMBB_Diagnosis_Deidentified_072020_Phecodes.txt.gz`](../data/PMBB_Exome/PMBB_Diagnosis_Deidentified_072020_Phecodes.txt.gz) — old PMBB v2 diagnosis-phecode mapping (Jul 2020 freeze)
- [`phecode_hl.txt.gz`](../data/PMBB_Exome/phecode_hl.txt.gz) — extracted phecode 389 column for HL
- [`cases_control.txt.gz`](../data/PMBB_Exome/cases_control.txt.gz) — final case/control assignment used in paper
- [`cases_control_HLcaseAudAndPhecode.txt.gz`](../data/PMBB_Exome/cases_control_HLcaseAudAndPhecode.txt.gz) — audiogram + phecode hybrid definition (the "leaky control" fix from Fig 1 of paper)
- [`cases_control_allowOnlyPhecode*.txt.gz`](../data/PMBB_Exome/) — phecode-only definition variants
- [`cases_control_degHL*.txt.gz`](../data/PMBB_Exome/) — 0–4 PTA-bin degree-of-HL phenotype (continuous)
- [`menieres_phecode_case_control.txt.gz`](../data/PMBB_Exome/menieres_phecode_case_control.txt.gz), [`tinnitus_phecode_case_control.txt.gz`](../data/PMBB_Exome/tinnitus_phecode_case_control.txt.gz) — sibling phenotypes (negative-control / specificity tests)
- [`PMBB_no3rdRels.txt.gz`](../data/PMBB_Exome/PMBB_no3rdRels.txt.gz) — 41,757 individuals with no 3rd-degree relatives (post-IBD QC; baseline cohort)
- [`tokeep_moreHLcases.txt.gz`](../data/PMBB_Exome/tokeep_moreHLcases.txt.gz), [`tokeep_moreHLcases_keep.txt.gz`](../data/PMBB_Exome/tokeep_moreHLcases_keep.txt.gz) — relaxed IBD filter that *keeps* HL cases even if related (paper-specific trick — see README.gz line 124-138)
- [`cases_all.txt.gz`](../data/PMBB_Exome/cases_all.txt.gz), [`cases_noRels.txt.gz`](../data/PMBB_Exome/cases_noRels.txt.gz), [`cases_withRels.txt.gz`](../data/PMBB_Exome/cases_withRels.txt.gz), [`kept_cases.txt.gz`](../data/PMBB_Exome/kept_cases.txt.gz), [`removed_cases.txt.gz`](../data/PMBB_Exome/removed_cases.txt.gz) — intermediate cohort files during IBD reconciliation
- [`hearingAid_*` files in DFNA/](../data/DFNA/) — hearing-aid prescription as alternative HL phenotype

**Covariates**
- [`covs.txt.gz`](../data/PMBB_Exome/covs.txt.gz) — base covariates (age, sex, etc.)
- [`covs_withAnc.txt.gz`](../data/PMBB_Exome/covs_withAnc.txt.gz), [`covs_withAnc_onlyEUR-AFR*.txt.gz`](../data/PMBB_Exome/) — with genetic ancestry; EUR/AFR splits
- [`covs_withAnc_onlyEUR-AFR_rmAncColumn_20PCs*.txt.gz`](../data/PMBB_Exome/) — 20-PC version for plink (final paper covariate set)
- [`anc.txt.gz`](../data/PMBB_Exome/anc.txt.gz), [`anc_only.txt.gz`](../data/PMBB_Exome/anc_only.txt.gz) — ancestry assignments (EUR / AFR / etc.)
- [`covs_rand.txt.gz`](../data/PMBB_Exome/covs_rand.txt.gz), [`covs_rand_withbinaryHL.txt.gz`](../data/PMBB_Exome/covs_rand_withbinaryHL.txt.gz) — randomized phenotype controls (paper inflation diagnostic, λ = 1.01)
- [`EUR.txt.gz`](../data/PMBB_Exome/EUR.txt.gz), [`AFR.txt.gz`](../data/PMBB_Exome/AFR.txt.gz), [`EUR_wHeader.txt.gz`](../data/PMBB_Exome/EUR_wHeader.txt.gz), [`AFR_wHeader.txt.gz`](../data/PMBB_Exome/AFR_wHeader.txt.gz) — keep-lists per ancestry

**Variant annotations** (PMBB v2 pVCF + REVEL + ClinVar)
- [`annot_genes_full.txt.gz`](../data/PMBB_Exome/annot_genes_full.txt.gz) — full annotation table joined with paper gene list (12 MB)
- [`annot_genes_full_funcToInclude.txt.gz`](../data/PMBB_Exome/annot_genes_full_funcToInclude.txt.gz) — filtered to: pLoF (frameshift / stopgain / splicing) + missense REVEL>0.6 (the paper's burden criteria)
- [`annot_genes_full_funcToInclude_extract.txt.gz`](../data/PMBB_Exome/annot_genes_full_funcToInclude_extract.txt.gz) — same, plink `--extract` format
- [`annot_genes_someMissing.txt.gz`](../data/PMBB_Exome/annot_genes_someMissing.txt.gz), [`annot_header.txt.gz`](../data/PMBB_Exome/annot_header.txt.gz), [`not_in_annot.txt.gz`](../data/PMBB_Exome/not_in_annot.txt.gz), [`not_in_annot_full.txt.gz`](../data/PMBB_Exome/not_in_annot_full.txt.gz) — QA on the annotation join (genes that were missing → manually checked in OMIM, see README.gz line 11-22)
- [`annot_Func.refGene_ExonicFunc.refGene_uniq.txt.gz`](../data/PMBB_Exome/annot_Func.refGene_ExonicFunc.refGene_uniq.txt.gz) — annotation function categories observed in the pVCF
- [`matched_snp_IDs_annot_pVCF*.txt.gz`](../data/PMBB_Exome/) — SNP IDs reconciled between annotation file and pVCF (some IDs differ; this is the lookup table)
- [`pLOFs_missenseREVEL.60.txt.gz`](../data/PMBB_Exome/pLOFs_missenseREVEL.60.txt.gz) — final variant set used in burden tests
- [`gene_list_regions.txt.gz`](../data/PMBB_Exome/gene_list_regions.txt.gz) — chr:start-end regions per gene for biobin
- [`exome_subset_for_hearing_loss_2021-06-03.csv.gz`](../data/PMBB_Exome/exome_subset_for_hearing_loss_2021-06-03.csv.gz) — snapshot of exome subset (Jun 2021)

**Other**
- [`README.gz`](../data/PMBB_Exome/README.gz) — **THE pipeline cookbook.** 1,262 lines of bash-history-style commands documenting every step from gene-list cleanup → annotation → plink extraction → biobin → ZNF175 deep-dive → manuscript figures. **READ THIS FIRST.**
- [`.RData.gz`](../data/PMBB_Exome/.RData.gz) — saved R workspace (6 MB)
- [`brant audbase 1.7.21.TXT.gz`](../data/PMBB_Exome/) — Brant's audiogram database (Jan 2021 version, 18 MB)
- [`pmbb_hearing_loss_45K_cohort_2021-08-26_PMBB_ID.txt.gz`](../data/PMBB_Exome/pmbb_hearing_loss_45K_cohort_2021-08-26_PMBB_ID.txt.gz) — 45K cohort ID map (Aug 2021)
- [`chr20_original_sites.vcf.gz`](../data/PMBB_Exome/chr20_original_sites.vcf.gz), [`chr21_original_sites.vcf.gz`](../data/PMBB_Exome/) — chromosome-level VCF snapshots (Oct 2021)
- [`highest_p_rand_gene_variables_OR1N2.txt.gz`](../data/PMBB_Exome/highest_p_rand_gene_variables_OR1N2.txt.gz) — randomization-phenotype top-hit diagnostic
- [`tmp.txt.gz`](../data/PMBB_Exome/tmp.txt.gz) — disposable

### Subdirectories

| Subdir | Files | Purpose |
|---|---:|---|
| [`scripts/`](../data/PMBB_Exome/scripts/) | 81 | **All Daniel's analysis scripts.** Python + R. See script catalog below. |
| [`genotypes/`](../data/PMBB_Exome/genotypes/) | 229 | Per-chromosome plink bed/bim/fam + VCF; both pre- and post-MAF-filter (`allIndvs_chr{1..22}*`). Inputs to biobin / burden tests. |
| [`allGenes/`](../data/PMBB_Exome/allGenes/) | 349 | Exome-wide rare-variant burden run (all genes, not just HL gene set). Per-chromosome biobin outputs, regression results, ClinVar-carrier files, meta-analysis with UKBB. |
| [`include_degHL1/`](../data/PMBB_Exome/include_degHL1/) | 349 | Same structure as `allGenes/` but using the *degree-of-HL* (0–4) continuous phenotype. Linear + binary versions. |
| [`addBack_multiallelic_stoploss/`](../data/PMBB_Exome/addBack_multiallelic_stoploss/) | 60 | Re-runs that added back multi-allelic + stoploss variants (initially excluded). Includes ESRRB sub-folder. |
| [`archive/`](../data/PMBB_Exome/archive/) | 200 | **Older versions of per-chromosome plink files** (`allIndvs_chr{N}_maf.001_noRels.*`) — pre-`keepHLcases` IBD trick. Two READMEs from Apr 2021 and Sep 2021. |
| [`vcf_SNP_IDs/`](../data/PMBB_Exome/vcf_SNP_IDs/) | 23 | Per-chr lookup tables: pVCF SNP IDs → chr:pos:ref:alt. Used to reconcile annotation IDs vs pVCF IDs. |
| [`ZNF175/`](../data/PMBB_Exome/ZNF175/) | 55 | **ZNF175 deep-dive.** Chr19 plink + VCF, with-relatedness and no-relatedness versions, plus Joe Park's analysis subfolder (`Joe_analyses/`: variants table + `51587727.txt.gz` — the chr19 position of the key pLoF variant). |
| [`TCOF1/`](../data/PMBB_Exome/TCOF1/) | 18 | TCOF1 deep-dive (replicated known HL gene). pLoF-only carriers, study-individual annotations. |
| [`ADGRV1/`](../data/PMBB_Exome/ADGRV1/) | 7 | ADGRV1 region files (gene list addition — was VLGR1/GPR98 in older nomenclature). |
| [`UKBB_analyses/`](../data/PMBB_Exome/UKBB_analyses/) | 28 | UKBB 200K exome replication. Phenotype = hearing-aid use + self-report. Subfolder `genos/`. |
| [`TWAS_files/`](../data/PMBB_Exome/TWAS_files/) | 2 | Transcriptome-wide association preparation (covs + genes); appears unfinished. |
| [`biobin/`](../data/PMBB_Exome/biobin/) | 11 | First biobin run on the merged-all-chrs cohort. `merged_maf.001_noRels_keepHLcases-bins.csv` (top hit was ESRRB — line 151 of README.gz). |
| [`renamed_pheno_files_for_looping/`](../data/PMBB_Exome/renamed_pheno_files_for_looping/) | 4 | Phenotype files renamed to loop over (audiogram-only / phecode-only / hybrid / no-NA versions). |
| [`audbase_feb252021/`](../data/PMBB_Exome/audbase_feb252021/) | 7 | Audiogram database snapshot Feb 2021. `RGC21_45k_aud_1.csv.gz` is the matched audiogram-PMBB-ID linkage. Includes an R Markdown phenotyping notebook. |
| [`with_email/`](../data/PMBB_Exome/with_email/) | 4 | Subset of ZNF175 carriers with email contact (for re-contact / chart-review follow-up). |
| [`with_email_living/`](../data/PMBB_Exome/with_email_living/) | 12 | Same but restricted to living individuals. Also contains ESRRB cases/controls. Has a shell script `get_ids.sh.gz`. |
| [`regeneron/`](../data/PMBB_Exome/regeneron/) | 25 | Regeneron meta-analysis prep (chr20/chr21 SNP files, log/linear regressions). |
| [`manuscript/`](../data/PMBB_Exome/manuscript/) | 1 | `known_hl_genes.txt.gz` — final gene list used in the published Fig 2. |
| [`mtg_aug8_2021/`](../data/PMBB_Exome/mtg_aug8_2021/) | 1 | One-off meeting prep file. |
| [`ppv/`](../data/PMBB_Exome/ppv/) | 1 | `hl_values.txt.gz` — positive predictive value calc. |

### Scripts catalog ([`PMBB_Exome/scripts/`](../data/PMBB_Exome/scripts/))

81 scripts, all gzipped. Decompress on read or `zcat`. Grouped by function:

**Annotation & SNP selection**
- `only_HL_genes.py.gz` — filter annotation file down to gene list
- `annot_IDs_vs_pVCF.py.gz` — reconcile SNP IDs between annotation and pVCF
- `only_func_cats_to_include.py.gz`, `only_func_cats_to_include_stoploss_newannot.py.gz`, `only_func_cats_to_include_UKBB.py.gz` — filter to pLoF + missense REVEL>0.6
- `missense_or_pLOF.py.gz`, `format_hl_genes_stable.py.gz` — variant-class subsets
- `make_region_file.py.gz` — produce chr:start-end region file for biobin
- `format_meta.py.gz`, `format_results_to_meta.py.gz` — prepare files for cross-cohort meta-analysis

**Phenotype & case/control**
- `case_control.py.gz`, `case_control_allowOnlyPhecode.py.gz`, `case_control_allowOnlyPhecode_rmAudNA.py.gz`, `case_control_degreeHL1.py.gz` — different case/control schemes
- `HL_case_aud_and_phecode.py.gz` — the audiogram+phecode hybrid scheme used in the paper
- `add_case_ctrl.py.gz`, `add_carrier_degHL.py.gz`, `add_degHL_to_pheno.py.gz` — splice case/control or degree-HL into existing tables
- `keep_HL_cases_IBD.py.gz` — the "keep HL cases even if related" trick (README line 127)
- `hearing_aid_pheno.py.gz` — hearing-aid usage as alternative phenotype
- `get_ctrls.py.gz`, `get_ctrls_4max.py.gz` — control matching (4 controls per case)
- `add_PMBB_to_AudBase.py.gz`, `add_PMBB_to_AudBase_matchDOB.py.gz` — join PMBB IDs to Brant's audiogram database
- `format_UKBB_phenos.py.gz` — UKBB phenotype prep

**Covariates**
- `make_covs.py.gz` — base covariate file
- `add_anc.py.gz` — genetic ancestry
- `add_PCs.py.gz` — principal components (20 PCs)
- `add_MRN.py.gz` — medical record numbers
- `add_pheno_covs_to_biobin.py.gz`, `add_pheno_covs_to_biobin_DegHL1.py.gz`, `add_pheno_covs_to_biobin_UKBB.py.gz` — pheno+covs into biobin format

**Burden & per-gene counts**
- `burden_bin.py.gz` — burden binning
- `sum_of_lof_variants_in_HL_genes.py.gz`, `sum_of_lof_variants_in_HL_genes_DFNX.py.gz`, `sum_of_lof_variants_in_HL_genes_excludeList.py.gz`, `sum_of_lof_variants_in_HL_genes_degHL1.py.gz`, `sum_of_lof_variants_in_HL_genes_UKBB.py.gz`, `sum_of_lof_variants_in_HL_genes_DFNX_UKBB.py.gz` — aggregate burden per individual
- `counts_per_gene.py.gz`, `counts_per_gene_cases.py.gz`, `counts_per_gene_cases_degHL1.py.gz`, `counts_per_gene_cases_degHL1_Ncarriers.py.gz`, `counts_per_gene_cases_rand.py.gz`, `counts_per_gene_cases_UKBB.py.gz`, `counts_per_gene_rand.py.gz` — per-gene carrier counts (cases / controls / randomized)
- `n_case_controls_biobin.py.gz` — N per group, biobin format
- `homozygous_DFNB.py.gz`, `hmz_DFNA.py.gz` — homozygous-genotype counts per inheritance pattern
- `two_by_two.py.gz` — 2×2 Fisher's tables
- `cat_biobin_chr_t_addChr.py.gz` — concatenate per-chr biobin outputs, transpose, add chromosome

**Gene-specific (ZNF175, TCOF1, ClinVar)**
- `ZNF175_carrier.py.gz`, `ZNF175_carrier_send.py.gz` — identify ZNF175 variant carriers from VCF
- `TCOF1_annots_carriers.py.gz` — TCOF1 carrier annotation
- `ClinVar_carrier.py.gz`, `ClinVar_carrier_addGene_HLstatus.py.gz` — ClinVar P/LP carrier identification
- `ClinVar_and_nonClinVar_counts.py.gz`, `ClinVar_and_nonClinVar_counts_UKBB.py.gz` — burden split by ClinVar status

**Regression (all R)**
- `run_regression.R.gz`, `run_regression_anc.R.gz` — base burden regression (logistic by default)
- `run_logistic_regression.R.gz`, `run_logistic_regression_UKBB.R.gz`, `run_logistic_regression_rand.R.gz` — explicit logistic
- `run_regression_degHL1_binary.R.gz`, `run_regression_degHL1_linear.R.gz` — degree-of-HL phenotype (binary vs continuous)
- `run_regression_rand.R.gz` — randomized phenotype null
- `run_poisson.R.gz`, `run_poisson_rand.R.gz`, `run_negbinom.R.gz`, `run_negbinom_rand.R.gz`, `run_quasipoisson.R.gz`, `run_quasipoisson_rand.R.gz` — count-model alternatives (allele counts as dep var)
- `meta.R.gz`, `parse_results_to_meta.R.gz`, `parse_results_to_meta_degHL.R.gz` — meta-analysis (PMBB × UKBB × Regeneron)

**Output formatting & viz prep**
- `format_regression_results.py.gz`, `format_regression_results_geneNames.py.gz`, `format_regression_results_geneNames_STable.py.gz`, `format_regression_results_geneNames_STable_degHL1.py.gz`, `format_regression_results_geneNames_logistic.py.gz` — clean regression output → paper-table form
- `manhattan_format.py.gz` — Manhattan-plot input prep

## `PMBB_Imputed/` — GWAS & PRS workspace

385 GB (mostly `genotypes/` per-chr imputed plink files — 16k+ files). Top-level files cluster around:

**Phenotypes & covariates** (same individuals as `PMBB_Exome/`, but on imputed GWAS data)
- `deghl.txt.gz`, `deghl_noheader.txt.gz`, `deghl_rand.txt.gz` — degree-of-HL phenotype
- `snhl.txt.gz`, `snhl_plus1.txt.gz`, `snhl_plus1_noheader.txt.gz`, `snhl_plus1_rand.txt.gz` — sensorineural HL (binary)
- `covs.txt.gz`, `covs_discrete.txt.gz`, `covs_quant.txt.gz`, `covs_noheader.txt.gz` — covariates
- `indvs.txt.gz`, `EUR_keep.txt.gz`, `has_aud.txt.gz` — sample lists
- `IID_Anc.txt.gz` — ancestry per individual

**GWAS summary statistics & PRS inputs**
- `GWAS_PRSice.txt.gz` (257 MB), `GWAS_hm3_hg19_formatted.txt.gz`, `GWAS_hm3_hg19_formatted_PRScs.txt.gz`, `GWAS_hm3_hg38_formatted.txt.gz` — formatted summary stats for PRSice / PRS-CS
- `hapmap3.txt.gz` (38 MB), `hm3_allchr_hg19.bim.gz` — HapMap3 LD reference
- `2247_1.v1.1.fastGWA.gz` (234 MB) — fastGWA output (UKBB or similar — needs check)
- `media-1 (2).txt.gz` — published supplement (looks like a paper supp-table download)
- `n_carriers_nVariants_knownHL_withDegHL1.txt.gz` — joint table of carrier + degree-HL

**Subdirectories**
- `genotypes/` (huge, per-chr imputed bed/bim/fam files)
- `paper_replication/` — chr-by-chr extract of just the paper's GWAS SNPs (for replication and reanalysis under dominant/recessive models — `merged_paperSNPs_{dominant,recessive}.assoc.linear.gz` exist)
- `prscs/`, `prsice/` — PRS-CS and PRSice outputs
- `gcta/` (16 KB index) — GCTA heritability runs
- `h2/` — heritability outputs (paper reports h² = 4.53%)
- `freqs/` — allele frequencies
- `bootstrap/` — bootstrap resampling
- `modeling/`, `results/` — downstream models and outputs
- `ldblk_ukbb_eur/` — LD blocks (UKBB EUR, for PRS-CS)
- `PLPPR5/` — deep-dive on the single GWAS hit (chr1:99058420:C:T)
- `backup/`
- `scripts/` — see below

### Scripts catalog ([`PMBB_Imputed/scripts/`](../data/PMBB_Imputed/scripts/))

16 scripts:
- `paper_replication.py.gz`, `paper_SNPs_extract.py.gz` — extract paper SNPs from current pVCF
- `regeneron_replication.py.gz` — Regeneron-cohort replication
- `make_pheno.py.gz` — phenotype file generation
- `format_ss_PRScs.py.gz`, `format_ss_freq.py.gz`, `format_ss_freq_alleles.py.gz` — format summary stats for PRS-CS
- `format_p_freq_for_qq.py.gz` — QQ-plot input
- `liftover_bim_hm3.py.gz` — hg19↔hg38 liftover on HapMap3 bim
- `impR2_gt.30.py.gz` — imputation R² > 0.3 filter
- `sum_scores.py.gz` — sum PRS scores per individual
- `bootstrap.R.gz` — bootstrap resampling
- `make_figs.R.gz`, `make_figs_binary.R.gz`, `make_figs_binary_rand.R.gz`, `make_figs_linear.R.gz` — final figure generation (R)

## `DFNA/` — Shadi's separate paper

This is a **sibling project** (Shadi's DFNA analysis), not ours. Kept here because the HL gene list was merged with Shadi's. Likely not needed for our pipeline beyond the gene list. 188 files, 910 MB.

Structure: top-level per-category VCFs (`cases_category{1,2,3}*`, `controls_category{1,2,3}*`), `archive/`, `revisions/` (paper revision response), `scripts/` (5 scripts: `count_carriers_per_gene.py.gz`, `count_perc_with_variants.py.gz`, `get_SNPs_to_extract.py.gz`, `only_HL_genes.py.gz`, `regions_per_gene.py.gz`).

## Findings & flags

### ⚠ "8 signal-driving cases" — discrepancy to confirm

Kickoff meeting and analysis plan refer to "8 signal-driving cases" / "8 cases driving the burden signal" for ZNF175. The runbook (now decompressed at [`analysis/daniel/runbook_hui2023.txt`](../analysis/daniel/runbook_hui2023.txt)) references **multiple** counts:

| Line | Definition | Carriers | Phecode-cases |
|---|---|---:|---:|
| 193-196 | All individuals at chr19:51587727 | 8 | 1 (PMBB2106731298975) |
| 197-200 | All individuals at chr19:51581437 | 90 | 2 (PMBB7501686571326, PMBB1245988577461) |
| 348-353 | Joe's named-individual email | 6 | 1 |
| (later) | ZNF175 pLoF aggregate from add-back stoploss + multiallelic run | varies | varies |

**Question for Daniel/Doug/Molly:** what does "8 signal-driving cases" mean? Most likely interpretation: the 8 carriers at chr19:51587727 (regardless of phecode status). But Joe's email also flags 6 individuals — partial overlap unclear. Needs resolution before Priority-1 deep-dive starts.

See [`pipeline_walkthrough.md`](pipeline_walkthrough.md) Phase 8 + 11 for the detailed cross-reference.

### ⚠ PMBB version mismatch

All Daniel's files are built against `PMBB-Release-2020-2.0` (v2). Our [`pmbb_v3/`](../data/pmbb_v3/) symlinks point at `PMBB-Release-2024-3.0` (v3). The replication will need to either:
- regenerate genotype/annotation tables on v3 (correct way — sample composition has changed)
- OR re-use Daniel's v2 cohort and call this a "documentation pass, not a re-derivation"

Plan currently implies the former. Need to confirm with Molly.

### Audiogram linkage

`audbase_feb252021/` contains `RGC21_45k_aud_1.csv.gz` (the Brant + RGC matched audiogram file with PMBB IDs, from Feb 2021). Kickoff notes 4K individuals with both audiograms and exomes — this file is the foundation for that subset. **Action:** confirm with ENT team whether a refreshed audiogram-PMBB ID linkage exists for v3, or whether the Feb 2021 snapshot is still authoritative.

### Where the original pipeline starts

If reading top-down: `data/PMBB_Exome/README.gz` line 1 → `data/PMBB_Exome/scripts/only_HL_genes.py.gz` → annotation table → plink extraction → biobin → regression. The README is essentially a runbook. **First pass at a pipeline walkthrough should follow the README line by line**, decompressing each referenced script as we go.

## Next steps

1. ✅ Decompress `data/PMBB_Exome/README.gz` → [`analysis/daniel/runbook_hui2023.txt`](../analysis/daniel/runbook_hui2023.txt)
2. ✅ Decompress 102 scripts → [`analysis/daniel/scripts/`](../analysis/daniel/scripts/)
3. ✅ Write [`pipeline_walkthrough.md`](pipeline_walkthrough.md) phase map covering all 19 phases of the runbook
4. ☐ Resolve the "8 cases" ambiguity (see Findings)
5. ☐ Pilot re-run of Phase 4 (first biobin run) on existing intermediates to confirm we can reproduce the published top hit (ESRRB)
6. ☐ Start v3 porting in [`analysis/01_phase1_exploration_pmbb_v3/`](../analysis/01_phase1_exploration_pmbb_v3/) — beginning with Phase 1 (gene-list curation)
