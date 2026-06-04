# Reference Paper Summary — Park et al. 2021 (Nature Medicine)

**Citation:** Park J, Lucas AM, Zhang X, Chaudhary K, Cho JH, Nadkarni G, Dobbyn A, Chittoor G, Josyula NS, Katz N, Breeyear JH, Ahmadmehrabi S, Drivas TG, Chavali VRM, Fasolino M, Sawada H, Daugherty A, Li Y, Zhang C, Bradford Y, Weaver J, Verma A, Judy RL, Kember RL, Overton JD, Reid JG, Ferreira MAR, Li AH, Baras A, Regeneron Genetics Center, LeMaire SA, Shen YH, Naji A, Kaestner KH, Vahedi G, Edwards TL, Chen J, Damrauer SM, Justice AE, Do R, Ritchie MD, Rader DJ. (2021). *Exome-wide evaluation of rare coding variants using electronic health records identifies new gene-phenotype associations.* **Nature Medicine 27(1):66–72.** doi: 10.1038/s41591-020-1133-8

**Local copy:** `docs/papers/nihms-1768013.pdf` (PMC author manuscript, PMC 2022 Jan 20)

---

## ⭐ Why this paper matters most to OUR project

**This is the published source of the project's ZNF175 priority.** ZNF175 → tinnitus is a robust, replicated finding in Table 2 of this paper (NOT in Hui et al. 2023 PLOS Genetics — that paper does not mention ZNF175). The Hall Lab / Epstein Lab ZNF175 follow-up *extends* this published Park 2021 finding; it does not start from scratch.

- **ZNF175 → Tinnitus:** discovery **p = 3.24 × 10⁻¹⁰**, **3 replications**, ✓ clinical/experimental evidence (Table 2, the strongest discovery p-value among the novel genes).
- **Hearing loss "barely missed" the significance threshold** (`p < E-06`) for ZNF175 — i.e. the binary-phecode HL signal sits just under the line while tinnitus clears it.
- Replicated in **BioMe, DiscovEHR, and UKB**.
- **Mouse ortholog *Zfp719*:** LOF mice are **profoundly deaf**, have abnormal Preyer reflex (auditory startle) and raised ABR thresholds (refs 20–21: Bowl 2017 Nat Commun; Ingham 2019 PLoS Biol).
- *Zfp719* is **expressed in inner and outer hair cells** of the mouse ear (ref 22: Liu 2014 J Neurosci).
- Human **ZNF175** has a suggested role in **neurotrophin production and neuronal survival** (ref 23).

This is the published anchor behind Doug Epstein's mouse biology motivation and Daniel Hui's gene-specific PMBB follow-up. See [[feedback-znf175-not-in-paper]] for the distinction vs. Hui 2023.

---

## TL;DR

A "genome-first" **exome-by-phenome-wide association study (ExoPheWAS)** in the Penn Medicine BioBank (PMBB). Collapses rare predicted loss-of-function (pLOF) variants into per-gene "gene burdens" and tests each against ~1,000 EHR-derived phecodes. Discovers 97 gene-phenotype associations, validates 26 as robust via the DiCE framework + multi-biobank replication. Five are positive controls; 21 are novel — spanning glaucoma, aortic ectasia, type 1 diabetes, muscular dystrophy, and **hearing/tinnitus (ZNF175)**.

This is the methodological parent of the rare-variant burden approach used throughout the project. (Joe Park, first author, is noted in project records as having recommended the burden methodology.)

---

## Cohort (discovery)

- **PMBB discovery:** 10,900 individuals with WES linked to EHR (40.7% female; median age 67)
  - Ancestry: EUR 8,198 (75.2%), AFR 2,172 (19.9%), AMR 304 (2.8%), EAS 79 (0.7%), SAS 114 (1.0%)
  - **Hearing loss prevalence: 579 (5.3%); Tinnitus is a separate phecode.**
- Note: much smaller and earlier than the ~40k cohort in Hui 2023, and earlier than the PMBB v3 the project is extending to.

## Replication cohorts
- **PMBB2:** 6,432 African American individuals (GRCh38)
- **BioMe (Mount Sinai):** 23,989 (6,470 AFR / 8,735 EUR / 8,784 Hispanic)
- **DiscovEHR (Geisinger):** 85,450 EUR (IDT + VCRome platforms, meta-analyzed)
- **UK Biobank:** 32,268 EUR (population-based — noted "healthy volunteer bias", lower disease prevalence)
- **BioVU (Vanderbilt):** 66,400 (genotype only, targeted single-variant replication)

---

## Methods (high-level)

### Discovery burden
- **pLOFs:** frameshift indels, stop gain/loss, canonical splice-site dinucleotide disruption (ANNOVAR, RefSeq)
- **MAF ≤ 0.1% in gnomAD v2**, ancestry-specific thresholds (NFE / AFR / Latino)
- Genes with **≥25 heterozygous pLOF carriers** → 1,518 genes tested
- Phecodes with **≥20 cases** → 1,000 phecodes (Phecode Map 1.2, R `PheWAS`; case = ≥2 ICD dates, control = never)
- Logistic regression adjusted for **age, age², sex, 10 PCs**; EUR and AFR run separately then **inverse-variance-weighted meta-analysis**
- Additive fixed-threshold collapsing; Firth penalized likelihood cross-check vs. exact logistic regression
- **Significance threshold: p < 10⁻⁶** (where observed QQ deviates from expected; λ∆95 = 1.558 overall)

### Robustness / replication
- Same-cohort: REVEL ≥ 0.5 missense burden + single variants (MAF > 0.1%), checked for non-overlapping carriers (no mutual-carrier driving)
- Cross-cohort: pLOF burden + missense burden + single variants in PMBB2, BioMe, DiscovEHR, UKB; targeted single variants in BioVU
- **DiCE (Diverse Convergent Evidence)** ranking: combines number of replications + clinical/experimental evidence; ≥2 total check marks → "robust"
- Chart review of carriers to adjudicate diagnoses; removed associations that lost significance or lacked a common etiology

### Note: UKB undercalling
- 3 genes (CES5A, CYP2D6, ZC3H3) overlap undercalled FE-pipeline regions; for these, variants with ≥65% call rate were retained (vs. <5% missingness elsewhere).

---

## Key results

### 26 robust genes (Table 2)

**5 positive controls:** CFTR (bronchiectasis / cystic fibrosis / pseudomonal pneumonia), TTN (cardiomyopathy / conduction / dysrhythmias), MYBPC3 (hypertrophic cardiomyopathy), BRCA2 (breast cancer), CYP2D6 (opiate adverse effects).

**21 novel** (ranked alphabetically below the line) — selected highlights:
| Gene | Phecode | Discovery P | Repl. (N) | Evid. |
|---|---|---|---|---|
| **ZNF175** | **Tinnitus** | **3.24E-10** | **3** | **✓** |
| DNAH6 | Lack of coordination | 7.93E-10 | 2 | |
| RGS12 | Type 1 diabetes | 6.48E-08 | 5 | ✓ |
| CILP | Aortic ectasia | 4.29E-08 | 3 | ✓ |
| BBS10 | Hypertrophic cardiomyopathy | 2.89E-08 | 1 | ✓ |
| EFCAB5 | Prolapse of vaginal walls | 3.19E-08 | 3 | |
| PPP1R13L | Primary open angle glaucoma | 7.29E-07 | 2 | ✓ |
| FER1L6 | Muscular wasting/disuse atrophy | 7.18E-07 | 3 | ✓ |
| MYCBP2 | Spasm of muscle | 2.08E-07 | 2 | ✓ |
| SCNN1D | Cardiac conduction disorders | 4.52E-07 | 5 | |
| CES5A | Abnormal coagulation profile | 8.10E-08 | 5 | |
| (also: ABCA10, CTC1, DNHD1, EPPK1, FLG2, RTKN2, TGM6, TRDN, WDR87, ZNF334) | | | | |

Functional validation generated for several novel genes (scRNA-seq, immunolocalization, GEO meta-analysis): PPP1R13L/glaucoma (iPSC-RGC oxidative stress), RGS12/T1D (islet macrophages), CILP/aorta (adventitial fibroblasts + TGF-β), MYCBP2/muscle (TMD expression).

### African-ancestry-specific signal
- 16 rare deleterious single variants, African-ancestry-specific, replicated their discovery gene-disease associations; none previously in the GWAS catalog or literature.

---

## Takeaways for project framing

1. **Medical biobanks (PMBB) > population biobanks (UKB)** for rare-variant disease discovery — sicker population, higher case prevalence (relevant to why ZNF175 HL signal behaves differently across cohorts).
2. The **burden methodology** here (pLOF + REVEL≥0.5 missense, MAF cutoffs, meta-analysis design) is the template the project's replication tracks follow.
3. **Don't expect a uniform QQ fit** for rare-pLOF burden p-values; validity comes as much from cross-cohort replication as from any single threshold.
4. ZNF175's tinnitus signal is strong while its HL signal sits just under threshold — consistent with the project's focus on quantitative-trait / audiogram analysis and the "modifier + second-hit" model rather than binary-phecode HL alone. See [[project-znf175-priority]].

---

## Relationship to the other reference paper

- **Park 2021 (this paper, Nat Med):** exome-wide ExoPheWAS, ~10.9k PMBB, source of **ZNF175 → tinnitus**. Methodological parent.
- **Hui 2023 (PLOS Genetics, `pgen.1010584.pdf`):** HL-focused gene-burden + GWAS in ~40k PMBB; novel HL candidates COL5A2/HMMR/NNT/RAPGEF3; **does NOT mention ZNF175**. See `paper_summary_hui2023.md`.

The project replicates Hui 2023's HL pipeline while the ZNF175 deep-dive traces its published origin to Park 2021.
