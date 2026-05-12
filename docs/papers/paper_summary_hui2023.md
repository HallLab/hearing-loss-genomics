# Reference Paper Summary — Hui et al. 2023 (PLOS Genetics)

**Citation:** Hui D, Mehrabi S, Quimby AE, Chen T, Chen S, Park J, Li B, Ruckenstein MJ, Rader DJ, Ritchie MD, Brant JA, Epstein DJ, Mathieson I. (2023). *Gene burden analysis identifies genes associated with increased risk of four congenital and adult-onset hearing loss phenotypes.* PLOS Genetics. doi: 10.1371/journal.pgen.1010584

**Local copy:** `docs/papers/pgen.1010584.pdf`

---

## TL;DR

A gene-burden study in the Penn Medicine BioBank (PMBB) that asked whether rare deleterious variants in **known congenital hearing loss (HL) genes** also raise risk of **adult-onset HL**. Answer: yes, with measurable effect. Four novel candidate genes also emerged. Common-variant GWAS found one new locus (PLPPR5). Polygenic and burden scores have modest predictive power.

---

## Cohort

- **PMBB total:** 40,627 individuals (51% M / 49% F; median age 58; ~75% European, ~25% African American ancestry)
- **Phecode-based phenotyping (phecode 389):**
  - ≥2 instances → case (3,304 potential)
  - 0 instances → control (34,704 potential)
  - 1 instance → missing/NA
- **Audiogram crosscheck (n=1,917 with audiograms):**
  - 65% of phecode-defined cases met audiogram HL (PTA > 25 dB)
  - 27% of phecode-defined controls also met audiogram HL → control set is "leaky"
- **Hybrid case definition adopted:** degree of HL on 0–4 scale (PTA-bin based) when audiogram available; phecode otherwise
- **Final analysis cohort:** 36,507 individuals (1,110 cases / 35,397 controls), 11,977,860 common variants, 542,910 rare variants

## Comparison cohorts
- **UK Biobank:** 148,970 (39,272 cases / 109,698 controls) — self-reported "hearing difficulty/problems" + "hearing difficulty/problems with background noise"
- **Regeneron meta-analysis:** 5-cohort meta-analysis, 125,749 common-variant controls / 469,497 rare-variant cases analysis

---

## Methods (high-level)

### Rare variant burden
- pLoFs: frameshift, stopgain, splicing
- Plus missense with REVEL > 0.6
- gnomAD MAF < 0.001, cohort MAF < 0.01
- Tests: gene-wise and total gene-burden (summed across known HL gene set)
- HL gene set: **173 known HL genes** (DFNA = dominant, DFNB = recessive)

### Common variant GWAS
- Additive model, MAF > 1%
- Phenotype randomization control for inflation (λ = 1.01 final)

### Other
- Heritability estimated
- Polygenic risk score via PRS-CS using UKBB GWAS summary statistics
- Replication: against Praveen et al. and Ivarsdottir et al. published results

---

## Key results — numbers to replicate

| Metric | Value |
|---|---|
| Heritability of adult HL in PMBB | **h² = 4.53%** |
| Cases / controls (final) | 1,110 / 35,397 |
| Controls carrying ≥1 deleterious variant in known HL gene | 72.8% |
| Cases carrying ≥1 deleterious variant in known HL gene | 74.0% (p = 0.51 by Fisher's) |
| Controls with ≥1 pathogenic/likely-pathogenic ClinVar variant | 6.80% |
| Cases with same | 7.93% (p = 0.147) |

Note: at the carrier-vs-non-carrier level the cases-vs-controls difference is small. The signal is in the *aggregated burden* and in specific genes, not in any single variant.

### Genes implicated (rare variant)
- **Replicated known HL genes with increased burden:** TCOF1, ESRRB (associated with HL in PMBB and replicated)
- **Novel candidate genes (HL not previously implicated):** COL5A2, HMMR, NNT, RAPGEF3
  - Plausible mechanisms inferred for 3 of 4 (NNT — mitochondrial / insulin / HL; RAPGEF3 — noise-induced inner ear pathology, beta-cell dysfunction; ZNF175 — the priority gene for our follow-up)

### GWAS hit
- **PLPPR5** (chr1:99058420:C:T, MAF=0.012, p=8.27×10⁻⁹)
- Single GW-significant locus
- Mouse cochlear scRNA-seq: Plppr5 selectively expressed in inner hair cells
- Did NOT replicate in Praveen et al. (p=0.759) — caveat
- Conversely: 9 of 45 Praveen et al. lead SNPs replicated at p<0.05 in PMBB (binomial p=0.0003)

### PRS / burden score
- PRS-CS scores constructed from UKBB GWAS
- Both PRS and rare-variant burden score have **low individual predictive power** — they explain risk in aggregate but don't reliably stratify single individuals

---

## Author summary takeaways (for project framing)

1. Loss-of-function variants in known Mendelian HL genes raise risk of adult-onset HL — supporting overlap between congenital and adult forms.
2. Many "Mendelian" HL variants are **incompletely penetrant** and may act cumulatively.
3. Hospital-recruited biobanks (like PMBB) can power HL genetics studies even when HL isn't the primary reason for contact with the health system.

---

## What's not in the paper (and why our project exists)

- The ZNF175 hit is mentioned but its biological mechanism in humans is not resolved.
- The 8 cases driving the original burden signal are not deep-dived.
- No analysis of *second-hit* variants in HL genes among ZNF175 carriers.
- No use of the ~4,000 PMBB participants with both audiograms AND exomes for quantitative-trait analysis.
- No UKBB / All of Us replication of the ZNF175 specific signal.

These are the gaps Phase 1+ of our project addresses.
