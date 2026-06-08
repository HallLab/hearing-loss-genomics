---
marp: true
title: ZNF175 → tinnitus across PMBB releases
paginate: true
---

# ZNF175 → Tinnitus: reproducing Park (2021) and tracking the signal from 11K to 44K

**Andre Rico** · Hall Lab · Chapter 2

*Why was the ZNF175–tinnitus burden signal strong in Park's 11K exomes but "lost" in the ~44K replication — artifact, or biology?*

---

## 1. PMBB data across releases

| Release | WES (exome) | Imputed (array) | Year | Used by |
|---|---|---|---|---|
| **Freeze One** | **11,451** | — | ~2019–20 (pre-release) | **Park (discovery)** |
| Release 2020 **v1.0** | — *(none)* | 21,263 | 2020 | — |
| Release 2020 **v2.0** | **~44,086** | yes | 2020 | **Hui** |
| Release 2024 **v3.0** | ~57,171 | 56,358 | 2024 | — |
| Release 2026 **v4.0** | ~70,925 | 70,493 | 2026 | — |

- **EHR phenotypes** (ICD-9/10, labs, vitals, meds) are linked for the cohort in every release.
- ⚠️ Key point: **"Park's v1" is the *Freeze One* exome freeze, NOT the public "Release v1.0"** (which is imputed-array only, no WES). WES first appears publicly at **v2.0**.

---

## 2. Park et al. 2021 (*Nature Medicine*) — the origin

- **Genome-first ExWAS**: collapse rare **pLOF** variants per gene → test against ~1,000 EHR phecodes, in **~11K PMBB exomes** (Freeze One).
- 97 gene–phenotype hits at *p* < 10⁻⁶ → **26 robust** (replicated across BioMe / DiscovEHR / UKB).
- ★ **ZNF175 → tinnitus**: ***p* = 3.24 × 10⁻¹⁰**, replicated in 3 biobanks.
  - Mouse ortholog *Zfp719* knockout is **profoundly deaf**; expressed in inner/outer **hair cells**.
- **This is the published origin of the project's ZNF175 priority** (it is *not* in Hui 2023).

---

## 3. Hui et al. 2023 (*PLOS Genetics*)

- **Hearing-loss-focused** gene-burden + GWAS in **~40K PMBB exomes (Release v2.0)**.
- Burden = **pLOF + missense REVEL > 0.6**; phenotype = **audiogram-hybrid hearing loss** (not just phecode).
- Findings: novel HL candidate genes; one common-variant GWAS hit (**PLPPR5**).
- ⚠️ **ZNF175 is *not* in this paper** — and the ZNF175–tinnitus signal had **decayed** by this larger freeze.

---

## 4. Park vs Hui — what differs

| | **Park 2021** | **Hui 2023** |
|---|---|---|
| Cohort | Freeze One, **~11K** WES | Release v2.0, **~44K** WES |
| Scope | **exome-wide** (~1,518 genes) | **HL-focused** (~173 known HL genes) |
| Phenotype | **tinnitus** (phecode) + others | **hearing loss** (audiogram-hybrid) |
| Variants | **pLOF** (discovery) | **pLOF + missense REVEL>0.6** |
| ZNF175 | ★ **lead hit** (tinnitus) | **absent** |

→ Different cohorts, scope, phenotype **and** variant set — so it was unclear whether the ZNF175 "loss" was **biology** or **methodology**.

---

## 5. Our analysis — Park's pipeline on BOTH cohorts

**Hold the pipeline fixed, vary only the cohort.** Apply Park's design identically to v1 (11K) and v2 (44K):
pLOF + MAF ≤ 0.1% → **carrier status**; **tinnitus** (ICD 388.3x/H93.1x, rule-of-2); logistic `~ carrier + age + age² + sex + PC1–10`, EUR/AFR-stratified + meta.

**Don't conflate three quantities:**

| | N | Tinnitus (total) | ZNF175 carriers | **Carrier AND tinnitus** |
|---|---|---|---|---|
| **v1 (~11K)** | 9,161 | 131 (1.4%) | 26 (0.28%) | **4** |
| **v2 (~44K)** | 42,614 | 841 (2.0%) | 69 (0.16%) | **4** |

Tinnitus is **common**; being a **carrier is rare**; the burden test asks if tinnitus is **enriched among carriers**.

---

## 6. Results — the signal reproduces at 11K and decays at 44K

![w:1000](results/07/fig_znf175_v1_v2.png)

| | Fisher 2×2 | Adjusted logistic |
|---|---|---|
| **v1 (~11K)** | OR **12.9**, *p* = 4.7×10⁻⁴ | OR **14.6**, ***p* = 4.8×10⁻⁶** |
| **v2 (~44K)** | OR 3.07, *p* = 0.048 | OR 3.5, *p* = 0.015 |

We **reproduced Park's discovery** at 11K (vs his published 3.24×10⁻¹⁰).

---

## 7. Conclusion

- With an **identical pipeline**, the signal is **strong at 11K and decays at 44K** → **not a pipeline artifact**.
- The **~4 carrier-cases anchoring the signal did not grow** with the cohort (carrier pool 26 → 69); new carriers are ultra-rare and phenotype-negative → **enrichment dilutes** → OR drops ~14 → ~3.5.
- This is the classic **winner's-curse** signature of a small-case-driven signal.
- Hui's **"loss"** = this effect-size regression **+** the **exome-wide bar** (*p* < 10⁻⁶ across ~1,500 genes), which the gene-focused v2 signal doesn't clear.
- **Bottom line:** the signal is **attenuated, not a clean null** — consistent with a real but modest, small-case-anchored effect.

---

## 8. Caveats & next steps

- **Small carrier-case counts (n = 4)** → estimates noisy; qualitative pattern (strong v1, weak v2) is the robust takeaway.
- Faithful to **Park's method**, not bit-identical to his exact code (VEP vs ANNOVAR; ICD vs phecode; standard vs Firth; GRCh38 vs GRCh37). A bit-for-bit replication is feasible (tools on LPC) — documented next step.
- Cross-cohort **per-individual matching not possible** (different ID schemes) — "4 in both, likely overlapping" is inferred.
- **Next:** (1) profile the 4 carrier-cases (variants, age, ancestry); (2) phecode via R `PheWAS` + Firth; (3) extend the burden to **v3/v4** (≥57K/70K) to see where the effect settles.

---

## Appendix — provenance

- Cohorts: Freeze One `/static/PMBB/PMBB_Freeze17/` (GRCh38, 11,451); Release 2020 v2.0 `data/pmbb_v2/` (GRCh38, ~44K).
- Pipeline: `analysis/chapter_2/scripts/0X/` (NB 00–07); outputs `analysis/chapter_2/results/0X/`.
- Reference papers: Park et al. 2021 *Nat Med* 27(1):66–72; Hui et al. 2023 *PLOS Genet*.
- Full write-up: `analysis/chapter_2/findings_znf175_11k_vs_44k.md`.
