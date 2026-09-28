# Q1 — preliminary result (PMBB v4, phecode 389 "hearing impairment")

**Date:** 2026-08-26 · **Author:** Andre Rico · **Status:** PRELIMINARY — charter not yet ratified
**Script:** [`../../analysis/q1_primary/run_q1_preliminary.py`](../../analysis/q1_primary/run_q1_preliminary.py)
**Compute:** none. The genes were already tested exome-wide by Elena; this applies the pre-specified
set and the multiple-testing discipline that were missing.

---

## 1. What was pre-declared (before looking)

| Choice | Value |
|---|---|
| Primary set | ClinGen HL GCEP **Definitive+Strong** ∩ phecode 389 → **61 genes** |
| Primary p-value | SAIGE-GENE+ **Cauchy omnibus** per (gene, mask, MAF, ancestry) |
| Multiple testing | Benjamini-Hochberg **within each (ancestry, mask, MAF) stratum** |
| MOI split | AD vs AR declared in advance; XL reported separately |
| Calibration | λ_GC computed exome-wide per stratum |
| Phenotype | `hearing impairment` (n=57,632; 6,712 cases) — **HL-only, not the combined HL/tinnitus** |

Primary-set composition: **AR 40 · AD 19 · XL 2**.

## 2. Calibration — clean

λ_GC across all **36 strata**: **0.85 – 1.06**. No inflation anywhere. Mild deflation in the
`pLOF_pDM` strata (0.84–0.92) is the usual behaviour of rare-variant tests with many low-MAC genes.
Full table: [`calibration_lambda_gc.csv`](calibration_lambda_gc.csv).

> Note: the technical controls declared in the charter (`BRCA1`/breast cancer, `TTN`/cardiomyopathy,
> `CFTR`/CF) could **not** be run — they require SAIGE on a different phenotype, which is compute we
> have not spent. λ_GC is the calibration evidence available today.

## 3. Result — 4 genes, all risk-increasing

22 gene×stratum rows pass FDR<0.05, but the strata are nested (MAF 0.001 ⊂ 0.01; pDM ⊂ pLOF_pDM;
EUR ⊂ combined). **The honest count is 4 distinct genes.**

| Gene | MOI | Best p | FDR | OR | MAC case / control | Ancestry | Mask |
|---|---|---|---|---|---|---|---|
| **COCH** | AD | 1.18×10⁻⁴ | 0.006 | **3.40** | 23 / 68 | EUR, combined | pDM, pLOF_pDM |
| **TMPRSS3** | AR | 2.68×10⁻⁴ | 0.015 | **4.73** | 14 / 34 | EUR, combined | pDM, pLOF_pDM |
| **GJB3** | AD | 2.17×10⁻⁴ | 0.012 | **3.14** | 16 / 55 | AFR | pDM, pLOF_pDM |
| **TPRN** | AR | 6.74×10⁻⁴ | 0.014 | **2.86** | 24 / 92 | EUR | pDM, pLOF_pDM |

All four have **positive beta** — increased risk. Full table: [`q1_primary_full.csv`](q1_primary_full.csv).

## 4. Reading it

**COCH is the headline, and it is the most a-priori plausible gene in the whole set.** COCH causes
DFNA9 — autosomal dominant, **late-onset progressive** sensorineural hearing loss with vestibular
involvement, typically presenting in the 3rd–5th decade. It is the one Mendelian HL gene whose natural
history actually matches "adult-onset hearing loss in a hospital biobank".

**The mask pattern is an internal consistency check that passes.** Every hit sits in `pDM`
(damaging missense) or `pLOF_pDM` — **none in pure `pLOF`**. For COCH this is the known mechanism:
DFNA9 is caused by *missense* variants (LCCL domain, dominant-negative / aggregation), not by
haploinsufficiency; COCH truncating carriers are generally unaffected. A noise hit would not be
expected to land on the biologically correct mask.

**The two AR genes are the more surprising and the weaker claim.** TMPRSS3 and TPRN are recessive;
heterozygous carriers showing OR 2.9–4.7 is exactly the "carrier state contributes to adult HL"
hypothesis. But it is equally consistent with a handful of true recessive (compound het) individuals
inside the case group. Distinguishing those requires genotype-level inspection, not the burden table.

**This partially contradicts the charter's stated expectation.** The charter says the expected outcome
is null. It is not null. That is a result worth stating plainly at review.

## 5. Caveats — read before quoting any of this

1. **Small counts.** MAC_case 14–29. Effect estimates at this scale are unstable and inflated at
   discovery — this is precisely the Cycle-1 winner's-curse lesson. The ORs will regress.
2. **No independent replication.** `combined` contains `EUR`, so those are not two observations.
   Only GJB3 (AFR) is ancestry-independent of the EUR signals.
3. **Wrong phenotype vs. the 2026-07-01 decision.** This is HL-only. The combined HL-and/or-tinnitus
   phenotype still does not exist.
4. **Masks are collapsed.** `pDM` = AlphaMissense OR REVEL, unioned. The meeting decided to keep them
   separate and to test REVEL 0.5 vs 0.6. Not done.
5. **3 primary genes never tested:** `POU3F4`, `SMPX` (X-linked — the pipeline runs chr1–22, exactly
   the silent drop flagged in advance) and `LRTOMT` (cause unknown; likely a symbol/alias mismatch —
   needs checking, it is on chr11 and should be present).
6. **Source-file defect.** `newmasks_combined_gene_results.csv` has 41 fields of which only the first
   16 are real, and one data row (gene `ABRA`) is glued into the header and therefore missing from
   every table built from this file. Recorded in [`manifest.json`](manifest.json). Elena should be
   told before anyone else reads that file.

## 6. Next

1. Tell Elena about the header defect (§5.6).
2. Chase `LRTOMT` — a missing primary gene is a pipeline question, not a biology one.
3. Genotype-level look at the TMPRSS3 / TPRN carriers: het carriers, or hidden compound hets?
4. Decide XL: declare the primary set autosomal (61→59) or extend the pipeline to chrX.
5. Build the combined HL/tinnitus phenotype and re-run — the one 2026-07-01 decision still unexecuted.
6. Once audiogram access lands: check whether COCH carriers show the expected progressive
   high-frequency pattern. That is the sharpest available test of whether this is real.
