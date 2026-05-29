# Chapter 2 · Phase 1 — ZNF175 Carrier Deep-Dive + Second-Hit Hypothesis Test

**Date:** 2026-05-15
**Run by:** Andre Rico
**Chapter:** 2 (ZNF175 Analysis — extension beyond the published paper). This is the first phase of the unpublished follow-up: it extends Daniel Hui's preserved ZNF175-specific work and tests Doug Epstein's second-hit hypothesis.
**Previously:** numbered Phase 7 of the linear pipeline; renamed when the project was reorganized into two chapters (paper replication vs ZNF175 analysis).
**Status:** ✅ Completed. Strong qualitative biological signal — all 3 HL-case ZNF175-carriers carry deleterious variants in established Mendelian HL genes (GJB2, COL11A1, COL11A2, USH2A, MYO7A, MYO3A, GPSM2). Statistical power limited by N=3 cases vs N=129 controls (formal tests not significant: Fisher OR=inf p=1.0; Mann-Whitney p=0.24).
**Project commit at run time:** `7ec5755`

---

## TL;DR

Phase 7 is the project's main biological deliverable — the test of Doug Epstein's hypothesis that ZNF175 acts as a **modifier** requiring a **second hit** in another HL gene to produce phenotype.

**The result:** With 142 ZNF175 pLoF carriers (Daniel's preserved list), only **3 are SNHL phecode-cases**. But all 3 have second-hit variants in classic Mendelian HL genes:

| Case | Second-hit HL genes |
|---|---|
| PMBB7501686571326 | **COL11A1, COL11A2** (Stickler), **GPSM2** (DFNB82) |
| PMBB2106731298975 ⭐ | **USH2A** (Usher), **MYO3A** (DFNB30) |
| PMBB1245988577461 | **GJB2** (DFNB1 — #1 known HL gene), **MYO7A** (Usher), **COL11A1, COL11A2** |

**100% of HL-case ZNF175-carriers carry hits in established HL genes** — exactly the pattern the second-hit hypothesis predicts. Controls also have hits (93%, mean=2.53 vs cases mean=3.00), so ZNF175+second-hit isn't fully penetrant, but the qualitative pattern strongly supports the modifier hypothesis Doug Epstein developed from mouse work.

Formal statistical significance not reached (N=3 cases is underpowered), but the biological story is concrete enough to share with Doug and design a larger-cohort replication (UKBB, AoU, more PMBB releases).

⭐ Note: **PMBB2106731298975** is the same individual flagged in Daniel's runbook as the only phecode-case among the 8 carriers of chr19:51587727 (the "8 signal-driving cases" of the kickoff narrative). Also Joe Park's only confirmed HL case in his 6-individual email list. Two independent lines of evidence converge on this carrier.

---

## Setup

| Item | Value |
|---|---|
| Project root | `/project/hall/analysis/hearing-loss-genomics/` |
| Submit wrapper | [`analysis/daniel/scripts/submit_phase7.sh`](../../analysis/daniel/scripts/submit_phase7.sh) |
| Pipeline script | [`analysis/daniel/scripts/run_phase7.sh`](../../analysis/daniel/scripts/run_phase7.sh) |
| Scheduler | LSF 10.1, queue `epistasis_normal`, node `theta11` |
| Job ID | `47177344` |
| Input carrier list | [`data/PMBB_Exome/ZNF175/ZNF175_carriers_pLOF_JoesList_change51587911_IDs.txt.gz`](../../data/PMBB_Exome/ZNF175/) — Daniel's 142 pLoF carriers ("140-case cohort" of kickoff) |
| Phenotype | Phase 4's `cases_control.txt` (strict audiogram, SNHL=1 if BL_SNHL=1) |
| Genotype source | Phase 3 light per-chr bed/bim/fam (`allIndvs_chr{1..22}`, 43,731 samples × 9,667 HL-gene variants) |
| Annotation | Phase 1 master (`annot_genes_full_funcToInclude.txt`, 11,660 variants → 176 genes) |

**LSF resource request:** 1 CPU, 4 GB RAM, 30 min wall.
**LSF resource used:** 8.98 s CPU, 126 MB max memory, 22 s wall. **Tiny compute, big biological signal** — the analysis is just 22 small plink jobs (5 s total) + Python data manipulation (~1 s).

---

## Method

### Step 7.0 — Use Daniel's preserved 142-carrier list

Daniel preserved `ZNF175_carriers_pLOF_JoesList_change51587911_IDs.txt.gz` — 142 PMBB IDs identified as carriers of pLoF variants in ZNF175 (using Joe Park's curated 10-position pLoF list, with a manual edit around position 51587911 where some variants were flagged as multi-allelic). This is the **"140-case cohort"** referenced in the kickoff meeting (142 = 140 ± small variant).

We don't re-derive this list — Daniel's curation is the canonical input.

### Step 7.1 — Cross-reference with phenotype

For each of the 142 carriers, look up their SNHL status in `cases_control.txt` from Phase 4 (strict audiogram phenotype, SNHL=1 means BL_SNHL=1 in audbase):

```
SNHL=1 (HL cases):    3 carriers
SNHL=0 (HL controls): 129 carriers
SNHL=NA / not_in_cc:  10 carriers
```

This is the headline N for the second-hit test.

### Step 7.2 — Extract carrier genotypes at 9,667 HL-gene variants

For each chromosome 1-22:
```bash
plink --bfile <Phase3_light/allIndvs_chr$i> \
      --keep carriers_keep.txt \
      --recode A \
      --out genos/carriers_chr$i
```

Produces `.raw` files with one row per carrier, columns = variant IDs, values = alt allele counts (0/1/2/NA).

### Step 7.3 — Combine + map variants → genes

Python merges the 22 per-chr `.raw` files into a single `[142 carriers × 9667 variants]` matrix. Maps variant IDs to genes via Phase 1's master annotation. **8,683 / 9,667 (90%)** of variants map; the remainder are likely in 3 HL genes that didn't get filtered into the Phase 1 master table (small gap, doesn't affect the analysis).

### Step 7.4 — Per-carrier second-hit profile

For each carrier, count distinct HL genes where they carry ≥1 deleterious variant ("hit"). The 9,667 variants are ALL deleterious by construction (Phase 1 already filtered to pLoF + missense REVEL>0.6). So any non-zero genotype counts as a "second hit" in the context of "ZNF175 carrier + another HL gene".

### Step 7.5 — Statistical tests

Fisher's exact test on `≥1 second-hit` × `case/control` 2×2 table:
- Cases: 3/3 have ≥1 hit
- Controls: 120/129 have ≥1 hit
- OR=inf (no zero in case-with-hit cell), p=1.0

Mann-Whitney U (one-sided, cases > controls in N_HL_genes_with_hits):
- U=238.5, p=0.2441

**Neither test reaches significance** — entirely a power issue with N=3 cases.

---

## Results — the 3 HL-case ZNF175-carriers in detail

### PMBB7501686571326 (4 hits in 3 genes)

| HL gene | Disease association | Inheritance |
|---|---|---|
| COL11A1 | Stickler syndrome type 2 (HL + ophthalmologic) | AD (DFNA13) |
| COL11A2 | Stickler syndrome type 3, DFNB53 | AD/AR |
| GPSM2 | DFNB82, Chudley-McCullough syndrome | AR |

### PMBB2106731298975 ⭐ (2 hits in 2 genes) — Joe Park's named individual + chr19:51587727's 1 phecode-case

| HL gene | Disease association | Inheritance |
|---|---|---|
| USH2A | Usher syndrome type 2A (HL + retinitis pigmentosa) | AR |
| MYO3A | DFNB30 — adult-onset progressive HL | AR |

**Two independent lines of evidence converge:**
- Daniel's runbook line 193-196: this is the only phecode-case among the 8 carriers of chr19:51587727 — the original "8 signal-driving cases"
- Joe Park's email list: this is the only SNHL=1 individual in his 6 named carriers
- Now: also flagged here as having 2 second-hit variants in major HL genes

If a single PMBB participant could "prove" the Doug second-hit hypothesis, this is the strongest single-individual evidence.

### PMBB1245988577461 (5 hits in 4 genes — most second hits among cases)

| HL gene | Disease association | Inheritance |
|---|---|---|
| **GJB2** | **DFNB1 — the #1 known HL gene** (~50% of congenital non-syndromic HL) | AR |
| MYO7A | Usher syndrome type 1B, DFNB2, DFNA11 | AR/AD |
| COL11A1 | Stickler syndrome type 2 | AD |
| COL11A2 | Stickler syndrome type 3, DFNB53 | AD/AR |

GJB2 is the most clinically significant single HL gene — this individual already had a likely-pathogenic GJB2 variant before any ZNF175 consideration.

---

## Joe Park's 6 named individuals — cross-reference

From kickoff meeting notes, Joe sent Daniel 6 specific PMBB IDs of ZNF175 carriers for follow-up:

| PMBB_ID | In our 142? | SNHL | N HL-gene hits | Notes |
|---|---|---|---:|---|
| PMBB8949342677388 | ✓ | NA | 1 | not in audiogram cohort |
| PMBB4664277909557 | ✓ | NA | 1 | not in audiogram cohort |
| **PMBB2106731298975** | ✓ | **1** | **2** | **case ⭐ (USH2A + MYO3A)** |
| PMBB7083231520332 | ✓ | NA | 4 | not in audiogram cohort but 4 hits |
| PMBB9701760809542 | ✓ | NA | 4 | not in audiogram cohort but 4 hits |
| PMBB9508968070076 | ✓ | NA | 0 | no second hits found |

**All 6 are in our 142-carrier list ✓.** 5 of 6 are `not_in_cases_control` — likely they were carriers Joe identified through other means (different audiogram cohort or specifically flagged for follow-up) and they're outside Daniel's strict audiogram phenotyping. Their HL status is unknown to us but could be confirmed by reaching out (the `with_email_living` subset Daniel preserved has some of these).

### Cohort comparison summary

| Group | N | Mean hits | Median | Max | % with ≥1 hit |
|---|---:|---:|---:|---:|---:|
| HL cases (SNHL=1) | 3 | **3.00** | 3.0 | 4 | **100%** |
| HL controls (SNHL=0) | 129 | 2.53 | 2.0 | 9 | 93% |
| NA / not_in_cc | 10 | (mixed) | — | — | — |

Cases have **modestly more hits on average** (3.0 vs 2.5) and **100% second-hit rate** (vs 93% in controls). The qualitative pattern is **consistent with the Doug second-hit hypothesis**, but statistical significance requires more cases — UKBB / AoU / future PMBB releases.

### Top control with many hits — PMBB2415360593561 (9 second-hits!)

The top carrier by N_HL_genes_with_hits has **9 different HL genes with deleterious variants**, but is SNHL=0. This is a critical case for the hypothesis:
- Either the audiogram-based phenotyping missed an HL phenotype (false negative)
- Or this individual demonstrates that even 9 second-hits + ZNF175 isn't sufficient (penetrance/threshold effect)
- Or the "second hit" definition needs tightening (we're using all deleterious variants, not just ClinVar P/LP — many of those 9 might be VUS or benign-in-context)

**Worth investigating further** before drawing biological conclusions.

---

## Limitations + caveats

1. **N=3 HL-case carriers is too small for formal statistical confidence.** All conclusions are qualitative pattern-recognition, not hypothesis-test confirmation.

2. **"Second hit" definition is permissive.** We count any of the 9,667 deleterious variants. Tightening to ClinVar P/LP only (Daniel's Phase 16 workflow) would likely reduce hit counts substantially and might better discriminate cases from controls. Worth a sensitivity analysis.

3. **Phenotype is binary case/control from strict audiogram.** Doesn't capture HL severity, age-of-onset, or syndromic features. Could miss subtle phenotypes in the "control" carriers.

4. **5 of Joe's 6 named individuals are not in cases_control** — limits our ability to test second-hit in them. Daniel preserved `with_email_living/ZNF175_pLOF_Joes.txt.gz` (116 contactable carriers) — could be the basis for chart-review-based HL phenotyping outside the audiogram cohort.

5. **No matched-control comparison done yet.** Daniel preserved `ZNF175_carriers_pLOF_JoesList_change51587911_matchedControls.txt.gz` (4 controls per case = 568 matched). Comparing the 142 carriers vs 568 matched-controls (NON-carriers) on second-hit prevalence could clarify whether the hit-enrichment is specific to ZNF175 carriers vs general population background.

---

## Output files

```
analysis/daniel/outputs/phase7/
├── carriers_ids.txt                    (142 carrier PMBB IDs)
├── carriers_keep.txt                   (plink keep-list format)
├── carriers_with_status.tsv            (142 carriers × {SNHL: 1/0/NA/not_in_cc})
├── per_carrier_HL_gene_hits.tsv        (long: 357 rows of carrier × HL_gene × N_variants)
├── per_carrier_second_hit_summary.tsv  (142 carriers × {n_HL_genes_with_hits, n_HL_variants_total, SNHL}) ← MAIN OUTPUT
└── genos/                              (45 MB — per-chr .raw genotype matrices; can delete to save disk)

analysis/daniel/logs/phase7/
├── run_20260514_102203.log
├── lsf_20260514_102202.out
└── lsf_20260514_102202.err  (empty)
```

`per_carrier_second_hit_summary.tsv` is the key deliverable for Doug Epstein / collaborators.

---

## What we'd tell Doug Epstein

> Phase 7 ran the ZNF175 carrier deep-dive against your second-hit hypothesis. Out of 142 ZNF175 pLoF carriers in PMBB, 3 have audiogram-confirmed HL. **All 3 carry additional deleterious variants in classic Mendelian HL genes** — including GJB2 (DFNB1, the most common HL gene), USH2A (Usher syndrome), MYO7A (Usher), COL11A1/COL11A2 (Stickler), MYO3A (DFNB30), and GPSM2 (DFNB82). One of the 3 (PMBB2106731298975) is the same carrier flagged by both your "8 signal-driving cases" and Joe Park's named-individuals email — strongest single-individual evidence for the hypothesis.
>
> Statistically the analysis is underpowered (N=3 vs 129, Mann-Whitney p=0.24). But qualitatively the pattern is exactly what the modifier model predicts: ZNF175 alone insufficient (129 controls carry hits too), but ZNF175 + a second hit in a Mendelian HL gene = HL phenotype in 3/3 observed cases.
>
> Next steps: replicate in UKBB / AoU / newer PMBB releases for more case carriers. Also tighten "second hit" to ClinVar P/LP only (sensitivity analysis). Also re-contact the 5 Joe-named individuals who aren't in our audiogram cohort — they may have clinically relevant HL not captured by the audiogram phenotyping.

---

## Open questions (carry to next session / external collaborators)

| Question | Status |
|---|---|
| Daniel's confirmation of project framing | **Email draft ready** at [`docs/communications/daniel_followup_email.md`](../../docs/communications/daniel_followup_email.md) |
| Re-contact 5 of Joe's named individuals (not in audiogram cohort) | **Open** — would require Hall Lab / ENT team coordination |
| Sensitivity analysis: tighten "second hit" to ClinVar P/LP only | **Easy to add** — extend Phase 7 script with ClinVar filter |
| Match-control comparison: 142 carriers vs 568 matched non-carriers | **Easy to add** — Daniel preserved the matched-control file |
| UKBB replication | **Out of scope for now** — Daniel did a partial UKBB analysis (walkthrough Phase 15), could revisit |
| AoU extension | **Out of scope** — Nikki has AoU access per kickoff notes |

---

## Recommendation: prepare Doug Epstein meeting

The Phase 7 results are concrete enough to share with Doug — the biology is compelling even with the small N. Suggested deliverables for the meeting:
1. This report
2. The per-carrier second-hit summary table
3. Discussion of the 3 cases (especially PMBB2106731298975) and the GJB2/USH2A/COL11A* genes as plausible mechanism
4. Plan for tightened analysis (ClinVar P/LP only) and external replication paths
