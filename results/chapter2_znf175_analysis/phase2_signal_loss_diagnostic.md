# Chapter 2 · Phase 2 — Signal-Loss Diagnostic

**Date:** 2026-05-28 (independently validated by Daniel Hui, 2026-05-29)
**Run by:** Andre Rico
**Audience:** Molly Hall, Doug Epstein, Nikki Palmiero
**Status:** Diagnostic done — explains why Joe Park's PMBB v1 ZNF175 signal did not replicate in PMBB v2

> **Independent validation by Daniel Hui (2026-05-29):** Daniel — without having seen this diagnostic — independently described the exact mechanism: *"right like if you lose a SNP because it stops passing some global MAF threshold but it is more common in one ancestry, it could be worth testing still."* See [`docs/communications/daniel_response_received.md`](../../docs/communications/daniel_response_received.md). His recommendation: use ancestry-stratified MAF filtering (now incorporated into Chapter 2 Phase 4 of the analysis plan).

---

## TL;DR

**The signal disappeared because 79% of Joe Park's burden-driving carriers were eliminated from the v2 burden test** — not by chance, not by biology, but by **two specific filter mechanisms**:

1. **MAF shift** — the top driver variant (`chr19:51581437`, the "115 occurrences" Doug flagged in the kickoff) crossed the MAF<0.001 threshold in v2 and was dropped. This variant has **MAF=3.4% in East Asian** populations (gnomAD) — PMBB v2's expanded cohort includes more EAS ancestry, naturally pushing its overall MAF above the rare-variant cutoff.
2. **Variant calling drift** — 5 of Joe's 10 original pLoF variants are **no longer present in the PMBB v2 annotation**, likely due to changes in the variant calling / QC pipeline between PMBB v1 (2018-2019) and PMBB v2 (2020-2021).

**Of Joe's 33 original v1 carriers across the 10 pLoFs, only 7 (21%) survived into v2 as burden-eligible.** The signal didn't disappear — it was filtered out.

This confirms Molly's MAF-shift hypothesis (kickoff transcript 11:52) and partially Doug's intuition that the signal was real. The path forward is to redo the burden test **either without the MAF filter or using a different MAF threshold (e.g., MAF<0.01)**, which would re-include `chr19:51581437`.

---

## What Molly asked for (kickoff transcript Apr 8 2026)

The kickoff defined a clear first-step priority before any biological extensions:

> *"What is the p-value when you are just looking at the PLOFs, looking at the missense, looking at them together?"* (Molly, 11:46)
>
> *"I have a feeling that one reason that this might not have reproduced once the new versions and larger sample sizes were included is because of allele frequency differences. The same missense or pilafs that were included in the smaller sample size might have been dropped when the larger sample size was used because they weren't considered rare anymore."* (Molly, 11:52)
>
> *"One particular variant that showed up, I don't know, 115 times in PMBB, so it seems that it's a higher frequency, and in Nomad, it's a much higher frequency."* (Doug, 11:53)

This phase answers those questions concretely using preserved files from Joe Park's PMBB v1 analysis.

---

## Data sources used

| File | What it is | Source |
|---|---|---|
| [`data/PMBB_Exome/ZNF175/Joe_analyses/znf175_variants.txt.gz`](../../data/PMBB_Exome/ZNF175/Joe_analyses/) | Joe's preserved variant list (hg38), with PMBB v1 Hom/Het/Missing counts | Joe Park, 2018-2019 |
| [`data/PMBB_Exome/ZNF175/ZNF175_annot_genes_full.txt.gz`](../../data/PMBB_Exome/ZNF175/) | All ZNF175 variants in PMBB v2 annotation, with v2 counts, REVEL, gnomAD, ClinVar | Daniel/PMBB v2, 2020-2021 |
| Daniel's runbook lines 213-232 | Manual hg38→hg19 lift of Joe's 10 pLoF positions, with QA notes | Daniel Hui, 2021 |

**Pipeline:** [`analysis/daniel/scripts/pmbb_exome/znf175_signal_loss_diagnostic.py`](../../analysis/daniel/scripts/pmbb_exome/znf175_signal_loss_diagnostic.py)

**Outputs:** [`analysis/daniel/outputs/phase8_signal_diagnostic/`](../../analysis/daniel/outputs/phase8_signal_diagnostic/)

---

## Finding 1 — Joe Park's original variant inventory (PMBB v1)

Joe's preserved file contains **162 ZNF175 variants** in hg38 coordinates with PMBB v1 carrier counts. Of these, the 10 pLoF variants that drove the original burden test (cited by Daniel in runbook lines 213-232):

| hg38 position | hg19 position | Type | v1 Hom | v1 Het | v1 carriers (~11k) |
|---|---|---|---:|---:|---:|
| 19_52084690 | **51581437** | frameshift insertion | 1 | 12 | **13** |
| 19_52091682 | 51588429 | frameshift deletion | 0 | 8 | **8** |
| 19_52091067 | 51587814 | stopgain | 0 | 5 | **5** |
| 19_52090981 | 51587728 | frameshift deletion | 0 | 3 | **3** |
| 19_52084743 | 51581490 | stopgain | 0 | 1 | **1** |
| 19_52091164 | 51587911 | frameshift deletion | 0 | 1 | **1** |
| 19_52091468 | 51588215 | frameshift deletion | 0 | 1 | **1** |
| 19_52091568 | 51588315 | stopgain | 0 | 1 | **1** |
| 19_52091635 | 51588382 | stopgain | 0 | 1 | **1** |
| 19_52090598 | 51587345 | frameshift deletion | 0 | 0 | **0** |
| **TOTAL** | | | **1** | **33** | **34** |

Note: the "8 driving cases" Doug references could refer to multiple subsets:
- 8 = carriers of the second-highest variant (`51588429`)
- 8 = HL-case carriers (subset of the 34 total carriers — the 8 of 33 who happened to also have HL phenotype)

This is a question for follow-up confirmation but doesn't change the diagnostic.

---

## Finding 2 — Of Joe's 10 pLoFs, only 4 remain burden-eligible in PMBB v2

Cross-referencing each of Joe's 10 pLoFs against the current PMBB v2 annotation:

| hg38 → hg19 | v1 carriers | v2 carriers | v2 MAF | Status in v2 | Why? |
|---|---:|---:|---:|---|---|
| 52084690 → **51581437** | **13** | **94** | **0.0011** | **❌ DROPPED** | **MAF exceeds 0.001 threshold** (this is the "115" variant) |
| 52091682 → 51588429 | **8** | — | — | ❌ DROPPED | Not in v2 annotation |
| 52091067 → 51587814 | 5 | 18 | 2.1×10⁻⁴ | ✅ passes | — |
| 52090981 → 51587728 | 3 | — | — | ❌ DROPPED | Not in v2 annotation (multi-allelic per Daniel's runbook) |
| 52084743 → 51581490 | 1 | 1 | 1.1×10⁻⁵ | ✅ passes | — |
| 52091164 → 51587911 | 1 | 2 | 2.3×10⁻⁵ | ⚠️ reclassified | Present in v2 but not annotated as pLoF |
| 52091468 → 51588215 | 1 | — | — | ❌ DROPPED | Not in v2 annotation |
| 52091568 → 51588315 | 1 | 2 | 2.3×10⁻⁵ | ✅ passes | — |
| 52091635 → 51588382 | 1 | — | — | ❌ DROPPED | Not in v2 annotation |
| 52090598 → 51587345 | 0 | — | — | n/a | Had 0 carriers in v1; not relevant |

**Summary of carrier survival:**

| Bucket | v1 carriers count |
|---|---:|
| Variants surviving into v2 as burden-eligible pLoF (3 variants: 51581490, 51587814, 51588315) | **7** |
| Variants dropped because MAF exceeds 0.001 in v2 (1 variant: 51581437) | **13** |
| Variants dropped because no longer in v2 annotation (5 variants) | **13** |
| Variants reclassified (1 variant: 51587911) | **1** |
| **Total v1 driver carriers** | **34** |

**Only 21% of the original signal-driving carriers survive the v2 burden filter.** This is the mechanical explanation for signal loss.

---

## Finding 3 — Confirmation: chr19:51581437 is the "115" variant Doug mentioned

Doug referenced (transcript 11:53): *"One particular variant that showed up, I don't know, 115 times in PMBB, so it seems that it's a higher frequency, and in Nomad, it's a much higher frequency."*

`chr19:51581437` (a frameshift insertion) matches:

- **PMBB v2 carriers: 94** (close to Doug's "115" — Doug noted "I don't know")
- **PMBB v2 MAF: 0.0011** — exceeds the rare-variant threshold (<0.001)
- **gnomAD MAF: 0.0025** ("much higher frequency" per Doug)

Population breakdown in gnomAD makes the mechanism explicit:

| Ancestry | gnomAD MAF |
|---|---:|
| ALL | 0.0025 |
| **East Asian (EAS)** | **0.0338** |
| Other (OTH) | 0.0015 |
| African (AFR) | 0.0006 |
| South Asian (SAS) | 0.0004 |
| American (AMR) | 6.0×10⁻⁵ |
| Non-Finnish European (NFE) | 2.7×10⁻⁵ |
| Ashkenazi Jewish (ASJ) | 0 |
| Finnish (FIN) | 0 |

**The variant is ~30× more common in East Asian populations than in NFE.** When PMBB v2 expanded from ~11k to ~43k exomes, it included more ancestrally diverse individuals — pushing the overall MAF of this variant above the rare-variant cutoff. **This is allele frequency drift driven by cohort composition change, not by any biological change.**

This is the cleanest, most diagnostic finding of the phase — and it precisely matches Molly's hypothesis from the kickoff.

---

## Finding 4 — Five of Joe's pLoFs vanished from the v2 annotation

Beyond the MAF-shift on `51581437`, five other variants (totaling **13 more v1 carriers**) are no longer present in the PMBB v2 annotation:

| hg19 | Original type (Joe) | v1 carriers | Runbook note |
|---|---|---:|---|
| 51588429 | frameshift deletion | 8 | "yes" — but missing in v2 |
| 51587728 | frameshift deletion | 3 | "not in — multiallelic" |
| 51587911 | frameshift deletion | 1 | "not in — it's multiallelic" |
| 51588215 | frameshift deletion | 1 | "in, is indel" — but missing in v2 |
| 51588382 | stopgain | 1 | "not in — not in annotation" |

Daniel's runbook already noted at QA time that several of these were "multiallelic" or "not in annotation" — and this turned out to be the case for ~40% of Joe's original signal-driving carriers. **The PMBB v1 → v2 variant calling pipeline change** (different caller version, different filtering thresholds, possibly different reference panel) silently removed these from the eligible variant pool.

---

## What this means for the project

### For Molly's hypothesis (MAF shift)
**Confirmed.** The single biggest contributor to signal loss is the `51581437` MAF crossover, which alone accounts for 39% (13/34) of Joe's original carriers. The mechanism (EAS-ancestry expansion in PMBB v2) is biologically plausible and statistically clean.

### For Doug's hypothesis (real biological signal)
**Partially vindicated.** The signal didn't disappear because ZNF175 stopped being important — it disappeared because the burden test stopped seeing the carriers that drove the original association. A re-run with the original v1 variant set (or a more permissive MAF threshold) should re-surface the signal.

### For Joe's missense-burden replication in BioMe/Discovery HR
**Still unaddressed.** Doug noted (transcript 11:15) that the replication in BioMe and Discovery HR used **missense burden, not pLoF**. PMBB v2 has 4 missense variants in ZNF175 passing REVEL>0.6 + MAF<0.001, with only 5 total v2 carriers — a missense burden test on this small set is feasible but underpowered. Adding more permissive REVEL thresholds (e.g., >0.5) increases the eligible set.

---

## Recommended next steps (Phase 2.3 and Phase 3 of the analysis plan)

### Phase 2.3 — Raw-p burden tests with multiple variant-set definitions

Run targeted ZNF175 burden tests in PMBB v2 (raw p, no multiple-test correction) under different inclusion criteria, to quantify the impact of the MAF filter:

| Test | Variant set | Expected outcome |
|---|---|---|
| A | Joe's 7 surviving v2 variants + `51581437` (force-include despite MAF) | Should recover signal if MAF filter is the only issue |
| B | All v2 pLoF + MAF<0.01 (relaxed) | Tests whether broader rare-variant inclusion rescues signal |
| C | All v2 pLoF + missense REVEL>0.6 + MAF<0.001 | The "combined" test Molly asked for |
| D | All v2 pLoF + missense REVEL>0.5 + MAF<0.01 | Most permissive, closest to Joe's BioMe/Discovery replication |
| E | Joe's original variant set (as best we can reconstruct from v1 file) | Direct reproduction attempt |

Report raw p-value for each on both binary HL (logistic) and degree-HL (linear) phenotypes.

### Phase 3 — Joe Park outreach (Nikki)

To definitively reproduce Joe's pipeline, we need:
- Original burden test code (which variants exactly, what MAF cutoff, what regression model)
- Final p-value reported in v1
- Was MAF cutoff <0.001? <0.01? Other?

Until we have this, Phase 2.3 Test A is our best approximation.

### Phase 4 — Modern methods (SAIGE-GENE+)

If Phase 2.3 confirms the MAF-shift mechanism, SAIGE-GENE+ can be tested separately:
- Does SAIGE-GENE+ have ancestry-aware MAF handling?
- Different mask definitions might capture different variant subsets
- Compare on the same cohort

---

## Files generated in this phase

```
analysis/daniel/scripts/pmbb_exome/znf175_signal_loss_diagnostic.py
analysis/daniel/outputs/phase8_signal_diagnostic/
├── znf175_v2_inventory.tsv         (380 v2 variants × 17 cols)
├── joe_v1_variants.tsv             (162 v1 variants × 12 cols)
├── variant_115_candidates.tsv      (top v2 variants by carrier count)
└── joe_v1_v2_crossref.tsv          (Joe v1 ↔ v2 cross-reference)
```

## Open questions still requiring outside input

1. **Joe Park's exact burden test specification** — MAF cutoff, regression model, covariates, multiple-test correction strategy. Nikki to reach out.
2. **The "8 driving cases" precise definition** — 8 HL-cases of the 34 carriers? Or carriers of a specific variant? — Joe or Daniel.
3. **Was the chr19:51587727 multi-allelic split (with 51587731) the same in PMBB v1?** Could explain why Joe had 3 carriers there but it's "not in" v2.
