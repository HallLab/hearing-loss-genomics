# Analysis Plan — Paper Replication + ZNF175 Signal Investigation

> Working document. Updated as we go. Reach out to Molly with questions before executing major changes.
>
> **Project status (2026-05-28):** Chapter 1 (paper replication, Phases 1-7) COMPLETE — Tables 3 & 4 of Hui et al. 2023 reproduced. Chapter 2 Phase 1 (ZNF175 second-hit on 142 carriers) DONE.
>
> **Reframe (2026-05-28, post Molly + Doug kickoff transcript review):** The kickoff meeting's first-step priority was specifically to **reproduce Joe Park's original PMBB v1 (~11k exomes) analysis** and **diagnose why the ZNF175 pLOF burden signal disappeared** when the cohort expanded to PMBB v2 (~43k exomes). Our Chapter 2 Phase 1 (second-hit hypothesis) addressed Doug's biological question but did NOT address Molly's statistical diagnostic ask. This document now reorganizes Chapter 2 to prioritize the diagnostic phase (P2) before continuing biological extensions. ZNF175 is **NOT in the published Hui et al. 2023 paper** — see [`docs/papers/paper_summary_hui2023.md`](papers/paper_summary_hui2023.md).

## Objectives

1. Reproduce the Hui et al. 2023 PMBB pipeline on the **same PMBB release Daniel used (v2, `PMBB-Release-2020-2.0`)** to confirm we can replicate the published numbers. This validates our infrastructure before changing data versions.
2. Re-run the ZNF175 burden test with **SAIGE-GENE+** as the modern rare-variant comparator alongside the legacy biobin pipeline. (Leaning: SAIGE — it's what most groups in the field are converging on.)
3. Deep-dive the 8 signal-driving exomes to identify candidate second-hit variants in other HL genes.
4. Document the pipeline thoroughly so it's reusable for downstream analyses (UKBB, AoU, audiogram subcohort).
5. **(Later phase)** Port the validated pipeline to PMBB v3 (`PMBB-Release-2024-3.0`) and v4 (`PMBB-Release-2026-4.0`) to extend the analysis to newer cohorts.

---

## Task list

### Setup & onboarding
- [x] Confirm Daniel Hui's scripts are complete and accessible at the project root
- [x] Read scripts end-to-end; produce `docs/pipeline_walkthrough.md` documenting each step's input/output/parameters
- [x] Reach out to Daniel Hui (daniel.hui.work@gmail.com) for any clarifications
- [ ] Reach out to Joe Park if methodology questions arise about the original burden test
- [x] Inventory `data/` — confirm we have: PMBB v2 genotype data, exome data, phecode 389 definitions, audiogram-PMBB ID linkage
- [x] Verify `venv` has required packages; create `requirements.txt` / `environment.yml` if missing

### Replication (on PMBB v2 — same release Daniel used)
- [ ] Recreate case/control definitions per Hui et al.: phecode 389 (≥2 instances = case, 0 = control, 1 = missing)
- [ ] Cross-check against audiograms (subset of ~1,917 in paper). Confirm 65% cases / 27% controls have audiogram-defined HL.
- [ ] Reproduce paper's Table S2 (sample counts after each QC filter): target 35,397 controls + 1,110 cases for the audiogram-hybrid cohort
- [ ] Re-run rare-variant gene burden across known HL genes (173 genes, paper Fig 2). Check pLoF + missense REVEL>0.6 thresholds.
- [ ] Re-run pLoF-only burden across all genes; confirm ZNF175 reappears at expected significance
- [ ] Recompute heritability estimate (paper reports h² = 4.53%)
- [ ] Note any deviations from published numbers; document expected vs observed in `results/chapter1_paper_replication/`

### Future phase — port to newer PMBB releases
- [ ] **(Later)** Port pipeline to PMBB v3 (`PMBB-Release-2024-3.0`) — re-derive case/control on the newer cohort, compare results to v2 baseline
- [ ] **(Later)** Port to PMBB v4 (`PMBB-Release-2026-4.0`) if release is finalized and stable
- [ ] Document differences in cohort size, ancestry composition, phecode definitions between releases

## Chapter 2 — ZNF175 Signal Investigation (Hall Lab × Epstein Lab — UNPUBLISHED)

**Background (from kickoff transcript Apr 8 2026):** Joe Park did the original ZNF175 burden analysis in **PMBB v1 (~11k exomes)** — found a pLOF burden association driven by **8 carrier individuals**. He replicated the **missense burden** (not pLOF) in BioMe and Discovery HR cohorts. When PMBB expanded to v2 (~43k exomes), Daniel re-ran Joe's pipeline → **signal lost**. Joe also re-ran it himself → also lost. The project is here to figure out why and to test Doug Epstein's digenic hypothesis (ZNF175 + second-hit modifier locus, based on mouse *Zfp719* biology).

**Cohort timeline (Doug's transcript + PMBB documentation):**

| PMBB version | ~Exomes | Era | Signal status |
|---|---:|---|---|
| v1 (Release 1.0) | ~11,000 | 2018-19 | Joe finds ZNF175 pLOF burden, 8 driving cases |
| v2 (Release 2020-2.0) | 43,731 | 2020-21 | Daniel + Joe re-run → signal lost |
| v3 (Release 2024-3.0) | ~50,000+ | 2024 | Not attempted |

**Molly's first-step ask (paraphrasing transcript 11:46-11:54):** before pursuing biological hypotheses, diagnose why signal disappeared. Hypotheses she explicitly named:
- **H1 — Allele frequency shift:** variants rare in 11k may have crossed MAF threshold in 43k and been dropped. Specifically Doug flagged one variant that *"showed up 115 times in PMBB, and in gnomAD it's much higher frequency"* (transcript 11:53).
- **H2 — Missense vs pLOF:** Joe replicated MISSENSE burden in BioMe/Discovery HR; PMBB v2 re-run may have only tested pLOF. Need to test pLOF, missense, and combined.
- **H3 — Multiple-test correction:** raw p-value on targeted ZNF175 test may be <0.05 even when not FDR-significant exome-wide. Need raw p, not adjusted.
- **H4 — Modern methods:** last 6 months of rare-variant analysis advances may surface signal that biobin misses (Molly references Nikki's experience with Charlene / Tom Coppola group).
- **Preliminary evidence we already have:** chr19:51587727 (the rare variant driving 8 cases) appears as **multi-allelic** in PMBB v2 pVCF (fused with adjacent variant at 51587731). Daniel's QA pipeline likely dropped it. This may be H1 + a multi-allelic mechanism.

### Chapter 2 · Phase 1 — ZNF175 carrier deep-dive ✅ DONE

**Completed 2026-05-15.** Tested Doug's second-hit hypothesis on the 142 ZNF175-pLoF carrier cohort. Findings:
- 3 of 142 carriers are audiogram-confirmed HL cases
- 100% of those 3 carry deleterious variants in classic Mendelian HL genes (GJB2, USH2A, MYO7A, COL11A1, COL11A2, MYO3A, GPSM2)
- Statistical power limited (N=3) — Fisher OR=inf p=1.0, Mann-Whitney p=0.24
- Strong qualitative biological signal consistent with digenic hypothesis

See [`results/chapter2_znf175_analysis/phase1_carrier_deep_dive.md`](../results/chapter2_znf175_analysis/phase1_carrier_deep_dive.md) and [`results/chapter2_znf175_analysis/phase1_carrier_numbers_explained.md`](../results/chapter2_znf175_analysis/phase1_carrier_numbers_explained.md). **This phase addressed Doug's biological hypothesis but did NOT address Molly's statistical diagnostic ask.**

### Chapter 2 · Phase 2 — Signal-loss diagnostic (HIGHEST PRIORITY — Molly's first-step ask)

**Goal:** answer with concrete numbers WHY the ZNF175 pLOF burden signal disappeared between PMBB v1 (~11k) and v2 (~43k). Doable without Joe Park's original code.

**P2.1 — Variant inventory.** List every ZNF175 variant in current PMBB v2 annotation:
- [ ] Pull from `variant-annotation-counts.txt` filtered to ZNF175 region (chr19:51,575,000-51,600,000)
- [ ] For each: MAF in PMBB v2, REVEL score (if missense), function class (pLOF / missense / synonymous / other), ClinVar annotation
- [ ] Flag multi-allelic status (overlap with adjacent variant positions)
- [ ] Cross-reference with Joe's 10-variant pLOF list (from Daniel's runbook lines 213-232)
- [ ] Tabulate: which of Joe's 10 still pass each filter (MAF<0.001, MAF<0.01, REVEL>0.6, not multi-allelic, in annotation)
- [ ] Output: `results/chapter2_znf175_analysis/phase2_variant_inventory.md` + supporting TSV

**P2.2 — The "115 occurrence" variant identification.**
- [ ] Search PMBB v2 for any ZNF175 variant with ~115 carriers — Doug's transcript flag (11:53)
- [ ] Compute its MAF in PMBB v2 (115 / 2 × 43,731 ≈ 0.00132)
- [ ] Check gnomAD MAF for the same position (Doug noted "in gnomAD it's much higher frequency")
- [ ] Confirm: would it pass MAF<0.001? Pass MAF<0.01? Was it in Joe's original list?
- [ ] If it's chr19:51581437 (the high-frequency one Daniel flagged as "should have been removed"), document that match

**P2.3 — Raw-p burden tests on ZNF175 alone (no multiple-test correction).**
- [ ] Run biobin or simple Fisher exact / logistic on ZNF175 carriers vs non-carriers in PMBB v2 (~36,507 hybrid cohort from Phase 7)
- [ ] Three variant sets:
  - pLOF only (Joe's list filtered by current annotation)
  - missense only (REVEL>0.6, REVEL>0.5, no threshold — three sub-tests)
  - pLOF + missense combined
- [ ] Report raw p-values for each — no FDR correction
- [ ] Also test against degHL (linear, 20 PCs) — same three variant sets, raw p
- [ ] Output: `results/chapter2_znf175_analysis/phase2_raw_p_burden.md`

**P2.4 — Diagnostic synthesis.** With outputs from P2.1–P2.3, write a one-page diagnostic explaining:
- Which variants from Joe's original list survive vs which were dropped, and by which filter
- Whether the signal-loss hypothesis (H1 / H2 / H3) is supported by data
- A clear recommendation: is the signal genuinely gone, or recoverable with a specific reanalysis?
- Output: `results/chapter2_znf175_analysis/phase2_signal_loss_diagnostic.md`

### Chapter 2 · Phase 3 — Original Joe Park pipeline reproduction (depends on Joe outreach)

**Goal:** if Phase 2 doesn't fully explain signal loss, reproduce Joe's exact PMBB v1 analysis on PMBB v2 to verify.

**Blocker:** need Joe's pipeline code + variant list. Per kickoff transcript (Doug 11:50-11:54), Joe is still at Penn (Sarah Tishkoff lab). Nikki has his contact.

- [ ] Coordinate with Nikki on Joe outreach. Initial ask:
  - Original pipeline / code repo (whatever format he has — bash, R, Python)
  - Final variant list used in PMBB v1 analysis (positions + filters applied)
  - Slides or write-up from when he originally presented the finding
  - Any preserved intermediate files
- [ ] When materials arrive: stand up Joe's exact pipeline locally
- [ ] Run on PMBB v2 cohort — verify if signal re-emerges or stays lost
- [ ] Output: `results/chapter2_znf175_analysis/phase3_joe_pipeline_reproduction.md`

### Chapter 2 · Phase 4 — Modern rare-variant methods (with Nikki)

**Goal:** apply current rare-variant methodology (post-2025 updates) to test if newer tools surface signal biobin misses.

- [ ] Coordinate with Nikki on modern method recommendations (Molly mentioned "last 6 months of advances" she's been doing with Charlene / Tom Coppola group)
- [ ] **SAIGE-GENE+** as the likely first candidate:
  - Install via fresh conda/mamba env or Apptainer container (stale LPC install at `/project/ritchie/env/modules/saige/1.5.0` likely won't work)
  - Decide configuration (group test type, MAF cutoffs, mask definitions)
- [ ] Run on ZNF175 specifically with the three variant sets from P2.3
- [ ] Run exome-wide as sanity check — confirm we get qualitatively similar paper results (TCOF1, ESRRB, etc.)
- [ ] Compare ZNF175 results: SAIGE-GENE+ vs biobin
- [ ] Output: `results/chapter2_znf175_analysis/phase4_modern_methods.md`

### Chapter 2 · Phase 5 — External replication (UKBB / AoU / PMBB v3)

**Goal:** if signal recovered in P2/P3/P4, test in external cohorts. Doug noted (transcript 11:38) that UKBB phenotype is self-reported (weaker than PMBB audiograms) but provides much larger N.

- [ ] PMBB v3 (~50k+ exomes): re-run the winning Phase 2/3/4 configuration. Daniel preserved relevant scripts; pipeline already walkthrough-documented.
- [ ] UKBB: out of scope for now (we don't have access — per user 2026-05). Document for handoff to Nikki if she has access.
- [ ] All of Us: Nikki has access. Coordinate scope and divide work.
- [ ] Output: `results/chapter2_znf175_analysis/phase5_external_replication.md`

### Chapter 2 · Phase 6 — Doug's 4K audiogram cohort exploration (longer term)

**Doug mentioned (transcript 11:43-11:44):** ~4,000 PMBB individuals have BOTH audiograms AND exomes, with potentially more available via the ENT team. This is a higher-quality subset than the hybrid cohort. Could be the basis for a quantitative-trait analysis specifically for the digenic hypothesis (carrier × second-hit, weighted by audiogram severity).

- [ ] Coordinate with Doug / ENT team on getting the linked audiogram-exome list
- [ ] Re-run Phase 7 (degree-HL linear) restricted to this 4K subset — does ZNF175 signal emerge in a cleaner subset?
- [ ] Re-run Phase 1 (second-hit) within the 4K subset for higher-quality phenotyping

### Reporting
- [ ] Draft results summary (markdown) in `results/chapter1_paper_replication/`
- [ ] Schedule review meeting with Molly + Doug once Priority 1 + 2 are complete
- [ ] Publish clean writeup to Confluence: https://halllab.atlassian.net/wiki/spaces/HLP/ (Andre - Notebook folder)

---

## Open questions (to raise with Daniel / Doug / Molly / Joe Park)

- **For Joe Park** (via Nikki — highest priority, blocks Phase 3):
  - Original PMBB v1 pipeline code / variant list
  - Slides or write-up from when the ZNF175 finding was presented
  - Any preserved intermediate files from the 11k analysis
- **For Daniel** (email draft ready at [`docs/communications/daniel_followup_email.md`](communications/daniel_followup_email.md) — needs revision now that we understand the v1 context better):
  - Was your PMBB v2 attempt to reproduce Joe's v1 result, or independent burden testing?
  - Do you still have the variant list / filter spec from Joe's original analysis?
  - The 140-case curated cohort — was that for the second-hit follow-up after signal-loss, or part of the original burden test?
- For Molly / Doug:
  - Phase 6 (4K audiogram cohort) — when to engage Doug for the ENT linkage?
  - Modern method install path: fresh conda / Apptainer / LPC IT request?

---

## Communication / status updates

- **Nikki:** weekly sync on Phase 2 progress; coordinate Joe Park outreach for Phase 3
- **Molly:** present Phase 2 diagnostic findings at next review meeting (before starting Phase 3 or 4)
- **Doug:** loop in once Phase 2 produces concrete signal-loss diagnostic; engage for Phase 6 audiogram linkage
- **Daniel:** update email draft with corrected v1/v2 context before sending

---

## Deliverables timeline (revised 2026-05-28)

| Phase | Output | Target |
|---|---|---|
| Ch1 | Paper replication + chapter1_summary | ✅ DONE |
| Ch2 P1 | Second-hit hypothesis on 142 carriers | ✅ DONE |
| Ch2 P2 | Signal-loss diagnostic (variant inventory + raw-p tests) | week 1-2 |
| Ch2 P3 | Joe Park pipeline reproduction | depends on Joe outreach |
| Ch2 P4 | SAIGE-GENE+ comparison | week 3-4 (parallel with P3) |
| Ch2 P5 | External replication (AoU / PMBB v3) | TBD after P2-P4 results |
| Ch2 P6 | 4K audiogram cohort analysis | longer-term |
| Review with Molly | Present P2 diagnostic | end of week 2 |
| Review with Doug | Present P2 + biological synthesis | after P2 done |

Adjust as we hit data and find what's already done vs what needs rebuilding.
