# Analysis Plan — Paper Replication + ZNF175 Extension

> Working document. Updated as we go. Reach out to Molly with questions before executing major changes.
>
> **Project status (2026-05-13):** Phases 1-6 of paper replication COMPLETE and validated (byte-equivalent on top hits). Phase 7+ pivots to the **unpublished ZNF175 extension** — the gene is **NOT in the published Hui et al. 2023 paper** but came from Doug Epstein's mouse biology brought to Daniel for PMBB investigation. See [`docs/papers/paper_summary_hui2023.md`](papers/paper_summary_hui2023.md) for the verified scope of the published paper.

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
- [ ] Note any deviations from published numbers; document expected vs observed in `results/phase1/replication_summary.md`

### Future phase — port to newer PMBB releases
- [ ] **(Later)** Port pipeline to PMBB v3 (`PMBB-Release-2024-3.0`) — re-derive case/control on the newer cohort, compare results to v2 baseline
- [ ] **(Later)** Port to PMBB v4 (`PMBB-Release-2026-4.0`) if release is finalized and stable
- [ ] Document differences in cohort size, ancestry composition, phecode definitions between releases

### Phase 7 — ZNF175 deep-dive (UNPUBLISHED extension — Doug Epstein → Daniel)

**Reframe (2026-05-13):** ZNF175 is **NOT in the published paper** and **does NOT reach FDR significance in any of Daniel's preserved burden tests** (4 phenotype variants + degree-of-HL meta — checked in `data/PMBB_Exome/allGenes/HL_*/meta_results/`). This phase is the biological deep-dive that Daniel started but didn't publish, motivated by Doug Epstein's mouse work on *Zfp719* (ZNF175's syntenic ortholog showing HL phenotype in knockout).

**Strategy: skip meta-analysis (won't surface ZNF175 anyway), go directly to gene-specific carrier deep-dive.**

- [ ] Identify carriers of ZNF175 pLoF variants (especially chr19:51587727 and chr19:51581437 — referenced in Daniel's runbook lines 193-196)
- [ ] Validate the "8 signal-driving cases" definition with Daniel (open question — see [`results/phase1/phase1_replication_report.md`](../results/phase1/phase1_replication_report.md) and email draft to Daniel at [`docs/communications/daniel_followup_email.md`](communications/daniel_followup_email.md))
- [ ] Extract their complete exomes
- [ ] Annotate variants in all 173 known HL genes for these carriers
- [ ] Filter for predicted-deleterious (pLoF + REVEL>0.6 + ClinVar pathogenic/likely-pathogenic)
- [ ] Tabulate: which HL genes have hits, allele counts, ClinVar status, inheritance pattern (DFNA/DFNB)
- [ ] Cross-reference against Daniel's curated 140-case cohort: do the same second-hit patterns appear in mutation carriers WITH HL but not in those WITHOUT?
- [ ] Statistical test (if N permits): rare-variant burden in HL genes, ZNF175-carriers-with-HL vs ZNF175-carriers-without-HL

### Updated methods comparison (SAIGE-GENE+)
- [ ] Install a current SAIGE-GENE+ build — LPC has a stale SAIGE 1.5.0 install at `/project/ritchie/env/modules/saige/1.5.0` (no `module load` file; conda env recipe from 2020). Current SAIGE-GENE+ (≥1.1+) needs a fresh conda/mamba env or an Apptainer container.
- [ ] Run SAIGE-GENE+ alongside the legacy biobin pipeline for direct comparison
- [ ] Document any signal differences (gained / lost associations)

### Reporting
- [ ] Draft results summary (markdown) in `results/phase1/`
- [ ] Schedule review meeting with Molly + Doug once Priority 1 + 2 are complete
- [ ] Publish clean writeup to Confluence: https://halllab.atlassian.net/wiki/spaces/HLP/ (Andre - Notebook folder)

---

## Open questions (to raise with Daniel / Doug / Molly)

- **For Daniel** (email draft ready at [`docs/communications/daniel_followup_email.md`](communications/daniel_followup_email.md)):
  - Was ZNF175 introduced by Doug from his mouse Zfp719 work, or did it emerge from your PMBB burden tests?
  - What exact filter defines the "8 signal-driving cases"? Carriers of chr19:51587727 specifically?
  - Is the 140-case curated cohort still canonical? Where is it documented?
- For Molly / Doug:
  - For the audiogram-only quantitative-trait analysis on the 4K cohort: in scope or deferred?
- ~~Which modern rare-variant method should be the comparator?~~ → **SAIGE-GENE+** (Andre's pick, 2026-05-12). Open sub-question: install via fresh conda env, Apptainer container, or LPC IT request for a new module build?

---

## Deliverables timeline (proposed)

| Milestone | Target |
|---|---|
| Pipeline walkthrough doc | week 1–2 |
| Replication of paper numbers | week 3–4 |
| 8-exome deep-dive results | week 5–6 |
| Updated methods comparison | week 7–8 |
| Phase 1 writeup + Confluence publish | week 9 |
| Review with Molly & Doug | end of week 9 |

Adjust as we hit the data and find what's already done vs what needs rebuilding.
