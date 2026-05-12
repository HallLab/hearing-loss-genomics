# Analysis Plan — Phase 1 (Replication & ZNF175 Deep-dive)

> Working document. Updated as we go. Checkbox items in priority order. Reach out to Molly with questions before executing major changes.

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

### ZNF175 deep-dive
- [ ] Identify the 8 signal-driving cases from the original burden test (their ZNF175 variants + PMBB IDs)
- [ ] Extract their complete exomes
- [ ] Annotate variants in all 173 known HL genes for these 8 individuals
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

## Open questions (to raise with Molly / Doug)

- What exact filtering criteria define the "8 signal-driving cases"? (pLoF only? or pLoF + missense? what allele-frequency cutoff?)
- Is the 140-case curated cohort still considered the working cohort, or should we re-curate with PMBB v3?
- For the audiogram-only quantitative-trait analysis on the 4K cohort: in scope for Phase 1 or Phase 2?
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
