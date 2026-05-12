# Analysis Plan — Phase 1 (Replication & ZNF175 Deep-dive)

> Working document. Updated as we go. Checkbox items in priority order. Reach out to Molly with questions before executing major changes.

## Objectives

1. Reproduce the Hui et al. 2023 PMBB pipeline on current PMBB release (v3) to confirm we can replicate the published numbers.
2. Re-run the ZNF175 burden test with updated rare-variant methods (Charlene's recent recommendations).
3. Deep-dive the 8 signal-driving exomes to identify candidate second-hit variants in other HL genes.
4. Document the pipeline thoroughly so it's reusable for downstream analyses (UKBB, AoU, audiogram subcohort).

---

## Task list

### Setup & onboarding
- [ ] Confirm Daniel Hui's scripts are complete and accessible at the project root
- [ ] Read scripts end-to-end; produce `docs/pipeline_walkthrough.md` documenting each step's input/output/parameters
- [ ] Reach out to Daniel Hui (still at Penn, Tishkoff lab — same email) for any clarifications
- [ ] Reach out to Joe Park if methodology questions arise about the original burden test
- [ ] Inventory `data/` — confirm we have: PMBB v3 genotype data, exome data, phecode 389 definitions, audiogram-PMBB ID linkage
- [ ] Verify `venv` has required packages; create `requirements.txt` / `environment.yml` if missing

### Replication
- [ ] Recreate case/control definitions per Hui et al.: phecode 389 (≥2 instances = case, 0 = control, 1 = missing)
- [ ] Cross-check against audiograms (subset of ~1,917 in paper). Confirm 65% cases / 27% controls have audiogram-defined HL.
- [ ] Reproduce paper's Table S2 (sample counts after each QC filter): target 35,397 controls + 1,110 cases for the audiogram-hybrid cohort
- [ ] Re-run rare-variant gene burden across known HL genes (173 genes, paper Fig 2). Check pLoF + missense REVEL>0.6 thresholds.
- [ ] Re-run pLoF-only burden across all genes; confirm ZNF175 reappears at expected significance
- [ ] Recompute heritability estimate (paper reports h² = 4.53%)
- [ ] Note any deviations from published numbers; document expected vs observed in `results/phase1/replication_summary.md`

### ZNF175 deep-dive
- [ ] Identify the 8 signal-driving cases from the original burden test (their ZNF175 variants + PMBB IDs)
- [ ] Extract their complete exomes
- [ ] Annotate variants in all 173 known HL genes for these 8 individuals
- [ ] Filter for predicted-deleterious (pLoF + REVEL>0.6 + ClinVar pathogenic/likely-pathogenic)
- [ ] Tabulate: which HL genes have hits, allele counts, ClinVar status, inheritance pattern (DFNA/DFNB)
- [ ] Cross-reference against Daniel's curated 140-case cohort: do the same second-hit patterns appear in mutation carriers WITH HL but not in those WITHOUT?
- [ ] Statistical test (if N permits): rare-variant burden in HL genes, ZNF175-carriers-with-HL vs ZNF175-carriers-without-HL

### Updated methods (Charlene's recommendations)
- [ ] Get the updated rare-variant pipeline specs from Charlene
- [ ] Run alongside the legacy pipeline for direct comparison
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
- Which updated methods is Charlene using — REGENIE? STAAR? SAIGE-GENE+? confirm before re-running.

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
