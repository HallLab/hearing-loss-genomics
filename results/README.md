# Results — Replication & ZNF175 Deep-Dive Reports

This folder contains the per-phase replication reports for our work on Daniel Hui's PMBB hearing-loss analysis.

The project has two halves:

- **Phases 1-6 — paper replication.** Reproduce Hui et al. 2023 (PLOS Genetics) end-to-end on PMBB v2, validating our infrastructure against Daniel's preserved intermediates.
- **Phase 7+ — unpublished ZNF175 follow-up.** Hall Lab × Epstein Lab extension. ZNF175 is **not** in the published paper; it came from Doug Epstein's mouse *Zfp719* biology. Phase 7 tests the second-hit hypothesis in PMBB.

See [`../docs/analysis_plan.md`](../docs/analysis_plan.md) for full project context.

---

## Phase reports

| Phase | Topic | Report | Status |
|---:|---|---|---|
| 1 | Gene list curation (173 known HL genes) | [`phase1/phase1_replication_report.md`](phase1/phase1_replication_report.md) | done |
| 2 | SNP ID reconciliation (annotation → extract list) | [`phase2/phase2_replication_report.md`](phase2/phase2_replication_report.md) | done |
| 3 | plink genotype extraction (byte-identical to Daniel) | [`phase3/phase3_replication_report.md`](phase3/phase3_replication_report.md) | done |
| 4 | Preparatory files (covariates, case/control, weights) | [`phase4/phase4_replication_report.md`](phase4/phase4_replication_report.md) | done |
| 5 | First burden test on HL genes (biobin + logistic) | [`phase5/phase5_replication_report.md`](phase5/phase5_replication_report.md) | done — ESRRB top hit, Fig 2 replicated |
| 5 | biobin technical reference (for BF4 development) | [`phase5/biobin_technical_reference.md`](phase5/biobin_technical_reference.md) | reference |
| 6 | Exome-wide all-genes burden | [`phase6/phase6_replication_report.md`](phase6/phase6_replication_report.md) | done — top 5 match Daniel exactly |
| 7 | ZNF175 carrier deep-dive — second-hit hypothesis | [`phase7/phase7_replication_report.md`](phase7/phase7_replication_report.md) | done (awaiting Daniel's response on framing) |
| 7 | ZNF175 carrier numbers, explained (team-facing) | [`phase7/phase7_carrier_numbers_explained.md`](phase7/phase7_carrier_numbers_explained.md) | done |

---

## Conventions

- Each `phaseN_replication_report.md` covers: inputs, methodology, validation against Daniel's preserved outputs, deviations (if any), and findings.
- Non-MD outputs (CSVs, plots, intermediate files) are gitignored — only the reports are committed.
- Cross-references use relative paths so links work locally and on GitHub.

## Pending

- **Awaiting Daniel Hui's response** on three framing questions before advancing to Phase 8.
- **Candidate Phase 8+ work:** sensitivity analysis (ClinVar P/LP only), matched-controls comparison, re-contact prep, UKBB / All of Us external replication.
