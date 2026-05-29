# Results — Reports

This folder contains the per-phase reports for our work on Daniel Hui's PMBB hearing-loss analysis, organized into two chapters.

- **Chapter 1 — Paper Replication.** Reproduce Hui et al. 2023 (PLOS Genetics) end-to-end on PMBB v2, validating our infrastructure against Daniel's preserved intermediates. **Done** — including the paper's headline result (Tables 3 & 4: TCOF1, ESRRB, COL5A1, HMMR, RAPGEF3, NNT).
- **Chapter 2 — ZNF175 Analysis.** Hall Lab × Epstein Lab unpublished follow-up. ZNF175 is **not** in the published paper; it came from Doug Epstein's mouse *Zfp719* biology. Phase 1 of Chapter 2 tests the second-hit hypothesis in PMBB.

See [`../docs/analysis_plan.md`](../docs/analysis_plan.md) for full project context.

---

## Chapter 1 — Paper Replication

Reports live in [`chapter1_paper_replication/`](chapter1_paper_replication/). For a team-facing summary of what was accomplished in Chapter 1, start with [`chapter1_paper_replication/chapter1_summary.md`](chapter1_paper_replication/chapter1_summary.md).

| Phase | Topic | Report | Status |
|---:|---|---|---|
| 1 | Gene list curation (173 known HL genes) | [`chapter1_paper_replication/phase1_replication_report.md`](chapter1_paper_replication/phase1_replication_report.md) | done |
| 2 | SNP ID reconciliation (annotation → extract list) | [`chapter1_paper_replication/phase2_replication_report.md`](chapter1_paper_replication/phase2_replication_report.md) | done |
| 3 | plink genotype extraction (byte-identical to Daniel) | [`chapter1_paper_replication/phase3_replication_report.md`](chapter1_paper_replication/phase3_replication_report.md) | done |
| 4 | Preparatory files (covariates, case/control, weights) | [`chapter1_paper_replication/phase4_replication_report.md`](chapter1_paper_replication/phase4_replication_report.md) | done |
| 5 | First burden test on HL genes (biobin + logistic) | [`chapter1_paper_replication/phase5_replication_report.md`](chapter1_paper_replication/phase5_replication_report.md) | done — ESRRB top hit, Fig 2 replicated |
| 5 | biobin technical reference (for BF4 development) | [`chapter1_paper_replication/phase5_biobin_technical_reference.md`](chapter1_paper_replication/phase5_biobin_technical_reference.md) | reference |
| 6 | Exome-wide all-genes burden (binary, logistic) | [`chapter1_paper_replication/phase6_replication_report.md`](chapter1_paper_replication/phase6_replication_report.md) | done — top 5 match Daniel exactly |
| 7 | Degree-HL linear burden (Paper Tables 3 & 4) | [`chapter1_paper_replication/phase7_degree_hl_burden.md`](chapter1_paper_replication/phase7_degree_hl_burden.md) | done — all 6 paper genes recovered at FDR<0.05 |
| — | **Authoritative p-values** (single source of truth) | [`chapter1_paper_replication/chapter1_authoritative_pvalues.md`](chapter1_paper_replication/chapter1_authoritative_pvalues.md) | reference — full precision, LOKI version annotated |

## Chapter 2 — ZNF175 Analysis

Reports live in [`chapter2_znf175_analysis/`](chapter2_znf175_analysis/).

| Phase | Topic | Report | Status |
|---:|---|---|---|
| 1 | ZNF175 carrier deep-dive + second-hit hypothesis test | [`chapter2_znf175_analysis/phase1_carrier_deep_dive.md`](chapter2_znf175_analysis/phase1_carrier_deep_dive.md) | done (awaiting Daniel's response on framing) |
| 1 | ZNF175 carrier numbers, explained (team-facing) | [`chapter2_znf175_analysis/phase1_carrier_numbers_explained.md`](chapter2_znf175_analysis/phase1_carrier_numbers_explained.md) | done |

---

## Conventions

- Each phase report covers: inputs, methodology, validation against Daniel's preserved outputs (Ch1) or against the second-hit hypothesis (Ch2), deviations (if any), and findings.
- Non-MD outputs (CSVs, plots, intermediate files) are gitignored — only the reports are committed.
- Cross-references use relative paths so links work locally and on GitHub.
- **Script naming:** the orchestration scripts under [`../analysis/daniel/scripts/`](../analysis/daniel/scripts/) still use linear numbering (`run_phase1.sh` … `run_phase7.sh`). `run_phase7.sh` produces Chapter 2 · Phase 1 outputs. This keeps execution-order history obvious without renaming a stable interface.

## Pending

- **Awaiting Daniel Hui's response** on three framing questions before advancing Chapter 2 (see [`../docs/communications/daniel_followup_email.md`](../docs/communications/daniel_followup_email.md)).
- **Candidate Chapter 2 phases:** sensitivity analysis (ClinVar P/LP only), matched-controls comparison, re-contact prep, UKBB / All of Us external replication.
