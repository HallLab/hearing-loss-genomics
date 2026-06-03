# Analysis Plan — Hearing Loss Genomics

> **Working document.** Updated as we go.
>
> **Last revised:** 2026-06-03 (post Daniel response + LOKI restore in progress)

---

## How this doc is organized

The project is split into **two chapters**, each with sequential **phases**:

| Chapter | What | Where it lives |
|---|---|---|
| **Chapter 1** — Paper Replication | Reproduce Hui et al. 2023 on PMBB v2 (~43k exomes) | [`results/chapter1_paper_replication/`](../results/chapter1_paper_replication/) |
| **Chapter 2** — ZNF175 Signal Investigation | Diagnose why Joe Park's PMBB v1 ZNF175 signal didn't replicate in v2; test digenic hypothesis | [`results/chapter2_znf175_analysis/`](../results/chapter2_znf175_analysis/) |

Orchestration scripts in `analysis/daniel/scripts/` use **linear** numbering by execution order.


---

## Current state

| Item | Status | Owner | Detail |
|---|---|---|---|
| **Chapter 1** — all 7 phases | ✅ DONE | Andre | All 6 paper genes recovered at FDR<0.05 |
| **Chapter 2 Phase 1** — Second-hit on 142 carriers | ✅ DONE | Andre | 3/3 HL cases carry second hits in Mendelian HL genes |
| **Chapter 2 Phase 2** — Signal-loss diagnostic | ✅ DONE | Andre | Independently validated by Daniel 2026-05-29 |
| **Chapter 2 Phase 3** — Joe Park pipeline reproduction | ⏸️ BLOCKED | Nikki (outreach) | Email draft ready, pending send |
| **Chapter 2 Phase 4** — Extended biobin burden tests | 📋 NEXT | ___ | Multiple variant sets + MAF strategies |
| Chapter 2 Phase 5 — External replication | 🔭 LATER | ___ | After P3/P4 give signal |
| Chapter 2 Phase 6 — 4K audiogram cohort | 🔭 LATER | ___ | Needs ENT linkage |
| **Infra:** LOKI 2022 restore | ⏸️ IN PROGRESS | Scott Dudek (ticket) | Should be recoverable |

Status icons: ✅ done · ⏸️ blocked / waiting · 📋 ready to start · 🔭 longer-term

---

## Chapter 1 — Paper Replication ✅ CLOSED

**Goal:** Reproduce Hui et al. 2023 on PMBB v2 (`PMBB-Release-2020-2.0`, 43,731 exomes) to confirm our infrastructure reproduces the published numbers before extending to ZNF175 work.

**Outcome:** All 6 headline paper genes recovered at FDR<0.05:
- **TCOF1, ESRRB** (Table 3 — known HL genes)
- **COL5A1, HMMR, RAPGEF3, NNT** (Table 4 — novel candidates)

### Phase-by-phase summary

The paper's 19-phase runbook (Daniel's `runbook_hui2023.txt`) was consolidated into 7 chapter-internal phases. Each addresses a specific scientific question and was validated against Daniel's preserved intermediates.

| Phase | Goal | Phenotype / model | Validation outcome |
|---|---|---|---|
| **1** — Gene list curation | Build the 173-gene known HL set + filter variant annotation to functional classes (pLoF + missense REVEL>0.6) | — (setup) | Set-equal to Daniel's `annot_genes_full_funcToInclude.txt`  (11,661 variants) and `Hearing_loss_genes.txt.gz` (Shadi list)|
| **2** — SNP ID reconciliation | Map annotation IDs ↔ pVCF IDs; produce plink `.extract` file with 9,667 SNPs | — (setup) | md5-identical to Daniel's `.extract` |
| **3** — plink genotype extraction | Per-chromosome `.bed`/`.bim`/`.fam` from raw pVCFs + IBD-based individual filter | — (data prep) | **chr21 byte-identical to Daniel** (md5 verified) |
| **4** — Preparatory files | Build `cases_control.txt` (binary HL phenotype), `covs.txt` (4 PCs), `gene_list_regions.txt` | Binary HL ([Phenotype A](../results/chapter1_paper_replication/chapter1_authoritative_pvalues.md#phenotype-a--binary-hl-hl_needaud--used-by-phase-5-and-phase-6)) | Set-equal to Daniel (Python 2→3 float repr difference on `AgeSq` only) |
| **5** — HL gene burden | biobin logistic on 173 HL genes, 4 PCs (binary HL) | Binary HL, logistic | **ESRRB p=8.6308×10⁻⁵ byte-identical to Daniel.** Top hit. Paper's Fig 2 reproduced. |
| **6** — Exome-wide burden | biobin logistic on all ~18,000 genes, 4 PCs (binary HL) | Binary HL, logistic | Top 5 non-LOC genes byte-equivalent to Daniel (DNAJC8, UPK3BL1, COL5A1, BOD1, ZNF670); ~30-50% LOKI version drift on other HL genes |
| **7** — Degree-HL linear burden | Light-mode replication of paper's primary test (Tables 3 & 4) using Daniel's preserved STable | [Degree-HL 0-4](../results/chapter1_paper_replication/chapter1_authoritative_pvalues.md#phenotype-b--degree-of-hl-0-4-deghl--used-by-phase-7), R `lm()`, 20 PCs | **All 6 paper genes recovered at FDR<0.05.** Daniel's STable β/p within 5-10% of published. |

### Key decisions made during Chapter 1

| Decision | Rationale |
|---|---|
| Use PMBB v2 (Daniel's release), not v3/v4 | Validates infra against published numbers before changing data |
| Default to "light mode" (use Daniel's preserved intermediates as inputs) | Validates downstream logic in hours vs days of compute; heavy mode reserved for chr21 spot-check |
| Phase 7 in light mode only (no heavy biobin → R `lm()` rerun) | Daniel preserved the final supplementary table; heavy rerun would take days for no scientific gain |
| Phenotype `HL_needAud` used in Phase 5/6 (not the 4 phenotype variants Daniel ran) | The paper's headline analyses use `HL_needAud`; the others are sensitivity checks documented for future use |

### Documentation entry points

| Doc | Use when |
|---|---|
| [`chapter1_summary.md`](../results/chapter1_paper_replication/chapter1_summary.md) | Team-facing narrative — 3 analyses, 6 genes, why each phase matters |
| [`chapter1_authoritative_pvalues.md`](../results/chapter1_paper_replication/chapter1_authoritative_pvalues.md) | **Single source of truth for p-values** — full precision, LOKI version annotated, phenotype definitions |
| Individual `phaseN_*.md` reports in `results/chapter1_paper_replication/` | Technical detail per phase (methodology, validation, gotchas) |
| Drift trilogy in [`docs/andre/`](andre/) | Plain-language explanations: Phase 5 byte-equivalence, Phase 6 LOKI drift, Phase 7 iteration drift |

### Caveats / open infrastructure items

- **LOKI version drift on Phase 6.** Our runs used `loki-20230816`; Daniel used `loki-20220926`. ~30-50% p-value drift on most HL genes (bounded, non-directional). Phase 5 and Phase 7 NOT affected. **LOKI 2022 restore in progress via Scott Dudek** (see Infrastructure side-quests below).
- **Phase 7 carrier counts differ slightly from paper** (e.g., NNT: 840 in STable vs 1,013 in paper). Cause: Daniel's STable is an intermediate iteration, not the final paper-submitted run. β within 5-10%, p within ~1 order of magnitude. Documented in [`chapter1_authoritative_pvalues.md`](../results/chapter1_paper_replication/chapter1_authoritative_pvalues.md).

---

## Chapter 2 — ZNF175 Signal Investigation

**Background (from kickoff transcript Apr 8 2026):** Joe Park did the original ZNF175 burden analysis in **PMBB v1 (~11k exomes)** in 2018-2019 — found a pLOF burden association driven by **8 carrier individuals**. He replicated the **missense burden** (not pLOF) in BioMe and Discovery HR cohorts. When PMBB expanded to v2 (~43k exomes), Daniel re-ran Joe's pipeline → **signal lost**. Joe also re-ran it himself → also lost. The project is here to figure out why and to test Doug Epstein's digenic hypothesis (ZNF175 + second-hit modifier locus, based on mouse *Zfp719* biology).

**ZNF175 is NOT in the published Hui et al. 2023 paper** (verified via PDF) — see [`docs/papers/paper_summary_hui2023.md`](papers/paper_summary_hui2023.md).

**Cohort timeline:**

| PMBB version | ~Exomes | Era | Signal status |
|---|---:|---|---|
| v1 (Release 1.0) | ~11,000 | 2018-19 | Joe finds ZNF175 pLOF burden, 8 driving cases |
| v2 (Release 2020-2.0) | 43,731 | 2020-21 | Daniel + Joe re-run → signal lost |
| v3 (Release 2024-3.0) | ~50,000+ | 2024 | Not attempted (would be Phase 5) |

**Molly's first-step hypotheses (kickoff transcript 11:46-11:54):**

| Hyp | Mechanism | Status after Phase 2 |
|---|---|---|
| **H1 — MAF shift** | Variants rare in 11k crossed MAF threshold in 43k → dropped | ✅ **Confirmed dominant cause** — chr19:51581437 (EAS MAF=3.4%, global 0.0025) excluded by MAF<0.001 filter |
| **H2 — Missense vs pLOF** | Joe replicated missense burden in BioMe/Discovery HR; PMBB v2 re-run may have only tested pLOF | ⏸️ To test in Phase 4 (multiple variant sets) |
| **H3 — Multiple-test correction** | Raw p-value on targeted ZNF175 may be <0.05 even when not FDR-significant exome-wide | ⏸️ To test in Phase 4 (raw-p reporting) |
| **H4 — Newer rare-variant methods** | Molly noted in kickoff that the rare-variant field has advanced in the last ~6 months | 🔭 Deferred — separate scoping with Nikki after Phase 4 |
| **H5 — Variant calling drift** | (new, from our Phase 2) — 5 of Joe's 10 pLOFs not in v2 annotation | ✅ Documented but not closable without Joe's exact variant list |

---

### Chapter 2 · Phase 1 — ZNF175 carrier deep-dive ✅ DONE

**Completed 2026-05-15.** Tested Doug's digenic hypothesis on the 142 ZNF175-pLoF carrier cohort.

**Findings:**
- 3 of 142 carriers are audiogram-confirmed HL cases
- **100%** of those 3 carry deleterious variants in classic Mendelian HL genes (GJB2, USH2A, MYO7A, COL11A1, COL11A2, MYO3A, GPSM2)
- Statistical power limited (N=3) — Fisher OR=inf p=1.0, Mann-Whitney p=0.24
- Strong qualitative biological signal consistent with digenic hypothesis

**Reports:**
- [`phase1_carrier_deep_dive.md`](../results/chapter2_znf175_analysis/phase1_carrier_deep_dive.md) (technical)
- [`phase1_carrier_numbers_explained.md`](../results/chapter2_znf175_analysis/phase1_carrier_numbers_explained.md) (team-facing)

**Caveat:** addresses Doug's biological question but did NOT address Molly's statistical diagnostic ask. Phase 2 closes that gap.

---

### Chapter 2 · Phase 2 — Signal-loss diagnostic ✅ DONE

**Completed 2026-05-28. Independently validated by Daniel Hui 2026-05-29.** Answered with concrete numbers why Joe's PMBB v1 ZNF175 signal didn't replicate in PMBB v2.

**Key findings:**
- **79% of Joe's original signal-driving carriers** (26 of 33) were eliminated from the v2 burden test
- **39%** lost via MAF shift on `chr19:51581437` (gnomAD MAF=3.4% in EAS, excluded by global MAF<0.001)
- **39%** lost via variant-calling drift (5 of Joe's 10 pLOFs missing from v2 annotation)
- Confirms Molly's MAF-shift hypothesis (H1)

**Sub-phase status:**

| Sub-phase | Output | Status |
|---|---|---|
| **P2.1** — ZNF175 variant inventory in PMBB v2 | [`phase2_signal_loss_diagnostic.md`](../results/chapter2_znf175_analysis/phase2_signal_loss_diagnostic.md) + 4 TSVs in `analysis/daniel/outputs/phase8_signal_diagnostic/` | ✅ done |
| **P2.2** — Identify the "115 occurrence" variant | Identified as `chr19:51581437` | ✅ done |
| **P2.3** — Raw-p burden tests (pLOF / missense / combined) | Spec exists; not yet executed | ⏸️ pending — see Phase 4 (subsumed into the extended biobin variant-set sweep) |
| **P2.4** — Diagnostic synthesis | Written into `phase2_signal_loss_diagnostic.md` | ✅ done |

**Why P2.3 is parked:** the variant inventory and "115" identification gave us enough to validate Molly's hypothesis qualitatively. Phase 4 will run biobin with multiple variant sets, multiple MAF thresholds, and per-ancestry subsets — that subsumes the original P2.3 design and answers Molly's hypotheses (H1, H2, H3) more comprehensively. **Decision: roll P2.3 spec into Phase 4 inputs.**

---

### Chapter 2 · Phase 3 — Joe Park pipeline reproduction ⏸️ BLOCKED on Joe outreach

**Goal:** if Phase 4 extended biobin tests still don't surface the v1 signal, reproduce Joe's exact PMBB v1 analysis on PMBB v2 to verify whether the v1 signal was biology, methodology, or an artifact specific to that cohort.

**Blocker:** need Joe's pipeline code + variant list. Daniel confirmed he no longer has access; Joe is the only path.

**Status:** outreach email drafted at [`docs/communications/joe_park_outreach_email.md`](communications/joe_park_outreach_email.md). Pending Nikki to send (per kickoff transcript designation).

**When materials arrive:**
- Stand up Joe's exact pipeline locally
- Run on PMBB v2 cohort — verify if signal re-emerges or stays lost
- Output: `results/chapter2_znf175_analysis/phase3_joe_pipeline_reproduction.md`

---

### Chapter 2 · Phase 4 — Extended biobin burden tests 📋 NEXT

**Goal:** systematically test the hypotheses from Molly's kickoff (H1-H5) by running biobin burden tests under different variant-set definitions and MAF strategies. Same pipeline as Chapter 1 Phase 5/6, varying inputs.

**Concrete actions (in priority order):**

1. **Variant-set definitions** (covers Molly's H2: missense vs pLOF):
   - [ ] pLOF only (Joe's 10 variants × current v2 annotation overlap)
   - [ ] pLOF only **including chr19:51581437** (force-include the MAF-excluded driver — directly tests Phase 2 finding)
   - [ ] Missense only (REVEL>0.6, REVEL>0.5 — two sub-tests)
   - [ ] pLOF + missense combined
   - [ ] All ZNF175 variants in annotation (no MAF filter — sanity check)

2. **MAF strategies** (covers Molly's H1: allele frequency shift):
   - [ ] **Multiple global MAF thresholds** — biobin run at MAF<0.001 (paper), <0.005, <0.01 — sweep to characterize threshold dependence
   - [ ] **Per-ancestry burden tests** — run biobin separately on EUR-only, AFR-only subsets. Variants like chr19:51581437 (rare in EUR, common in EAS) get tested within ancestries where they're still rare. Daniel preserved ancestry-stratified covariates files (`EUR.txt.gz`, `AFR.txt.gz`) — reusable.

3. **Statistical reporting** (covers Molly's H3: multiple-test correction):
   - [ ] Report raw p-values (no FDR correction) — Molly's kickoff ask 11:47: *"it's possible that the multiple test correction was very different... now we're not doing the EWAS"*
   - [ ] Also report BH-FDR as reference

4. **Scope:**
   - [ ] Run on ZNF175 specifically (primary goal)
   - [ ] Sanity check: confirm Chapter 1 Phase 5 results still reproduce on this configuration
   - [ ] Build comparison table: variant set × MAF strategy × phenotype (binary HL vs degree-HL)

**Output:** `results/chapter2_znf175_analysis/phase4_extended_burden.md`

**Estimated effort:** 1 week. Pipeline already validated in Chapter 1 — Phase 4 is mostly running biobin with different inputs and consolidating.

**Methodology note for later discussion with team:** Daniel suggested in his email (2026-05-29) that newer rare-variant approaches could be worth exploring as a follow-up after Phase 4 (preserved verbatim in [`daniel_response_received.md`](communications/daniel_response_received.md)). Whether/when to scope that depends on Phase 4 results and team alignment.

---

### Chapter 2 · Phase 5 — External replication 🔭 LATER

**Goal:** if signal recovered in P2/P4 (or via P3 Joe pipeline), test in external cohorts. Doug noted (transcript 11:38) UKBB phenotype is self-reported (weaker than PMBB audiograms) but provides much larger N.

- [ ] **PMBB v3** (~50k+ exomes): re-run winning Phase 2/3/4 configuration. Daniel preserved relevant scripts; pipeline already walkthrough-documented.
- [ ] **UKBB:** out of scope (no access — per Andre 2026-05). Document for handoff to Nikki if she has access.
- [ ] **All of Us:** Nikki has access. Coordinate scope and divide work.

**Output:** `results/chapter2_znf175_analysis/phase5_external_replication.md`

---

### Chapter 2 · Phase 6 — Doug's 4K audiogram cohort 🔭 LATER

**Doug mentioned (transcript 11:43-11:44):** ~4,000 PMBB individuals have BOTH audiograms AND exomes, with potentially more available via the ENT team. This is a higher-quality subset than the hybrid cohort. Could be the basis for a quantitative-trait analysis specifically for the digenic hypothesis (carrier × second-hit, weighted by audiogram severity).

- [ ] Coordinate with Doug / ENT team on getting the linked audiogram-exome list
- [ ] Re-run Phase 7-equivalent (degree-HL linear) restricted to this 4K subset — does ZNF175 signal emerge in a cleaner subset?
- [ ] Re-run Chapter 2 Phase 1 (second-hit) within the 4K subset for higher-quality phenotyping

---

## Infrastructure & documentation side-quests

Work that doesn't fit a specific scientific phase but is real and tracked.

### LOKI 2022 restore ⏸️ IN PROGRESS (Scott Dudek)

**Context:** Daniel ran biobin in 2021 with `loki-20220926.db`; we ran with `loki-20230816.db`. The 2022 file was deleted in a Ritchie Lab data cleanup. Symlink at `/project/ritchie/datasets/loki/loki.db -> loki-20220926.db` confirms the file existed.

**Status:** Scott Dudek confirmed he'll submit a restore ticket (Slack DM 2026-05-29).

**Impact:** if recovered, would enable byte-equivalent Phase 6 replication. If not, current ~30-50% drift on Phase 6 remains the documented caveat; Phase 5 and Phase 7 unaffected.

**Plan B if restore fails:** documented in current Slack draft to Molly + Nikki — filter our 2023 LOKI removing new LOC/LINC pseudogenes, re-run Phase 6 with "trimmed" version. Won't reach byte-equivalent (boundary refinements remain), but should reduce drift.

### Documentation deliverables ✅ ongoing

- [x] `welcome.md` — onboarding doc (Elena Ch1 replication + continuation context)
- [x] `REPRODUCTION_GUIDE.md` — step-by-step per phase
- [x] `scripts/setup_env.sh` — automated environment validation
- [x] `chapter1_summary.md` + `chapter1_authoritative_pvalues.md` — single source of truth for Ch1
- [x] Drift trilogy in `docs/andre/` — plain-language Phase 5/6/7 drift explanations
- [x] All-genes p-values CSV in `analysis/daniel/outputs/phase8/comparison/`
- [ ] Confluence writeup (deferred — once Phase 4 has results, write canonical project summary for Hall Lab Confluence: https://halllab.atlassian.net/wiki/spaces/HLP/)

---

## External dependencies / communications

| Person | Status | Outstanding ask |
|---|---|---|
| **Daniel Hui** | ✅ Responded 2026-05-29 — see [`daniel_response_received.md`](communications/daniel_response_received.md) | Reply draft ready at [`daniel_reply_2026-05-29.md`](communications/daniel_reply_2026-05-29.md) — Andre to send |
| **Joe Park** | ⏸️ Outreach drafted at [`joe_park_outreach_email.md`](communications/joe_park_outreach_email.md) | Pending Nikki to send (per kickoff designation) |
| **Scott Dudek** (Ritchie Lab) | ⏸️ Ticket open for LOKI 2022 restore | Wait for restore confirmation |
| **Nikki Palmiero** | 📋 Active sync | Joe outreach + Phase 4 variant-set + MAF strategy review |
| **Molly Hall** | 📋 Active sync | Slack message drafted with Phase 2 findings + Plan B for LOKI |
| **Doug Epstein** | 🔭 Loop in after Phase 4 | Phase 6 ENT linkage when ready |
| **Elena** (interim, replicating Ch1) | 📋 Active | Welcome.md + REPRODUCTION_GUIDE ready; ping if blocked |

---

## Open questions (truly open — others moved to relevant phase)

1. **Was Joe Park's MAF cutoff in v1 exactly <0.001?** If yes, the MAF shift on `chr19:51581437` fully explains signal loss. If <0.01 or some other value, the picture is messier. → For Joe.
2. **Why are 5 of Joe's pLOFs no longer in PMBB v2 annotation?** Multi-allelic re-calling, different caller, different filter thresholds — we don't know the exact cause. → For Joe or Daniel.
3. **What's the right ancestry stratification for MAF filtering?** gnomAD's predefined subpops, PMBB's internal cluster labels, or something else? → For Nikki + Andre (decide before Phase 4 execution).
4. **Whether to extend to newer rare-variant methods after Phase 4** — depends on Phase 4 results and team alignment. Scope separately with Nikki when relevant.
5. **Doug's 4K audiogram-exome cohort** — when is the right time to engage Doug + ENT for linkage? → After Phase 4 results, suggest at Molly sync.
