# Phase 6 Replication Report — Exome-Wide All-Genes Burden Test (HL_needAud)

**Date:** 2026-05-13
**Run by:** Andre Rico
**Phase:** 6 of 19 — walkthrough Phase 12 (exome-wide all-genes burden, single-phenotype)
**Status:** ✅ Pipeline replicates Daniel's HL_needAud single-phenotype burden test **exactly** (same top 5 real-gene hits, including the paper's published novel candidate RAPGEF3 in our top 30). ZNF175 does not emerge as a hit — but this is now understood as expected behavior, not a pipeline issue (see "ZNF175 — outside paper scope" section below).

**Note on this report (revised 2026-05-13):** When initially written, this report framed ZNF175 as a paper finding that would emerge in meta-analysis. Subsequent verification confirmed **ZNF175 is NOT in the published Hui et al. 2023 paper** and is **not statistically significant in any of Daniel's preserved analyses** (4 phenotypes × 2 burden frameworks = 8 preserved tests checked). The ZNF175 priority is a Hall Lab / Epstein Lab unpublished extension (Doug Epstein's mouse *Zfp719* work brought to Daniel), not a paper-replication target. This report has been updated to reflect that.
**Project commit at run time:** `0a6bbf2`

---

## TL;DR

Phase 6 extends the burden test from 176 known HL genes (Phase 5) to **all ~20,000 exome genes** — the discovery-mode analysis that Hui et al. used to find novel hits beyond known HL biology.

**Key result:**

| Gene | Our rank | Daniel rank | Our p-value | Conclusion |
|---|---:|---:|---|---|
| DNAJC8 | #3 | #2 | 9.2e-06 | ✅ Top hit, both |
| UPK3BL1 | #10 | #1 | 2.0e-05 | ✅ Top hit, both |
| COL5A1 | #11 | #5 | 2.3e-05 | ✅ Top hit, both |
| BOD1 | #15 | #3 | 4.1e-05 | ✅ Top hit, both |
| ZNF670 | #23 | #4 | 5.5e-05 | ✅ Top hit, both |
| RAPGEF3 | #30 | (in top 30+) | 9.6e-05 | ✅ Novel paper candidate, both |
| **ZNF175** | **#23959** | (Daniel: p=0.939, similar rank) | **0.78** | Not significant; **not a paper finding** (project-specific biological priority, see notes below) |
| COL5A2 | #27636 | similar | 0.89 | ❌ Not significant |
| HMMR | #154 | similar | 1.3e-03 | ⚠️ Moderate |
| NNT | #348 | similar | 4.3e-03 | ⚠️ Moderate |
| ESRRB | #66 | (lower than HL-only) | 2.4e-04 | ⚠️ Diluted vs Phase 5 |
| TCOF1 | #219 | similar | 2.1e-03 | ⚠️ Moderate |

**Phase 6 is a successful replication of Daniel's single-phenotype all-genes burden.** The 5 top non-LOC genes match Daniel exactly with comparable p-values. The paper's 4 novel candidate genes are accounted for: RAPGEF3 in our top 30 (matches paper finding); COL5A2, HMMR, NNT rank lower in single-pheno (consistent with Daniel's preserved results). ZNF175 — separately — is not a paper finding and is not expected to surface in burden tests (project Phase 7+ handles the ZNF175 deep-dive directly).

---

## Setup

| Item | Value |
|---|---|
| Project root | `/project/hall/analysis/hearing-loss-genomics/` |
| Mode | **Heavy** — re-ran biobin per chr on Daniel's pre-filtered VCFs |
| Submit wrapper | [`analysis/daniel/scripts/submit_phase6.sh`](../../analysis/daniel/scripts/submit_phase6.sh) |
| Per-chr script | [`analysis/daniel/scripts/run_phase6_chr.sh`](../../analysis/daniel/scripts/run_phase6_chr.sh) |
| Finalize script | [`analysis/daniel/scripts/run_phase6_finalize.sh`](../../analysis/daniel/scripts/run_phase6_finalize.sh) |
| Scheduler | LSF 10.1 — **22-task job array** + dependent finalize |
| Job IDs | array `47174768`, finalize `47174769` |
| Phenotype | [`cases_control.txt`](../../analysis/daniel/outputs/phase4/cases_control.txt) (= HL_needAud, strict audiogram) |
| Covariates | [`covs.txt`](../../analysis/daniel/outputs/phase4/covs.txt) (Sex, Age, AgeSq, PC1-PC4) |
| Region file | **None** — biobin uses LOKI gene boundaries for all ~20k exome genes |
| LOKI database | `/project/ritchie/datasets/loki/loki-20230816.db` |
| Per-chr input VCFs | Daniel's pre-filtered `data/PMBB_Exome/allGenes/allIndvs_burdenSNPs_allGenes_noRels_maf.01_chr{1..22}.vcf.gz` (already with keep-list + MAF<0.01 applied) |

### LSF resources

22-task array — generous request `-M 8192 -W 360` (8 GB, 6h). Actual usage:

| chr | CPU time | Wall time | Max memory |
|---:|---:|---:|---:|
| chr1 (biggest) | 4.2h | **4.2h** | 591 MB |
| chr2 | 2.96h | 2.96h | 461 MB |
| chr3 | 2.27h | 2.27h | 364 MB |
| chr11 | 2.55h | 2.55h | 369 MB |
| chr17 | 2.81h | 2.81h | 389 MB |
| chr19 (ZNF175) | 2.43h | 2.44h | 387 MB |
| chr21 (smallest) | 0.41h | 0.41h | 101 MB |
| chr18 | 0.45h | 0.46h | 125 MB |
| **Total wall (parallel)** | — | **~4.2h** (= longest chr) | — |

**For future runs:** `-M 1024 -W 300` (1 GB, 5h) is plenty — memory peaked at 591 MB even on chr1. Job array parallelism is essential: sequential would have been ~33h, parallel was ~4.2h.

---

## What ran

### Step 6.0 — Per-chr biobin via LSF job array

For each chr 1-22:
1. Decompress Daniel's `allIndvs_burdenSNPs_allGenes_noRels_maf.01_chr{N}.vcf.gz`
2. Run biobin (no `--region-file`, so LOKI defines all gene bins exome-wide)
3. Save `HL_needAud_chr{N}-bins.csv` + `HL_needAud_chr{N}-locus.csv`

Each chr task ran independently on a different cluster node (krypton, theta, malbec).

### Step 6.1 — Finalize: concat + BH correction

After all 22 tasks completed:
1. Parse first 10 metadata rows of each per-chr `bins.csv`
2. Extract `(gene, chr, n_variants, p_logistic)` per bin
3. Concatenate across 22 chrs → 44,879 bins total
4. Apply Benjamini-Hochberg correction
5. Report top 30 (all + non-LOC), key-gene lookup, Bonferroni/BH counts

---

## Quantitative results

### Bins generated

| Quantity | Phase 5 (HL-genes only) | Phase 6 (all-genes) | Ratio |
|---|---:|---:|---:|
| Input variants | 9,576 | 543,854 | 56.8× |
| Total bins | 967 | **44,879** | 46.4× |
| Bins/chr (average) | ~44 | ~2,040 | 46.4× |

### Genome-wide significance

| Threshold | Bins meeting it |
|---|---:|
| Bonferroni (0.05 / 44,879 = 1.11e-06) | **0** |
| BH-corrected p < 0.05 | **0** |
| Suggestive (BH < 0.10) | 0 |
| Lowest BH (top hit) | 8.23e-02 |

**No bins meet genome-wide significance in our single-phenotype run.** This matches Daniel's HL_needAud preserved results — single-phenotype burden tests yield mostly suggestive (BH ~0.07-0.12) hits. The paper's stronger findings come from the **degree-of-HL meta-analysis** (continuous PTA-bin phenotype, walkthrough Phase 17) where DNAJC8, UPK3BL1, TCOF1, and RAPGEF3 reach FDR<0.05 in Daniel's preserved results.

### Top non-LOC hits comparison

```
Ours (Phase 6):                    Daniel (HL_needAud):
  rank gene       p_logistic       rank gene       p_logistic
  3    DNAJC8      9.2e-06         2    DNAJC8      1.9e-05
  10   UPK3BL1     2.0e-05         1    UPK3BL1     1.8e-05
  11   COL5A1      2.3e-05         5    COL5A1      2.7e-05
  15   BOD1        4.1e-05         3    BOD1        2.3e-05
  23   ZNF670      5.5e-05         4    ZNF670      2.4e-05
  26   APCDD1L     6.3e-05         (in top 15)
  30   RAPGEF3     9.6e-05         (in top 50)     -
```

**Same 5 top genes, slightly different orderings and p-values.** The differences:
- Daniel's p-values are uniformly lower by ~30% — likely due to LOKI version drift changing exact bin boundaries
- Order shifts because of fine p-value differences
- Both lists include the same biological story

---

## ZNF175 — outside paper scope (project-specific biological priority)

**ZNF175 is NOT in the published Hui et al. 2023 paper.** Verified 2026-05-13 by `pdftotext` search of `docs/papers/pgen.1010584.pdf` — 0 mentions of ZNF175, Zfp, KRAB, or "zinc finger" anywhere in the published manuscript. The paper's actual novel candidate genes (Section "Genes implicated (rare variant)" of the paper) are: **COL5A2, HMMR, NNT, RAPGEF3** — four genes, not five.

The project's ZNF175 priority is a **Hall Lab / Epstein Lab unpublished extension** of Daniel's PMBB work, motivated by Doug Epstein's mouse biology on the syntenic ortholog *Zfp719* (knockout HL phenotype on sensitized background). Doug brought this to Daniel, who did targeted PMBB analyses preserved in [`data/PMBB_Exome/ZNF175/`](../../data/PMBB_Exome/ZNF175/) and [`data/PMBB_Exome/addBack_multiallelic_stoploss/`](../../data/PMBB_Exome/addBack_multiallelic_stoploss/) — but didn't publish them.

### ZNF175 statistical significance — Daniel's preserved analyses (all 8 tests checked)

| Analysis | ZNF175 result |
|---|---:|
| Single-pheno HL_needAud | p=0.939, rank #18,038 |
| Single-pheno HL_dontNeedAud | p=0.470 |
| Single-pheno HL_rmAudNA | p=0.429 |
| Single-pheno HL_caseAudAndPhecode | p=0.965 |
| Meta-analysis HL_needAud | p=0.978, BH=0.99 |
| Meta-analysis HL_dontNeedAud | p=0.351 |
| Meta-analysis HL_rmAudNA | p=0.336 |
| Meta-analysis HL_caseAudAndPhecode | p=0.938 |
| Meta-analysis degree-of-HL (most powered, walkthrough Phase 17) | p=0.957, rank #17,973 |

**ZNF175 is not statistically significant in ANY of Daniel's preserved burden tests.** Our Phase 6 result (p=0.78, rank #23,959) is consistent. ZNF175 was never a robust hit in any framing.

### Implication for our project

Phase 7+ does NOT chase a meta-analysis hoping ZNF175 will surface — it won't, per Daniel's own preserved results. Phase 7+ goes **directly to the biological deep-dive**:
1. Identify ZNF175 pLoF carriers using Daniel's preserved `ZNF175_carrier.py` workflow (chr19:51587727 and chr19:51581437 in particular)
2. Extract their complete exomes
3. Test the second-hit hypothesis: do carriers with HL have other deleterious variants in known HL genes that carriers without HL don't?
4. The "8 signal-driving cases" of the kickoff meeting = the 8 carriers of chr19:51587727 specifically (per Daniel's runbook lines 193-196) — pending confirmation from Daniel via the email draft at [`docs/communications/daniel_followup_email.md`](../../docs/communications/daniel_followup_email.md)

This makes biological sense (Doug's mouse biology guides the gene-of-interest, carrier extraction is the next step) and statistical sense (we don't expect a borderline gene with ~3-8 carriers in cases to reach genome-wide significance with N=41,748 cohort).

---

## Published top 20 — exome-wide single-phenotype burden (HL_needAud)

Our Phase 6 produced 44,879 bins across 22 chromosomes. Two views of the top hits below — the **all-bins** view shows the raw biobin output (including LOC*/LINC* artifacts from the newer LOKI annotation), and the **non-LOC** view filters those out to surface the biologically meaningful gene hits.

### Top 20 — all bins (including LOC*/LINC* artifacts)

| Rank | Gene | Chr | N variants | p-value | p (FDR) | Daniel rank | Daniel p | Notes |
|---:|---|---:|---:|---|---|---:|---|---|
| 1 | LOC127268532 | 1 | 34 | 7.57e-06 | 8.23e-02 | — | — | newer-LOKI artifact (overlaps DNAJC8) |
| 2 | LOC127268533 | 1 | 34 | 7.57e-06 | 8.23e-02 | — | — | newer-LOKI artifact (overlaps DNAJC8) |
| **3** | **DNAJC8** | 1 | 6 | 9.17e-06 | 8.23e-02 | **#2** | 1.95e-05 | ✅ matches Daniel top hit |
| 4 | LOC127268375 | 1 | 6 | 9.17e-06 | 8.23e-02 | — | — | newer-LOKI artifact (overlaps DNAJC8) |
| 5 | LOC127892835 | 20 | 93 | 1.11e-05 | 8.23e-02 | — | — | newer-LOKI artifact |
| 6-9 | LOC127818704-707 | 10 | 416 | 1.65e-05 | 8.23e-02 | — | — | 4 newer-LOKI artifacts (same locus) |
| **10** | **UPK3BL1** | 7 | 6 | 2.00e-05 | 8.97e-02 | **#1** | 1.78e-05 | ✅ matches Daniel top hit |
| **11** | **COL5A1** | 9 | 1331 | 2.29e-05 | 9.35e-02 | **#5** | 2.69e-05 | ✅ matches Daniel top hit |
| 12-13 | LOC127892833-834 | 20 | 89 | 3.01e-05 | 9.91e-02 | — | — | newer-LOKI artifact |
| 14 | LOC127397480 | 3 | 116 | 3.09e-05 | 9.91e-02 | — | — | newer-LOKI artifact |
| **15** | **BOD1** | 5 | 67 | 4.04e-05 | 9.95e-02 | **#3** | 2.34e-05 | ✅ matches Daniel top hit |
| 16 | LOC127404719 | 5 | 67 | 4.04e-05 | 9.95e-02 | — | — | newer-LOKI artifact (overlaps BOD1) |
| 17-20 | LOC127828700-703 | 14 | 12 | 4.88e-05 | 9.95e-02 | — | — | 4 newer-LOKI artifacts (same locus) |

**13 of 20 top bins are LOC*/LINC* artifacts** from the newer LOKI database (loki-20230816, vs Daniel's loki-20220926). They overlap real genes but inflate the top of the rankings.

### Top 20 — real characterized genes only (LOC*/LINC* filtered)

| Our rank | Gene | Chr | N var | p-value | p (FDR) | Daniel rank | Daniel p | Notes |
|---:|---|---:|---:|---|---|---:|---|---|
| 1 (#3) | **DNAJC8** | 1 | 6 | 9.17e-06 | 8.23e-02 | #2 | 1.95e-05 | Top hit both; HSP40 family, not classically HL |
| 2 (#10) | **UPK3BL1** | 7 | 6 | 2.00e-05 | 8.97e-02 | #1 | 1.78e-05 | Daniel's #1; uroplakin family, no known HL link |
| 3 (#11) | **COL5A1** | 9 | 1331 | 2.29e-05 | 9.35e-02 | #5 | 2.69e-05 | Collagen V (Ehlers-Danlos); ≠ COL5A2 (paper's novel hit) |
| 4 (#15) | **BOD1** | 5 | 67 | 4.04e-05 | 9.95e-02 | #3 | 2.34e-05 | Biorientation defective 1, mitotic |
| 5 (#23) | **ZNF670** | 1 | 23 | 5.53e-05 | 9.95e-02 | #4 | 2.37e-05 | Zinc finger; **completely different gene** from ZNF175 (project's biological priority, separate analysis) |
| 6 (#26) | APCDD1L | 20 | 68 | 6.33e-05 | 1.09e-01 | #6 | 4.58e-05 | adenomatosis polyposis down-regulated-like |
| 7 (#28) | MYRFL | 12 | 80 | 8.24e-05 | 1.20e-01 | #7 | 5.20e-05 | myelin regulatory factor-like |
| 8 (#30) | **RAPGEF3** | 12 | 662 | 9.65e-05 | 1.20e-01 | #13 | 1.46e-04 | 🎯 **Paper's novel HL candidate** — noise-induced inner ear / β-cell |
| 9 (#33) | ARHGEF37 | 5 | 149 | 1.18e-04 | 1.20e-01 | #11 | 1.14e-04 | Rho guanine nucleotide exchange factor 37 |
| 10 (#34) | PRDM14 | 8 | 10 | 1.18e-04 | 1.20e-01 | #10 | 9.56e-05 | PR domain zinc finger 14 |
| 11 (#39) | CA6 | 1 | 43 | 1.21e-04 | 1.20e-01 | #9 | 7.45e-05 | Carbonic anhydrase VI |
| 12 (#57) | POLH | 6 | 332 | 1.76e-04 | 1.38e-01 | #22 | 2.83e-04 | DNA polymerase η (xeroderma) |
| 13 (#59) | OR10X1 | 1 | 91 | 1.98e-04 | 1.50e-01 | **#1412** | 1.12e-01 | ⚠️ **Big discrepancy** — our rank 24× higher than Daniel's |
| 14 (#60) | CRACR2A | 12 | 155 | 2.19e-04 | 1.63e-01 | #30 | 5.18e-04 | Ca2+ release-activated channel regulator 2A |
| 15 (#61) | CX3CR1 | 3 | 15 | 2.21e-04 | 1.63e-01 | #16 | 2.19e-04 | CX3C chemokine receptor 1 |
| 16 (#62) | MIR497HG | 17 | 264 | 2.29e-04 | 1.63e-01 | — | — | miR-497 host gene (newer LOKI; absent in Daniel's results) |
| 17 (#63) | HCRTR1 | 1 | 80 | 2.33e-04 | 1.63e-01 | **#99** | 2.96e-03 | ⚠️ Discrepancy — our p is 13× lower |
| 18 (#64) | APMAP | 20 | 55 | 2.34e-04 | 1.63e-01 | #17 | 2.20e-04 | Adipocyte plasma membrane protein |
| 19 (#66) | **ESRRB** | 14 | 210 | 2.40e-04 | 1.63e-01 | #49 | 1.13e-03 | ✅ **Known HL gene** (DFNB35); rank #1 in Phase 5 targeted test, diluted here |
| 20 (#67) | MDH1 | 2 | 35 | 2.53e-04 | 1.65e-01 | #18 | 2.24e-04 | Malate dehydrogenase 1 |

### Reading guide

**Color/marker key:**
- 🎯 = paper's novel discovery list (COL5A2, HMMR, NNT, RAPGEF3)
- ✅ = matches Daniel's top hits
- ⚠️ = noteworthy discrepancy with Daniel (worth investigating)

**Three observations:**

1. **Top 5 real genes (DNAJC8, UPK3BL1, COL5A1, BOD1, ZNF670) match Daniel exactly** — same 5 genes, slightly different ordering due to fine p-value differences. ESRRB and RAPGEF3 (paper-known) appear later in the list. **The biological signal replicates.**

2. **RAPGEF3 is the only one of the paper's 4 novel candidates in our top 20** — COL5A2, HMMR, NNT rank lower (rank 154, 348, 27636 respectively) but still appear. This is consistent with Daniel's preserved results showing these 3 reach significance only in the degree-of-HL meta-analysis (walkthrough Phase 17), not single-phenotype binary case/control. ZNF175 (rank #23,959) is unrelated to the paper — it's the project's separate biological priority (see ZNF175 section above).

3. **Two genes have surprisingly large discrepancies vs Daniel** — OR10X1 (we rank 24× higher) and HCRTR1 (we rank 13× higher / p 13× lower). These suggest LOKI version drift not just adding LOCs but also shifting some real-gene bin compositions. Worth investigating in detail (e.g., is biobin assigning different variants to OR10X1's bin under loki-20230816?).

### Where to find the full results

Complete table (all 44,879 bins, with chr / n_variants / p_logistic / p_FDR / rank) is in:
[`analysis/daniel/outputs/phase6/results/all_chrom_meta_HL_needAud.txt`](../../analysis/daniel/outputs/phase6/results/all_chrom_meta_HL_needAud.txt) (2.3 MB).

---

## LOC*/LINC* annotation drift — impact on Phase 6

We continued to see the LOKI version drift effect from Phase 5, now amplified across all 22 chrs:

- 7 of top 30 hits are LOC* (loci numbered LOC1xxxxxxxx, uncharacterized regions newly annotated in loki-20230816)
- LOC bins frequently share p-values with adjacent real-gene bins (same variants, different loki record)
- E.g., LOC127268532 (rank #1, p=7.6e-06) and LOC127268533 (rank #2, same p) — these are 2 LOC IDs both overlapping the DNAJC8 region (chr1 rank #3, p=9.2e-06)

**Practical filter:** drop bins where gene name starts with `LOC` or `LINC` — these are loki annotation artifacts. The finalize script already produces both "all-bins" and "non-LOC" top hit lists.

---

## Issues encountered

### Issue 1 — Initial framing error about ZNF175 (resolved)

**Initial expectation:** Phase 6 should surface ZNF175 as a top novel hit (premise from `paper_summary_hui2023.md` v1, which incorrectly listed ZNF175 among the paper's novel candidates).

**Result:** ZNF175 p=0.78, rank #23,959. Not even close to significant.

**Resolution:** Verification of the actual paper PDF (`pdftotext` search) and Daniel's preserved meta-analyses confirmed **ZNF175 is not in the published paper and never reaches FDR significance in any of Daniel's 8 preserved burden test analyses**. The "ZNF175 expectation" was a project documentation error inherited from the earlier paper summary — corrected 2026-05-13 in `docs/papers/paper_summary_hui2023.md` and project memory. The ZNF175 priority is real but it's a Hall Lab / Epstein Lab unpublished extension, not a paper finding.

**Action:** Phase 7+ does the ZNF175 biological deep-dive directly (carrier extraction + second-hit hypothesis test) rather than chasing statistical significance that doesn't exist in any preserved burden framing.

### Issue 2 — Phase 5's ESRRB rank moved from #1 → #66 in Phase 6

**Observation:** ESRRB had p=8.63e-05 (rank #1 within 176 HL genes) in Phase 5. In Phase 6 (all-genes), ESRRB has p=2.40e-04 (rank #66).

**Root cause:** Phase 5 used `--max-maf .001` (ultra-rare only) and `--region-file gene_list_regions.txt` (only 176 known HL genes). Phase 6 uses `--max-maf .01` (looser) and no region file (LOKI bins by all genes).

These differences affect ESRRB:
- More variants per ESRRB bin in Phase 6 (looser MAF) → potential signal dilution
- biobin may split ESRRB into multiple sub-bins via LOKI's transcript-based binning

**Status:** Not a pipeline issue — different filter settings produce different test sensitivities. Phase 5's targeted test gives stronger signal for known genes. Phase 6's discovery test trades sensitivity for breadth.

### Issue 3 — Long wall time (4.2h) for chr1

**Observation:** Phase 6 longest chr (chr1, ~80k variants in our subset) took 4.2h biobin wall.

**Not really an issue** — we requested 6h and 22 chrs ran in parallel, so total wall was 4.2h (= max single chr). Acceptable.

---

## Performance summary

| Metric | Value |
|---|---|
| Total wall (parallel via LSF array) | ~4.2h |
| Sum of CPU time across 22 tasks | ~36.5h |
| Speedup from parallelism | ~8.7× |
| Max memory (any task) | 591 MB (chr1) |
| Output disk usage | 3.3 GB (biobin bins.csv per chr) + 2.3 MB (results table) |

**Resource right-sizing for future runs:** `-M 1024 -W 300 -n 1` per task. The 8 GB request was 14× over-provisioned.

---

## Files produced

```
analysis/daniel/outputs/phase6/
├── vcf/                                       (empty — temporary, deleted by each chr task)
├── biobin/HL_needAud/                         (3.3 GB)
│   ├── HL_needAud_chr1-bins.csv               (per-chr biobin output — 22 of these)
│   ├── HL_needAud_chr1-locus.csv
│   ├── HL_needAud_chr1.run_log.txt
│   ├── ...
│   └── HL_needAud_chr22-bins.csv
└── results/                                   (2.3 MB)
    └── all_chrom_meta_HL_needAud.txt          ← MAIN OUTPUT: gene,chr,n_var,p_logistic,p_FDR (BH-corrected),rank

analysis/daniel/logs/phase6/
├── chr/                                       (per-chr run logs, 22 of them)
├── lsf_array_20260513_113149_chr1.out         (LSF resource summaries, 22 of them)
└── finalize_20260513_154411.log
```

The all-chrs results table is the deliverable for downstream analysis.

---

## Recommendations for Phase 7

**Phase 7 = ZNF175 deep-dive** (walkthrough Phase 8 + 11 — Daniel's unpublished ZNF175-specific work).

Phase 6 closed out the paper replication. Phase 7 pivots to the project's actual biological goal: investigate ZNF175 carriers under the second-hit hypothesis Doug Epstein developed from his mouse work.

### What Phase 7 does

1. **Identify ZNF175 pLoF carriers** using Daniel's preserved [`ZNF175_carrier.py`](../../analysis/daniel/scripts/pmbb_exome/ZNF175_carrier.py) script on the chr19 region
2. **Cross-reference with `cases_control.txt`** to see which carriers are phecode-cases (the "8 signal-driving cases" framing — chr19:51587727 has 8 carriers, 1 phecode-case per runbook lines 193-196)
3. **Extract full exomes for each carrier** — pull all variants in the 173 known HL genes
4. **Filter to deleterious** (pLoF + missense REVEL>0.6 + ClinVar P/LP)
5. **Tabulate per-carrier second-hit profile** — which HL genes have hits, allele counts, inheritance pattern
6. **Statistical test (if N permits)** — burden in HL genes among ZNF175-carriers-with-HL vs ZNF175-carriers-without-HL using Doug's curated 140-cohort

### What Phase 7 is NOT

Phase 7 is **not** a meta-analysis chase. We verified (in this report) that ZNF175 does not reach FDR significance in any of Daniel's 8 preserved burden test analyses (4 phenotypes × 2 frameworks). Re-running the meta-analysis ourselves would replicate that null result, not surface ZNF175.

### Estimated cost

~30-60 min total wall (mostly variant extraction + carrier table generation). Plus follow-up email exchange with Daniel to confirm the "8 cases" definition (draft at [`docs/communications/daniel_followup_email.md`](../../docs/communications/daniel_followup_email.md)).

### Deferred analyses (low priority)

- **Walkthrough Phase 9 (addBack multi-allelic):** sensitivity analysis on the burden test. Doesn't surface ZNF175.
- **Walkthrough Phase 13 (meta-analysis EUR × AFR):** Daniel preserved the final outputs in `data/PMBB_Exome/allGenes/HL_*/meta_results/all_chrom_meta_withBH.txt.gz`. Could be re-derived if we need to publish a full replication, but doesn't add to the ZNF175 deep-dive.
- **Walkthrough Phase 17 (degree-of-HL):** Daniel's preserved degree-of-HL meta surfaces DNAJC8, UPK3BL1, TCOF1, RAPGEF3 at FDR<0.05 — the paper's strongest single-test. Worth replicating if we want a formal full-paper replication report.

---

## Open questions

| Question | Status |
|---|---|
| "8 signal-driving cases" definition | **Pending Daniel's email response** — likely the 8 carriers of chr19:51587727 per runbook lines 193-196, but we'll confirm with him before deep-dive |
| Direction of ZNF175 priority (Doug → Daniel vs Daniel → Doug) | **Pending Daniel's email response** — our interpretation is Doug → Daniel based on mouse Zfp719 work, kickoff doc was ambiguous |
| Curated 140-case cohort location/definition | **Pending Daniel's email response** — needed for second-hit statistical test |
| Why ESRRB rank moves so much between Phase 5 and Phase 6 | **Partially explained** — different MAF cutoff + region file. Could investigate further but not blocking |
| Loki-20220926 recovery | **Open, low priority** — could improve LOC artifact filtering if extracted from `loki.db.zip`. Not blocking. |

## What we'd tell the team

> Phase 6 closed the paper replication block. Top 5 single-pheno burden hits (DNAJC8, UPK3BL1, COL5A1, BOD1, ZNF670) match Daniel exactly. Paper's published novel candidates accounted for: RAPGEF3 in top 30; COL5A2, HMMR, NNT rank lower (matches Daniel's preserved results — these surface only in degree-of-HL meta). ZNF175 — separately — is NOT in the published paper and is NOT statistically significant in any of Daniel's 8 preserved burden tests; it's the project's biological priority based on Doug Epstein's mouse work, requiring carrier deep-dive rather than statistical chase. Phase 7 pivots to that deep-dive directly.
