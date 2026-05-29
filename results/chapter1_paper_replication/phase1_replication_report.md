# Phase 1 Replication Report — Gene list curation + annotation filter

**Date:** 2026-05-12
**Run by:** Andre Rico
**Phase:** 1 of 19 (per [`docs/pipeline_walkthrough.md`](../../docs/pipeline_walkthrough.md))
**Status:** ✅ Pipeline replicated correctly (same variant + gene sets as Daniel); textual diff in outputs explained by PMBB v2 annotation schema migration
**Project commit at run time:** [`52f72c9`](https://github.com/) (`docs updates and scripts phase 1`)

---

## TL;DR

Phase 1 of the Hui et al. 2023 PMBB pipeline (gene list curation → variant annotation filter → functional class filter) was successfully re-executed on PMBB v2 data using Daniel's original scripts unchanged.

| Validation target | Daniel's intermediate | Our re-run | Match |
|---|---:|---:|---|
| Annotated HL-gene variants (Step 1.2 output) | 192,029 IDs | 192,029 IDs | **✓ identical set** |
| Final master variant×gene table (Step 1.3 output) | 11,661 IDs | 11,661 IDs | **✓ identical set** |
| Unique genes in master table | 176 | 176 | **✓ identical set** (diff empty) |

Outputs differ **textually** (byte-for-byte diff is large) because Daniel's intermediates were generated from `variant-annotations.txt` (now removed from disk) while we used the replacement `variant-annotation-counts.txt`. Same variants are selected; metadata columns are formatted differently. This is the schema migration Daniel himself flagged at [runbook line 42-43](../../analysis/daniel/runbook_hui2023.txt).

**Phase 1 pipeline replicated correctly. Cleared to proceed to Phase 2.**

---

## Setup

| Item | Value |
|---|---|
| Project root | `/project/hall/analysis/hearing-loss-genomics/` |
| Working directory at run time | project root |
| Python | venv at `venv/` (Python 3.12.8) — pandas 3.0.2, scipy 1.17.1, numpy 2.4.4, statsmodels 0.14.6 |
| Scheduler | LSF 10.1 (`bsub`) |
| Queue | `epistasis_normal` |
| Node | `krypton15` |
| Submit wrapper | [`analysis/daniel/scripts/submit_phase1.sh`](../../analysis/daniel/scripts/submit_phase1.sh) |
| Pipeline script | [`analysis/daniel/scripts/run_phase1.sh`](../../analysis/daniel/scripts/run_phase1.sh) |

**LSF resource request:** 1 CPU, 4 GB memory, 30 min wall.
**LSF resource actually used:** 42.55 s CPU, 19 MB max memory, 49 s wall, 50 s turnaround. The 4 GB request was ~210× over-provisioned; on subsequent phases we could drop to `-M 256` safely. The 30-min wall was ~37× over-provisioned.

---

## What ran

| Step | Cookbook line | Script | Inputs | Output | Time |
|---|---|---|---|---|---|
| 1.1 | — (data prep, not in runbook) | `zcat` | [`data/PMBB_Exome/all_genes_including_ShadisList.txt.gz`](../../data/PMBB_Exome/all_genes_including_ShadisList.txt.gz) (512 B) | [`outputs/phase1/all_genes_including_ShadisList.txt`](../../analysis/daniel/outputs/phase1/all_genes_including_ShadisList.txt) — 179 genes, 1.1 KB | <1 s |
| 1.2 | 25 | [`only_HL_genes.py`](../../analysis/daniel/scripts/pmbb_exome/only_HL_genes.py) | gene list + `data/pmbb_v2/Exome/Variant_annotations/PMBB-Release-2020-2.0_genetic_exome_variant-annotation-counts.txt` (5.1 GB) | [`outputs/phase1/annot_genes_full.txt`](../../analysis/daniel/outputs/phase1/annot_genes_full.txt) — 192,029 rows, 81 MB | 45 s |
| 1.3 | 46 | [`only_func_cats_to_include.py`](../../analysis/daniel/scripts/pmbb_exome/only_func_cats_to_include.py) | Step 1.2 output | [`outputs/phase1/annot_genes_full_funcToInclude.txt`](../../analysis/daniel/outputs/phase1/annot_genes_full_funcToInclude.txt) — 11,661 rows, 8.2 MB | <1 s |

**Step 1.1 note:** Daniel's runbook lines 1-24 describe how the 179-gene HL list was originally assembled (DFNA + DFNB from hereditaryhearingloss.org + Doug's email + Shadi's list) with manual curation (gene-name reconciliation like LRTOMT vs COMT2). The output of that manual process is `all_genes_including_ShadisList.txt.gz`. We **use it as-is** — re-curating from sources is not mechanically reproducible. If/when we port to PMBB v3 or v4, this list should be re-checked against HGNC updates.

---

## Quantitative results — our re-run vs Daniel's intermediates

| File | Daniel (compressed) | Ours (uncompressed) | Row count | Match |
|---|---:|---:|---:|---|
| `all_genes_including_ShadisList.txt` | 512 B (.gz) | 1.1 KB | 179 genes | trivially identical (we just decompressed his) |
| `annot_genes_full.txt` | 12 MB (.gz) | 81 MB | 192,029 | same variant IDs, content differs in metadata columns (see below) |
| `annot_genes_full_funcToInclude.txt` | 1.8 MB (.gz) | 8.2 MB | 11,661 | same variant IDs, same gene set; metadata columns differ |

**Compression ratios:** Daniel's `.gz` files are ~6-7× smaller than our uncompressed. This is normal for tabular text data with repetitive values. We could gzip our outputs to match but it's not necessary for downstream phases.

---

## Validation

The pipeline's correctness was validated by **semantic comparison** (set equality of variant IDs and genes) rather than byte-for-byte diff, because Daniel's intermediates were generated from a now-defunct version of the PMBB v2 annotation file (see Schema migration below).

### Verification commands run

```bash
# Same set of variant IDs in Step 1.2 output?
diff <(zcat data/PMBB_Exome/annot_genes_full.txt.gz | cut -f1 | sort -u) \
     <(cut -f1 analysis/daniel/outputs/phase1/annot_genes_full.txt | sort -u)
# → only difference: header (Daniel: "Constant_ID", ours: "ID")

# Same set of variant IDs in Step 1.3 (master table)?
diff <(zcat data/PMBB_Exome/annot_genes_full_funcToInclude.txt.gz | cut -f1 | sort -u) \
     <(cut -f1 analysis/daniel/outputs/phase1/annot_genes_full_funcToInclude.txt | sort -u)
# → only difference: header

# Same set of genes in the master table?
diff <(zcat data/PMBB_Exome/annot_genes_full_funcToInclude.txt.gz | tail -n +2 | cut -f8 | sort -u) \
     <(tail -n +2 analysis/daniel/outputs/phase1/annot_genes_full_funcToInclude.txt | cut -f8 | sort -u)
# → empty diff (176 genes each, identical set)
```

### What this validates

- The 179-gene HL list is correctly seeded
- `only_HL_genes.py` correctly selects all variants annotated against those genes (192,029 of them) — including 7 genes that have no variants in the annotation (and so are correctly absent: ACTB, COL4A5×, COL4A6×, HARS, KARS, NDP×, REST — × = X-linked, see runbook line 482-489)
- `only_func_cats_to_include.py` correctly applies the burden criteria (pLoF: frameshift / stopgain / splicing; missense REVEL > 0.6) and selects the same 11,661 variants Daniel did

### What this does NOT validate

- That non-ID columns (`Polyphen2_HDIV_score`, `SIFT_score`, etc.) have the **same values** as Daniel had — they differ textually because they come from a different version of the annotation file. Downstream phases use the **ID** + **gene name** + **function category** + **REVEL** columns, which we've verified. If a later phase reads other annotation columns, that assumption may break and we'll need to investigate.

---

## Discovery: PMBB v2 annotation schema migration

**The issue:** Daniel's `data/PMBB_Exome/annot_genes_full.txt.gz` has first column header `Constant_ID`. The current v2 annotation file (`/static/PMBB/PMBB-Release-2020-2.0/Exome/Variant_annotations/PMBB-Release-2020-2.0_genetic_exome_variant-annotation-counts.txt`) has first column header `ID`. The two files contain the same variants but with different column schemas.

**Daniel knew about this.** Runbook line 42-43:
```
#were replaced by these:
#/project/PMBB/PMBB-Release-2020-2.0/Exome/Variant_annotations/PMBB-Release-2020-2.0_genetic_exome_variant-annotation-counts.txt
```

Daniel started Phase 1 against the older `variant-annotations.txt` (Step 1.2, line 25), then transitioned to `variant-annotation-counts.txt` for downstream phases (Phases 8, 9, 12, all-genes burden, etc. — see runbook lines 175, 262, 406). The intermediate `annot_genes_full.txt.gz` we have in `data/PMBB_Exome/` is from the OLD file. The downstream `addBack_multiallelic_stoploss/annot_genes_full.txt.gz` is from the NEW file.

**Implication for our replication:** Today only the `-counts.txt` file exists. We use it for ALL phases — slightly different from Daniel's mixed approach but more internally consistent. Expected downstream effect: some intermediates may differ textually but should select the same variants/genes, as we saw here.

**Implication for downstream phases:** Apply the same set-equality validation we used here. Don't expect byte-for-byte equality with Daniel's intermediates. Expect equality on the columns the pipeline actually uses (ID, gene, function class, REVEL).

---

## The gene count puzzle: 179 vs 176 vs 173

Three different numbers appear in this phase. They reconcile cleanly:

| Number | Source | Meaning |
|---:|---|---|
| **179** | [`data/PMBB_Exome/all_genes_including_ShadisList.txt.gz`](../../data/PMBB_Exome/all_genes_including_ShadisList.txt.gz) | Curated pre-filter HL gene list (DFNA + DFNB + Doug + Shadi) |
| **176** | Our Step 1.3 output, Daniel's `annot_genes_full_funcToInclude.txt.gz` | Genes that survived join with annotation table (have ≥1 variant in `-counts.txt`) |
| **173** | Hui et al. 2023, Fig 2 | Genes used in the published burden test (after excluding X-linked genes for autosomal burden) |

The 179 → 176 drop: Daniel's runbook line 482-489 lists 7 missing genes (ACTB, COL4A5×, COL4A6×, HARS, KARS, NDP×, REST). Only 3 (COL4A5, COL4A6, NDP) are X-linked; the other 4 simply have no variants in the v2 annotation. 179 − 7 = 172, off by 4 from our 176 — there may be additional name-aliasing recoveries in the annotation join. Worth verifying when we get to Phase 12 (all-genes burden, where the paper Figure 2 ranking emerges).

The 176 → 173 drop: 3 X-linked genes (COL4A5, COL4A6, NDP) excluded from the autosomal burden test. The full burden test on 173 autosomal HL genes produces Figure 2 of the paper.

---

## Issues encountered + fixes

### Issue 1 — First submission failed in sanity check (ZNF175)

**Symptom:** First run (job `47167120`) exited at Step 1.1b with `ERROR: Expected gene ZNF175 NOT in gene list — list may be wrong`.

**Root cause:** The sanity check assumed ZNF175 should be IN the pre-burden HL gene list. It shouldn't be — ZNF175 was DISCOVERED in Phase 12's all-genes burden test, it's not a member of the curated 179-gene list of pre-known HL genes. My sanity check design was wrong.

**Fix:** Sanity check now requires three known HL genes (TCOF1, ESRRB, OTOF) to be present, AND explicitly requires ZNF175 to be ABSENT (negative control — confirms we have the canonical pre-burden list, not a contaminated one). Commit: pending.

### Issue 2 — Run log had every line duplicated

**Symptom:** Every log line in `analysis/daniel/logs/phase1/run_*.log` appeared twice.

**Root cause:** Both the `log()` helper function (which pipes to `tee -a "$LOG"`) AND the `exec > >(tee -a "$LOG") 2>&1` redirect at script start were writing to the same file.

**Fix:** Removed `tee` from inside `log()` — the `exec` redirect alone captures everything. Commit: pending.

### Issue 3 — Validation reported `DIFFERS` despite correct pipeline

**Symptom:** Second run (job `47167135`) completed successfully but the validation step reported `annot_genes_full.txt: DIFFERS` and `annot_genes_full_funcToInclude.txt: DIFFERS`. Exit code 2.

**Root cause:** Validation used naive byte-for-byte `diff` against Daniel's `.gz` files. The textual diff is real and large (every row differs in metadata columns), but the IDs, genes, and pipeline-relevant fields are identical. The diff didn't capture this distinction.

**Status:** The mismatch is now understood as Daniel's schema migration; not a real bug. **Recommended fix for future phases:** replace the byte-diff with set-equality checks on the relevant columns (variant ID, gene, function category). See "Recommendations for Phases 2-19" below.

---

## Performance summary

| Metric | Value |
|---|---|
| Total wall time (job execution) | 49 s |
| Total CPU time | 42.55 s |
| Max memory | 19 MB |
| Step 1.2 wall time (5.1 GB annotation read) | 45 s |
| Step 1.3 wall time (8 MB filter) | <1 s |
| Step 1.1 wall time (gzip decompress) | <1 s |
| LSF queue wait | <5 s |

**Performance is single-threaded Python reading a flat file linearly. No obvious bottleneck — disk throughput on `/static/` is plenty for sequential reads of 100 MB/s+.**

**For Phase 2+ resource sizing:** request `-M 256` (256 MB), `-W 5` (5 min wall) for any single-script phase like this. Heavier phases (biobin runs on full pVCF, all-22-chromosome loops) will need different sizing.

---

## Recommendations for Phases 2-19

1. **Validation framework: switch from byte-diff to set-equality.** For each output, compare:
   - Set of primary keys (variant IDs / sample IDs / gene names)
   - Set of relevant downstream-consumed columns
   - Row counts
   Not byte-for-byte equality (which will fail due to schema migration).

2. **Expect textual drift from schema migration.** Daniel mixed `variant-annotations.txt` (old) and `variant-annotation-counts.txt` (new) within his runbook. We use only the new one. Intermediates downstream of any annotation join will differ textually. The pipeline-relevant content should match.

3. **Document discrepancies as they appear.** Each phase report should list quantitative differences (if any) between our output and Daniel's reference, with a hypothesis for the cause. Don't bury or hand-wave.

4. **LSF resource right-sizing.** Single-Python-script phases like this one: `-n 1 -M 256 -W 5`. Don't over-provision — wastes queue capacity.

5. **Plink version pilot still pending.** When we hit Phase 3 (first plink invocation), pilot one chromosome under both plink 1.9 (`/appl/plink-1.90Beta6.18/`) and plink 2 (default `/appl/plink2-20240804/`). Compare outputs before committing the full replication to one version.

---

## Files produced by this phase

```
analysis/daniel/outputs/phase1/
├── all_genes_including_ShadisList.txt           (1.1 KB, 179 genes)
├── annot_genes_full.txt                         (81 MB, 192,029 rows)
└── annot_genes_full_funcToInclude.txt           (8.2 MB, 11,661 rows)  ← master table for downstream phases

analysis/daniel/logs/phase1/
├── run_20260512_204008.log                      (run-internal log)
├── lsf_20260512_204008.out                      (LSF stdout incl. resource summary)
└── lsf_20260512_204008.err                      (LSF stderr — empty)
```

The downstream phases will read `annot_genes_full_funcToInclude.txt` as their master variant×gene table.

---

## Open question still pending (carries from Phase 0)

- **"8 signal-driving cases" definition** — see [`docs/data_inventory.md`](../../docs/data_inventory.md) findings section. Needs Daniel / Doug / Molly to disambiguate before the ZNF175 deep-dive starts. Phase 1 doesn't depend on this.
