# Reproduction Guide — Chapter 1 (Hui et al. 2023 Replication)

**Audience:** Elena and anyone else replicating the published paper-replication phases.
**Scope:** Chapter 1 only. Chapter 2 (ZNF175 follow-up) is not validated and is excluded from this guide.
**Prereq:** completed `bash scripts/setup_env.sh` with no failures.

This guide tells you, for each Chapter 1 phase, what the inputs are, what command to run, how long it takes, and how to validate. Detailed methodology lives in each phase's report in `results/chapter1_paper_replication/`.

---

## Quick reference: Chapter 1 scripts

| Phase | Script to run | Output dir | Estimated time |
|---|---|---|---:|
| Ch1 P1 — gene list | `run_phase1.sh` | `analysis/daniel/outputs/phase1/` | < 1 min |
| Ch1 P2 — SNP IDs | `run_phase2.sh` | `analysis/daniel/outputs/phase2/` | ~20 s (light) / 15-20 min (heavy) |
| Ch1 P3 — plink extraction | `run_phase3.sh` | `analysis/daniel/outputs/phase3/` | ~2 s (light) / ~5 min/chr (heavy) |
| Ch1 P4 — prep files | `run_phase4.sh` | `analysis/daniel/outputs/phase4/` | < 1 min |
| Ch1 P5 — HL burden test | `run_phase5.sh` | `analysis/daniel/outputs/phase5/` | ~46 min (LSF) |
| Ch1 P6 — exome-wide burden | `submit_phase6.sh` then `run_phase6_finalize.sh` | `analysis/daniel/outputs/phase6/` | ~2 h (LSF, 22-task array) |
| **Ch1 P7 — degree-HL burden (light)** | **`run_phase8.sh`** ← phase8 = Ch1 P7 | `analysis/daniel/outputs/phase8/` | < 1 s |

> ⚠️ **Do NOT run `run_phase7.sh`** — that's a Chapter 2 script (ZNF175 deep-dive), not Chapter 1. Chapter 2 is in-progress and not validated. Use `run_phase8.sh` for Chapter 1 Phase 7.

Each script has a `submit_*.sh` partner that wraps it for LSF job submission. Use `submit_*.sh` for heavy phases (5, 6); plain `bash run_*.sh` for light phases.

---

## Reproducing Chapter 1 — paper replication

### Ch1 P1 — Gene list curation

**Goal:** produce `annot_genes_full_funcToInclude.txt` (master variant×gene table — input to every downstream phase).

**Run:**
```bash
source venv/bin/activate
bash analysis/daniel/scripts/run_phase1.sh
```

**Outputs to check:**
- `analysis/daniel/outputs/phase1/annot_genes_full_funcToInclude.txt` — should have 11,661 rows
- `analysis/daniel/logs/phase1/run_*.log` — should end with "Phase 1 complete"

**Validate:** compare to Daniel's preserved version
```bash
diff <(sort analysis/daniel/outputs/phase1/annot_genes_full_funcToInclude.txt) \
     <(zcat data/PMBB_Exome/annot_genes_full_funcToInclude.txt.gz | sort) | wc -l
```
Expected: a small number (header drift OK; not zero is fine — see Phase 1 report for known schema differences).

**Light mode:** there is no light mode — this is the gene list construction, must run from scratch. Fast though (< 1 min).

**Detailed methodology:** [`results/chapter1_paper_replication/phase1_replication_report.md`](results/chapter1_paper_replication/phase1_replication_report.md)

---

### Ch1 P2 — SNP ID reconciliation

**Goal:** produce `matched_snp_IDs_annot_pVCF_noNA_noMultiallelic.extract` (9,667 SNP IDs for plink `--extract`).

**Run (light mode — recommended):**
```bash
bash analysis/daniel/scripts/run_phase2.sh
```
Uses Daniel's preserved `vcf_SNP_IDs_allchr.txt.gz`. Skips the 22-task LSF extraction from raw pVCFs.

**Run (heavy mode):**
```bash
bash analysis/daniel/scripts/submit_phase2.sh   # submits LSF job array
```

**Validate:**
```bash
md5sum analysis/daniel/outputs/phase2/matched_snp_IDs_annot_pVCF_noNA_noMultiallelic.extract
# Expected: 5e80ebc0faa5e68277cfeb948af8b1da (when sorted unique)
```

**Detailed methodology:** [`results/chapter1_paper_replication/phase2_replication_report.md`](results/chapter1_paper_replication/phase2_replication_report.md)

---

### Ch1 P3 — plink genotype extraction

**Goal:** per-chr `allIndvs_chr{1..22}.{bed,bim,fam}` from Phase 2's `.extract` file.

**Run (light mode):**
```bash
bash analysis/daniel/scripts/run_phase3.sh
```
Decompresses Daniel's preserved files. No plink invocation.

**Run (heavy pilot on chr21):**
```bash
PHASE3_MODE=heavy CHR=21 bash analysis/daniel/scripts/run_phase3.sh
```

**Validate (heavy chr21 only):**
```bash
md5sum analysis/daniel/outputs/phase3/allIndvs_chr21.bed
# Expected: 9a0ccc62fc3fdb10cd48d9e9061ad7a6 (matches Daniel byte-for-byte)
```

**Heavy full (all 22 chrs):** requires LSF — see Phase 3 report.

**Detailed methodology:** [`results/chapter1_paper_replication/phase3_replication_report.md`](results/chapter1_paper_replication/phase3_replication_report.md)

---

### Ch1 P4 — Preparatory files

**Goal:** produce `cases_control.txt`, `covs.txt`, `gene_list_regions.txt` (inputs to Phase 5's biobin).

**Run:**
```bash
bash analysis/daniel/scripts/run_phase4.sh
```

**Validate:** semantic set-equality with Daniel
```bash
python3 analysis/daniel/scripts/validate_phase4.py
```

**Known issue (Python 2 vs 3 float repr):** `AgeSq` column will differ textually from Daniel's but be numerically identical (round to 4 decimals to compare).

**Detailed methodology:** [`results/chapter1_paper_replication/phase4_replication_report.md`](results/chapter1_paper_replication/phase4_replication_report.md)

---

### Ch1 P5 — First burden test (HL genes only, binary logistic)

**Goal:** biobin on 173 known HL genes → `merged_maf.001_noRels_keepHLcases-bins.csv`. Top hit should be ESRRB.

**Run (LSF required):**
```bash
bash analysis/daniel/scripts/submit_phase5.sh
```
Reserves 1 CPU, 1 GB RAM, 60 min wall (very over-provisioned — actual usage ~177 MB). Watch with `bjobs`.

**Run (interactive, no LSF):**
```bash
bash analysis/daniel/scripts/run_phase5.sh
```
Same script — biobin is single-threaded, so LSF doesn't speed it up. Just blocks your shell for ~46 min.

**Validate — ESRRB p-value byte-equivalent to Daniel:**
```bash
python3 -c "
import csv
with open('analysis/daniel/outputs/phase5/biobin/merged_maf.001_noRels_keepHLcases-bins.csv') as f:
    r = list(csv.reader(f))
genes = r[0][2:]
ps = r[8][2:]
for g, p in zip(genes, ps):
    if g == 'ESRRB':
        print(f'ESRRB p = {p}')
        # Expected: 8.6308e-05 (byte-equivalent to Daniel)
        break
"
```

**Detailed methodology:** [`results/chapter1_paper_replication/phase5_replication_report.md`](results/chapter1_paper_replication/phase5_replication_report.md)
**biobin reference (if you need it):** [`results/chapter1_paper_replication/phase5_biobin_technical_reference.md`](results/chapter1_paper_replication/phase5_biobin_technical_reference.md)

---

### Ch1 P6 — Exome-wide all-genes burden

**Goal:** biobin chr1-22, exome-wide. Top 5 should match Daniel's preserved list.

**Run (LSF job array, recommended):**
```bash
bash analysis/daniel/scripts/submit_phase6.sh
```
Submits a 22-task array (1 task per chr). Total wall ~2 h. After all chrs finish:
```bash
bash analysis/daniel/scripts/run_phase6_finalize.sh   # concat + BH-FDR
```

**Validate — top 5 non-LOC genes:**
```bash
head -n 1 analysis/daniel/outputs/phase6/results/all_chrom_meta_HL_needAud.txt
awk 'NR>1 && $1 !~ /^LOC|^LINC/' analysis/daniel/outputs/phase6/results/all_chrom_meta_HL_needAud.txt | head -5
# Expected top 5: DNAJC8, UPK3BL1, COL5A1, BOD1, ZNF670 (matches Daniel)
```

**Known issue:** LOC*/LINC* artifacts dominate the raw top 20 due to newer LOKI database. Filter `^LOC|^LINC` to see real genes. See Phase 6 report.

**Detailed methodology:** [`results/chapter1_paper_replication/phase6_replication_report.md`](results/chapter1_paper_replication/phase6_replication_report.md)

---

### Ch1 P7 — Degree-HL linear burden (light mode only)

**Goal:** reproduce Tables 3 & 4 of the paper (TCOF1, ESRRB, COL5A1, HMMR, RAPGEF3, NNT at FDR<0.05).

**Run:**
```bash
bash analysis/daniel/scripts/run_phase8.sh    # YES — script is run_phase8.sh, output is Ch1 P7
```

**What it does (light mode):** loads Daniel's preserved STable `allChr_STable_degHL.txt.gz`, applies paper filters (HL genes ∩ case_carriers>0 → Table 3; non-HL ∩ case_carriers>25 → Table 4), computes BH-FDR. Takes <1 second.

**Outputs:**
- `analysis/daniel/outputs/phase8/light_mode/table3_known_hl_genes.tsv`
- `analysis/daniel/outputs/phase8/light_mode/table4_novel_genes.tsv`

**Validate — all 6 paper genes recovered at FDR<0.05:**
```bash
awk -F'\t' 'NR==1 || /^(TCOF1|ESRRB)/' analysis/daniel/outputs/phase8/light_mode/table3_known_hl_genes.tsv
awk -F'\t' 'NR==1 || /^(COL5A1|HMMR|RAPGEF3|NNT)/' analysis/daniel/outputs/phase8/light_mode/table4_novel_genes.tsv
# Expected: all 6 genes present, FDR column < 0.05
```

**No heavy mode is provided here** — the heavy recipe (biobin → R `lm()`) is documented in the Phase 7 report but intentionally not scripted. We use Daniel's preserved STable directly.

**Detailed methodology:** [`results/chapter1_paper_replication/phase7_degree_hl_burden.md`](results/chapter1_paper_replication/phase7_degree_hl_burden.md)

---

## Light mode vs heavy mode (when to choose which)

**Light mode** = use Daniel's preserved intermediates from `data/PMBB_Exome/` as inputs, skip the expensive compute. Validates downstream logic against known-good inputs. Recommended for first-time reproduction and most validation work.

**Heavy mode** = re-derive intermediates from raw PMBB v2 data (`data/pmbb_v2/...`). Validates the full pipeline including raw-data ingestion. Required if porting to PMBB v3/v4 or if you suspect schema drift.

| Phase | Light works? | Heavy required if... |
|---|---|---|
| Ch1 P1 | — (just runs from raw, fast) | always run as-is |
| Ch1 P2 | ✅ default | porting to v3 (different pVCFs) |
| Ch1 P3 | ✅ default | porting to v3, OR validating plink invocation |
| Ch1 P4 | ✅ default | source CSVs change |
| Ch1 P5 | ✅ ish (light = use Daniel's already-extracted geno + run biobin) | actually re-running biobin is the point |
| Ch1 P6 | n/a (always heavy — that's the test) | — |
| Ch1 P7 | ✅ only | heavy not scripted (intentional) |

When in doubt, default to **light mode** — it validates 80% of correctness for 5% of the time.

---

## If something fails

1. **Re-run `bash scripts/setup_env.sh`** — confirm environment is still good
2. **Read the script header** — every `run_*.sh` has a docstring listing inputs/outputs
3. **Check the log** — `analysis/daniel/logs/phaseN/run_*.log` has timestamped error messages
4. **Read the phase report** — known gotchas are documented in each report's "Issues encountered" section
5. **For biobin issues** — see [`results/chapter1_paper_replication/phase5_biobin_technical_reference.md`](results/chapter1_paper_replication/phase5_biobin_technical_reference.md)
6. **For LOKI / variant assignment** — see "LOKI database drift" in [`chapter1_summary.md`](results/chapter1_paper_replication/chapter1_summary.md)
