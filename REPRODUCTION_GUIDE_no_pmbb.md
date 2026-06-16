# Reproduction Guide — Chapter 1, **No-PMBB-Access Mode** (validate-only)

**Audience:** Elena — and anyone replicating Chapter 1 **without access to the raw PMBB release** (`/static/PMBB/`) or the `ritchie` group.
**Scope:** Chapter 1 only. Chapter 2 (ZNF175 follow-up) is excluded — not validated.
**Prereq:** membership in the `hall` group (only). Activate the shared project venv (`source venv/bin/activate`) — do **not** run `scripts/setup_env.sh` (see Setup). **No `ritchie` group, no `/static/PMBB`, no biobin, no plink needed.**

> **What this guide is.** A version of [`REPRODUCTION_GUIDE.md`](REPRODUCTION_GUIDE.md) that runs **entirely on the preserved intermediates** in `data/PMBB_Exome/`. It never reads `/static/PMBB/PMBB-Release-2020-2.0/`. Every phase either (a) runs a light script that only touches intermediates, or (b) is **validate-only**: you inspect Daniel's preserved output directly and confirm the published numbers, instead of re-deriving it from raw data.
>
> **What you give up.** The three phases that re-derive a file from raw PMBB (P1 gene-list construction, P4 phenotype build, P5/P6 biobin compute) are not re-run here — you validate their **preserved outputs** instead. The scientific result is identical; you just don't re-execute the heavy compute. If you later get `ritchie` + `/static` access, use the full [`REPRODUCTION_GUIDE.md`](REPRODUCTION_GUIDE.md) to re-run those from scratch.

---

## Why this mode exists

The raw PMBB v2 release lives at `/static/PMBB/PMBB-Release-2020-2.0/` and requires the `ritchie` LPC group. If you don't have that access yet, you can **still validate the entire Chapter 1 result** because every intermediate Daniel produced is preserved (gzipped) inside `data/PMBB_Exome/`. Those preserved files are the inputs and outputs of every phase.

| Resource | Full mode | This mode (no-PMBB) |
|---|---|---|
| `ritchie` group | required | **not needed** |
| `/static/PMBB/...` (raw) | read | **never touched** |
| `data/pmbb_v2/` symlink | followed | **not followed** |
| biobin / plink / LOKI (`/project/ritchie/...`) | required | **not needed** |
| `data/PMBB_Exome/` intermediates | inputs | **the only data source** |

> ⚠️ **Governance note.** The intermediates in `data/PMBB_Exome/` are still **PMBB-derived individual-level data** (case/control by PMBB_ID, carrier lists, audiogram linkage). This mode removes the *filesystem* dependency on `/static/PMBB/`, not the requirement to be authorized to view PMBB data. You must still be on the project (`hall` group). If you are not authorized to view any PMBB-derived data at all, this mode does **not** apply — talk to Nikki or Andre first.

---

## Your validation target: **our results**, not the paper

Our re-run does **not** reproduce the paper's published p-values exactly — and that's expected, for two well-characterized reasons, **neither of which is an error**. So the number you should match is **ours**, not the paper's. The full-precision source of truth is [`results/chapter1_paper_replication/chapter1_authoritative_pvalues.md`](results/chapter1_paper_replication/chapter1_authoritative_pvalues.md) — when any other doc disagrees with it, that doc wins.

| Phase | Reproduces paper? | Your target | Why it differs from the paper |
|---|---|---|---|
| P5 (HL-only, biobin) | **yes, byte-exact** | ESRRB p = 8.6308×10⁻⁵ | — (LOKI not consulted; region-file pins genes) |
| P6 (exome-wide, biobin) | no — ~30%–5× off | **our `loki-20230816` values** (P6 table below) | **LOKI version drift** |
| P7 (degree-HL, lm) | close — within ~1 order of magnitude | **our STable values** (Tables 3 & 4) | **iteration drift** |

**Two drift sources — keep them separate, they are unrelated:**
- **LOKI drift (affects P6 only):** biobin's `--bin-regions` mode asks LOKI which gene each variant belongs to. Our LOKI (`loki-20230816`) differs from the paper's (`loki-20220926`, no longer on LPC) → different gene bins → shifted p-values. All 6 genes still rank top-tier.
- **Iteration drift (affects P7 only):** the preserved degree-HL STable is one of Daniel's *intermediate* runs (its carrier counts differ from the paper's), not the final submitted run. Both sides used the **same** LOKI — so this is **not** a LOKI issue.

→ **What this means for you:** a phase "passes" when your output matches **our** numbers (within tolerance), not the paper's. Chasing the paper's exact values would be chasing known, documented drift.

---

## Setup (no-PMBB)

You work in **your own folder inside the shared project** — no clone, no copy. The preserved data is already at `data/PMBB_Exome/`.

```bash
cd /project/hall/analysis/hearing-loss-genomics
mkdir -p analysis/elena          # your workspace: notes, summary, scratch
source venv/bin/activate         # reuse the shared venv — do NOT run setup_env.sh
```

- **Do not run `scripts/setup_env.sh`** — it assumes a personal project root and would tell you to *move* the shared venv. Just activate the existing one.
- **Write only inside `analysis/elena/`.** The 4 validate-only phases are pure inspection (redirect output there if you want to keep it). The 3 light run-scripts (P2/P3/P7) pin the project root and write to the shared `analysis/daniel/outputs/` — that's safe, they're deterministic; copy + repoint them into `analysis/elena/scripts/` if you want isolated outputs.
- **Never modify the shared `data/` symlinks** — that's what broke everyone's paths last time.

You need only: the `hall` group (ask Nikki or Andre if missing), the shared venv, and standard tools (`gunzip`/`zcat`, `awk`, `md5sum`). **Nothing under `/project/ritchie/`** — no `ritchie` group, no biobin, no plink.

---

## Quick reference

| Phase | No-PMBB action | How | Time |
|---|---|---|---:|
| Ch1 P1 — gene list | **validate-only** (preserved) | inspect `annot_genes_full_funcToInclude.txt.gz` | < 1 s |
| Ch1 P2 — SNP IDs | **run** (light, intermediates) | `run_phase2.sh` | ~20 s |
| Ch1 P3 — plink extract | **run** (light, intermediates) | `run_phase3.sh` (light only) | ~2 s |
| Ch1 P4 — prep files | **validate-only** (preserved) | inspect `cases_control` / `covs` / `regions` | < 1 s |
| Ch1 P5 — HL burden | **validate-only** (preserved) | inspect biobin `-bins.csv.gz` → ESRRB | < 1 s |
| Ch1 P6 — exome-wide burden | **validate-only** (preserved) | inspect `all_chrom_meta_withBH.txt.gz` | < 1 s |
| **Ch1 P7 — degree-HL burden** | **run** (light, intermediates) | `run_phase8.sh` (NOT `run_phase7.sh`) | < 1 s |

> ⚠️ **Do NOT run** `run_phase1.sh`, `run_phase4.sh`, `run_phase5.sh`, `submit_phase6.sh`, or the chr21 heavy pilot of `run_phase3.sh` — those read `/static/PMBB/` (or biobin) and will fail without `ritchie`. In this mode you validate their preserved outputs instead (steps P1, P4, P5, P6 below).
> ⚠️ **Do NOT run `run_phase7.sh`** — that's a Chapter 2 script. Chapter 1 Phase 7 = `run_phase8.sh`.

---

## Ch1 P1 — Gene list / variant×gene table (validate-only)

Full mode rebuilds `annot_genes_full_funcToInclude.txt` by filtering the 5.4 GB raw annotation. That raw file is under `/static/PMBB/`. In no-PMBB mode, validate the **preserved output** directly.

```bash
zcat data/PMBB_Exome/annot_genes_full_funcToInclude.txt.gz | wc -l
# Expected: 11,661 rows (the master variant×gene table; input to every downstream phase)
```

**Detailed methodology:** [`results/chapter1_paper_replication/phase1_replication_report.md`](results/chapter1_paper_replication/phase1_replication_report.md)

---

## Ch1 P2 — SNP ID reconciliation (run, light)

Already intermediates-only — uses Daniel's preserved `vcf_SNP_IDs_allchr.txt.gz`. Never touches raw pVCFs.

```bash
source venv/bin/activate
bash analysis/daniel/scripts/run_phase2.sh
```

**Validate:**
```bash
sort -u analysis/daniel/outputs/phase2/matched_snp_IDs_annot_pVCF_noNA_noMultiallelic.extract | md5sum
# Expected: 5e80ebc0faa5e68277cfeb948af8b1da
```

**Detailed methodology:** [`results/chapter1_paper_replication/phase2_replication_report.md`](results/chapter1_paper_replication/phase2_replication_report.md)

---

## Ch1 P3 — plink genotype extraction (run, light only)

Light mode **decompresses Daniel's preserved per-chr `.bed/.bim/.fam`** — no plink, no raw pVCF.

```bash
bash analysis/daniel/scripts/run_phase3.sh
```

> ⚠️ Do **not** set `PHASE3_MODE=heavy` — that runs plink on the raw chr21 pVCF and needs `/static/PMBB/`. Plain `run_phase3.sh` is light and is all you need.

**Validate:** confirm 22 per-chr filesets landed in `analysis/daniel/outputs/phase3/light/`:
```bash
ls analysis/daniel/outputs/phase3/light/allIndvs_chr*.bed | wc -l   # Expected: 22
```

**Detailed methodology:** [`results/chapter1_paper_replication/phase3_replication_report.md`](results/chapter1_paper_replication/phase3_replication_report.md)

---

## Ch1 P4 — Preparatory files (validate-only)

Full mode rebuilds `covs.txt` from the raw `Phenotype/2.0/...covariates.txt` under `/static/PMBB/`. In no-PMBB mode, validate the **preserved outputs**.

```bash
zcat data/PMBB_Exome/cases_control.txt.gz   | wc -l   # case/control assignments
zcat data/PMBB_Exome/covs.txt.gz            | wc -l   # covariates (age, sex, PCs)
zcat data/PMBB_Exome/gene_list_regions.txt.gz | head  # chr:start-end per gene (biobin regions)
```

These three are the exact inputs Phase 5's burden test consumes. Confirm they decompress cleanly and have the expected columns (see report).

**Detailed methodology:** [`results/chapter1_paper_replication/phase4_replication_report.md`](results/chapter1_paper_replication/phase4_replication_report.md)

---

## Ch1 P5 — First burden test, HL genes (validate-only)

Full mode re-runs biobin (needs the `ritchie` biobin tool) and re-derives the relatedness keep-list from raw IBD files under `/static/PMBB/`. In no-PMBB mode, validate the **preserved biobin result** — the headline is that ESRRB is the top hit.

```bash
python3 -c "
import csv, gzip
with gzip.open('data/PMBB_Exome/biobin/merged_maf.001_noRels_keepHLcases-bins.csv.gz','rt') as f:
    r = list(csv.reader(f))
genes, ps = r[0][2:], r[8][2:]
for g,p in zip(genes,ps):
    if g=='ESRRB':
        print(f'ESRRB p = {p}')   # Expected: 8.6308e-05
        break
"
```

**Detailed methodology:** [`results/chapter1_paper_replication/phase5_replication_report.md`](results/chapter1_paper_replication/phase5_replication_report.md)

---

## Ch1 P6 — Exome-wide all-genes burden (validate-only)

Full mode runs a 22-task biobin array (~2 h, needs biobin) then concatenates + BH-FDR. In no-PMBB mode, validate the **preserved meta result** (already concatenated, with BH-corrected p):

```bash
F=data/PMBB_Exome/allGenes/HL_needAud/meta_results/all_chrom_meta_withBH.txt.gz
zcat "$F" | head -1                                  # Gene  Beta  SE  p  BH-corrected_p
zcat "$F" | awk -F'\t' 'NR>1 && $1 !~ /^LOC|^LINC/' | sort -t$'\t' -gk4 | head -8
# Top non-LOC genes include the paper's novel candidates: UPK3BL1, BOD1, ZNF670, COL5A1, DNAJC8
```

**Our target p-values for the 6 headline genes** (this is what you target, *not* the paper):

| Gene | **Our target** (`loki-20230816`) | Paper / Daniel (`loki-20220926`) |
|---|---:|---:|
| TCOF1   | 2.10×10⁻³ | 1.60×10⁻³ |
| ESRRB   | 2.40×10⁻⁴ | 1.13×10⁻³ |
| COL5A1  | 2.29×10⁻⁵ | 2.69×10⁻⁵ |
| HMMR    | 1.30×10⁻³ | 4.20×10⁻³ |
| RAPGEF3 | 9.65×10⁻⁵ | 1.46×10⁻⁴ |
| NNT     | 4.30×10⁻³ | 1.60×10⁻³ |

> ⚠️ **No-PMBB caveat for P6:** our exact `loki-20230816` numbers come from biobin, which is `ritchie`-gated — so in no-PMBB mode you **cannot regenerate them yourself**. The preserved file you inspect above is *Daniel's* (`loki-20220926`), so its p-values match the right-hand column. In no-PMBB mode, P6 validation = **the 6 genes rank top-tier** in the preserved result; the precise our-vs-paper comparison is recorded in [`chapter1_authoritative_pvalues.md`](results/chapter1_paper_replication/chapter1_authoritative_pvalues.md) and becomes reproducible only once biobin access is available.

> **Known issue:** `LOC*`/`LINC*` pseudogenes dominate the raw top hits due to a newer LOKI database. Filter `^LOC|^LINC` to see the real genes. See the Phase 6 report.

**Detailed methodology:** [`results/chapter1_paper_replication/phase6_replication_report.md`](results/chapter1_paper_replication/phase6_replication_report.md)

---

## Ch1 P7 — Degree-HL linear burden (run, light)

Already intermediates-only — loads the preserved supplementary table `allChr_STable_degHL.txt.gz` and applies the paper filters. No raw, no biobin.

```bash
bash analysis/daniel/scripts/run_phase8.sh    # YES — script is run_phase8.sh, output is Ch1 P7
```

**Validate — all 6 paper genes at FDR<0.05:**
```bash
awk -F'\t' 'NR==1 || /^(TCOF1|ESRRB)/'              analysis/daniel/outputs/phase8/light_mode/table3_known_hl_genes.tsv
awk -F'\t' 'NR==1 || /^(COL5A1|HMMR|RAPGEF3|NNT)/'  analysis/daniel/outputs/phase8/light_mode/table4_novel_genes.tsv
```

**Your target = our STable-derived values, not the paper.** Phase 7 runs on Daniel's preserved degree-HL STable (our Phase 7 source). The bar is **all 6 genes recovered at FDR<0.05**. The exact β/p differ from the paper by **iteration drift** (within ~1 order of magnitude; HMMR is the worst at ~6×) — expected, documented, not a failure. Full-precision our-vs-paper values: [`chapter1_authoritative_pvalues.md`](results/chapter1_paper_replication/chapter1_authoritative_pvalues.md) (Phase 7 table).

**Detailed methodology:** [`results/chapter1_paper_replication/phase7_degree_hl_burden.md`](results/chapter1_paper_replication/phase7_degree_hl_burden.md)

---

## Summary: what you validated without ever touching `/static/PMBB/`

- **P1** — 11,661-row variant×gene table (preserved)
- **P2** — 9,667 reconciled SNP IDs (re-run, md5 match)
- **P3** — 22 per-chr genotype filesets (re-run, decompressed)
- **P4** — case/control + covariates + regions (preserved)
- **P5** — ESRRB top hit, p = 8.6308e-05 (preserved)
- **P6** — exome-wide meta with BH-FDR; novel candidates in the top (preserved)
- **P7** — Tables 3 & 4: TCOF1, ESRRB, COL5A1, HMMR, RAPGEF3, NNT at FDR<0.05 (re-run)

That is the complete Chapter 1 result. When/if `ritchie` + `/static` access comes through, re-run P1/P4/P5/P6 from raw via the full [`REPRODUCTION_GUIDE.md`](REPRODUCTION_GUIDE.md) to also validate raw-data ingestion.

## If something fails

1. Confirm you're in the `hall` group (`groups | grep hall`) and the shared venv is active (`source venv/bin/activate`). Don't run `setup_env.sh`.
2. Read the script header — every `run_*.sh` has a docstring listing inputs/outputs.
3. Check the log — `analysis/daniel/logs/phaseN/run_*.log`.
4. Read the phase report — known gotchas under "Issues encountered".
