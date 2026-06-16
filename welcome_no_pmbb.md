# Welcome, Elena — Chapter 1 Replication, **No-PMBB-Access Mode**

Welcome to the project. This is the version of [`welcome.md`](welcome.md) for when you **don't yet have access to the raw PMBB release** (`/static/PMBB/`) or the `ritchie` LPC group. You can still validate the **entire Chapter 1 result** using the preserved intermediates that already live in the project.

**Estimated read: 12 minutes.**

**Maintained by:** Andre Rico
**Hall Lab × Epstein Lab project**

---

## 1. Your scope, in 90 seconds

You're replicating **Chapter 1** — Daniel Hui's PMBB hearing-loss analysis as published in [Hui et al. 2023 PLOS Genetics](docs/papers/pgen.1010584.pdf) (PMBB v2, 43,731 exomes). The headline result is 6 genes at FDR<0.05: **TCOF1, ESRRB** (known HL genes — Table 3) and **COL5A1, HMMR, RAPGEF3, NNT** (novel — Table 4).

**The twist for you:** you'll validate all 7 phases **without reading the raw PMBB data**. Every intermediate Daniel produced is preserved (gzipped) in `data/PMBB_Exome/`, so you can re-run the light phases and validate the preserved outputs of the heavy ones — getting the full published result without `ritchie` access.

**What success looks like:**
- All 7 Chapter 1 phases validated in no-PMBB mode (3 re-run, 4 inspect-preserved),
- Your numbers match **our results** — the authoritative reference [`chapter1_authoritative_pvalues.md`](results/chapter1_paper_replication/chapter1_authoritative_pvalues.md), **not** the paper. Our re-run differs from the published paper for two documented reasons (LOKI drift in P6, iteration drift in P7); matching the paper exactly is **not** the bar. See "Your validation target" in [`REPRODUCTION_GUIDE_no_pmbb.md`](REPRODUCTION_GUIDE_no_pmbb.md).
- A short replication summary report.

**Out of scope:** **Chapter 2** (ZNF175 follow-up) — in progress, not validated. Don't run/modify anything under `results/chapter2_znf175_analysis/`. If you see "Phase 2.3", "signal-loss diagnostic", or "Joe Park outreach", those are open items — not yours.

---

## 2. Why no-PMBB mode works

The raw PMBB v2 release (`/static/PMBB/PMBB-Release-2020-2.0/`) needs the `ritchie` group. You may not have it yet. But you don't need it for Chapter 1 validation, because:

- The **light phases** (P2 SNP-IDs, P3 genotype-decompress, P7 degree-HL) already run only on preserved intermediates.
- The **heavy phases** (P1 gene-list, P4 phenotype, P5 & P6 biobin) each produced an output that is **preserved in `data/PMBB_Exome/`** — you validate those directly instead of re-deriving them from raw.

So the only data source is `data/PMBB_Exome/`. No `/static/PMBB/`, no biobin, no plink.

> ⚠️ **Governance note.** The preserved intermediates are still **PMBB-derived individual-level data**. No-PMBB mode removes the *filesystem* need for `/static/PMBB/` — it does **not** remove the need to be authorized to view PMBB data. You must be on the project (`hall` group). If you're not authorized to view PMBB-derived data at all, stop and talk to Nikki or Andre — this mode won't apply and we'd use synthetic data instead.

---

## 3. First-day setup

### a. Filesystem access — you only need `hall`

```bash
groups | tr ' ' '\n' | grep -E 'hall|ritchie'
```
- **`hall`** — required (project files). If missing, ask Nikki or Andre.
- **`ritchie`** — **not needed** in no-PMBB mode. Don't block on it.

### b. Create your workspace **inside the project**

You work in **your own folder inside the shared project** — not a separate clone. This is the fast path: the preserved data is already right there at `data/PMBB_Exome/`, so there is nothing to clone, copy, or symlink.

```bash
cd /project/hall/analysis/hearing-loss-genomics
mkdir -p analysis/elena          # your workspace: notes, summary, scratch
```

- **Data** is already at the project root (`data/PMBB_Exome/`) — read it directly; you have access through the `hall` group. Nothing to set up.
- **Write only inside `analysis/elena/`** — your summary report, notes, and any validation output you redirect there (e.g. `... > analysis/elena/p5_esrrb.txt`).
- **The 4 validate-only phases** (P1, P4, P5, P6) are pure inspection — run the command, redirect anything you want to keep into your folder.
- **The 3 light run-scripts** (P2, P3, P7) write to the shared `analysis/daniel/outputs/` (they pin the project root internally). That's safe — they're deterministic: identical preserved inputs produce identical outputs. If you want fully separate outputs, copy those three scripts into `analysis/elena/scripts/` and edit `PROJECT_ROOT` / `OUT_DIR` at the top.

> ⚠️ **The one rule that keeps the project safe:** only ever **write inside `analysis/elena/`**, and **never modify the shared `data/` symlinks** (or anyone else's folder). Reading shared data is fine — changing shared state is what broke everyone's paths last time.

### c. Activate the environment

The shared project already has a configured Python venv. Just activate it — **do not run `scripts/setup_env.sh`**: that script assumes a personal project root and would tell you to *move* the shared venv (which you must not do).

```bash
cd /project/hall/analysis/hearing-loss-genomics
source venv/bin/activate
```

Validate-only mode needs **nothing under `/project/ritchie/`** — no `ritchie` group, no biobin, no plink. Just the shared venv plus standard tools (`zcat`, `awk`, `md5sum`). If you ever need a Python package that's missing, make your **own** venv in `analysis/elena/venv` rather than touching the shared one.

### d. Light reading (in order)

1. This document
2. [`results/chapter1_paper_replication/chapter1_summary.md`](results/chapter1_paper_replication/chapter1_summary.md) — the replication narrative
3. [`REPRODUCTION_GUIDE_no_pmbb.md`](REPRODUCTION_GUIDE_no_pmbb.md) — **your step-by-step**, no-PMBB version
4. [`docs/data_inventory.md`](docs/data_inventory.md) — data paths

Skip `STATUS_SNAPSHOT.md`, `docs/communications/`, and `results/chapter2_znf175_analysis/`.

---

## 4. Repository map (only what you need)

You work inside the shared project, in your own `analysis/elena/` folder:

```
/project/hall/analysis/hearing-loss-genomics/   ← the shared project
├── welcome_no_pmbb.md            ← you are here
├── REPRODUCTION_GUIDE_no_pmbb.md ← your step-by-step (no-PMBB)
├── data/PMBB_Exome/              ← preserved intermediates — already here, read directly
│       (you do NOT touch data/pmbb_v2/ or /static/PMBB/ in this mode)
├── analysis/
│   ├── elena/                    ← YOUR folder: summary, notes, scratch (write only here)
│   └── daniel/scripts/run_phase2.sh, run_phase3.sh, run_phase8.sh   ← the 3 you run
├── venv/                         ← shared Python env — activate, don't rebuild
└── results/chapter1_paper_replication/   ← 7 phase reports + authoritative p-values
```

---

## 5. The 7 phases in no-PMBB mode

| Phase | Mode | What you do |
|---|---|---|
| Ch1 P1 — gene list | **validate-only** | inspect preserved `annot_genes_full_funcToInclude.txt.gz` (11,661 rows) |
| Ch1 P2 — SNP IDs | **run** | `run_phase2.sh` (light) |
| Ch1 P3 — plink extract | **run** | `run_phase3.sh` (light only — never `PHASE3_MODE=heavy`) |
| Ch1 P4 — prep files | **validate-only** | inspect preserved `cases_control` / `covs` / `regions` |
| Ch1 P5 — HL burden | **validate-only** | inspect preserved biobin → ESRRB p = 8.6308e-05 |
| Ch1 P6 — exome-wide | **validate-only** | inspect preserved `all_chrom_meta_withBH.txt.gz` |
| **Ch1 P7 — degree-HL** | **run** | **`run_phase8.sh`** (NOT `run_phase7.sh`) |

Exact commands and validation criteria: [`REPRODUCTION_GUIDE_no_pmbb.md`](REPRODUCTION_GUIDE_no_pmbb.md).

> ⚠️ **Do NOT run** `run_phase1.sh`, `run_phase4.sh`, `run_phase5.sh`, `submit_phase6.sh`, or `run_phase3.sh` in heavy mode — they read `/static/PMBB/` (or biobin) and will fail without `ritchie`. Validate their preserved outputs instead.
> ⚠️ **`run_phase7.sh` is Chapter 2, not Chapter 1 Phase 7.** Use `run_phase8.sh`.

---

## 6. Suggested plan (≈1–2 days)

Because the heavy compute is replaced by inspection, no-PMBB mode is fast.

**Day 1 — orientation + validate**
- Create your folder (`mkdir -p analysis/elena`) and activate the shared venv: `source venv/bin/activate` (no `setup_env.sh` — see §3.c).
- Read [`chapter1_summary.md`](results/chapter1_paper_replication/chapter1_summary.md) and skim each phase report.
- Open the [Hui et al. 2023 PDF](docs/papers/pgen.1010584.pdf) — Abstract, Figures 2–4, Tables 3 & 4.
- Walk through P1→P7 in [`REPRODUCTION_GUIDE_no_pmbb.md`](REPRODUCTION_GUIDE_no_pmbb.md). The 3 "run" phases take seconds; the 4 "validate-only" phases are one command each.

**Day 2 — write up + sync**
- Write a short replication summary (what you validated, any deviations) into `analysis/elena/`.
- Sync with Nikki or Andre. If you want to also re-run the heavy phases from raw, that's the trigger to request `ritchie` access and switch to the full [`REPRODUCTION_GUIDE.md`](REPRODUCTION_GUIDE.md).

---

## 7. Common gotchas

1. **You don't run `setup_env.sh`** in this mode — activate the shared venv instead (§3.c). Validate-only needs only the `hall` group; no `ritchie`, biobin, or plink.
2. **`run_phase7.sh` is NOT Chapter 1 Phase 7.** Use `run_phase8.sh`.
3. **Never set `PHASE3_MODE=heavy`** — that needs the raw chr21 pVCF. Plain `run_phase3.sh` is light.
4. **Don't run P1/P4/P5/P6 scripts** — they read raw PMBB / biobin. Validate the preserved outputs (steps in the guide).
5. **`LOC*`/`LINC*` dominate P6's top hits** — newer LOKI database artifact. Filter them to see real genes.
6. **Some `.txt.gz` lack a header line** — common in Daniel's preserved outputs; phase reports note these.
7. **Write only inside `analysis/elena/`** — never repoint the shared `data/` symlinks. Reading shared data is fine; changing shared state breaks everyone's paths.

---

## 8. If you later get `ritchie` + `/static` access

Nothing here is wasted — you've validated the full result. To *also* re-run the heavy phases from raw data (validates raw-data ingestion, needed for any future v3/v4 port), switch to [`welcome.md`](welcome.md) + [`REPRODUCTION_GUIDE.md`](REPRODUCTION_GUIDE.md) and run the heavy/full modes of P1, P4, P5, P6.
