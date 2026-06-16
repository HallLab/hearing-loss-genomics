# Welcome, Elena — Chapter 1 Replication Guide

Welcome to the project. This document is written specifically for you, with a clear scope: **replicate Chapter 1 (Hui et al. 2023 paper on PMBB v2)**. You're not picking up the open analytical work — that's still in development and unvalidated. Your work focuses on the published, validated paper-replication phases.

**Estimated read: 15 minutes.**

**Maintained by:** Andre Rico (on vacation — see section 9 if you need to reach me)
**Hall Lab × Epstein Lab project**

---

## 1. Your scope, in 90 seconds


You're replicating **Chapter 1** of the project, which reproduces Daniel Hui's PMBB hearing-loss analysis as published in [Hui et al. 2023 PLOS Genetics](docs/papers/pgen.1010584.pdf). The paper used PMBB v2 (`PMBB-Release-2020-2.0`, 43,731 exomes). Your goal is to re-run the 7 phases of Chapter 1 and confirm you get the same numerical results as the published paper and as my prior runs.

**What success looks like:**

- All 7 phases of Chapter 1 execute cleanly in your environment,
- Your outputs match either Daniel's preserved intermediates or my prior re-run outputs (validation criteria are spelled out in each phase)
- You produce a short summary report of your replication run

**What's out of scope for you:**

- **Chapter 2** (ZNF175 follow-up) — work in progress, not validated. Files exist in `results/chapter2_znf175_analysis/` but please **don't run, modify, or rely on those for replication**. If you see references to "Phase 2.3", "signal-loss diagnostic", "Joe Park outreach" anywhere — those are open work items.

If you finish Chapter 1 ahead of schedule and want more to do, talk to Molly or Nikki

---

## 2. Chapter 1, in one paragraph

Hui et al. 2023 ran rare-variant gene-burden tests across the exome in PMBB v2 to find genes whose deleterious variants correlate with adult hearing loss. The paper's headline result is 6 genes at FDR<0.05: **TCOF1, ESRRB** (already known HL genes — Table 3) and **COL5A1, HMMR, RAPGEF3, NNT** (novel candidates — Table 4). Replicating the paper means walking through 7 phases: building the gene list, reconciling variant IDs, extracting genotypes, preparing covariates, running the burden test on known HL genes, then exome-wide, and finally the degree-of-HL linear regression that produces Tables 3 & 4.

For the full narrative — what each phase does, what genes show up where, why the methodology matters — read [`results/chapter1_paper_replication/chapter1_summary.md`](results/chapter1_paper_replication/chapter1_summary.md). It's the single best document for understanding what you're replicating.

---

## 3. First-day setup

### a. Get filesystem access

You need membership in two LPC groups:
- `hall` — for project files at `/project/hall/analysis/hearing-loss-genomics/`
- `ritchie` — for raw PMBB data at `/static/PMBB/PMBB-Release-2020-2.0/`

Check with:
```bash
groups | tr ' ' '\n' | grep -E 'hall|ritchie'
```
If either is missing, ask Molly to request access via LPC IT.

### b. Run the setup script

```bash
cd /project/hall/analysis/hearing-loss-genomics
bash scripts/setup_env.sh
```

This validates LPC tool paths (biobin, plink, R, Python), creates your own `venv/`, and installs Python dependencies. **All 6 checks should pass.** If any fail, fix before continuing — running phases with a broken environment wastes time.

### c. Light reading (in this order)

1. This document (you're here)
2. [`results/chapter1_paper_replication/chapter1_summary.md`](results/chapter1_paper_replication/chapter1_summary.md) — the paper-replication narrative
3. [`REPRODUCTION_GUIDE.md`](REPRODUCTION_GUIDE.md) — step-by-step commands for each phase
4. ['docs/data_inventory.md'](data_inventory.md) — Data Inventory used and paths

You can safely skip `STATUS_SNAPSHOT.md`, `docs/communications/`, and anything under `results/chapter2_znf175_analysis/` — those are for the project-continuation.

---

## 4. Repository map (only what you need)

```
hearing-loss-genomics/
├── HANDOFF.md                 ← you are here
├── README.md                  ← public-facing description
├── REPRODUCTION_GUIDE.md      ← step-by-step for each Ch1 phase
├── scripts/
│   └── setup_env.sh           ← run on first day
├── data/                      ← GITIGNORED, lives on disk only
│   ├── PMBB_Exome/            ← Daniel's preserved intermediates (5.6 GB)
│   ├── pmbb_v2/ → /static/PMBB/PMBB-Release-2020-2.0/   (raw PMBB v2)
│   └── ... (other dirs not relevant for you)
├── analysis/daniel/
│   ├── runbook_hui2023.txt    ← Daniel's decompressed cookbook (reference only)
│   ├── scripts/
│   │   └── run_phase1.sh ... run_phase8.sh   ← run these (see section 5)
│   ├── outputs/phaseN/        ← your runs will write here
│   └── logs/phaseN/           ← run logs
├── results/
│   ├── README.md              ← results index
│   ├── chapter1_paper_replication/   ← your target — 7 phase reports + chapter1_summary
│   └── chapter2_znf175_analysis/     ← ignore — not validated, in-progress
└── venv/                      ← your Python env (created by setup_env.sh)
```

---

## 5. Script numbering — IMPORTANT to read carefully

The orchestration scripts in [`analysis/daniel/scripts/`](analysis/daniel/scripts/) are numbered linearly by execution order, **not** by chapter/phase. So one of the scripts has a confusing name. For Chapter 1, use these scripts in this order:

| Phase | Script to run | What it produces |
|---|---|---|
| Ch1 P1 — gene list | `run_phase1.sh` | 173-gene HL set + filtered annotation |
| Ch1 P2 — SNP IDs | `run_phase2.sh` | 9,667-SNP `.extract` file |
| Ch1 P3 — plink extract | `run_phase3.sh` | per-chr `.bed`/`.bim`/`.fam` |
| Ch1 P4 — prep files | `run_phase4.sh` | case/control + covariates + region file |
| Ch1 P5 — HL burden | `run_phase5.sh` | biobin on 173 HL genes → ESRRB top hit |
| Ch1 P6 — exome-wide burden | `submit_phase6.sh` + `run_phase6_finalize.sh` | biobin on all genes |
| **Ch1 P7 — degree-HL burden** | **`run_phase8.sh`** ← **YES, phase8 produces Ch1 P7** | Tables 3 & 4 of the paper |

**Why `run_phase8.sh` for Chapter 1 Phase 7?** Because `run_phase7.sh` is a Chapter 2 script (ZNF175 deep-dive) that was written before the Chapter 1 degree-HL gap was identified. Renaming would touch many references and risk breakage, so we deferred it. **Do NOT run `run_phase7.sh` — that's Chapter 2, not your scope.**

Every script has a docstring at the top stating exactly what it produces. When in doubt, read the script header before running.

---

## 6. Suggested day-by-day plan

This is a flexible suggestion — adjust based on how fast you move and what you want to learn.

### Day 1 — setup + orientation

- Run `bash scripts/setup_env.sh`
- Read [`chapter1_summary.md`](results/chapter1_paper_replication/chapter1_summary.md) (15 min)
- Skim each individual Ch1 phase report in `results/chapter1_paper_replication/`
- Open the [Hui et al. 2023 paper PDF](docs/papers/pgen.1010584.pdf) — read at least the Abstract, Figures 2-4, Tables 3 & 4
- Run Phase 1 (light mode, < 1 min) as a sanity check that your environment works end-to-end

### Day 2-3 — light-mode replication

Run each phase in **light mode** (uses Daniel's preserved intermediates as inputs — fast, validates downstream logic). Commands and validation criteria are in [`REPRODUCTION_GUIDE.md`](REPRODUCTION_GUIDE.md).

- Phase 1 — gene list (< 1 min)
- Phase 2 — SNP IDs (~20 s in light mode)
- Phase 3 — plink extraction (~2 s in light mode)
- Phase 4 — preparatory files (< 1 min)
- Phase 5 — HL gene burden test (**~46 min** because biobin is single-threaded — start this and read paper sections while it runs)

After each phase, run the validation commands documented in the report and in `REPRODUCTION_GUIDE.md`. Note any deviations from expected values.

### Day 4 — Phase 6 (exome-wide)

This is the longest phase. Submit via LSF and let it run in the background:
```bash
bash analysis/daniel/scripts/submit_phase6.sh
```
~2 h on 22-task array. After all chrs finish, run `run_phase6_finalize.sh` to concat + apply BH-FDR.

While Phase 6 runs, you can also do Phase 7 (it's < 1 second, light mode only):
```bash
bash analysis/daniel/scripts/run_phase8.sh   # YES — phase8 = Ch1 P7
```
Validate that all 6 paper genes (TCOF1, ESRRB, COL5A1, HMMR, RAPGEF3, NNT) appear at FDR<0.05.

### Day 5 — write up + sync with Molly

Write a short replication summary (template below in section 8). Sync with Molly to walk her through your results.

**If you finish in 3 days instead of 5, talk to Molly about what to do next.** Do not extend into Chapter 2 unilaterally.

---

## 7. Common gotchas (saves you hours)

1. **`run_phase7.sh` is NOT Chapter 1 Phase 7.** It's Chapter 2 Phase 1. Use `run_phase8.sh` for Chapter 1 Phase 7. See section 5.

2. **biobin needs a liblzma shim on RHEL 9.** All `run_phaseN.sh` scripts already handle this with `LD_LIBRARY_PATH`. Don't try to run biobin standalone without the shim.

3. **plink 1.9 vs plink 2.0.** Default `plink` on PATH is plink 2.0 which over-reserves RAM and gets killed by LSF. Use `/appl/plink-1.90Beta6.18/plink` explicitly. All scripts already do this.

4. **R 3.6.3 is the default but broken** (libreadline issue). Use `/appl/R-4.4/bin/Rscript`. All scripts already do this.

5. **Daniel's original scripts reference `/project/PMBB/...` paths.** That path no longer exists. We use `data/pmbb_v2/...` (symlinked). Don't try to run Daniel's original scripts directly — use `analysis/daniel/scripts/run_phaseN.sh` instead, which has the corrected paths.

6. **Phase 5 takes ~46 min.** biobin is single-threaded — there's no way to speed this up. Plan for it.

7. **Phase 6's top 20 will be dominated by `LOC*`/`LINC*` artifacts.** This is from a newer LOKI database that added pseudogenes. Filter them out to see the real top genes (DNAJC8, UPK3BL1, COL5A1, BOD1, ZNF670). The Phase 6 report documents this.

8. **Phase 7 only has light mode.** Daniel's actual degree-HL pipeline used R `lm()` separately, not biobin's built-in linear test. We use his preserved supplementary table directly. Don't try to write a heavy-mode version — it would take days and isn't needed.

9. **Some `.txt.gz` files lack a header line.** Common in Daniel's preserved outputs. The phase reports note these where relevant.

---

