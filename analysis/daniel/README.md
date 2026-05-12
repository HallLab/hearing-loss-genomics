# `analysis/daniel/` — Hui et al. 2023 pipeline replication (PMBB v2)

Workspace for **re-executing Daniel Hui's original pipeline** on PMBB v2 data, to confirm we reproduce the published numbers from [Hui et al. 2023, PLOS Genetics](../../docs/papers/pgen.1010584.pdf) (doi:10.1371/journal.pgen.1010584).

This is the **active Phase 1 workspace.** We use **PMBB v2** (`PMBB-Release-2020-2.0`) — the exact release Daniel used — to verify replicability before porting to newer releases. Porting to v3 (`PMBB-Release-2024-3.0`) and v4 (`PMBB-Release-2026-4.0`) is a separate, later phase.

Data access:
- Daniel's v2 intermediates: [`data/PMBB_Exome/`](../../data/PMBB_Exome/), [`data/PMBB_Imputed/`](../../data/PMBB_Imputed/), [`data/DFNA/`](../../data/DFNA/) — plink/biobin outputs from the 2021 runs
- Raw v2 release: [`data/pmbb_v2/`](../../data/pmbb_v2/) symlink → `/static/PMBB/PMBB-Release-2020-2.0/` — original pVCFs, annotations, IBD, PCA, phenotypes

## Layout

```
analysis/daniel/
├── README.md                 # this file
├── runbook_hui2023.txt       # decompressed verbatim from data/PMBB_Exome/README.gz
│                             # — Daniel's original 1,262-line bash cookbook
├── scripts/
│   ├── pmbb_exome/           # 81 .py + .R scripts (was data/PMBB_Exome/scripts/*.gz)
│   ├── pmbb_imputed/         # 16 scripts (GWAS, PRS) (was data/PMBB_Imputed/scripts/*.gz)
│   └── dfna/                 # 5 scripts (was data/DFNA/scripts/*.gz)
├── configs/                  # parameter files for our re-runs (empty — to be populated)
├── outputs/                  # outputs from our re-runs (empty — to be populated)
└── logs/                     # job logs (LSF / SLURM / interactive)
```

## What "replicate" means here

Daniel's pipeline produced intermediate files that are still in [`data/PMBB_Exome/`](../../data/PMBB_Exome/), [`data/PMBB_Imputed/`](../../data/PMBB_Imputed/), and [`data/DFNA/`](../../data/DFNA/). We have **two complementary modes** of replication:

| Mode | What it does | Cost | What it validates |
|---|---|---|---|
| **Re-run from intermediates** | Re-execute the burden/regression scripts on the existing plink/biobin intermediates | Hours | Output reproducibility — we get the same numbers from the same intermediates |
| **Re-run from raw pVCF** | Re-execute the full pipeline: annotation → plink extraction → biobin → regression → meta | Days, needs cluster | End-to-end reproducibility — our infrastructure regenerates the published numbers |

Start with **re-run from intermediates** (faster, validates we read inputs/run scripts correctly). Only move to **re-run from raw pVCF** if needed for v3 porting.

## How to use the runbook

[`runbook_hui2023.txt`](runbook_hui2023.txt) is Daniel's working bash log — not a clean, parameterized script. Each block:
1. Defines inputs (often inline)
2. Runs a command (sometimes via `bsub` — IBM LSF; LPC may use SLURM now)
3. Has notes with cohort counts, sanity checks, top hits

To re-run a section: read the cookbook lines, identify the scripts and inputs, run interactively or wrap into a SLURM job. The decompressed scripts live in `scripts/pmbb_exome/`, `scripts/pmbb_imputed/`, `scripts/dfna/`.

The curated phase-by-phase guide is in [`docs/pipeline_walkthrough.md`](../../docs/pipeline_walkthrough.md). Read that first — it maps the cookbook into ~14 phases with line ranges, inputs, outputs, and gotchas.

## Tool requirements

All confirmed available on LPC (`superman` node) as of 2026-05-12:

| Tool | Path / how to load | Notes |
|---|---|---|
| **LSF scheduler** | `bsub`, `bjobs` already on PATH (`/lsf/10.1/...`) | Same as Daniel used — no SLURM translation needed |
| **biobin** | `module load rlsoftware/latest` (puts `biobin` on PATH) or `module load biobin/r4221` | Core burden tool |
| **LOKI database** | `/project/ritchie/datasets/loki/loki.db` → `loki-20220926.db` (default symlink) | Daniel's `biobin -D ~/group/datasets/loki/loki.db` resolves here via symlink |
| **plink 1.9** | `/appl/plink-1.90Beta6.18/plink` (explicit path) | Default `plink` on PATH is now **plink2** — use the full path to lock the 1.9 version for strict replication |
| **plink2** | `/appl/plink2-20240804/plink` (default on PATH as `plink`) | Most plink1.9 commands work in plink2 but pilot-test first |
| **Python** | `source venv/bin/activate` (Python 3.12) | Project venv — ⚠ missing `pandas scipy statsmodels numpy`, install before running scripts |
| **R** | `/appl/R-4.4/bin/Rscript` (explicit path) | ⚠ Default `R-3.6.3` is **broken** on this RHEL 9 node (`libreadline.so.6` missing); R-3.6 through R-4.3 all broken. **Use R-4.4 or R-4.5.** Daniel's scripts use base R, runs fine on 4.x |
| **regenie** | `module load regenie/3.2.1` (and older) | **Likely Charlene's "updated methods" tool** for the side-by-side comparison |
| **VEP / ANNOVAR** | `module load variant_effect_predictor/92` or `annovar/20191024` | If re-annotating from raw VCFs |

**Quick env setup** for a fresh shell on `superman`:
```bash
cd /project/hall/analysis/hearing-loss-genomics
source venv/bin/activate
module load rlsoftware/latest                                            # biobin + Ritchie Lab tools on PATH
export LD_LIBRARY_PATH=$PWD/analysis/daniel/configs/lib-shims:$LD_LIBRARY_PATH   # liblzma.so.0 shim — biobin needs it on RHEL 9
export PATH=/appl/plink-1.90Beta6.18:$PATH                               # lock plink 1.9 ahead of plink2
alias Rscript=/appl/R-4.4/bin/Rscript                                    # use working R 4.4 (default R-3.6.3 is broken)
```

**About the lib shim:** [`configs/lib-shims/liblzma.so.0`](configs/lib-shims/liblzma.so.0) → `/usr/lib64/liblzma.so.5`. `biobin` was compiled against XZ 5.0-era `liblzma.so.0`, but RHEL 9 ships only `liblzma.so.5`. The two versions are ABI-compatible enough for biobin's read-only VCF/bz2 usage — the symlink lets the dynamic loader find a working library. Verified: `biobin --help` runs cleanly with the shim. If a future biobin operation fails with an `lzma_*` error, this assumption breaks and we need to rebuild biobin against current libs.

**One-time Python install (already done 2026-05-12):**
```bash
pip install pandas scipy statsmodels numpy
```
Snapshot in [`requirements.txt`](../../requirements.txt) (110 packages, full venv).

**Python dependency audit:** Daniel's 102 scripts use only `sys` and `gzip` (stdlib). The pandas/scipy install is insurance, not strictly needed for the legacy pipeline.

**R dependency audit:** Daniel's R scripts use `library(MASS)`, `library(boot)`, `library(meta)`, `library(qqman)`. All four already installed in R-4.4's library — no additional installs needed.

## Path mapping: old `/project/PMBB/...` → new `/static/PMBB/...`

Daniel's runbook references PMBB v2 at `/project/PMBB/PMBB-Release-2020-2.0/...`. That path no longer exists — v2 has been moved to `/static/PMBB/PMBB-Release-2020-2.0/`. To replay Daniel's commands:

| Daniel's path | Use instead |
|---|---|
| `/project/PMBB/PMBB-Release-2020-2.0/` | `/static/PMBB/PMBB-Release-2020-2.0/` (or `data/pmbb_v2/`) |
| `/project/pmbb_all/PMBB-Release-2020-2.0/` | same → `/static/PMBB/PMBB-Release-2020-2.0/` |
| `/project/ritchie07/personal/daniel/HearingLoss/PMBB_Exome/` | `data/PMBB_Exome/` (his intermediates copied here) |

The cleanest pattern in our re-runs is to use [`data/pmbb_v2/`](../../data/pmbb_v2/) (the symlink) so paths stay short and portable. When porting to v3 later, this becomes [`data/pmbb_v3/`](../../data/pmbb_v3/).

## Open questions to resolve before starting

1. **"8 signal-driving cases" definition.** The runbook references multiple ZNF175 carrier counts:
   - chr19:51587727 → 8 carriers, 1 phecode-case (line 193-196)
   - chr19:51581437 → 90 carriers, 2 phecode-cases (line 197-200)
   - Joe's email list of 6 named individuals (lines 346-353)
   - Aggregate ZNF175 pLoF carriers via [`scripts/ZNF175_carrier.py`](scripts/pmbb_exome/ZNF175_carrier.py) on the multiallelic+stoploss variant set
   Confirm with Doug/Daniel which set is the "8 driving cases" for the deep-dive.

2. **plink version pilot.** Daniel used plink 1.9; default on LPC is now plink2. Run one Phase 3 chromosome under both 1.9 and 2.0 and compare outputs — confirm no divergence before committing the full replication to plink2.

3. **R session reproducibility.** The runbook does inline `d <- read.table(...)` in interactive R sessions. We should script these into `.R` files so the regression results are reproducible without manual REPL work.
