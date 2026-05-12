# `analysis/daniel/` — Hui et al. 2023 pipeline replication (PMBB v2)

Workspace for **re-executing Daniel Hui's original pipeline** on PMBB v2 data, to confirm we reproduce the published numbers from [Hui et al. 2023, PLOS Genetics](../../docs/papers/pgen.1010584.pdf) (doi:10.1371/journal.pgen.1010584).

This is one of two parallel analysis tracks (see `../README.md` and [`docs/CLAUDE.md`](../../CLAUDE.md)):

- **Track 1 (here):** Daniel's v2 pipeline → confirm replicability of published results
- **Track 2:** [`analysis/01_phase1_exploration_pmbb_v3/`](../01_phase1_exploration_pmbb_v3/) — Andre re-derives on PMBB v3 (`PMBB-Release-2024-3.0`)

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

Identified from the runbook:

| Tool | Used for | Notes |
|---|---|---|
| `plink` (v1.9) | genotype extraction, MAF/IBD filtering, format conversion | older plink — check version compatibility |
| `biobin` | gene-burden binning + logistic regression in one step | `biobin -D loki.db ...`; needs `loki.db` (Ritchie Lab's gene/region DB at `~/group/datasets/loki/loki.db`) |
| Python 3 | scripts in `scripts/*/` | uses stdlib mostly; check `sys.argv` patterns |
| R | regression (`run_*.R`), meta-analysis (`meta.R`), figures (`make_figs*.R`) | scripts read tab/comma files; check for required packages |
| LSF (`bsub`) | job submission | LPC may use SLURM now — will need translation |
| `bcftools`, `vcftools` | VCF manipulation | seen in some sections |

**Action item:** verify our `venv/` has Python 3 + R bindings + any required pip packages; check `loki.db` is accessible from this account; decide LSF→SLURM translation strategy.

## Open questions to resolve before starting

1. **"8 signal-driving cases" definition.** The runbook references multiple ZNF175 carrier counts:
   - chr19:51587727 → 8 carriers, 1 phecode-case (line 193-196)
   - chr19:51581437 → 90 carriers, 2 phecode-cases (line 197-200)
   - Joe's email list of 6 named individuals (lines 346-353)
   - Aggregate ZNF175 pLoF carriers via [`scripts/ZNF175_carrier.py`](scripts/pmbb_exome/ZNF175_carrier.py) on the multiallelic+stoploss variant set
   Confirm with Doug/Daniel which set is the "8 driving cases" for the deep-dive.

2. **biobin availability.** The Ritchie Lab's `biobin` tool with `loki.db` is the core of every burden test. Verify it still runs on LPC.

3. **PMBB v2 file path stability.** All `bsub` commands reference `/project/PMBB/PMBB-Release-2020-2.0/...`. Confirm these paths still resolve. If v2 has been archived, we need an alternate path.

4. **LSF → SLURM.** Daniel used `bsub`. Confirm LPC's current scheduler.

5. **R session reproducibility.** The runbook does inline `d <- read.table(...)` in interactive R sessions. We should script these into `.R` files so the regression results are reproducible without manual REPL work.
