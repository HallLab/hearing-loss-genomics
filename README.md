# Hearing Loss Genomics — PMBB

Genomic analysis of age-related hearing loss in the **Penn Medicine BioBank (PMBB)**, extending Hui et al., PLOS Genetics 2023 ([doi:10.1371/journal.pgen.1010584](https://doi.org/10.1371/journal.pgen.1010584)).

> **First time here?**
> - **Replicating Chapter 1 (paper replication)?** Start with [`HANDOFF.md`](HANDOFF.md) — the Chapter 1 guide written for Elena (the interim person doing this work). Walks you through setup, day-by-day plan, and validation. Then [`REPRODUCTION_GUIDE.md`](REPRODUCTION_GUIDE.md) has the per-phase commands.
> - **Picking up Chapter 2 (ZNF175 follow-up)?** Start with [`STATUS_SNAPSHOT.md`](STATUS_SNAPSHOT.md) — captures the in-progress state of Chapter 2 work, blockers, pending tasks. Chapter 2 is not validated yet.

## Collaboration

- **Hall Lab** (Penn) — Molly Hall (PI), Nikki Palmiero, Andre Rico
- **Epstein Lab** (Penn) — Douglas Epstein

## Current focus

Deep-dive replication and follow-up of the **ZNF175** novel HL gene signal in PMBB. See `docs/analysis_plan.md` for the full task list.

## Repository structure

```
.
├── README.md           # This file
├── data/               # Raw PMBB data, reference papers (gitignored)
├── docs/               # Plans, meeting notes, methodology references
│   ├── analysis_plan.md
|   ├── meetings
│       └── kickoff_meeting_summary.md
|   ├── papers
│       └── paper_summary_hui2023.md
├── notebooks/          # Jupyter notebooks for exploration
├── results/            # Analysis outputs (selective gitignore)
└── venv/               # Python virtual env (gitignored)
```

## Environment

Run on **UPenn LPC** (`superman`), Red Hat Enterprise Linux 9.4.

```bash
source /project/hall/analysis/hearing-loss-genomics/venv/bin/activate
```

## Documentation

Detailed project docs live in `docs/`. Project-wide notes and writeups go to the Hall Lab Confluence:
**https://halllab.atlassian.net/wiki/spaces/HLP/overview?homepageId=718504536** (Andre - Notebook folder)

## Data handling

PMBB data is PHI-adjacent. **Never** commit anything from `data/`, paste participant IDs into chat, or copy raw data outside `/project/hall/`.
