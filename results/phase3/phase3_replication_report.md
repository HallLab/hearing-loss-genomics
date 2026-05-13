# Phase 3 Replication Report — plink genotype extraction

**Date:** 2026-05-13
**Run by:** Andre Rico
**Phase:** 3 of 19 (per [`docs/pipeline_walkthrough.md`](../../docs/pipeline_walkthrough.md))
**Status:** ✅ Light mode PASSED (all 22 chrs validated); ✅ chr21 heavy pilot showed plink 1.9 reproduces Daniel **byte-for-byte**; ⚠ plink 2.0 OOM-killed by LSF (cosmetic — `--memory` flag needed for cluster use)
**Project commit at run time:** `8c66825`

---

## TL;DR

Phase 3 introduced **per-individual genotype data** into the pipeline for the first time. We ran two modes simultaneously:

1. **Light mode** (all 22 chrs): decompressed Daniel's per-chr plink files into our workspace. All counts and consistency checks PASSED.
2. **Heavy pilot** (chr21 only): re-ran `plink --vcf … --extract … --make-bed` ourselves on the raw 8.7 GB chr21 pVCF, with both **plink 1.9** and **plink 2.0** for a side-by-side validation.

**Headline finding:** plink 1.9 (v1.90b6.18) reproduces Daniel's chr21 outputs **byte-for-byte** — `.bim` and `.fam` diff = 0 lines, `.bed` md5sum identical. plink 1.9 is the version to use for all future heavy plink work.

plink 2.0 was killed by LSF after ~57 s for exceeding the 4 GB memory limit — plink 2.0's default is to reserve ~50% of node RAM (515 GB on a 1 TB node) for its workspace, regardless of how few variants it'll actually process. Not a bug, just a default that needs the `--memory 2048` flag in cluster environments. Comparison v1.9 vs v2.0 deferred (not needed since v1.9 = Daniel exactly).

**Phase 3 pipeline replicated correctly. Cleared to proceed to Phase 4 — but using plink 1.9 from now on for any heavy plink work.**

---

## Setup

| Item | Value |
|---|---|
| Project root | `/project/hall/analysis/hearing-loss-genomics/` |
| Mode | **Light** for 22 chrs + **Heavy pilot** on chr21 |
| Submit wrapper | [`analysis/daniel/scripts/submit_phase3.sh`](../../analysis/daniel/scripts/submit_phase3.sh) |
| Pipeline script | [`analysis/daniel/scripts/run_phase3.sh`](../../analysis/daniel/scripts/run_phase3.sh) |
| Scheduler | LSF 10.1, queue `epistasis_normal`, node `krypton04` |
| Job ID | `47172811` |
| plink 1.9 binary | `/appl/plink-1.90Beta6.18/plink` — `PLINK v1.90b6.18 64-bit (16 Jun 2020)` |
| plink 2.0 binary | `/appl/plink2-20240804/plink` — `PLINK v2.00a6LM 64-bit Intel (4 Aug 2024)` |
| Phase 2 input | [`analysis/daniel/outputs/phase2/matched_snp_IDs_annot_pVCF_noNA_noMultiallelic.extract`](../../analysis/daniel/outputs/phase2/matched_snp_IDs_annot_pVCF_noNA_noMultiallelic.extract) (9,667 IDs) |

**LSF resource request:** 1 CPU, 4 GB memory, 60 min wall.
**LSF resource used:** 358.81 s CPU, 4 GB memory (= hit ceiling exactly when plink 2.0 was OOM-killed), 362 s wall. Job exited with status EXIT due to the OOM kill at the end of the pilot — but light mode and plink 1.9 completed cleanly within the 4 GB envelope before plink 2.0 was attempted.

---

## Part A — Light mode (all 22 chrs)

### What ran

Single step: decompress Daniel's `data/PMBB_Exome/genotypes/allIndvs_chr{1..22}.{bed,bim,fam}.gz` (66 files) into [`analysis/daniel/outputs/phase3/light/`](../../analysis/daniel/outputs/phase3/light/).

Daniel keeps these as `.gz` (PMBB intermediates are gzipped in his archive). plink reads `.bed/.bim/.fam` directly — so we just decompress and they're ready for Phase 4's biobin runs.

### Performance

- **Wall time:** 2 s for all 66 files
- **Output size:** 153 MB decompressed (down from ~30 MB compressed)

### Validation (light mode)

Three integrity checks, all PASSED ✓:

#### (a) Per-chr variant counts sum to Phase 2 .extract count

| chr | variants | chr | variants | chr | variants | chr | variants |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 932 | 7 | 468 | 13 | 255 | 19 | 328 |
| 2 | 881 | 8 | 224 | 14 | 116 | 20 | 38 |
| 3 | 570 | 9 | 236 | 15 | 152 | 21 | 173 |
| 4 | 358 | 10 | 628 | 16 | 229 | 22 | 258 |
| 5 | 729 | 11 | 1065 | 17 | 487 | | |
| 6 | 717 | 12 | 689 | 18 | 134 | | |

**Total: 9,667 variants across 22 chromosomes ≡ Phase 2 .extract count.** ✓

#### (b) Sample count consistent across chromosomes

All 22 `.fam` files report **43,731 samples** (identical sets). ✓

This is the PMBB v2 cohort size after IBD filtering at the pVCF level — same N that Daniel's runbook flagged at lines 65-66 (`N=43731 in chr22 pVCF`).

#### (c) Variant subset check

All 9,667 unique variant IDs across the 22 `.bim` files are present in our Phase 2 `.extract`. ✓

### Output (ready for Phase 4)

```
analysis/daniel/outputs/phase3/light/
├── allIndvs_chr1.{bed,bim,fam}
├── allIndvs_chr2.{bed,bim,fam}
├── ...
├── allIndvs_chr21.{bed,bim,fam}   ← 173 variants, 43,731 samples
└── allIndvs_chr22.{bed,bim,fam}
```

These 66 files are the canonical Step 3.1 outputs that Phases 4-8+ will consume.

---

## Part B — Heavy pilot on chr21 (plink 1.9 vs plink 2.0)

### Goal

Two questions at once:
1. **Can plink actually run on this LPC against the v2 pVCFs?** (Smoke test — we hadn't invoked plink in this environment yet.)
2. **Do plink 1.9 and plink 2.0 produce the same output?** (Open question from Phase 1 report — see [`results/phase1/phase1_replication_report.md`](../phase1/phase1_replication_report.md) recommendation 5.)

Picked chr21 because it's the smallest autosomal pVCF (8.7 GB). Daniel's reference for chr21 has 173 variants in 43,731 samples.

### What ran

```bash
# plink 1.9
/appl/plink-1.90Beta6.18/plink \
    --vcf data/pmbb_v2/Exome/pVCF/GL_by_chrom/PMBB-Release-2020-2.0_genetic_exome_chr21_GL.vcf.gz \
    --vcf-half-call m \
    --extract analysis/daniel/outputs/phase2/matched_snp_IDs_annot_pVCF_noNA_noMultiallelic.extract \
    --make-bed \
    --out analysis/daniel/outputs/phase3/pilot_chr21/v19/allIndvs_chr21

# plink 2.0 — same args
/appl/plink2-20240804/plink \
    --vcf data/pmbb_v2/Exome/pVCF/GL_by_chrom/PMBB-Release-2020-2.0_genetic_exome_chr21_GL.vcf.gz \
    --vcf-half-call m \
    --extract analysis/daniel/outputs/phase2/matched_snp_IDs_annot_pVCF_noNA_noMultiallelic.extract \
    --make-bed \
    --out analysis/daniel/outputs/phase3/pilot_chr21/v20/allIndvs_chr21
```

### Result B.1 — plink 1.9: ✅ byte-for-byte match to Daniel

| Metric | Daniel's reference | Our plink 1.9 | Match |
|---|---|---|---|
| `.bim` (variants) | 173 rows | 173 rows | ✓ diff 0 lines |
| `.fam` (samples) | 43,731 rows | 43,731 rows | ✓ diff 0 lines |
| `.bed` md5sum | `9a0ccc62fc3fdb10cd48d9e9061ad7a6` | `9a0ccc62fc3fdb10cd48d9e9061ad7a6` | ✓ identical |

**plink 1.9 reproduces Daniel's chr21 output exactly.** This validates our pipeline + environment for heavy plink work. Wall time: 301 s (5 min) to read the 8.7 GB pVCF and extract 173 variants.

### Result B.2 — plink 2.0: ❌ killed by LSF (OOM)

`TERM_MEMLIMIT: job killed after reaching LSF memory usage limit. Exited with exit code 1.`

**Root cause** (from plink 2.0's stdout before it died):
```
1031007 MiB RAM detected, ~1020756 available; reserving 515503 MiB for main workspace.
```

plink 2.0 detected ~1 TB of RAM on the node and tried to reserve ~515 GB for its own workspace, even though it would only be processing a small subset of variants. This is plink 2.0's default behavior (designed for biobank-scale workloads on large nodes); it doesn't know about LSF's memory limit. LSF killed it after 57 s when it tried to allocate beyond our 4 GB request.

**Not a bug** — plink 2.0 needs the `--memory <MB>` flag to behave in cluster environments:
```bash
plink2 ... --memory 2048    # cap workspace at 2 GB
```

**Status of the v1.9 vs v2.0 comparison:** **deferred.** Since plink 1.9 reproduces Daniel byte-for-byte, the comparison isn't load-bearing for this project. If we ever need plink 2.0 (e.g., for a feature not available in 1.9), we'll add `--memory 2048` and re-run. For now: **use plink 1.9 for all heavy plink work in Phase 3+.**

### Performance — pilot chr21

| Step | Wall time |
|---|---|
| plink 1.9 reads 8.7 GB pVCF + extracts 173 vars + writes bed/bim/fam | **301 s (5 min)** |
| plink 2.0 (killed at ~57 s) | n/a |

Extrapolation for heavy mode on all 22 chrs (sequential): chr21 was the smallest pVCF; chr1 at 79 GB would be ~9× slower (~45 min). Heavy mode total (sequential) ≈ **~3-6 hours**. With 22-task LSF job array, wall time would be ~max single-chr time ≈ **~45 min**.

---

## Comparison vs Daniel's reference (chr21 specifics)

| File | Daniel (compressed) | Ours light mode | Ours plink 1.9 | Match |
|---|---:|---:|---:|---|
| `allIndvs_chr21.bim` | 1.4 KB (.gz) | identical (decompressed) | 6.1 KB | ✓ diff 0 |
| `allIndvs_chr21.fam` | 485 KB (.gz) | identical (decompressed) | ~700 KB | ✓ diff 0 |
| `allIndvs_chr21.bed` | 14 KB (.gz) | identical (decompressed) | ~937 KB | ✓ md5 match |

---

## Issues encountered + fixes

### Issue 1 — plink 2.0 OOM-killed by LSF

**Symptom:** Job exited at status EXIT after 6 minutes; `TERM_MEMLIMIT` in LSF summary. Light mode + plink 1.9 had already completed cleanly.

**Root cause:** plink 2.0 by default reserves ~50% of detected RAM for its main workspace, regardless of actual workload size. On a node with 1 TB RAM, that's 515 GB — far exceeding our 4 GB LSF request.

**Fix (for any future plink 2.0 use):** add `--memory 2048` (or whatever cap matches the LSF request, in MB). plink 1.9 doesn't have this issue (more conservative default), so for now we just use plink 1.9.

**Status:** documented in [`docs/pipeline_walkthrough.md`](../../docs/pipeline_walkthrough.md) Phase 3 section and in memory (`project_lpc_environment.md`). The Phase 3 script remains unchanged — the comparison isn't blocking and we don't need plink 2.0 for the replication.

---

## Performance summary

| Metric | Value |
|---|---|
| Total wall time | 362 s (~6 min) |
| CPU time | 358.81 s |
| Max memory | 4096 MB (hit ceiling at plink 2.0 attempt) |
| Light mode (66 files decompressed) | 2 s |
| plink 1.9 on chr21 (8.7 GB pVCF) | 301 s = 5 min |
| plink 2.0 on chr21 (killed) | ~57 s of 5+ min that 2.0 would have taken |

**For future heavy plink runs:** request `-M 4096 -W 30` per chr (or per task in a job array). Memory ceiling could likely go to 2 GB safely with plink 1.9 — chr21 used <4 GB even with a lot of slack. Chr1 (79 GB pVCF) may want a bit more — recommend `-M 4096 -W 60` for the largest chrs.

---

## Decision: use plink 1.9 for all heavy plink work

Confirmed by byte-for-byte match against Daniel's reference. **Path:** `/appl/plink-1.90Beta6.18/plink`. This was set up in our shell env from Phase 1 — see [`analysis/daniel/README.md`](../../analysis/daniel/README.md) "Quick env setup".

Recorded in memory (`project_lpc_environment.md`) for future sessions.

---

## Recommendations for Phase 4+

1. **Phase 4+ inputs come from `analysis/daniel/outputs/phase3/light/`** — the 22 per-chr `allIndvs_chr*.{bed,bim,fam}` files. Note these are pre-IBD-filter (43,731 samples each, no MAF cap). Phase 4 will need Phase 6's `keepHLcases` filter applied — or use Daniel's already-Phase-6-applied files at `data/PMBB_Exome/genotypes/allIndvs_chr*_maf.001_noRels_keepHLcases.{bed,bim,fam}.gz`.

2. **When Phase 6 / heavy runs are needed:** use plink 1.9 (`/appl/plink-1.90Beta6.18/plink`) explicitly. The default `plink` on PATH is plink 2.0 and will misbehave on LSF without `--memory`.

3. **Heavy-mode planning for v3 port:** when we eventually replicate Phase 3 against v3 pVCFs, use a 22-task LSF job array (`bsub -J "name[1-22]"`) with each task processing one chr — total wall ~45 min vs ~6 hours sequential.

4. **plink 2.0 still worth re-testing once** (low priority): re-run the same command with `--memory 2048` to confirm it produces a bed file matching plink 1.9. Not blocking, but if it differs we'd want to know why before any future "modernization" pass.

5. **biobin (Phase 4) is the next plink-adjacent tool** — already validated environment-wise (shim for `liblzma.so.0`, see `analysis/daniel/configs/lib-shims/`). It reads VCFs, not bed/bim/fam, so we'll feed it the VCF version of the per-chr files (Phase 3 Step 3.4 in Daniel's runbook produces these — but Daniel didn't preserve them, only the keepHLcases versions).

---

## Files produced by this phase

```
analysis/daniel/outputs/phase3/
├── light/                                       (153 MB, ready for Phase 4+)
│   ├── allIndvs_chr1.{bed,bim,fam}              ← 932 variants
│   ├── ...
│   ├── allIndvs_chr19.{bed,bim,fam}             ← 328 variants (ZNF175 region)
│   ├── allIndvs_chr21.{bed,bim,fam}             ← 173 variants
│   └── allIndvs_chr22.{bed,bim,fam}
└── pilot_chr21/
    ├── v19/                                     (6.3 MB — plink 1.9 chr21 outputs)
    │   ├── allIndvs_chr21.{bed,bim,fam}         ← byte-identical to Daniel
    │   ├── allIndvs_chr21.log
    │   └── plink_{stdout,stderr}.txt
    └── v20/                                     (1.3 MB — plink 2.0 partial before OOM)
        ├── allIndvs_chr21.log                   ← shows the 515 GB reservation attempt
        └── plink_{stdout,stderr}.txt

analysis/daniel/logs/phase3/
├── run_20260513_075014.log                      (pipeline log)
├── lsf_20260513_075013.out                      (LSF stdout — shows TERM_MEMLIMIT)
└── lsf_20260513_075013.err                      (empty)
```

---

## Open questions

| Question | Status |
|---|---|
| plink 1.9 vs plink 2.0 byte equivalence | **Deferred** — plink 1.9 = Daniel exactly, so 2.0 comparison isn't blocking |
| "8 signal-driving cases" definition | **Still pending** — needs Daniel/Doug/Molly (Phase 0 carry-over) |
| Phase 4+ heavy mode strategy | **Resolved** — use plink 1.9, job array per chr when needed |

## What Phase 4 will do (preview)

Phase 4 takes the per-chr genotype files (currently in `outputs/phase3/light/`) and runs the **first end-to-end biobin burden test** of the HL gene panel. New ingredients enter:
- **Phenotype:** `cases_control.txt` from `data/PMBB_Exome/` (case = audiogram BL_SNHL=TRUE, control = audiogram or phecode FALSE)
- **Covariates:** age, sex, PCs, ancestry — `covs.txt` and variants
- **biobin** invocation: `biobin -D loki.db -V <merged_vcf> -p cases_control.txt --covariates covs.txt --bin-regions Y --region-file gene_list_regions.txt -G 38 --test logistic`

Expected output: a per-gene burden test result. Top hit should be **ESRRB** (Daniel notes this at runbook line 151: *"it's ESRRB, just like the paper"*). That's our Phase 4 success criterion.
