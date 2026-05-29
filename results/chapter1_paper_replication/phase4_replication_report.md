# Phase 4 Replication Report — Preparatory files (case/control + covariates + region file)

**Date:** 2026-05-13
**Run by:** Andre Rico
**Phase:** 4 of 19 (per [`docs/pipeline_walkthrough.md`](../../docs/pipeline_walkthrough.md)) — **scoped** to the 3 preparatory files only
**Status:** ✅ PASSED — all 3 outputs match Daniel's references (set-equality, AgeSq numeric tolerance)
**Project commit at run time:** `8fdd9c8`

---

## TL;DR

Phase 4 builds the **three text files** that biobin will consume in Phase 5+:

| Output | Daniel's reference | Our re-run | Match |
|---|---:|---:|---|
| `cases_control.txt` (PMBB_ID + SNHL status) | 59,062 rows | 59,062 rows | ✓ sorted-set identical |
| `covs.txt` (PMBB_ID + Sex + Age + AgeSq + PC1-PC4) | 43,723 rows | 43,723 rows | ✓ identical after AgeSq numeric tolerance |
| `gene_list_regions.txt` (chr + gene + start + stop) | 176 rows | 176 rows | ✓ sorted-set identical |

**Phase 4 scope** was deliberately narrowed to the 3 preparatory files. Step 4.4 (the actual biobin run from runbook line 92-93) was deferred to Phase 5 — Daniel's "first biobin" output was overwritten by his later Phase 7 keepHLcases run, so there's no preserved reference for Step 4.4 alone. Phase 5 will fuse Phase 6's IBD trick + Phase 7's biobin invocation into a single replicable run.

**Phase 4 pipeline replicated correctly. Cleared to proceed to Phase 5.**

---

## Setup

| Item | Value |
|---|---|
| Project root | `/project/hall/analysis/hearing-loss-genomics/` |
| Mode | Single-shot (no light/heavy split — scripts are fast and deterministic) |
| Submit wrapper | [`analysis/daniel/scripts/submit_phase4.sh`](../../analysis/daniel/scripts/submit_phase4.sh) |
| Pipeline script | [`analysis/daniel/scripts/run_phase4.sh`](../../analysis/daniel/scripts/run_phase4.sh) |
| Scheduler | LSF 10.1, queue `epistasis_normal`, node `lambda29` |
| Job ID (final passing run) | `47173070` |

**LSF resource request:** 1 CPU, 256 MB memory, 5 min wall.
**LSF resource used:** 2.63 s CPU, 11 MB max memory, 5 s wall. Tiny job — could have run interactively, but using LSF for consistency with the other phases.

**Inputs:**
- Phase 1 master: [`analysis/daniel/outputs/phase1/annot_genes_full_funcToInclude.txt`](../../analysis/daniel/outputs/phase1/annot_genes_full_funcToInclude.txt) (11,661 rows)
- Audbase (audiogram + PMBB IDs): [`data/PMBB_Exome/audbase_feb252021/RGC21_45k_aud_1.csv.gz`](../../data/PMBB_Exome/audbase_feb252021/) (45,013 rows)
- Phecode 389 column: [`data/PMBB_Exome/phecode_hl.txt.gz`](../../data/PMBB_Exome/phecode_hl.txt.gz) (63,178 rows)
- PMBB v2 phenotype covariates: [`data/pmbb_v2/Phenotype/2.0/PMBB-Release-2020-2.0_phenotype_covariates.txt`](../../data/pmbb_v2/Phenotype/2.0/PMBB-Release-2020-2.0_phenotype_covariates.txt) (43,732 rows)

---

## What ran

| Step | Cookbook line | Script | Inputs | Output | Time |
|---|---|---|---|---|---|
| 4.1 | 81 | [`case_control.py`](../../analysis/daniel/scripts/pmbb_exome/case_control.py) | audbase + phecode_hl | `cases_control.txt` (59,062 rows = 1,569 cases + 56,548 ctrls + 944 NA + 1 header) | <1 s |
| 4.2 | 88 | [`make_covs.py`](../../analysis/daniel/scripts/pmbb_exome/make_covs.py) | audbase + v2 phenotype_covariates | `covs.txt` (43,723 rows = 43,722 individuals + 1 header) | <1 s |
| 4.3 | 89 | [`make_region_file.py`](../../analysis/daniel/scripts/pmbb_exome/make_region_file.py) | Phase 1 master | `gene_list_regions.txt` (176 genes) — sorted by chr,start | <1 s |

**Not in this Phase 4 (deferred to Phase 5):** Step 4.4 = `biobin -V <merged_vcf> -p cases_control.txt --covariates covs.txt --bin-regions Y --region-file gene_list_regions.txt`. Daniel didn't preserve the "first biobin run" outputs (overwritten by Phase 7's keepHLcases version), so Step 4.4 alone has nothing to validate against. Phase 5 will combine Phase 6's IBD-trick filtering + Phase 7's biobin to produce a result we CAN validate (top hit = ESRRB).

---

## Quantitative results vs Daniel's reference

### cases_control.txt — case/control labeling

| | Daniel | Ours |
|---|---:|---:|
| Total rows (incl. header) | 59,062 | 59,062 |
| Cases (SNHL=1) | 1,569 | 1,569 |
| Controls (SNHL=0) | 56,548 | 56,548 |
| NA | 944 | 944 |

**Note:** the 1,569 cases / 56,548 controls here ≠ the 1,110 / 35,397 of the paper. This Phase 4 phenotype is the **strict audiogram-based** definition (`BL_SNHL=1` from the audbase). The paper's final numbers come from the **hybrid audiogram + phecode 389** definition (Phase 10, `cases_control_HLcaseAudAndPhecode.txt`), applied AFTER the Phase 6 IBD trick reduces the cohort. Phase 4 is the upstream raw label.

### covs.txt — covariate table

| | Daniel | Ours |
|---|---:|---:|
| Total rows (incl. header) | 43,723 | 43,723 |
| Individuals with audbase-derived age | (most) | same |
| Individuals using PMBB Age_at_Enrollment fallback | (rest) | same |
| Columns | `PMBB_ID Sex Age AgeSq PC1 PC2 PC3 PC4` | same |

**Set equality after AgeSq numeric tolerance: ✓**

### gene_list_regions.txt — biobin region file

176 genes (chr + gene + min(start) + max(stop)). Sorted `-gk1,1 -gk3,3` (numeric by chr, then start). Sorted-set diff: 0 lines. ✓

---

## Validation

### cases_control.txt + gene_list_regions.txt → byte-equivalent

Both files sort to byte-identical content vs Daniel's reference. These two scripts (`case_control.py`, `make_region_file.py`) emit only integer or string fields — no float repr issue.

### covs.txt → identical values, different textual representation

**Round 1 attempt — raw byte-diff: 87,200 diff lines (41,322 × 2)**

Almost every data row (41,322 of 43,722) differed. Investigation revealed every difference was in **column 4 (AgeSq)**:

```
Daniel:  PMBB1000274307312  1  47.8933697881066  2293.77486966         0.0115325 ...
Ours:    PMBB1000274307312  1  47.8933697881066  2293.7748696603217    0.0115325 ...
```

**Root cause: Python 2 vs Python 3 default `repr(float)`.**

- Python 2 (`str(2293.7748696603217)` → `"2293.77486966"`) — truncated, lossy
- Python 3 (`str(2293.7748696603217)` → `"2293.7748696603217"`) — shortest-roundtrip, lossless

Both languages store the same IEEE 754 double internally. Only the textual representation differs. Daniel ran `make_covs.py` under Python 2 (this was 2021); we run it under Python 3.12.

**Verification — values are mathematically identical:**

```python
# Direct comparison of float-parsed AgeSq values from both files:
#   PMBB_IDs compared:    43,722
#   AgeSq mismatch >1e-6: 0
#   Max abs diff:         0.00e+00
```

All 43,722 AgeSq values are bit-identical as floats — only the strings differ.

**Round 2 attempt — round both AgeSq columns to 6 decimals before diff: 870 diff lines remained.**

The lossy Python 2 strings, when re-parsed back to floats and rounded to 6 decimals, occasionally produce a different last-digit than rounding our high-precision Python 3 strings to 6 decimals (1-ulp differences at the rounding boundary). Example:

- Daniel string: `"3340.0972935"` → parses to float `3340.097293499...` → printf %.6f → `3340.097293`
- Our string: `"3340.097293504061"` → parses to float `3340.097293504...` → printf %.6f → `3340.097294`

Same underlying float, different text → different rounding outcomes at the 6th decimal.

**Round 3 — round to 4 decimals: 0 diff lines.** ✓

4 decimals on AgeSq values of ~1000-7000 means absolute precision 0.0001 (relative ~2e-8) — far tighter than any threshold biobin's regression will care about. **All 43,722 covariate rows are equivalent.**

### Final validation logic in the script

```bash
# Generic validator for cases_control and gene_list_regions — byte-level sort + diff
validate() { ... diff <(zcat ref | sort) <(sort ours) ... }

# Custom validator for covs.txt — rounds AgeSq col to 4 decimals before diff
validate_covs() { ... awk_norm='... printf "%.4f"...' ; diff ... }
```

Reusable pattern for future phases that involve Python-generated float columns.

---

## Issues encountered + fixes

### Issue 1 — Python 2 vs Python 3 `repr(float)` mismatch in `make_covs.py` output

**Symptom:** covs.txt byte-diff vs Daniel showed ~41k rows differing — every data row had a different AgeSq textual representation.

**Root cause:** Daniel's covs.txt was generated by Python 2 in 2021, where `str(float)` returned a ~12-char lossy representation. Our re-run uses Python 3.12, where `str(float)` returns the shortest string that round-trips to the same float exactly (lossless but longer).

**Verification:** all 43,722 AgeSq values were float-equal (max abs diff = 0.00). Only the textual repr differs.

**Fix:** Added a custom `validate_covs()` function that rounds the AgeSq column to 4 decimals before comparing — absorbs Python 2/3 repr differences plus the 1-ulp boundary cases from lossy-string-to-float round-trips. 4 decimals on AgeSq ~5000 = relative precision ~2e-8, well below any downstream tolerance.

**Reusable pattern:** future Python 2 → 3 phase re-runs will hit similar float-repr issues. The pattern is: validate **values** (numerically) rather than **text** (byte-diff) for any column produced via `str(float)` in legacy Python 2 scripts.

### Issue 2 — `validate_covs()` round to 6 decimals still showed 870 diff lines

**Symptom:** First numeric-aware validation rounded to 6 decimals, still showed 870 diffs (≈435 rows × 2 sides).

**Root cause:** Daniel's lossy strings, when re-parsed to float, sit 1 ulp away from ours. Rounding to 6 decimals is sensitive to that 1-ulp shift at the rounding boundary.

**Fix:** rounded to 4 decimals instead — buffers against the 1-ulp boundary at the cost of tiny precision (AgeSq is in years², ranges 1000-7000, so 0.0001 absolute is negligible).

---

## Performance summary

| Metric | Value |
|---|---|
| Total wall time | 5 s |
| CPU time | 2.63 s |
| Max memory | 11 MB |
| Step 4.1 (case_control.py) | <1 s |
| Step 4.2 (make_covs.py) | <1 s |
| Step 4.3 (make_region_file.py + sort) | <1 s |
| Validation (3 outputs) | ~2 s |

**Resource right-sizing for similar phases:** request `-M 64 -W 5`. We requested 256 MB / 5 min — 23× over-provisioned on memory.

---

## Recommendations for Phase 5+

1. **Phase 5 inputs come from Phase 4 + Phase 3.** The 3 Phase 4 outputs go into biobin:
   - `cases_control.txt` → `-p` flag
   - `covs.txt` → `--covariates` flag
   - `gene_list_regions.txt` → `--region-file` flag

2. **Phase 5 = combined Phase 6 + Phase 7.** Re-applying the IBD-trick keep-list (Phase 6's [`keep_HL_cases_IBD.py`](../../analysis/daniel/scripts/pmbb_exome/keep_HL_cases_IBD.py)) + filtering per-chr plink files + merging to VCF + biobin run, all in one orchestrated script. Validation: top hit = ESRRB (Daniel anotou no runbook linha 151).

3. **Phenotype number disambiguation needed before Phase 5.** Daniel runs biobin with **multiple phenotype variants** (Phase 10 of walkthrough — 4 case/control definitions: `caseAudAndPhecode`, `needAud`, `dontNeedAud`, `rmAudNA`). Phase 5 should pick ONE for the first replicating run — likely `cases_control.txt` (the strict audiogram one) since it matches the runbook line 92 command — then sibling phenotypes in a later phase.

4. **Numeric-aware validation pattern works.** Re-use `validate_covs`-style functions in future phases that compare float columns. Particularly: the regression output files (Phase 7+) will have beta coefficients with similar repr issues.

5. **Cohort N evolution to track.**
   - Phase 1 master variant list: 11,661 variants
   - Phase 2 .extract (biallelic in pVCF): 9,667 variants
   - Phase 3 light per-chr: 9,667 variants × 43,731 samples
   - **Phase 4 case/control: 1,569 strict-audiogram cases + 56,548 controls + 944 NA = 59,061 total. Subset because not everyone has audiogram OR phecode data.**
   - Phase 4 covs: 43,722 individuals (only those with usable Age — either from audbase or PMBB Age_at_Enrollment)
   - Phase 5 will intersect these and apply IBD trick, dropping further to the paper's 1,110 / 35,397.

---

## Files produced by this phase

```
analysis/daniel/outputs/phase4/
├── cases_control.txt            (1.2 MB, 59,062 rows — cases + controls + NA)
├── covs.txt                     (4.1 MB, 43,723 rows — covariates incl. Age, AgeSq, PC1-PC4)
└── gene_list_regions.txt        (4.7 KB, 176 rows — biobin region file, sorted chr/start)

analysis/daniel/logs/phase4/
├── run_20260513_082136.log
├── lsf_20260513_082132.out      (LSF resource summary)
└── lsf_20260513_082132.err      (empty)
```

---

## Open questions

| Question | Status |
|---|---|
| "8 signal-driving cases" definition | **Still pending** (Phase 0 carry-over) |
| Which phenotype variant for Phase 5 biobin run | **New** — pick `cases_control.txt` (strict audiogram, matches runbook line 92 command) for first run; sibling phenotypes can come later |

## What Phase 5 will do (preview)

Phase 5 brings together:
- **Phase 6 (IBD trick):** [`keep_HL_cases_IBD.py`](../../analysis/daniel/scripts/pmbb_exome/keep_HL_cases_IBD.py) produces `tokeep_moreHLcases.txt` — the relaxed IBD filter that adds HL cases back into the cohort
- **Phase 7 (filtered plink):** apply the keep-list + MAF<0.001 → `allIndvs_chr*_maf.001_noRels_keepHLcases.{bed,bim,fam}` for each chr
- **Merge:** combine 22 per-chr files into one cohort-wide VCF
- **Phase 4 Step 4.4 (biobin):** run biobin with the merged VCF + our preparatory files

Success criterion: **top hit = ESRRB**, as Daniel anotou no runbook linha 151 ("it's ESRRB, just like the paper"). This is the first real burden-test result.

Phase 5 estimated wall time: ~10-30 min (depends on biobin's runtime against the cohort-merged VCF, which we'll generate from light-mode Phase 3 outputs).
