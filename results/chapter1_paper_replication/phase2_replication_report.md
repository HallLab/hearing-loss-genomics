# Phase 2 Replication Report — SNP ID reconciliation

**Date:** 2026-05-12
**Run by:** Andre Rico
**Phase:** 2 of 19 (per [`docs/pipeline_walkthrough.md`](../../docs/pipeline_walkthrough.md))
**Status:** ✅ PASSED — both validations green, sorted unique extract is bit-identical to Daniel's (md5sum match)
**Project commit at run time:** `52f72c9`

---

## TL;DR

Phase 2 reconciled the master variant list from Phase 1 (11,661 variants, annotation-IDs) against the PMBB v2 pVCF SNP IDs, producing the **9,667-ID `.extract` file** that plink will consume in Phase 3 to pull genotypes.

| Validation target | Daniel's reference | Our re-run | Match |
|---|---:|---:|---|
| `matched_clean` VCF_ID set (after NA + multi-allelic filter) | 9,668 rows | 9,668 rows | ✓ identical |
| `.extract` line set (final plink input) | 9,667 lines | 9,667 lines | ✓ identical (md5sum-equal when sorted unique) |

Ran in **light mode** — skipped Step 2.1 (the `zgrep` over 781 GB of pVCFs) by reusing Daniel's pre-extracted `vcf_SNP_IDs_allchr.txt.gz`. This was safe: that file is a pure mechanical extract (`zgrep -v '#' | cut -f 1-5`) of the same pVCFs we still have access to today, with no schema-level transformation. Total wall time: **20 s** vs ~1 h that heavy mode would have cost.

**Phase 2 pipeline replicated correctly. Cleared to proceed to Phase 3.**

---

## Setup

| Item | Value |
|---|---|
| Project root | `/project/hall/analysis/hearing-loss-genomics/` |
| Mode | **Light** — Step 2.1 reused Daniel's pre-extracted pVCF IDs |
| Submit wrapper | [`analysis/daniel/scripts/submit_phase2.sh`](../../analysis/daniel/scripts/submit_phase2.sh) |
| Pipeline script | [`analysis/daniel/scripts/run_phase2.sh`](../../analysis/daniel/scripts/run_phase2.sh) |
| Scheduler | LSF 10.1, queue `epistasis_normal`, node `krypton11` |
| Job ID (successful run) | `47167714` |

**LSF resource request:** 1 CPU, 1 GB memory, 15 min wall.
**LSF resource actually used:** 10.01 s CPU, 8 MB max memory, 20 s wall, 12 s turnaround. The 1 GB memory request was ~128× over-provisioned; even 64 MB would have been plenty.

**Why light mode:** heavy mode would re-execute Daniel's parallel 22-chromosome `zgrep -v '#' | cut -f 1-5` over the raw pVCFs at `data/pmbb_v2/Exome/pVCF/GL_by_chrom/`. Those files total **781 GB compressed** (e.g. chr1 = 79 GB, chr19 = 56 GB, chr21 = 8.7 GB) — even with 22 parallel LSF tasks, we'd be paying for ~10-20 min of cluster I/O for output we can validate against Daniel's existing files directly. Light mode reuses Daniel's pre-extracted `vcf_SNP_IDs_allchr.txt.gz` (77 MB compressed → 369 MB decompressed) which is a pure mechanical extract of the same pVCFs. No schema risk because cols 1-5 of a VCF (CHROM, POS, ID, REF, ALT) are part of the VCF spec, not Daniel-specific.

---

## What ran

| Step | Cookbook line | Action | Time |
|---|---|---|---|
| 2.1 (light) | — (Daniel's offline cut -f 1-5 on 22 pVCFs) | `zcat data/PMBB_Exome/vcf_SNP_IDs/vcf_SNP_IDs_allchr.txt.gz` → 369 MB, 11,337,300 rows | 3 s |
| 2.3 | 53 | [`annot_IDs_vs_pVCF.py`](../../analysis/daniel/scripts/pmbb_exome/annot_IDs_vs_pVCF.py): join master list (11,661 rows) with pVCF IDs (11.3 M rows) by `(chr, pos)` | 7 s |
| 2.4 | 55 | `grep -v NA \| grep -v ";"` — drop rows not in pVCF + drop multi-allelic | <1 s |
| 2.5 | 57-58 | Sanity check: `awk '$4 != $7'` (VCF ref vs annot ref) — expect 0 mismatches | <1 s |
| 2.6 | 60 | `tail -n +2 \| cut -f 3` — strip header, extract VCF_ID column | <1 s |

**Step 2.2 (`cat per-chr → allchr`) was implicit** — Daniel's `vcf_SNP_IDs_allchr.txt.gz` IS already the concatenation, so we use it directly.

---

## Quantitative results — variant attrition through the phase

```
Phase 1 master list:     11,661 variants (in 176 HL genes, pLoF + missense REVEL>0.6)
                              │
                              │ join by (chr, pos) — dedupes 425 multi-allelic-in-annotation entries
                              ▼
matched_snp_IDs_annot_pVCF.txt:   11,236 unique chr:pos positions
                              │
                              │ drop 395 NA (in annotation but not in pVCF — not sequenced)
                              │ drop 1,173 multi-allelic (";" in VCF_ID — plink can't handle)
                              ▼
matched_clean.txt:               9,668 rows (incl. header)
                              │
                              │ tail -n +2 | cut -f 3
                              ▼
.extract (plink input):          9,667 IDs ← THIS feeds Phase 3
```

**Total attrition Phase 1 → Phase 3 input:** 11,661 → 9,667 = **17% loss**. Breakdown:
- **3.6%** lost to multi-allelic-in-annotation dedup (same chr:pos with multiple ALT alleles collapsed) — 425 vars
- **3.4%** lost because annotation has the variant but the pVCF doesn't (not sequenced / failed QC) — 395 vars
- **10.1%** lost because the pVCF reports them as multi-allelic — 1,173 vars

**Daniel's note on multi-allelic loss:** these are added back later in Phase 9 (`addBack_multiallelic_stoploss/` workspace — see runbook line 261 onward). Phase 3 + downstream initially proceed on the 9,667 biallelic-only subset, then a parallel "+multi-allelic" arm runs separately for sensitivity.

## Comparison vs Daniel's reference

| File | Daniel (compressed) | Ours (uncompressed) | Row count | Match |
|---|---:|---:|---:|---|
| `vcf_SNP_IDs_allchr.txt` | 77 MB (.gz) | 369 MB | 11,337,300 | identical (we just decompressed his) |
| `matched_snp_IDs_annot_pVCF.txt` | 133 KB (.gz) | 595 KB | 11,236 | ✓ same variant IDs |
| `matched_..._noNA_noMultiallelic.txt` | 104 KB (.gz) | 488 KB | 9,668 (with header) | ✓ same VCF_ID set |
| `matched_..._noNA_noMultiallelic.extract` | 40 KB (.gz) | 150 KB | 9,667 (no header) | ✓ **md5sum-identical** when sorted unique |

The `.extract` files are bit-identical (`md5sum 5e80ebc0faa5e68277cfeb948af8b1da` on both, after sort + uniq), confirming the pipeline produces the same plink input that Daniel had.

---

## Validation

We used **semantic set-equality** rather than byte-for-byte diff (same pattern established in [`results/phase1/phase1_replication_report.md`](../phase1/phase1_replication_report.md)).

### Verification commands

```bash
# (a) matched_clean VCF_ID set equivalence
diff <(zcat data/PMBB_Exome/matched_snp_IDs_annot_pVCF_noNA_noMultiallelic.txt.gz \
        | tail -n +2 | cut -f3 | sort -u) \
     <(tail -n +2 analysis/daniel/outputs/phase2/matched_snp_IDs_annot_pVCF_noNA_noMultiallelic.txt \
        | cut -f3 | sort -u)
# → empty (0 diff lines)

# (b) .extract line set equivalence
diff <(zcat data/PMBB_Exome/matched_snp_IDs_annot_pVCF_noNA_noMultiallelic.extract.gz | sort -u) \
     <(sort -u analysis/daniel/outputs/phase2/matched_snp_IDs_annot_pVCF_noNA_noMultiallelic.extract)
# → empty (0 diff lines)

# Bit-identity check (sorted-unique stream md5sum)
zcat data/PMBB_Exome/matched_snp_IDs_annot_pVCF_noNA_noMultiallelic.extract.gz | sort -u | md5sum
sort -u analysis/daniel/outputs/phase2/matched_snp_IDs_annot_pVCF_noNA_noMultiallelic.extract | md5sum
# → 5e80ebc0faa5e68277cfeb948af8b1da on both
```

### What this validates

- The Phase 1 master list correctly seeds Phase 2 (no drift)
- The Python join logic in [`annot_IDs_vs_pVCF.py`](../../analysis/daniel/scripts/pmbb_exome/annot_IDs_vs_pVCF.py) produces the same chr:pos matches Daniel had
- The NA + multi-allelic filter applies identically
- The header-stripping + col-3 extract produces the same 9,667 IDs Daniel passed to plink

### What this does NOT validate

- That heavy mode (re-running `zgrep | cut -f 1-5` on the raw v2 pVCFs) produces the same per-chr files Daniel did. Light mode skipped this. If the v2 pVCFs were silently re-released between 2021 and now, we'd miss it. Low likelihood — those files are 5 years old and stable.
- That the 1,568 dropped variants (NA + multi-allelic) recover identically in Phase 9's add-back run.

---

## Issues encountered + fixes

### Issue 1 — `.extract` had an unwanted header line

**Symptom:** First successful-looking submission produced `.extract` with 9,668 lines instead of Daniel's 9,667. First line was `VCF_ID` (the column header).

**Root cause:** Daniel's runbook says `cut -f 3 matched_clean > extract`, which retains the header. But Daniel's actual `.extract.gz` on disk has 9,667 lines (no header) — he must have stripped it elsewhere (not documented in the runbook). **Critical:** `plink --extract` expects one SNP ID per line **with no header**. Keeping `VCF_ID` as the first line would make plink look for a SNP literally named "VCF_ID", fail to find it, and warn while silently dropping the search.

**Fix:** Changed Step 2.6 from `cut -f 3` to `tail -n +2 | cut -f 3` — explicit header skip. Comment in script flags this as undocumented behavior in the runbook.

### Issue 2 — `set -e` killed the validation block

**Symptom:** First submission's run log ended mid-validation. Job exited with status EXIT (not DONE), no error message in LSF stderr. The matched_clean validation logged success, but the `.extract` validation never ran.

**Root cause:** The validation pattern `diff_count=$( diff <(...) <(...) | wc -l )` has a subtle interaction with `set -euo pipefail`:
- When `diff` finds differences, it returns exit code 1
- `pipefail` makes the pipeline's exit code the leftmost non-zero — so the subshell returns 1
- `set -e` aborts the script when an unhandled non-zero exit occurs
- The script dies between computing `diff_count` and the `if [[ "$diff_count" -eq 0 ]]` check that would have decided what to log

In the first submission, the matched_clean check PASSED (diff returned 0, no issue), but the `.extract` check would have FAILED (because of the header issue from #1), making `diff` return 1 and killing the script.

**Fix:**
- Wrapped the validation block in `set +e` / `set -e` toggle, so non-zero `diff` exits don't abort the script
- Used `{ diff ... || true; } | wc -l` pattern to be doubly defensive
- The validation now correctly **reports** mismatches via the if/else rather than crashing

### Issue 3 — `du -h $file | cut -f1` reported "0" for sub-MB files

**Symptom:** Log messages like `(0, 11236 rows, 8s)` showed file size as "0" even though `ls -la` confirmed the file was 595 KB.

**Root cause:** Cosmetic only — `du -h` and `cut` interaction edge case for files smaller than du's default block reporting threshold. Doesn't affect correctness.

**Status:** Not fixed; cosmetic. The row count and elapsed time are the load-bearing metrics in logs.

---

## Performance summary

| Metric | Value |
|---|---|
| Total wall time | 20 s |
| CPU time | 10.01 s |
| Max memory | 8 MB |
| Step 2.1 (decompress 77 MB → 369 MB) | 3 s |
| Step 2.3 (Python join, 11.3 M rows × 11.6 k dict lookups) | 7 s |
| Steps 2.4-2.6 + validation | <1 s combined |

**For Phase 3+ resource sizing:** request `-M 64` (64 MB) is plenty for this kind of single-script Python-streaming work. We over-provisioned by ~128× in this run.

**If we ever need heavy mode (e.g. for v3 port):** expect ~10-20 min wall time with 22 parallel LSF tasks (one per chr) — bottleneck is gzip-decompress of 56-79 GB files, ~150 MB/s typical. Each task wants `-M 512 -W 30 -n 1`.

---

## Recommendations for Phases 3+

1. **The `.extract` file from Phase 2 is the canonical input for Phase 3's plink runs.** Path: [`analysis/daniel/outputs/phase2/matched_snp_IDs_annot_pVCF_noNA_noMultiallelic.extract`](../../analysis/daniel/outputs/phase2/matched_snp_IDs_annot_pVCF_noNA_noMultiallelic.extract). Pass to plink as `--extract <path>`.

2. **`set -e` / `pipefail` + `diff` validation needs explicit handling.** Phase 1's validation worked because the diff happened to return 0; Phase 2 caught the trap. **Pattern for all future phase scripts:** wrap validation blocks in `set +e` / `set -e` toggles, OR use `{ cmd || true; } | downstream` patterns.

3. **Plink 1.9 vs plink2 pilot is still pending.** Phase 3 is the natural place to do it — pick one chromosome, run `plink --vcf ... --extract --make-bed` with both versions, confirm the bed/bim/fam outputs are identical.

4. **Expect the 1,568 dropped variants to come back in Phase 9.** Daniel deliberately starts strict (no multi-allelic, no NA) for the first pass, then adds them back in `addBack_multiallelic_stoploss/`. When we replicate Phase 9, the master list there will be larger than 9,667.

5. **Heavy mode for v3 port:** when we eventually port Phase 2 to PMBB v3, we'll need to re-do Step 2.1 against the v3 pVCFs (Daniel's `vcf_SNP_IDs_allchr.txt.gz` is v2-only). That's the LSF job-array case — 22 parallel tasks, est. 10-20 min wall.

---

## Files produced by this phase

```
analysis/daniel/outputs/phase2/
├── vcf_SNP_IDs_allchr.txt                                    (369 MB, 11,337,300 rows — pVCF ID dump, all chrs)
├── matched_snp_IDs_annot_pVCF.txt                            (595 KB, 11,236 rows — incl. NA, multi-allelic)
├── matched_snp_IDs_annot_pVCF_noNA_noMultiallelic.txt        (488 KB, 9,668 rows incl. header — biallelic in-pVCF subset)
└── matched_snp_IDs_annot_pVCF_noNA_noMultiallelic.extract    (150 KB, 9,667 lines — plink --extract input)  ← Phase 3 reads this

analysis/daniel/logs/phase2/
├── run_20260512_213947.log
├── lsf_20260512_213946.out
└── lsf_20260512_213946.err  (empty — clean run)
```

---

## Open questions still pending (carries from Phase 0/1)

- **"8 signal-driving cases" definition** — unchanged. Phase 2 doesn't depend on it.
- **plink 1.9 vs plink2 pilot** — now actually relevant: Phase 3 is the first plink invocation.

## What Phase 3 will do (preview)

Phase 3 takes our 9,667-ID `.extract` file and runs `plink --vcf <chr_pVCF> --extract <our_extract> --make-bed --out genotypes/allIndvs_chr<N>` for each chromosome 1-22. The output is plink bed/bim/fam — a binary matrix of `[variant × person]` genotypes. **This is the first time per-individual genotype data enters the pipeline.** Once we have the matrix, all downstream burden tests (Phase 4+) work on those genotypes.

Estimated heavy mode (per chr): ~5-15 min plink wall per chr, 22 parallel = ~15-20 min total via job array.
Estimated light mode: skip — Daniel's per-chr plink files at `data/PMBB_Exome/genotypes/allIndvs_chr*_maf.001_noRels.{bed,bim,fam}.gz` exist and we'd reuse them.
