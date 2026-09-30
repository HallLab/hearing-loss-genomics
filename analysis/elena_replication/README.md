# Elena Replication — independent verification of the PMBB v4 rare-variant pipeline

**Status:** OPEN — scope fixed, no runs yet
**Author:** Andre Rico · **Opened:** 2026-09-28
**Target:** Elena's PMBB v4 SAIGE-GENE+ pipeline (`analysis/elena/`)

---

## 0. Why this exists

The PMBB v4 rare-variant gene-burden pipeline produces one summary table that downstream analyses
read as fact:

```
analysis/elena/HL_only_rarevariant/newmasks_combined_gene_results.csv
```

A review of that pipeline on 2026-09-28 found three defects in it (§3). One of them silently drops
~5,400 gene results; another silently drops eleven whole chromosome × MAF strata; the third is a
covariate choice that does not match the analysis decision on record. None was reported by the
pipeline, which exited successfully in every case.

That is enough to say the table cannot be used as input without being re-derived. This directory
re-derives it, from the institutional release, independently.

This is **not** an audit of Elena's work as a person or a judgement on it. It is the ordinary
practice of re-deriving a number before building on it. The pipeline is large and mostly sound; the
defects found are the kind that any long compute chain accumulates when nothing re-reads the output.

---

## 1. Scope — where the replication starts, and why not earlier

Elena's work splits into two lines. The question was whether the v4 pipeline inherits anything from
the Chapter 1 replication, i.e. whether replicating the v4 work requires replicating the v2 work
first. It does not. Checked 2026-09-28:

| | Line 1 — `Andre_Ch1Replication/` | Line 2 — `rarevariantExWAS/` |
|---|---|---|
| PMBB release | `PMBB-Release-2020-2.0` (v2) | `PMBB-Release-2026-4.0` (v4) |
| Method | BioBin (37 refs) + PLINK (61 refs) | SAIGE-GENE+, no BioBin anywhere |
| Inputs | Daniel's preserved intermediates in `data/PMBB_Exome/` | the institutional release only |
| Goal | reproduce Hui et al. 2023 | new discovery in v4 |
| Cross-references | — | none to Line 1 |

**Decision: the replication starts at Line 2.** The two lines share no data release, no method and no
inputs. Line 1 is a separate question (does the published paper reproduce?) and is out of scope here.

### In scope

1. **`rarevariantExWAS/`** — the backbone (steps 1–11): phenotyping, sample QC, variant annotation,
   ancestry/PCs, SAIGE step 1, and the covariate files everything downstream consumes.
2. **`HL_only_rarevariant/`** — the run that produced the summary table named in §0.

### Out of scope, for now

`rarevariant_geneburden/`, `tinnitus_only_rarevariant/`, both `*_ExWAS/` single-variant directories,
`ZNF175CarrierStudy/`, and Line 1. They answer different questions or feed nothing that is currently
being built on. `tinnitus_only_rarevariant/` carries the identical defect described in §3.1 and will
need the same treatment before it is used.

---

## 2. What "replicated" means here

Re-deriving every step from the institutional release is the strong form, and it is what this is.
The two directories in scope hold ~1.5 TB and the filesystem has 8.4 TB free, so the cost is real
but not prohibitive.

**The replication runs in pipeline order — Phase 1 through Phase 5 — not in order of cost.** Each
phase consumes the previous phase's output, so checking them out of order produces findings that
cannot be attributed: a disagreement at Phase 4 could originate in Phase 1, and there would be no
way to tell which. Running in dependency order means every phase is checked against a foundation
that has already been verified, so a disagreement belongs to the phase where it appears.

The phases are those of the pipeline itself, described in
[`pipeline_plan.md`](pipeline_plan.md) §2.

| # | Phase | What the replication does | Cost |
|---|---|---|---|
| 1 | Who is in the study | Rebuild cases, controls and the sample list from the release; diff against the pipeline's. | hours |
| 2 | Which variants count | Diff the per-gene variant masks against the group files the release publishes itself. No annotation is re-run. | minutes |
| 3 | What is adjusted away | Rebuild ancestry PCs and covariate files; diff against the pipeline's. | hours |
| 4 | The statistical test | Re-execute SAIGE step 1 and step 2 for a sample of strata. | days, ~1 TB |
| 5 | Reading the result | Rebuild the summary table from the per-stratum outputs. | minutes |

### Success criteria

- **Phase 1 passes** if the re-derived case, control and sample sets match the pipeline's
  person-for-person, or if every difference is explained by a documented rule.
- **Phase 2 passes** if the masks and the release group files agree on which variants belong to each
  gene and on the functional call for each, or if every disagreement is explained by a documented
  difference in threshold or transcript choice.
- **Phase 3 passes** if re-derived covariates match row-for-row on IID, PHENO, AGE, SEX and Batch,
  and the PC discrepancy in §3.3 is either justified by the scree plots or corrected.
- **Phase 4 passes** if re-run SAIGE p-values agree with the originals to within numerical tolerance
  for the tested strata.
- **Phase 5 passes** if the rebuilt summary table reproduces every row of the per-stratum outputs,
  including the ~5,400 currently lost (§3.1).

A phase that **fails** is a result, not a setback. It tells us which numbers taken from this
pipeline have to be restated — and, because the phases it depends on were verified first, it tells
us where the problem is rather than only that there is one.

### One thing that is not part of the replication

The summary table on disk is missing ~5,400 rows (§3.1), so nobody can presently read what the
existing results even say. Re-parsing the raw per-stratum outputs to recover them takes minutes and
is worth doing whenever convenient. It is **triage, not replication**: it makes existing output
legible and provides no independent evidence about whether that output is correct.

---

## 3. Defects this replication must resolve

Found 2026-09-28. These are the specific things a successful replication has to either reproduce or fix.

### 3.1 Three SAIGE outputs written without a header line

Of 781 result files under `HL_only_rarevariant/saige_results_newmasks/`, three have no header row:

```
AFR/pDM/AFR_chr8_pDM_maf0.01.txt      1,815 lines
AFR/pDM/AFR_chr8_pDM_maf0.001.txt     1,803 lines
AFR/pDM/AFR_chr8_pDM_maf0.0001.txt    1,776 lines
```

The merge script [`step3_2_merge.bsub`](../elena/HL_only_rarevariant/step3_2_merge.bsub) is correct;
`pandas.read_csv` took each file's first data row as column names, and `pd.concat` unioned the
result. That is why the merged CSV has 41 fields with columns literally named `ABRA`, `pDM`,
`1e-04`. Only the first 16 are real.

**Consequence:** ~5,400 gene rows — the whole AFR × pDM × chr8 stratum at all three MAF thresholds —
are absent from every table built from that file. Among them is **GRHL2**, a Definitive
autosomal-dominant hearing-loss gene (DFNA28), which was therefore never tested in that stratum.

**Not random.** The identical three files failed in the tinnitus run
(`tinnitus_only_rarevariant/saige_results_newmasks/`). Two independent executions, the same chunk.
The replication must determine what makes AFR/chr8/pDM reproducibly emit a headerless file.

### 3.2 Eleven SAIGE jobs lost to OOM without retry

`combined/ALL` holds 55 of the expected 66 files (22 chromosomes × 3 MAF). Missing:

```
chr1:0.01  chr1:0.001  chr2:0.01  chr11:0.01  chr12:0.01
chr16:0.01 chr16:0.001 chr17:0.01 chr17:0.001 chr19:0.01 chr19:0.001
```

All in the unrestricted mask, largest ancestry group, most permissive MAF. The killed job is
visible in the Nextflow work directory: `step2_SPAtests.R --chrom=1 --maxMAF_in_groupTest=0.01`
ended in `Killed`. Nextflow did not retry with more memory and the pipeline reported success.

The damage is contained to one mask: `pDM`, `pLOF` and `pLOF_pDM` are complete at 66/66 in all three
ancestry groups, and `EUR/ALL` and `AFR/ALL` are complete too. Only `combined/ALL` is short. Any
λ_GC or exome-wide summary computed for `combined/ALL` at MAF 0.01 or 0.001 was computed on partial
data.

### 3.3 PC counts diverge from the 2026-07-01 decision

[`step2_newcovariate.bsub`](../elena/HL_only_rarevariant/step2_newcovariate.bsub) consumes
`covariates_combined_5PCs`, `covariates_EUR_9PCs` and `covariates_AFR_10PCs`. The meeting decided
5–6 PCs. Per-ancestry scree plots exist (`rarevariantExWAS/{EUR,AFR}_ScreePlot.png` and the matching
`*_variance_explained.tsv`, 2026-07-14), so there is likely an empirical justification — but it was
never carried back to the recorded decision, and the scree plots have not been audited.

Within-ancestry PC counts differing from the combined-cohort count is defensible on its face. The
replication should either produce that justification or correct the covariates.

### Defects found during the replication itself

The three above are what the 2026-09-28 review found before any replication ran. Defects the
replication finds are recorded in the phase that found them, and listed in the status log (§7).

---

## 4. Isolation rules

- **`analysis/elena/` is read-only.** Nothing in this replication writes, moves or deletes anything
  there. Elena's directory is the reference being checked, and a reference that gets edited stops
  being one.
- **All outputs land here**, under `analysis/elena_replication/`.
- **Inputs come from the institutional release**, `/static/PMBB/PMBB-Release-2026-4.0/`, not from
  Elena's intermediates — except where a stage exists specifically to diff against them.
- **Every comparison is recorded**, including the ones that agree. A replication that only documents
  its disagreements cannot be distinguished from one that did not look.

---

## 5. Structure

```
analysis/elena_replication/
  README.md          this file — scope, defects, success criteria
  pipeline_plan.md   what the v4 pipeline does and what the replication does to each phase
  docs/              pages staged for Confluence — decision requests and phase reference records
  phase_1/           who is in the study
  phase_2/           which variants count
  phase_3/           what is adjusted away
  phase_4/           the statistical test
  phase_5/           reading the result
```

Each phase directory holds its own work, so a phase can be read and audited on its own:

```
phase_N/
  scripts/   checks, numbered in the order they run
  results/   FINDINGS.md plus the tables and diffs behind it
  data/      re-derived intermediates (gitignored)
  logs/      LSF job output
```

Phase directories are created as the replication reaches them, not up front.

Two conventions apply throughout: numbers are generated by scripts and never transcribed by hand,
and every finding names the script that produced it.

---

## 6. Open questions

| # | Question | Blocks |
|---|---|---|
| 1 | What makes AFR/chr8/pDM reproducibly emit a headerless SAIGE file? | §3.1, Phase 4 |
| 2 | Do the scree plots justify 9 and 10 PCs, or do the covariates need rebuilding? | §3.3, Phase 3 |
| 3 | `LRTOMT`, a known hearing-loss gene on chr11, is absent from the results for reasons unrelated to §3.1 — symbol/alias mismatch, or a real gap? | Phase 2 |
| 4 | ~~Which sample frame is right?~~ **Resolved, Phase 1; mechanism corrected 2026-09-30.** The cohort was built by merging the exome phenotype against the **imputed genotype `.fam`** (70,493 samples; 70,493 − 85 without exome = the 70,408 sample list). The 517 dropped have exome data and complete exome PCs but no imputed genotypes at all. An intermediate version of this row blamed a PC completeness check and rejected the `.fam` explanation — that rejection tested the *exome* LD-pruned `.fam` rather than the *imputed* one, and was itself the error. Whether the exclusion was deliberate goes to Nikki and Elena. See [`phase_1/results/FINDINGS.md`](phase_1/results/FINDINGS.md). | — |
| 5 | Do the pipeline's masks agree with the group files the release publishes? If not, is a third annotator needed to adjudicate? | Phase 2 |
| 6 | Variant QC appears to have run on the full 70,925 cohort while the association test ran on 57,632. Defensible, arguably better, but undocumented. Noticed during Phase 1; **to be established by Phase 2's own check**, not assumed. | Phase 2 |
| 7 | Why was the phenotype file overwritten 36 minutes after the SAIGE covariates were built from it, and did anyone know? 556 people the rules exclude were analysed as controls as a result. | Nikki, Elena — [`phase_1/results/FINDINGS.md`](phase_1/results/FINDINGS.md) Finding 3 |

---

## 7. Status log

| Date | Event |
|---|---|
| 2026-09-28 | Review of `analysis/elena/` completed; three defects recorded (§3). Scope fixed to Line 2 after confirming the two lines are independent. No replication runs yet. |
| 2026-09-29 | Phase 2 (variant annotation) brought into scope, after finding that the release publishes its own SAIGE group files — the check costs a file diff, not an annotation run. Rationale in [`pipeline_plan.md`](pipeline_plan.md) §7. |
| 2026-09-29 | Reorganised around the pipeline's own phases. An earlier draft ordered the work as Stages A–D by cost; that vocabulary was dropped because it mapped almost one-to-one onto the phases while reordering them, which made findings hard to attribute. |
| 2026-09-30 | **Phase 1 mechanism corrected.** Elena's published documents (`elena_publishes/`) log the merge that built the sample list: `FAM samples: 70493 | Phenotype samples: 70925 | Matched samples: 70408`. The 517 were dropped by merging against the imputed genotype `.fam`, not by a PC completeness check. The consequence is unchanged — 517 exome-sequenced people, 40 of them cases, excluded for lacking array data — but the stated mechanism was wrong, and so was the earlier "correction" that rejected the `.fam` explanation after testing the exome `.fam` instead of the imputed one. |
| 2026-09-29 | **Phase 1 Finding 3.** The phenotype file on disk is not the one the analysis used — it was overwritten 36 minutes after the SAIGE covariates were built from it. The run consumed 57,632 people, not 57,080: Finding 1's 427 are still absent, but 556 people the phenotype rules exclude were analysed as controls. Two defects in opposite directions. Also removed a premature claim about Phase 2 from the Phase 1 reference page. |
| 2026-09-29 | **Phase 1 complete — split verdict.** Phenotype definition reproduces exactly: 70,925/70,925 person-for-person, cases 6,752. Sample frame does not: 427 analysable participants (40 cases, 387 controls) are dropped for missing imputed PCs the analysis does not use — a defect found by the replication, not by the review. Open question 4 resolved and its stated mechanism corrected. [`phase_1/results/FINDINGS.md`](phase_1/results/FINDINGS.md) |
