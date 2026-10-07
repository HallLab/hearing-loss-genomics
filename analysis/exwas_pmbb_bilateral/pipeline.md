# Pipeline

Five phases, from diagnosis codes to a ranked gene list. Two are inherited intact; three are run
here. Each stage names its inputs, its outputs and the premise that governs it.

```
   RELEASE                  PHASE 1          PHASE 2         PHASE 3
   PMBB v4            who is in the      which variants   what is adjusted
                        study              count            away
┌──────────────┐      ┌────────────┐     ┌───────────┐    ┌─────────────┐
│ phecode_x    │─────▶│  cases     │     │   VEP     │    │  ancestry   │
│ observation  │      │   3,164    │     │  masks    │    │    PCA      │
│ exome VEP    │─────────────────────────▶│ 3 masks  │    │  5 / 4 / 3  │
│ exome PLINK  │      │  controls  │     └─────┬─────┘    └──────┬──────┘
└──────────────┘      │   50,755   │           │                 │
                      └─────┬──────┘           │          ┌──────▼──────┐
                            │                  │          │ covariates  │
                            └──────────────────┼─────────▶│ per cohort  │
                                               │          └──────┬──────┘
                                        PHASE 4│                 │
                                     ┌─────────▼─────────────────▼───┐
                                     │  SAIGE step 1  null models ×3 │
                                     │  SAIGE step 2  594 tasks      │
                                     └───────────────┬───────────────┘
                                             PHASE 5 │
                                     ┌───────────────▼───────────────┐
                                     │  merge · Manhattan · QQ       │
                                     │  tables · ours vs replication │
                                     └───────────────────────────────┘
```

---

## Phase 1 — who is in the study · **built here**

The only place this analysis differs from the replication it inherits.

| | |
|---|---|
| in | `conditions_phecode_x`, OMOP `observation`, exome ancestry list |
| out | `phase_1/results/cohort.tsv` — 3,164 cases, 50,755 controls |
| premise | **P1** cases, **P2** controls, **P3** exome frame, **P4** observation table |

**Step 1.1 — cache the release tables.** `00_extract_release_tables.sh` pulls every `SO_39x`
diagnosis row and every tinnitus row from `observation`. Both source files are multi-GB and are
scanned once.

**Step 1.2 — assign case, control or excluded.** `01_phenotype_bilateral.py`. A case needs
`SO_396.2` and `SO_396.8` on the same date, twice. A control needs no ear-family evidence at all.
Everyone else is excluded, and the reasons are counted rather than lumped: 7,595 have hearing loss
that is not bilateral sensorineural, 9,411 have other ear evidence only.

Nobody moves between groups. Narrowing the cases does not widen the controls.

## Phase 2 — which variants count · **inherited, not re-run**

| | |
|---|---|
| in | release VEP (already consumed upstream) |
| out | `phase_2/results/masks_v2/{pLOF,pDM,pLOF_pDM}.txt` |
| premise | **P5** rebuilt masks, **P6** no `ALL` |

Masks assign each variant to a gene and a damage class. They depend on the genome annotation and
not at all on who is a case, so the rebuild done in the replication carries over unchanged. The
scripts are here and will reproduce the files; the files are here so nothing has to be rerun.

`ALL.txt` is present but unused — the decision is reversible without rebuilding anything.

```
pLOF       17,841 genes     462,144 entries
pDM        17,623 genes     890,332 entries
pLOF_pDM   17,945 genes   1,340,936 entries
```

## Phase 3 — what is adjusted away · **PCA inherited, covariates built here**

| | |
|---|---|
| in | Phase 1 cohort, copied PCA, release covariates |
| out | `phase_3/results/covariates/{combined,EUR,AFR}.txt` |
| premise | **P7** PC counts, **P8** categorical batch |

**Step 3.1 — ancestry PCA.** *Inherited.* Ran over everyone with exome data per ancestry (51,867
EUR, 14,927 AFR), not over the analysis cohort, which is exactly why a new phenotype does not
require a new PCA. The LD-pruned sets double as the GRMs for step 1.

**Step 3.2 — covariate files.** *Built here.* `IID PHENO AGE AGE2 SEX Batch PC1..PCn`, one per
cohort, with the new phenotype. People without a recorded age drop out — age is a covariate the
model consumes, so excluding them is correct, unlike the imputed-PC cut P3 reverses.

## Phase 4 — the test · **built here**

| | |
|---|---|
| in | covariates, masks, GRMs, exome genotypes |
| out | 594 result files |
| premise | **P6** three masks, **P8** categorical batch |

**Step 4.1 — null model.** Three, one per cohort. Fits the expectation of who has hearing loss from
age, sex, batch, ancestry and relatedness, without looking at any gene. Relatedness is the expensive
part and does not depend on the gene, which is why it is computed once and reused.

**Step 4.2 — association.** 3 cohorts × 22 chromosomes × 3 masks × 3 max-MAF cutoffs = **594 tasks**,
as an LSF array. An array element that fails does not touch its siblings — the pipeline this
replaces used `errorStrategy = 'terminate'`, where one out-of-memory task cancelled ten healthy ones
and the merge never noticed.

Each gene yields three p-values: Burden (the sum), SKAT (the dispersion) and SKAT-O (the mixture).

## Phase 5 — reading the result · **built here**

| | |
|---|---|
| in | the 594 result files |
| out | merged tables, Manhattan, QQ, top hits |
| premise | **P9** both corrections, both levels |

**Step 5.1 — merge.** Refuses rather than globs: all 594 cells present, every file carrying the
expected header, no file still being written, no blank gene, no non-numeric p-value.

**Step 5.2 — calibration.** QQ plots and λ. Establishes whether a null result is a real absence or a
deflated test, which a p-value table alone cannot distinguish.

**Step 5.3 — tables.** Top hits with a fragility flag, because a p-value alone ranks a five-allele
gene alongside a four-hundred-allele one.

**Step 5.4 — against the replication.** The same genes under the broad `SO_396` phenotype and under
this one. That comparison is the point of running this at all: it measures what the phenotype
decision costs and buys, with everything else held fixed.

---

## What is not in this pipeline, and why

**Audiograms.** The strategic direction, and the clinical lead's preference —
*"we certainly think the audiograms are going to be a safer bet"*. They are quantitative and remove
the guesswork that P1 is trying to patch with diagnosis codes. They need a REDCap-to-PMBB ID bridge
that does not exist yet. When it does, this pipeline changes only in Phase 1.

**The single-variant ExWAS.** A different analysis on the same cohort: common variants, MAF ≥ 0.01,
genome-wide threshold 5 × 10⁻⁸. Complementary to this one, not comparable — the two test disjoint
sets of variants.

**Replication in another biobank.** All of Us or UK Biobank. The third sense of "replication" the
clinical lead named, and the one neither this folder nor the replication addresses.
