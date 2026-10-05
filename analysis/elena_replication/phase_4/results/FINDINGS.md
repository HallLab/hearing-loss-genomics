# Phase 4 — the association test

**Status:** step 1 complete · step 2 running (job 49232320, 594 tasks) · **Arm:** corrected only
**Scripts:** `01_split_masks_by_chrom.sh` · `02_saige_step1.bsub` · `03_saige_step2.bsub` · `04_merge_results.py`
**Design and deviations:** [`../PLAN.md`](../PLAN.md)

Phase 4 covers the pipeline's Step 8. Findings 1–4 below are defects in the results the pipeline
delivered, found while building ours against them. Our own results are not in yet.

---

## Step 1 reproduces, on different inputs

| cohort | N | PCs | our variance ratio | hers | our tau₂ |
|---|---:|---:|---|---|---:|
| AFR | 11,334 | 3 | 0.999999998377 | 0.999999998107 | 0 |
| EUR | 42,779 | 4 | 0.985991635 | 0.985951387 | 0.0872 |
| combined | 57,498 | 5 | 0.982009077 | 0.979613817 | 0.1182 |

Agreement to the third or fourth decimal, with **our** GRM, roughly half the PCs (3 vs 10, 4 vs 9)
and `Batch` declared categorical. The null model is robust to those choices, which is what we
expected: they move power, not calibration.

`tau₂ = 0` in AFR means the polygenic variance component converged to zero, so there SAIGE reduces
to logistic regression on the covariates. Not a defect of ours — her AFR variance ratio is also 1 to
nine decimals, so hers did the same. It is the sample: 11,334 people, 1,285 cases.

---

## Finding 1 — 11 of her 792 step-2 cells produced no output

All in the combined cohort's `ALL` mask, at the two looser MAF cutoffs, on the gene-densest
chromosomes:

```
combined / ALL / maf0.001 : chr1, chr16, chr17, chr19
combined / ALL / maf0.01  : chr1, chr2, chr11, chr12, chr16, chr17, chr19
```

The cause is in her Nextflow work directory, and it is a chain rather than a single failure:

1. **`combined / chr1 / ALL / maf0.01` exceeded the 64 GB limit.** LSF killed it —
   `TERM_MEMLIMIT`, exit 137, `.command.log` ending in `Killed`.
2. **The other 10 were cancelled, not failed.** Her config sets `errorStrategy = 'terminate'` and
   `maxRetries = 0`, so when the chr1 task died Nextflow killed every task still in flight. They
   carry exit 130, SIGINT.
3. **Nothing downstream noticed** — see Finding 3.

So one task hit a real resource ceiling and took ten healthy ones with it.

**This is why our step 2 is an LSF job array rather than a Nextflow pipeline.** An array element
that fails does not touch its siblings; we lose that element and resubmit it. It is also the mask we
are not running: `ALL` carries 15,951,289 mask entries against 1,340,936 for `pLOF_pDM`, so the cell
that OOM'd is about twelve times heavier than our heaviest remaining one.

---

## Finding 2 — three of her result files lost their header line

`AFR / pDM / chr8`, at all three MAF cutoffs. The first data row sits where the header should be, so
anything reading them takes gene `ABRA` and its numbers as the column names:

```
ABRA | pDM | 1e-04 | 0.0289903099999999 | ...     <- line 1 of AFR_chr8_pDM_maf0.0001.txt
```

778 of her 781 files carry the correct 13-column header; these three do not.

Two consequences, and the second is worse than the first:

- **Gene `ABRA`'s result is gone** in all three — consumed as the header.
- **5,391 rows of her merged table are unusable.** `pd.concat` aligns by column name, so those
  files' rows went into new columns named `ABRA`, `pDM`, `1e-04` and so on, leaving `Region`,
  `Group` and `Pvalue` blank for every one of them. The real values are present, parked in phantom
  columns. Any filter on `Pvalue` drops these rows silently; any `groupby("Region")` collapses them
  into a single nameless gene.

---

## Finding 3 — the merge checks neither completeness nor schema

`step3_2_merge.bsub` calls `rglob("*.txt")`, concatenates what comes back, prints `value_counts`
for ancestry, mask and MAF, and writes the CSV. It never asks whether the expected cells are all
there, nor whether the files agree on their columns. Findings 1 and 2 are therefore invisible in its
output: `newmasks_combined_gene_results.csv` has 1,830,387 rows and looks complete.

This is the finding that generalises. The two upstream defects are ordinary operational accidents —
a job too big for its memory limit, three files written without a header. What turns them into
wrong numbers is that the step which could have caught both was written to glob rather than to
check.

Our `04_merge_results.py` refuses instead. Before it writes anything it asserts that all 594 cells
exist with no duplicates and nothing unexpected, that every file carries the 13 expected columns,
that each file's `max_MAF` matches its filename, and that no row has a blank or unnamed gene or a
non-numeric p-value. Any failure stops the merge and names the files.

---

## Finding 4 — a genome-spanning pseudo-gene was tested

Her mask files carry a gene called `-`, holding every variant VEP could not assign to a symbol.
SAIGE ran burden tests on it like any other gene: **1,591 rows** across her results, pooling
variants from across the whole genome into one "gene".

Kept in proportion: the smallest p-value the pseudo-gene reached is 1.44 × 10⁻⁴, far from
significant after correction, so it did not contaminate her top hits. It is a test that cannot mean
anything rather than a false positive.

Our rebuilt masks dropped it in Phase 2 — `-` is not `protein_coding`, so it fell out with the 1,101
non-coding genes Nikki raised. Our smoke-test output confirms it: 186 rows, 186 distinct genes, no
blank or `-` gene, and no `Cauchy` row either, since we pass one annotation per mask.
