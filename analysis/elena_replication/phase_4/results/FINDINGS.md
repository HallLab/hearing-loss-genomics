# Phase 4 — the association test

**Status:** complete · **Arm:** corrected only
**Scripts:** `01_split_masks_by_chrom.sh` · `02_saige_step1.bsub` · `03_saige_step2.bsub` · `04_merge_results.py` · `05_calibration_and_consistency.py`
**Design and deviations:** [`../PLAN.md`](../PLAN.md)

Phase 4 covers the pipeline's Step 8. The results are below; Findings 1–4 that follow them are
defects in the results the pipeline delivered, found while building ours against them.

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

## Results — 594 of 594 cells, 472,453 tests

| cohort | tests | genes | min p | 0.05/tests | 0.05/genes |
|---|---:|---:|---|---|---|
| combined | 159,960 | 17,943 | 4.21 × 10⁻⁶ | 3.13 × 10⁻⁷ | 2.79 × 10⁻⁶ |
| EUR | 159,197 | 17,928 | 1.54 × 10⁻⁵ | 3.14 × 10⁻⁷ | 2.79 × 10⁻⁶ |
| AFR | 153,296 | 17,789 | 7.23 × 10⁻⁶ | 3.26 × 10⁻⁷ | 2.81 × 10⁻⁶ |

**No gene reaches exome-wide significance in any cohort**, under either defensible denominator —
one test per mask × MAF combination, or one per gene. The corrected arm finds nothing. A third,
still more permissive denominator is examined below and changes no answer.

Nor does the lenient correction. Phase 5 ran Benjamini-Hochberg at two levels, and found **zero
genes at q < 0.05** in every cohort, with the smallest q between 0.26 and 0.79 — see
[Phase 5, Finding 0](../../phase_5/results/FINDINGS.md). That matters because FDR is where a weak
but real effect shows up first, so "nothing under Bonferroni" and "nothing under FDR either" are
different strengths of claim.

**Neither does hers**, on either basis, and each is checked against *her own* gene and test counts
rather than ours:

| | her best p | gene | 0.05/her genes | 0.05/her tests |
|---|---|---|---|---|
| on the same three masks | 3.33 × 10⁻⁶ | `TMC3-AS1` (combined) | 2.63 × 10⁻⁶ | 9.09 × 10⁻⁸ |
| | 9.05 × 10⁻⁶ | `TMC3-AS1` (EUR) | 2.64 × 10⁻⁶ | 9.15 × 10⁻⁸ |
| | 4.95 × 10⁻⁶ | `IGSF9` (AFR) | 2.68 × 10⁻⁶ | 9.54 × 10⁻⁸ |
| including her `ALL` mask | 2.07 × 10⁻⁶ | `PRIMPOL` (EUR) | 2.04 × 10⁻⁶ | 8.06 × 10⁻⁸ |
| | 2.96 × 10⁻⁶ | `AARS1` (AFR) | 2.04 × 10⁻⁶ | 8.36 × 10⁻⁸ |

Nothing passes. Her closest call is `PRIMPOL` at 2.07 × 10⁻⁶ against a per-gene bar of
2.04 × 10⁻⁶ — short by 1.5%, which is as near as either arm gets to a finding.

Her bars are *stricter* than ours, and the mask defect is why: her gene count is about 24,500
against our 17,943, because her masks carry the non-coding genes. The defect inflates her
denominator at the same time as it populates her top of list.

So the headline conclusion is concordant: this cohort, at this size, does not support a
hearing-impairment gene burden finding.

That concordance is worth stating plainly, because it is the opposite of what a reader might expect
from Phases 2 and 3. Four mask defects and surplus PCs did not manufacture a false positive. What
they did instead is reshuffle the ranking underneath a null result.

### The ranking is materially different

Top 50 genes per cohort, by best p-value, our arm against hers on the same three masks:

| cohort | shared in top 50 |
|---|---:|
| combined | 24 / 50 |
| EUR | 17 / 50 |
| AFR | 15 / 50 |

Half to two thirds of each top-50 list changes. If this cohort were larger — or if these lists were
used to pick genes for follow-up, which is what top-gene tables are for — the defects would decide
which genes got looked at.

### Her rank-1 gene is a lncRNA, and it disappears

**`TMC3-AS1` is the top gene in her combined and EUR results**, at p = 3.33 × 10⁻⁶ and
9.05 × 10⁻⁶. It is annotated `lncRNA`. A pLOF burden test asks whether losing the protein's function
associates with the phenotype, and this gene has no protein to lose. The test rests on
`Number_rare = 2`: two rare variants carry the whole result.

It is absent from our masks entirely — it fell out in Phase 2 with the 1,101 non-coding genes Nikki
raised. In AFR, where it survived to rank 6,707 in her results, it never mattered.

The name is part of why this is worth flagging rather than filing. `TMC3-AS1` is antisense to
`TMC3`, and `TMC1` is an established deafness gene. A reader skimming a top-gene table sees "TMC"
and reads plausibility into it. An artefact that looks like a finding is more dangerous than one
that looks like noise.

Three others in her top 50 are not protein-coding either, with their release biotype:

| gene | cohort, rank | biotype |
|---|---|---|
| `SLCO1B7` | EUR, rank 8 | `transcribed_unprocessed_pseudogene` |
| `NUP153-AS1` | combined, rank 48 | `lncRNA` |
| `TRGV9` | EUR, rank 47 | `TR_V_gene` (T-cell receptor variable segment) |

AFR's top 50 has none.

### Known deafness genes move both ways

This is where the result refuses to be tidy, so it is recorded as it is rather than as it would be
convenient:

| gene | combined, ours | combined, hers | EUR, ours | EUR, hers |
|---|---|---|---|---|
| `SIX1` (DFNA23 / BOR) | **rank 7**, p 5.24 × 10⁻⁵ | rank 38, p 5.22 × 10⁻⁴ | rank 6 | rank 4 |
| `COCH` (DFNA9) | rank 22, p 3.88 × 10⁻⁴ | **rank 12**, p 1.18 × 10⁻⁴ | rank 373 | **rank 68** |

`SIX1` improves markedly in the corrected arm. `COCH` gets worse, in both cohorts. Two established
genes moving in opposite directions is **not** evidence that the corrected arm recovers known
biology better — it is evidence that the lists differ, which we already knew. Neither gene is
anywhere near significant in either arm, so neither movement should be read as a result.

What can be said without overreach: one specific artefact at rank 1 is gone, and the rest of the
reshuffling has no established direction. A systematic enrichment test against the ClinGen hearing
loss gene list would settle it, and is deliberately not done here — that list lives in `cycle_2`,
and this folder's isolation rule keeps `cycle_2` out. It is a public endpoint, so it can be fetched
independently into Phase 5 if the lab wants the test.

### The absence is real, not a deflated analysis

"No gene passes the threshold" and "there is nothing here" are different claims, and a conservative
test produces the first without supporting the second. Under the null, the fraction of tests below
*p* should be *p*. Measured at one MAF cutoff, so the same gene is not counted three times:

| cohort | p<0.05 | p<0.01 | p<0.001 |
|---|---|---|---|
| combined | 0.92 | 0.95 | 0.86 |
| EUR | 0.89 | 0.93 | 1.15 |
| AFR | 0.93 | 0.99 | 0.93 |

Well calibrated, slightly conservative. Neither inflated, which would make the p-values
meaningless, nor deflated, which would hide real signal.

**And nothing replicates across ancestries.** EUR and AFR are independent samples:

| | observed | expected | ratio | p |
|---|---:|---:|---:|---|
| min-p < 0.01 in both | 22 | 18.8 | 1.17 | 0.26 |
| min-p < 0.005 in both | 3 | 5.3 | 0.57 | 0.90 |
| min-p < 0.001 in both | 0 | 0.3 | — | 1 |

**The expectation here must use the measured marginal rates, not the threshold.** A first pass used
0.01 and produced an apparent 12-fold excess that does not exist: per gene we take the minimum p
over 9 correlated tests, so the marginal rate of min-p < 0.01 is 0.031 in EUR and 0.034 in AFR, not
0.01. The warning is written into the script, because that arithmetic manufactures findings.

### The three MAF cutoffs are largely the same test

`--maxMAF_in_groupTest` decides which variants enter the gene's set: up to 1%, 0.1% or 0.01% in
frequency. The cutoffs are nested, and in exome data a pLOF or damaging-missense variant is almost
always very rare, so they often select the same variants. Over 156,946 (cohort, mask, gene) triples:

| | share |
|---|---:|
| all three cutoffs give an **identical** p | 31.7% |
| only 0.001 and 0.01 identical | 42.3% |
| all three distinct | **26.0%** |

The variant counts say why. The average gene's set holds 27.9 variants at the strictest cutoff and
30.1 at the loosest — loosening from 0.01% to 1% adds about two variants.

**This makes our reported Bonferroni bar conservative, and the conclusion survives it anyway:**

| cohort | min p | 0.05/tests | collapsing identical MAF | 0.05/genes | passes? |
|---|---|---|---|---|---|
| combined | 4.21 × 10⁻⁶ | 3.13 × 10⁻⁷ | 4.92 × 10⁻⁷ | 2.79 × 10⁻⁶ | no |
| EUR | 1.54 × 10⁻⁵ | 3.14 × 10⁻⁷ | 5.24 × 10⁻⁷ | 2.79 × 10⁻⁶ | no |
| AFR | 7.23 × 10⁻⁶ | 3.26 × 10⁻⁷ | 4.92 × 10⁻⁷ | 2.81 × 10⁻⁶ | no |

The masks are correlated too — `pLOF_pDM` contains `pLOF` and `pDM` — so the effective number of
independent tests is smaller still. It does not matter: nothing clears even the per-gene bar, which
is the most permissive one anyone would defend.

### Where the MAF cutoff does change the answer, and GJB3

68 tests improve more than 100-fold at the strictest cutoff. That pattern — signal concentrated in
the rarest variants and diluting as slightly commoner ones enter — is how a real gene would behave,
since a variant that genuinely destroys function is kept rare by selection.

| gene | cohort | p at 10⁻⁴ | p at 10⁻² | fold |
|---|---|---|---|---:|
| `MLLT6` | AFR | 8.5 × 10⁻⁶ | 0.053 | 6,263 |
| `ENTREP1` | AFR | 6.2 × 10⁻⁵ | 0.037 | 594 |
| `ACLY` | AFR | 1.3 × 10⁻⁵ | 0.0037 | 292 |
| `DHCR7` | EUR | 1.6 × 10⁻⁴ | 0.032 | 203 |
| `PTPN23` | combined | 4.0 × 10⁻⁵ | 0.0074 | 185 |
| **`GJB3`** | combined | 2.4 × 10⁻⁴ | 0.024 | 100 |

**`GJB3` is worth a note for Phase 5.** Connexin 31, DFNA2B — an established non-syndromic deafness
gene, at rank 14 of 17,943 in the combined cohort. And it shows the pattern cleanly: at the
strictest cutoff `Pvalue_Burden` is 2.4 × 10⁻⁴ with all variants collapsed; loosening to 0.1% lets
nine commoner variants in and the burden p falls to 0.093.

This is a flag, not a result. p = 2.4 × 10⁻⁴ against a bar of 2.8 × 10⁻⁶ is two orders of magnitude
short, and with 17,943 genes tested, a known deafness gene landing at rank 14 by chance is not
surprising. It is recorded because it is the one top-ranked gene whose identity and behaviour both
point the same way.

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
