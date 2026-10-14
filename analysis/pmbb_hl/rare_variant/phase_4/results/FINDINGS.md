# Phase 4 — the association test

**Run 2026-10-14.** Step 1 job `49686445` (6/6), step 2 job `49687057` (132 tasks) plus
`49687xxx` re-running one cell. SAIGE-GENE+ 1.5.0.

## Result: nothing is significant, under either bar, in any cell of the 2 × 3 design

| phenotype | cohort | genes | min omnibus p | best gene | Bonferroni bar | hits | min q | FDR hits |
|---|---|---:|---:|---|---:|---:|---:|---:|
| broad | combined | 17,943 | 2.88e-05 | `BPTF` | 2.787e-06 | 0 | 0.1984 | 0 |
| broad | EUR | 17,928 | 3.31e-05 | `GYS2` | 2.789e-06 | 0 | 0.4770 | 0 |
| broad | AFR | 17,789 | 1.40e-05 | `CFH` | 2.811e-06 | 0 | 0.2482 | 0 |
| bilateral | combined | 17,941 | 2.06e-05 | `PCBD1` | 2.787e-06 | 0 | 0.2394 | 0 |
| bilateral | EUR | 17,921 | **6.33e-06** | **`FBXO25`** | 2.790e-06 | 0 | **0.1134** | 0 |
| bilateral | AFR | 17,778 | 1.46e-05 | `IGSF9` | 2.812e-06 | 0 | 0.2593 | 0 |

**The closest anything comes is `FBXO25` in bilateral EUR**, at q = 0.1134 against the
0.10 bar — short by 13%. It is 2.3× short of the Bonferroni bar. It was also the top
bilateral EUR gene in the first pass, at 1.37e-05; with the release's 15 PCs it is
6.33e-06, so it strengthened rather than weakened under better ancestry correction.
That is worth watching and is not a result.

## The three levels of multiplicity, and who pays for each

| level | corrected by |
|---|---|
| SKAT-O's ρ sweep | SAIGE, internally |
| the 9 cells of one gene | **Cauchy/ACAT**, in the merge |
| the ~17,900 genes | **Bonferroni and FDR**, in the merge |

The nine cells are not nine hypotheses. They are nine ways of asking one question about
one gene, and they are strongly correlated — the MAF cutoffs are nested, so everything
in 0.0001 is also in 0.001 and in 0.01. Bonferroni across them would charge nine times
for one question. ACAT combines dependent p-values without needing to know the
dependence structure.

### The two bars

**Bonferroni, `0.05 / genes`** — fixed, known before the run, the same number for every
gene. It controls the chance of **even one** false positive in the whole analysis.

**FDR (Benjamini-Hochberg), q < 0.10** — asked for by Nikki Palmiero on 2026-10-14,
*"keep the Bonferroni, just also add FDR P of less than 0.1"*. It has no fixed p-value
cutoff: BH finds the largest k with p(k) ≤ (k/n)·q, so the requirement loosens as you go
down the list. In the broad combined cohort the first gene would need p ≤ 5.57e-06, the
tenth p ≤ 5.57e-05, the hundredth p ≤ 5.57e-04. The effective cutoff is a result, not a
parameter.

Both are applied to the **omnibus**, one value per gene, not to the cells.

FDR is reported for the gene burden only. *"You do not do FDR for ExWAS"* — and the
split is principled: FDR is a threshold choice over a count of tests we performed, while
the ExWAS is judged against the 5e-8 convention, which counts nothing we did.

## Deduplication, before the combination

Raised by Nikki Palmiero on 2026-10-14:

> *"Let's say these two have the exact same rows […] all of the numbers are the same.
> What that means is there's no pDM. It's doing the test again. So we would drop this
> because it's a duplicate. SAIGE doesn't notice that it's a duplicate."*

Duplicate cells are detected on identity of `(max_MAF, MAC, Pvalue)` within a gene and
removed before the ACAT:

| phenotype | cohort | cells | dropped | `pLOF == pLOF;pDM` | `pDM == pLOF;pDM` |
|---|---|---:|---:|---:|---:|
| broad | combined | 159,948 | 1,542 | 1,100 | 442 |
| broad | EUR | 159,179 | 2,164 | 1,366 | 798 |
| broad | AFR | 153,251 | **6,512** | 2,575 | 3,937 |
| bilateral | combined | 159,859 | 1,614 | 1,138 | 475 |
| bilateral | EUR | 159,007 | 2,271 | 1,409 | 861 |
| bilateral | AFR | 152,749 | **6,911** | 2,717 | 4,194 |

**AFR carries four times the duplicates of the combined cohort.** With 11,334 people
rather than 57,498, far more genes have no qualifying variant of one class at a given
MAF cutoff, so `pLOF:pDM` collapses onto whichever class remains. The smaller the
cohort, the more of the grid is redundant.

### Two mechanisms produce them

**Mask-level.** A gene with no pDM variants at all has `pLOF:pDM` identical to `pLOF`.
445 of 17,945 genes, measured directly on the phase 2 masks.

**MAF-level.** After a cutoff is applied, the surviving variants can coincide even when
the full masks differ — a gene whose pLOF variants are ultra-rare and whose pDM variants
are not will have identical sets at `maxMAF = 0.0001`. Measured on the first pass: 1,606
duplicate cells against the 1,335 the mask-level count predicts, so 271 came from this.

### Why it usually does not matter, and when it does

ACAT is a **mean** of transformed p-values. If every cell is duplicated exactly once —
the common case, three masks collapsing to two — the mean is unchanged and the omnibus
is identical. Measured on the first pass: 550 genes affected, **median ratio 1.000**,
range 0.69 to 1.73, no verdict changed.

The effect appears only when duplication is **uneven**. Here it changes the omnibus in
**424 to 875 genes per arm**, 2.4% to 4.9%:

| | genes where our omnibus differs from SAIGE's Cauchy |
|---|---:|
| broad combined | 715 of 17,943 |
| broad EUR | 481 of 17,928 |
| broad AFR | 875 of 17,789 |
| bilateral combined | 525 of 17,941 |
| bilateral EUR | 424 of 17,921 |
| bilateral AFR | 787 of 17,778 |

It is corrected regardless of size, because the direction of the error is the dangerous
one: if the duplicated cell is a gene's most extreme, counting it twice pulls the
omnibus **toward** significance. Harmless while nothing is near a bar; not harmless when
something is.

### SAIGE's own Cauchy row is discarded

SAIGE emits one `Cauchy` row per gene, combining all nine cells **including the
duplicates**. Verified on the first pass's chr21: it equals an ACAT over all cells in 36
of 36 genes that contain duplicates, and over the deduplicated cells in only 35. The
merge recomputes the omnibus itself.

## One definitional choice, recorded

11,540 variants are both pLOF and pDM. The SAIGE group file carries one annotation per
variant and gives pLOF precedence, so those variants enter the pLOF test and are
**absent from the pDM test**. In this analysis **`pDM` means "probably damaging and not
truncating"**, not "probably damaging".

Kept deliberately. A variant that truncates the protein is already tested by the arm
that exists for truncation; counting it again under pDM measures the same thing twice —
the same objection as the duplicate-cell problem, one level down. Reversing it means a
second group file with pDM precedence and a separate run of the pDM arm: three times the
step 2 tasks, and overlapping masks that raise the correlation between the cells ACAT
combines.

## What the design was

| | |
|---|---|
| step 1 | 6 null models, 2 phenotypes × 3 cohorts, PCs from the release at 80% variance (11/15/16) |
| step 2 | 132 tasks, 2 × 3 × 22 chromosomes |
| per call | the whole grid — `--maxMAF_in_groupTest=0.0001,0.001,0.01` and `--annotation_in_groupTest="pLOF,pDM,pLOF:pDM"` |

One call produces all nine cells per gene. The first pass's replication arm learned this
the expensive way, running 594 tasks where 66 would do.

## Three defects the checks caught before anything ran

| | |
|---|---|
| the chromosome split dropped every annotation line | the mask has two lines per gene and the chromosome is readable only on the first. `REFUSING: 17945 lines lost in the split` — SAIGE would have received genes with no annotations |
| `--outputPrefix_varRatio` would have doubled the suffix | SAIGE appends `.varianceRatio.txt` to whatever prefix it is given; the first pass produced `EUR.varianceRatio.txt.varianceRatio.txt`. Caught before the six models were written |
| the merge checked only the main array | with a retry task running under a different job name, it would have merged while chr17 was still being written |

One cell, `broad_AFR` chr17, exited before SAIGE started with an empty stderr — the same
transient node failure seen twice before in this project. Re-run cleanly.

## Files

| | |
|---|---|
| `results/omnibus_<phenotype>_<cohort>.tsv` | one row per gene: omnibus p, best cell, both bars, both verdicts |
| `results/cells_<phenotype>_<cohort>.tsv` | the cells, before deduplication |
| `results/04_merge.json` | counts, collisions, and the comparison against SAIGE's Cauchy |
| `results/step1/`, `results/step2/` | SAIGE's own output |
