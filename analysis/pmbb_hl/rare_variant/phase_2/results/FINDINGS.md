# Phase 2 — the burden masks

**Run 2026-10-14**, job `49676518`, 22/22 tasks. Built from the release's own VEP
annotation (`Exome/vep_annotations/`, 106 GB), one chromosome per task.

## Result

| mask | genes | entries |
|---|---:|---:|
| pLOF | 17,841 | 462,144 |
| pDM | 17,623 | 890,332 |
| pLOF_pDM | 17,945 | 1,340,936 |

139,846,060 annotation rows fall inside the cohort's genotypes.

**One build serves both phenotypes.** A mask depends on the annotation and on the
genotype set, and both arms read the same `chr*_deduplicated` files, so the masks are
identical for broad and bilateral. Building them twice would only create a chance for
them to differ.

`ALL` is not built. An unfiltered burden test answers no question this project asks, and
it was dropped from the analysis in the first pass, so building it would leave a mask on
disk that nothing tests. Per-gene variant counts are emitted separately, which is the
only thing it was serving.

## It reproduces the first pass exactly

Not merely the same counts — the same files:

| mask | genes | entries | content |
|---|---|---|---|
| pLOF | 17,841 = 17,841 | 462,144 = 462,144 | md5 identical |
| pDM | 17,623 = 17,623 | 890,332 = 890,332 | md5 identical |
| pLOF_pDM | 17,945 = 17,945 | 1,340,936 = 1,340,936 | md5 identical |

The first pass read the 106 GB serially; this reads it as 22 parallel tasks and merges.
Getting byte-identical output from a restructured implementation is the evidence that
the restructuring introduced nothing.

## The three fixes, measured

Each corrects a defect the first pass established in the original run. All three are
measurably active:

| | count | what it fixes |
|---|---:|---|
| rows dropped, not `protein_coding` | 71,443,642 | the original's gene set carries 651 antisense, 88 LINC and 14 microRNA genes. `pLOF` means predicted loss of function; a gene with no protein has no function to lose. Its top burden gene in two of three cohorts was `TMC3-AS1`, an antisense lncRNA. |
| splice annotations rejected by the SpliceAI gate | 2,780,725 | the original matched `"splice"` as a substring anywhere in the consequence string, admitting three IMPACT=LOW terms unconditionally. Roughly 64% of its pLOF mask fell outside its own documented definition. The `SpliceAI >= 0.2` gate its analysis plan specifies was never applied. |
| pDM assignments from REVEL | 1,366,320 | the original coerced VEP's multi-valued REVEL entries to NaN, silently dropping roughly three quarters of pDM. Found by Nikki Palmiero on 2026-10-05, after the first pass had missed it. |

For scale: the gate rejects 2,780,725 splice candidates and admits 426,418 — it is
removing about 87% of what substring matching would have let through.

## The rules, in full

```
pLOF      exact term match on {frameshift_variant, stop_gained, start_lost,
          stop_lost, splice_acceptor_variant, splice_donor_variant,
          transcript_ablation}
          OR any other splice annotation WITH SpliceAI >= 0.2

pDM       am_class in {pathogenic, likely_pathogenic}
          OR max REVEL across transcripts >= 0.5

pLOF_pDM  the union, pLOF taking precedence in the annotation column

all three AND BIOTYPE == protein_coding
```

Terms are matched against the consequence list and **not** against VEP's `IMPACT`
column, which is per row: a row reading `stop_gained,splice_polypyrimidine_tract_variant`
is HIGH as a whole, so filtering on IMPACT would readmit exactly what is being excluded.

REVEL is read as the per-transcript **list** VEP writes, and the maximum taken. That is
the whole of the second fix.

## What the merge refuses on

`02_merge_masks.py` will not write a partial mask, because a short mask does not fail
downstream — SAIGE tests whatever genes it is handed and the result looks finished.

| check | why |
|---|---|
| LSF asked, not mtimes trusted | a task still running has already created its output file |
| alternating `var`/`anno` lines of equal length | a truncated write leaves a gene with more variants than annotations |
| entry counts match each chromosome's manifest | catches a partial read |
| no gene on two chromosomes | a symbol appearing twice is a name collision or an annotation problem, and merging the two would pool unrelated variants into one burden test |
| `pLOF_pDM` is exactly the union of `pLOF` and `pDM` | cheap, and if the union were built wrong every `pLOF_pDM` test would be subtly incorrect with nothing else noticing |

All passed.

## Files

| | |
|---|---|
| `results/masks/{pLOF,pDM,pLOF_pDM}.txt` | the masks, SAIGE group-file format |
| `results/02_masks_manifest.json` | counts and the fix statistics |
| `data/by_chrom/` | per-chromosome partials and manifests, regenerable |
| `data/by_chrom/gene_variant_counts.chr*.tsv` | per-gene variant totals, for QC |
