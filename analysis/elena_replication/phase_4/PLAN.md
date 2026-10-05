# Phase 4 — the association test

**Status:** step 1 submitted 2026-10-05 (job 49231968) · **Arm:** corrected only
**Scripts:** `01_split_masks_by_chrom.sh` · `02_saige_step1.bsub` · `03_saige_step2.bsub`

Phase 4 covers the pipeline's Step 8: SAIGE-GENE+ step 1 (null model) and step 2 (gene burden
association). It is the first phase that produces p-values rather than inputs.

---

## 1. One arm, not two

Phases 1–3 maintained two arms in parallel: a **reproduction** arm (the cohort and masks the
pipeline actually delivered, so divergences could be attributed) and a **corrected** arm (the
same pipeline with its four mask defects and three cohort defects repaired).

Phase 4 runs **only the corrected arm**: our cohort of 57,507 analysable, 57,498 with age.

The reproduction arm's inputs are all still on disk — `phase_1/results/` holds the 57,632 cohort,
`phase_2/results/masks/` the original masks — so it stays runnable if the lab wants the
side-by-side later. What we do not do is spend 594 more jobs reproducing numbers we already know
to be built on a mask where 63.8% of pLOF entries fail the pipeline's own SpliceAI gate.

---

## 2. What goes in

| | combined | EUR | AFR |
|---|---:|---:|---:|
| samples | 57,498 | 42,779 | 11,334 |
| cases | 6,752 | 5,183 | 1,285 |
| controls | 50,746 | 37,596 | 10,049 |
| PCs | 5 | 4 | 3 |

PC counts are ours, from the scree, not the pipeline's 5/9/10 — see
[`../docs/10_fase_3_resumo.pt.md`](../docs/10_fase_3_resumo.pt.md).

**GRM genotypes (step 1).** EUR and AFR use the LD-pruned sets we built in Phase 3
(53,920 and 114,694 markers). The combined cohort uses Elena's
`exome_ldpruned/ALL.common_ldpruned`, because Phase 3 only re-ran the two within-ancestry PCAs and
we never built a combined set. Borrowing hers is safe on the point that matters: it was made with
`--keep HL_TIN_PMBBv4_keep.txt`, all 70,925 people, so it does **not** inherit the imputed-`.fam`
defect that Phase 1 found in Step 1, and it covers all 57,498 of our combined cohort with zero gaps.

**Burden genotypes (step 2).** Elena's `rarevariant_geneburden/plink_deduplicated/chr<N>_deduplicated`,
reused rather than rebuilt. All 70,925 are present, so all three cohorts are covered with zero gaps,
and all 47,818 chr1 pLOF variant IDs from our rebuilt masks resolve against the chr1 `.bim` — the
`chr:pos:ref:alt` convention matches ours exactly.

**Masks.** Three of the four rebuilt in Phase 2, `phase_2/results/masks_v2/`, split per chromosome
by script 01. Gene counts survive the split exactly: pLOF 17,841 · pDM 17,623 · pLOF_pDM 17,945,
and every gene is autosomal, so chromosomes 1–22 lose nothing.

**`ALL` is not run.** Andre's call, and the reasoning holds: a burden test over every variant in the
gene pools synonymous and intronic variants that carry no mechanism in with the ones that do, so it
dilutes rather than tests anything. It is also by far the most expensive mask — 15,951,289 entries
against 1,340,936 for `pLOF_pDM`, which is where most of the 2.1 GB and most of the CPU would have
gone. The group files stay on disk, in `phase_2/results/masks_v2/ALL.txt` and
`phase_4/data/masks_by_chrom/ALL/`, so the decision is reversible without rebuilding anything.
Consequence for Phase 5: Elena ran `ALL`, and we will have nothing to set against those results.

---

## 3. Parameters

Taken from Elena's `saige_hearing_pipeline.config` unless marked.

```
step 1   --traitType=binary  --nThreads=8
         --covarColList=AGE,AGE2,SEX,Batch,PC1..PCn
         --qCovarColList=SEX,Batch                        <- deviation, see 4.1
step 2   --minMAF=0  --minMAC=1  --LOCO=FALSE
         --maxMAF_in_groupTest = 0.0001 | 0.001 | 0.01
         --annotation_in_groupTest = one per mask          <- deviation, see 4.2
         --is_Firth_beta=TRUE  --is_single_in_groupTest=FALSE
         --is_output_moreDetails=TRUE  --is_no_weight_in_groupTest=TRUE
```

SAIGE 1.5.0, `/project/hall/tools/saige/1.5.0/saige_1.5.0.sif`, run under
`apptainer/1.4.1`. Elena ran the same container under `DEV/singularity` via Nextflow; we run it as
two LSF job arrays instead, which is the same calls in a form that can be audited line by line.

Scale: 3 cohorts × 22 chromosomes × 3 masks × 3 MAF cutoffs = **594 step-2 tasks**. The index decode
in script 03 was checked to generate all 594 combinations exactly once.

---

## 4. Deviations from the pipeline, deliberate

### 4.1 Batch declared categorical

Her `--covarColList` includes `Batch` but she passes no `--qCovarColList`, so SAIGE reads the
values 1/2/3 as a continuous variable — it assumes batch 2 sits exactly halfway between batch 1 and
batch 3. Batch is a label, not a quantity. We declare it categorical so SAIGE dummy-codes it.

Batch is not evenly spread, which is what makes the linearity assumption bite: in AFR the three
batches hold 8,492 / 686 / 2,156 people.

`SEX` is listed alongside it. With two levels this changes nothing numerically, but it is
categorical and saying so costs nothing.

### 4.2 One annotation per mask

Her step 2 passes `--annotation_in_groupTest='ALL,pLOF,pDM,pLOF:pDM'` for every mask, against a
group file that contains one annotation label. Two consequences: SAIGE looks for categories the
file does not hold, and the `pLOF_pDM` run re-emits pLOF-only and pDM-only tests that the
standalone `pLOF` and `pDM` runs already produce. We pass the single category each file holds, so
each mask yields exactly one test per gene.

The duplicate tests are not quite identical, which is why this is worth stating rather than
silently fixing: 11,540 variants are both pLOF and pDM, and in `pLOF_pDM.txt` each carries one
label. So pLOF-within-pLOF_pDM and standalone pLOF differ for those variants.

### 4.3 Group files split by chromosome

Mechanical. SAIGE step 2 runs one chromosome but reads the whole group file, so each mask would be
re-parsed 198 times. No gene line in any mask spans two chromosomes, so the split cannot change
which genes are tested. The split also converts tab to space, matching the
format Elena's files use and that this container has demonstrably accepted.

---

## 5. Carried to Phase 5

**The combined GRM is thin.** 38,833 markers, against 53,920 for our EUR and 114,694 for our AFR.
The cause is in her log: `--hwe 1e-6` applied across ancestries removes 176,135 variants where the
same filter within EUR removes 18,691. Hardy–Weinberg departures across a structured sample are
expected and are not evidence of genotyping error, so the filter is doing the wrong job here. Below
the ~100k markers usually wanted for a stable variance ratio. Flagged, not fixed — rebuilding it
would mean deciding how to handle HWE in a multi-ancestry sample, which is a question for the lab.

**No X chromosome.** Both arms are autosomes only. Our Phase 2 rebuild covered chromosomes 1–22
because the release VEP sweep did, and Elena's config lists `chromosomes = 1..22`. Consistent
between arms, but it means a hearing-loss study has no coverage of the X — worth raising, since
X-linked deafness genes exist.
