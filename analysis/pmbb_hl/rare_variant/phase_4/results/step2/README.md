# Every cell of the burden arm, in one place

Built by [`../../scripts/05_cell_table.py`](../../scripts/05_cell_table.py). The two
files here are a reshape of `cells_*.tsv` and `omnibus_*.tsv`, not a new analysis —
nothing is recomputed except the two columns marked below, and the script **refuses to
write** if its numbers do not reproduce the shipped ones exactly.

| file | rows | what it is |
|---|---:|---|
| `all_cells.tsv` | 1,415,668 | one row per (phenotype, cohort, gene, mask, max_MAF) |
| `genes.tsv` | 160,940 | one row per (phenotype, cohort, gene) |

Both are sorted by `omnibus_rank` within each arm, so the genes worth looking at are at
the top. `*.tsv` is gitignored — these live on disk, not in the repository.

**The nine subdirectories here are SAIGE's own raw output**, one per arm, 22 chromosome
files each (`<arm>_chr<N>.txt`, its `.index`, and `.singleAssoc.txt`). They are the
input `04_merge_results.py` read; the two tables above are the readable form of the
same thing. 1.1 GB in total, so prefer the tables unless you need a specific chromosome
file.

## Columns

| | |
|---|---|
| `phenotype` · `cohort` | one of broad / bilateral / tinnitus × combined / EUR / AFR |
| `gene` · `chr` | the gene symbol as it appears in the mask |
| `mask` | `pLOF`, `pDM`, or `pLOF;pDM` |
| `max_maf` | 0.0001, 0.001, 0.01 — the frequency ceiling for that cell |
| `p_cell` | SAIGE-GENE+ SKAT-O p for that cell, after its own ρ sweep |
| `p_burden` · `p_skat` | the two component tests, for diagnosis |
| `beta_burden` · `se_burden` | log-odds and its error, from the burden component |
| `mac` · `mac_case` · `mac_control` | cumulative minor allele counts in that cell |
| `n_rare` · `n_ultra_rare` | variants above and below SAIGE's ultra-rare collapsing threshold |
| `in_omnibus` | did this cell enter the combination — see below |
| `n_cells_identical` | how many combined cells carry this exact result |
| `is_representative` | the tightest cutoff of each identical group |
| `p_omnibus` | **the gene's verdict.** ACAT over the `in_omnibus` cells |
| `p_omnibus_distinct` | recomputed: ACAT over the `is_representative` cells only |
| `omnibus_rank` | rank by `p_omnibus` within the arm, 1 = best |
| `bonferroni_bar` | `0.05 / genes tested in that arm` |
| `fdr_q` · `fdr_q_distinct` | Benjamini-Hochberg q for each omnibus |

`genes.tsv` adds `best_mask`, `best_max_maf`, `best_cell_p` (where the gene's strongest
cell was — ties broken toward the **tightest** cutoff), `cells`, `cells_used`,
`distinct_tests` and `mac_widest`.

## Two kinds of repeated cell, and they are different

A gene can have the same result reported more than once for two unrelated reasons. Both
are flagged, and confusing them is easy.

**1. Same cutoff, different mask** — `in_omnibus = False`.
When a gene has no pDM variants, `pLOF` and `pLOF;pDM` score the identical variant set
and return the identical p. This is the case Nikki raised on 2026-10-14. Phase 4 already
removes it before combining; the flag here just shows which rows those were. **This is
applied in everything the project reports.**

**2. Same mask, looser cutoff** — `is_representative = False`.
When a gene has no variants between two cutoffs, the looser cell contains exactly the
variants of the tighter one and returns the identical p. `FBXO25` is the clear example:
16 pLOF variants, all with MAC ≤ 8, so its three pLOF cells are one test counted three
times.

This second kind is **not** removed in what the project reports. It is not a bug — ACAT
admits arbitrary non-negative weights, so `p_omnibus` is a valid p-value either way. But
the implicit weighting is *"how many cutoffs did this test survive"*, which was never a
stated choice and systematically favours genes whose signal is ultra-rare.
`p_omnibus_distinct` is the same combination with each distinct test counted once.

### What the difference is worth

| arm | min q, as reported | min q, distinct tests |
|---|---:|---:|
| broad combined | 0.198 | 0.198 |
| broad EUR | 0.477 | 0.446 |
| broad AFR | 0.248 | 0.248 |
| bilateral combined | 0.239 | 0.247 |
| **bilateral EUR** | **0.113** | **0.188** |
| bilateral AFR | 0.259 | 0.344 |
| tinnitus combined | 0.343 | 0.411 |
| tinnitus EUR | 0.712 | 0.801 |
| tinnitus AFR | 0.496 | 0.608 |

**No arm crosses any bar under either, in any cohort** — the choice changes no
conclusion. It changes one sentence: `FBXO25` is "13% short of the FDR bar" as reported
and "88% short" with distinct tests. Since that gene is the closest the project gets,
the sentence is worth getting right, and the decision is open — see the open question in
[`../FINDINGS.md`](../FINDINGS.md).

## Worked example — the question this table was built to answer

> *Where was `FBXO25` found — which MAF?*

```
awk -F'\t' 'NR==1 || ($1=="bilateral" && $2=="EUR" && $3=="FBXO25")' all_cells.tsv
```

| mask | max_maf | p_cell | mac | case | control | n_rare | n_ultra_rare | in_omnibus | repr. |
|---|---:|---:|---:|---:|---:|---:|---:|:--:|:--:|
| pLOF | 0.0001 | **2.11e-06** | 18 | 7 | 11 | 0 | 16 | ✓ | ✓ |
| pLOF | 0.001 | 2.11e-06 | 18 | 7 | 11 | 0 | 16 | ✓ | |
| pLOF | 0.01 | 2.11e-06 | 18 | 7 | 11 | 0 | 16 | ✓ | |
| pLOF;pDM | 0.0001 | 4.19e-04 | 38 | 8 | 30 | 0 | 30 | ✓ | ✓ |
| pLOF;pDM | 0.001 | 3.09e-02 | 146 | 19 | 127 | 2 | 31 | ✓ | ✓ |
| pLOF;pDM | 0.01 | 3.09e-02 | 146 | 19 | 127 | 2 | 31 | ✓ | |
| pDM | 0.0001 | 0.829 | 20 | 1 | 19 | 0 | 14 | ✓ | ✓ |
| pDM | 0.001 | 0.868 | 128 | 12 | 116 | 2 | 15 | ✓ | ✓ |
| pDM | 0.01 | 0.868 | 128 | 12 | 116 | 2 | 15 | ✓ | |

The answer: **`pLOF`, in the ultra-rare band, and the MAF cutoff does not matter** —
`n_rare = 0` means no pLOF variant sits above SAIGE's ultra-rare threshold, so raising
the ceiling adds nothing. `pDM` alone is null (p = 0.83) and adding it to `pLOF` only
dilutes, 4.19e-04 then 0.031.

And the thing the table makes impossible to miss: **the whole result is 7 case alleles
against 11 control alleles**, across 16 near-singleton variants in 2,547 cases and
37,596 controls. `beta_burden` of 3.71 is an odds ratio of 41; a Fisher test on the same
counts gives 9.4. That gap between the two is what an effect estimated on seven alleles
looks like.
