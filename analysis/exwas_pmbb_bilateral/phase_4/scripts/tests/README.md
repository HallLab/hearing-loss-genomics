# Design test — one call over the grid, against nine calls over its cells

Run 2026-10-07, before committing this analysis to the one-call design (premise **P10**).
Script: `native_grid_chr21.bsub`. Output: `../results/test_native/`.

Deliberately run against the **replication's** null model and chr21 results, so the comparison
isolates the job structure from the phenotype change. The question is about SAIGE's behaviour, not
about this cohort.

## What came back

| comparison | genes | identical |
|---|---:|---|
| `pLOF:pDM` — same group file either way | 623 | **623 (100%)** |
| `pLOF` — read from pLOF_pDM vs its own file | 620 | **620 (100%)** |
| `pDM` — read from pLOF_pDM vs its own file | 522 | 381 (73%) |

One call produced 1,973 rows against 1,765 from the nine — the extra 208 being **one Cauchy row per
gene**, for all 208 genes on chr21. Matched against the ACAT computed by hand over the nine cells:
165 of 208 identical, median ratio 1.000.

## Why only pDM moves

Of the 11,517 variants our masks label both pLOF and pDM, **11,516 are written as `pLOF`** in the
combined file. Label precedence, not a bug. So `pDM` read out of `pLOF_pDM` means "damaging missense
that is not also loss-of-function", while the standalone file means "damaging missense, including
those".

These are variants that truncate the protein on one transcript and are missense on another. Either
reading is defensible, and the native one is arguably cleaner — a variant that truncates is not a
missense story. What matters is that a `pDM` number from this analysis and a `pDM` number from the
replication are answering slightly different questions, and should not be set side by side without
saying so.

## Verdict

Adopt the one-call design. 66 tasks instead of 594, the headline p-value comes from the tool rather
than from a calculation that runs afterwards, and the only thing that moves is a mask whose
definition becomes the sharper of the two.
