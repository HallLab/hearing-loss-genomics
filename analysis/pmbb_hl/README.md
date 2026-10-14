# pmbb_hl — hearing loss in PMBB v4, rebuilt

**Opened 2026-10-14.** A clean rebuild of the analysis, carrying forward what the first
pass established and fixing what it got wrong.

```
pmbb_hl/
├── rare_variant/   SAIGE-GENE+, MAF <= 0.01 / 0.001 / 0.0001
├── exwas/          SAIGE single variant, MAF >= 0.01
└── synthesis/      both arms, both phenotypes, one place
```

Each analysis runs **two phenotypes** from the same code:

| | rule |
|---|---|
| **broad** | `SO_396` or any child, on ≥ 2 distinct dates. The phenotype the original run used. |
| **bilateral** | `SO_396.2` and `SO_396.8` on the **same date**, on ≥ 2 distinct dates. |

---

## What this rebuild changes, and why

### 1. PCs come from the release, not from a PCA of our own

The first pass re-ran a within-ancestry PCA. The stated reason was that the original
pipeline's PCA had a coverage gap. **That reason was wrong**, and it was wrong in a
specific, checkable way:

- The original run computed **no PCA at all**. Its PCs are byte-identical to the
  release's `exome_PC*` columns — `max|difference| = 0.00e+00` across all 24 components
  in the three cohorts.
- What its `step7_ldprune.bsub` builds is the **GRM**, not a PCA. There is no `--pca`
  anywhere in its scripts and no eigenvector file in its output.
- The release's by-ancestry exome PCA covers **100%** of the corrected cohort — zero
  people missing, in either ancestry.

So a GRM-building step was read as a PCA step, and a PCA was recomputed to close a gap
that did not exist there. That recomputation is the origin of every later difficulty
with PCs: it created the only real methodological divergence from the original run, made
the eigenvalue profiles incomparable, and is the most likely source of the MHC inflation
the first pass found (λ = 1.34 in EUR against the original's 0.97).

**This rebuild uses `exome_PC*` from the release**, as the original did.

### 2. The GRM is still rebuilt, because that gap is real

196 EUR and 34 AFR members of the corrected cohort are absent from the original's
LD-pruned files, which were built over a keep list that excluded them. That is a genuine
gap and rebuilding is the right response. It is a separate thing from the PCs, and the
first pass conflated the two.

### 3. Both phenotypes are built by one script

So the control set is shared by construction rather than by coincidence, and the nesting
the later phases rely on — every bilateral case is also a broad case — is asserted in
code rather than assumed.

### 4. Analysable and tested are reported separately from the start

Phase 1 defines who is *analysable*. Phase 3 drops anyone missing a covariate. In the
first pass those two numbers (57,507 and 57,498) drifted into different documents and
had to be reconciled afterwards. Both are emitted here, at every phase.

---

## What carries over unchanged

These were checked in the first pass and do not need redoing:

| | |
|---|---|
| tinnitus lives in OMOP `observation` | the phecode files hold ~3,016 of 28,110 events; a control set built without that source contains tinnitus patients |
| the bilateral rule needs same-date pairing | the loose reading gives 4,029 cases against 3,164 |
| masks exclude non-coding genes | the original's set carries 651 antisense, 88 LINC and 14 microRNA genes; `pLOF` presupposes a protein |
| both bars, fixed before running | Bonferroni and, on the burden arm, FDR at q < 0.1 — asked for on 2026-10-14 |
| the poster's variants fail its own QC | call rate 0.20 / 0.96 / 0.95 against a stated floor of 0.99 — an open question for the group, not for this rebuild |

The first pass is kept at `analysis/gene_burden_*` and `analysis/exwas_*` and is not
deleted: its findings are what this rebuild is built on.
