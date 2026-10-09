# Pre-registered premises — Common-variant ExWAS — broad hearing loss (replication)

**Written before step 2 was run.** Nothing in this file was chosen after seeing a p-value.
Last edited: 2026-10-09.

This arm is the **single-variant, common-variant** half of Elena's pipeline:

```
ExWAS                     MAF >= 0.01                 -> SAIGE          <- THIS ARM
Rare Variant Gene Burden  MAF <= 0.01, 0.001, 0.0001  -> SAIGE-GENE     <- ../gene_burden_elena_replication/
```

The two arms test **disjoint variant sets**. A gene can appear in one and not the other
without any inconsistency.

---

## Why this arm is being run at all

Douglas Epstein is on record in the 2026-10-02 meeting deprioritising it:

> *"just with respect to the common variants, I wasn't excited about any of those genes."*
>
> *"if the variants that associated were not pathogenic, then I didn't understand the value
> of that [...] from a biological standpoint, it's unclear what that would mean."*
>
> *"that was just contributing noise. But then when she did the rare variant, there was at
> least some potential signal there that made sense to me from a biological standpoint."*

It is run anyway for one reason: **it is the only arm with a declared target.** Elena's
summer poster reported `KNL1` and `ZFPM1` from this analysis at a 5e-8 threshold. Every
other arm we have run had no prior hypothesis to confirm or refute. This one can actually
fail, which makes it the only arm where replication means something falsifiable.

The biological reservation above stands regardless of the outcome, and belongs in any
write-up of the result.

---

## P1 — Cohort: inherited, not recomputed

Cohort inherited from the replication arm: the broad hearing-loss phenotype,
after the two cohort defects found in round 1 were corrected.

| cohort | samples | cases | PCs |
|---|---:|---:|---:|
| combined | 57,498 | 6,752 | 5 |
| EUR | 42,779 | 5,183 | 4 |
| AFR | 11,334 | 1,285 | 3 |

Source: `../gene_burden_elena_replication/phase_3/results/covariates/`
Not rebuilt here. If that folder changes, this arm is stale.

## P2 — Covariates: inherited, including our deviation

`AGE, AGE2, SEX, Batch, PC1..PCn`, with `SEX` and `Batch` declared categorical
(`--qCovarColList`).

PC counts are **ours** (scree-supported: 5/4/3), not Elena's 5/9/10. This is the same
deviation already documented for the burden arm, carried forward deliberately so the two
arms remain comparable to each other.

## P3 — Step 1 null model: reused, not refitted

`../gene_burden_elena_replication/phase_4/results/step1/{combined,EUR,AFR}/`

The null model does not depend on the variant filter — same phenotype, same covariates,
same GRM. Refitting it would produce the same numbers at the cost of the expensive
relatedness computation.

**Elena did the same thing**, which is confirmation rather than assumption:
`saige_hearing_pipeline.config` sets
`step1_outdir = "${base}/HL_only_rarevariant/saige_step1"`.

Note the filename as written by our step 1: `{GROUP}.varianceRatio.txt.varianceRatio.txt`
(doubled suffix — a prefix was passed where SAIGE wanted a stem). Harmless, but scripts
must use the literal path.

## P4 — Genotypes: the same files, unchanged

`analysis/elena/rarevariant_geneburden/plink_deduplicated/chr{1..22}_deduplicated`

These are the files our burden arms already read, and the files Elena's own ExWAS read.
Genotype call rate > 0.99 was applied when they were built. Nothing is rebuilt here, so
there is no opportunity for this arm to diverge from the burden arm on genotypes.

## P5 — Variant filter

| | |
|---|---|
| `--minMAF` | `0.01` — **inclusive**, MAF >= 0.01 |
| `--minMAC` | `1` (moot: MAF >= 0.01 in ~50k samples implies MAC ~1000) |
| `--LOCO` | `FALSE` |

Mirrors `min_maf = 0.01`, `min_mac = 1`, `loco = false` from her config. Her inline comment
reads *"only variants with MAF >= min_maf are tested"* — so the `MAF > 0.01` on the pipeline
diagram is a rounding of `>=`.

## P6 — Significance threshold: fixed now, two bars

**Primary: `5e-8`.** This is a replication, so it is judged at the original's bar. It is also
the conservative choice of the two.

**Secondary: `0.05 / n_variants_tested`**, counted per cohort from the phase 2 census.
For exome-only testing this lands near 1e-7, i.e. **less strict** than 5e-8. Reported
alongside, never instead.

Fixing both before running is the point. The choice decides whether `KNL1` and `ZFPM1`
"replicate", and that decision must not be made after seeing where their p-values land.

## P7 — The declared targets

`KNL1` and `ZFPM1`, from the summer poster. Their p-values in this arm will be reported
**whatever they are**, in their own named section, not left to be found in a top-50 table.
A target that fails to replicate is a result.

## P8 — What this arm does not have

No masks. No functional annotation. No Cauchy combination. No MAF nesting. The whole
cell structure of the burden arm's phase 4 collapses to one test per variant.

One test per variant per cohort, 3 cohorts x 22 chromosomes = **66 step 2 tasks**.

## P9 — Functional annotation of whatever clears a bar

Registered **before** step 2 runs, so that annotating the hits is a planned step and
not a response to the hits being unconvincing.

Every variant clearing either bar (P6) gets annotated and reported next to its p-value:
VEP consequence, affected gene and transcript, and AlphaMissense score via ANNOVAR
`dbnsfp45a_canonical`. The same tooling the burden arms used — see
`reference_vep_annovar_lpc` and Biofilter 4 on LPC.

**Why this is not optional here.** In the burden arm, annotation *builds* the test:
SAIGE-GENE tests genes, and a gene test needs a rule for which variants enter it, which
is what the `pLOF` and `pDM` masks are. This arm has no grouping, so annotation decides
nothing about what is tested — every variant is tested individually whatever its
consequence class.

That is precisely the gap the 2026-10-02 meeting identified. Douglas Epstein:

> *"if the variants that associated were not pathogenic, then I just... I didn't
> understand the value of that [...] from a biological standpoint, it's unclear what
> that would mean."*

Nikki Palmiero, on the original ExWAS results:

> *"these all mostly focused on common variants. And for the most part, we don't really
> see much of the damaging there anyway."*

A common-variant ExWAS produces statistical association without saying whether the
variant does anything. Annotating the hits turns that from an unanswerable question in
a meeting into a column in a table. If the associated variants are not damaging, that
is the finding, and it is reported as one.

**Pre-registered now so it is not a choice made later:** the annotation is reported for
every variant that clears a bar, including ones where it weakens the story.

---

## Phase map

| phase | status |
|---|---|
| 1 — cohort | **inherited** from `../gene_burden_elena_replication/phase_1/` |
| 2 — variants | MAF census only; no build |
| 3 — covariates | **inherited** from `../gene_burden_elena_replication/phase_3/` |
| 4 — test | SAIGE single-variant, step 1 reused |
| 5 — results | this arm, annotation of hits (P9), four-way comparison |
