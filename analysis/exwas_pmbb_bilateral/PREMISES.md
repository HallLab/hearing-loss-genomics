# Pre-registered premises — Common-variant ExWAS — bilateral sensorineural hearing loss

**Written before step 2 was run.** Nothing in this file was chosen after seeing a p-value.
Last edited: 2026-10-09.

This arm is the **single-variant, common-variant** half of Elena's pipeline:

```
ExWAS                     MAF >= 0.01                 -> SAIGE          <- THIS ARM
Rare Variant Gene Burden  MAF <= 0.01, 0.001, 0.0001  -> SAIGE-GENE     <- ../gene_burden_pmbb_bilateral/
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

Cohort inherited from the bilateral arm: sensorineural AND bilateral recorded
on the same date, on two distinct dates.

| cohort | samples | cases | PCs |
|---|---:|---:|---:|
| combined | 53,910 | 3,164 | 5 |
| EUR | 40,143 | 2,547 | 4 |
| AFR | 10,539 | 490 | 3 |

Source: `../gene_burden_pmbb_bilateral/phase_3/results/covariates/`
Not rebuilt here. If that folder changes, this arm is stale.

## P2 — Covariates: inherited, including our deviation

`AGE, AGE2, SEX, Batch, PC1..PCn`, with `SEX` and `Batch` declared categorical
(`--qCovarColList`).

PC counts are **ours** (scree-supported: 5/4/3), not Elena's 5/9/10. This is the same
deviation already documented for the burden arm, carried forward deliberately so the two
arms remain comparable to each other.

## P3 — Step 1 null model: reused, not refitted

`../gene_burden_pmbb_bilateral/phase_4/results/step1/{combined,EUR,AFR}/`

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

**Primary: `5e-8`.** This arm is *not* a replication — the bilateral phenotype is new, so
there is no original p-value to be judged against. 5e-8 is adopted anyway, for two reasons:
it is the convention, and it keeps this arm numerically comparable to
`../exwas_elena_replication/`, which is the comparison that matters.

**Secondary: `0.05 / n_variants_tested`**, counted per cohort from the phase 2 census.
For exome-only testing this lands near 1e-7, i.e. **less strict** than 5e-8. Reported
alongside, never instead.

Fixing both before running is the point, and it is fixed here at the same values as the
replication arm so that a difference between the two arms is a difference in phenotype
rather than in bookkeeping.

## P7 — The declared targets

`KNL1` and `ZFPM1`, from the summer poster — but with a caveat this arm must state
plainly: they were found on the **broad** phenotype, so this arm cannot replicate them.
Restricting to bilateral sensorineural drops 3,588 cases, and a target that weakens here
is confounded between "the association was phenotype-specific" and "we lost the power to
see it". The power control in `../gene_burden_pmbb_bilateral/phase_5/power_control/`
already showed that confound is real for the burden arm.

Their p-values are reported **whatever they are**, in their own named section, with that
caveat attached. The arm that can actually replicate them is
`../exwas_elena_replication/`.

## P8 — What this arm does not have

No masks. No functional annotation. No Cauchy combination. No MAF nesting. The whole
cell structure of the burden arm's phase 4 collapses to one test per variant.

One test per variant per cohort, 3 cohorts x 22 chromosomes = **66 step 2 tasks**.

---

## Phase map

| phase | status |
|---|---|
| 1 — cohort | **inherited** from `../gene_burden_pmbb_bilateral/phase_1/` |
| 2 — variants | MAF census only; no build |
| 3 — covariates | **inherited** from `../gene_burden_pmbb_bilateral/phase_3/` |
| 4 — test | SAIGE single-variant, step 1 reused |
| 5 — results | this arm, plus the four-way comparison |
