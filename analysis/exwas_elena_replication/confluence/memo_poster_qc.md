# ExWAS replication — one QC question before the next meeting

**Andre Rico · 9 October 2026 · for Nikki**

---

## Short version

I finished replicating the common-variant ExWAS. **Our results and the files you sent
agree almost exactly** — r = 0.98–0.99 per variant, 99.9% of the same variants tested,
λ within 0.01 of each other. Neither run finds anything significant, and your files'
own `Significance` column already says so.

One thing I could not reconcile, and I don't think it's an analysis difference: the
three variants behind the poster's six genome-wide significant results **fail the
quality filter the poster itself describes.** I'd like to check with you which variant
set produced that table before we present any of this.

---

## 1. The replication

| cohort | your variants | ours | shared | Pearson r | Spearman ρ | λ yours | λ ours |
|---|---:|---:|---:|---:|---:|---:|---:|
| combined | 70,429 | 70,366 | 70,333 (99.9%) | 0.989 | 0.978 | 0.992 | 0.989 |
| EUR | 74,528 | 74,536 | 74,483 (99.9%) | 0.979 | 0.962 | 0.991 | 0.999 |
| AFR | 117,943 | 117,953 | 117,723 (99.8%) | 0.990 | 0.981 | 1.012 | 1.021 |

Minimum p-values: yours 2.57e-06 / 1.99e-05 / 3.70e-06, ours 3.02e-06 / 2.12e-05 /
4.11e-06. Nothing within two orders of magnitude of 5e-8 in either run.

The few hundred variants that appear on only one side come from the cohort correction —
ours has 40 more cases and 174 fewer controls, which moves a handful across the 1% MAF
line in both directions. Every per-variant disagreement larger than one log unit is in
the MHC: 14 of 14 in EUR, 16 of 16 in AFR.

**So the pipeline reproduces.** That part is settled.

## 2. The poster's six results

From `Horizontal Draft LANS Poster - SUIP-PGSIP.pptx`, under *"Genome Wide Significant
Results Total: 6"* — three unique variants:

| gene | variant | cohort | p |
|---|---|---|---|
| KNL1 | chr15:40611486:A:G | Combined | 3.51e-12 |
| KNL1 | chr15:40611486:A:G | EUR | 7.21e-12 |
| ZFPM1 | chr16:88533289:A:C | EUR | 2.46e-10 |
| ZFPM1 | chr16:88533290:G:C | EUR | 2.52e-10 |
| ZFPM1 | chr16:88533290:G:C | Combined | 4.54e-09 |
| ZFPM1 | chr16:88533289:A:C | Combined | 6.18e-09 |

None appears in your delivered files, in the raw SAIGE output under
`HL_only_ExWAS/saige_results/`, or in ours.

**They do exist in the PMBB v4 release.** Measured there directly with plink2 — this
number comes from the institutional release, not from either of our pipelines:

| variant | call rate | MAF |
|---|---:|---:|
| chr15:40611486:A:G | **0.2013** | 0.284 |
| chr16:88533289:A:C | 0.9617 | 0.295 |
| chr16:88533290:G:C | 0.9473 | 0.283 |

The poster's methods section reads: *"variants with a genotype call rate < 0.99 and a
Minor Allele Frequency (MAF) < 0.01 were removed."* All three are below 0.99. The
variant carrying the 3.51e-12 has a genotype call in **20% of participants** — 56,646 of
70,925 are missing.

For scale: the call rate of everything this arm tested runs from 0.9534 to 1.0, median
0.9999.

## 3. Where they drop out — and one thing that may be worth a look

| stage | chr15 | chr16 |
|---|---|---|
| PMBB v4 release | present (546,760) | present (780,404) |
| `rarevariantExWAS/variant_qc/chrN_callrate.bim` | gone (476,227) | **present** (670,730) |
| `plink_deduplicated/chrN_allele_ids.bim` | gone (469,424) | gone (662,670) |
| `chrN_deduplicated.bim` | gone (469,422) | gone (662,667) |

Applying a true `--geno 0.01` filter to the release leaves **475,376** variants on
chr16. `chr16_callrate.bim` holds **670,730**. So the stage named for call rate appears
not to be filtering at 0.99 — which is how two variants at 0.96 and 0.95 passed through
it. They were eventually dropped at the `allele_ids` step, for what looks like a
different reason.

It did catch the chr15 variant at 0.20, so some filter is running there; it just seems
to be looser than the one described.

## 4. Something of ours that is worse than yours

While splitting λ by region I found a cost of our PC choice that I should flag rather
than leave for someone else to find:

| cohort | λ outside MHC (yours / ours) | λ inside MHC (yours / ours) |
|---|---|---|
| combined | 0.990 / 0.989 | **1.145** / 1.041 |
| EUR | 0.992 / 0.995 | 0.965 / **1.339** |
| AFR | 1.012 / 1.020 | 0.998 / **1.092** |

Outside the MHC we are indistinguishable. Inside it, in EUR and AFR, **our 4 and 3 PCs
leave inflation that your 9 and 10 do not.** We chose ours from scree plots across the
genome and I still think that is right genome-wide — outside the MHC it costs nothing —
but in the window where population structure is strongest, your choice is better. In
combined the direction reverses.

No conclusion depends on it (λ = 1.34 on 1,032 variants does not produce a 5e-8 hit),
but it should not be inherited unexamined by anything that cares about the MHC.

## 5. What I'd like to check with you

**Which variant set produced the poster's results table?** It is not
`rarevariant_geneburden/plink_deduplicated/`, which is what the delivered pipeline
reads.

My guess is chronology rather than anything else: the poster was prepared over the
summer, the QC was tightened afterwards, and the files you sent on 2 October are the
later run. If that's right, the sequence is normal and the only thing needed is that we
don't carry the earlier numbers into the next presentation — I'm a co-author on that
poster, so this is my name on it too.

Two things I can do once you confirm, whichever is more useful:

- re-run the ExWAS with the call-rate filter applied strictly at 0.99 from the release,
  so there is one clean number to quote
- leave it and simply note in the write-up that the earlier table is superseded

## Where everything is

| | |
|---|---|
| results | `analysis/exwas_elena_replication/phase_4/results/` |
| the audit, re-runnable | `phase_5/scripts/02_poster_audit.py` |
| figures | `analysis/exwas_elena_replication/figures/` |
| full write-up | `phase_5/results/FINDINGS.md` |

The audit script reads the claims out of the `.pptx` rather than from anything I typed,
measures each variant on the release, and refuses to produce a report if a variant can't
be measured — so it can be re-run against the poster file directly.
