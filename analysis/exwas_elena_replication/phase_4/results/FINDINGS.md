# Phase 4 — single-variant ExWAS, results

**Run 2026-10-09** · job `49569375` (66 tasks) + `49571242` (1 retry) · 66/66 cells merged.
No step 1 was fitted: the null models are inherited from the burden arm (premise P3),
which is what Elena's own config did.

## Nothing is significant, in any cohort, at either bar

| cohort | tested | secondary bar | min p | hits at 5e-8 | hits at secondary |
|---|---:|---:|---:|---:|---:|
| combined | 70,366 | 7.106e-07 | 3.02e-06 | 0 | 0 |
| EUR | 74,536 | 6.708e-07 | 2.12e-05 | 0 | 0 |
| AFR | 117,953 | 4.239e-07 | 4.11e-06 | 0 | 0 |

The best result anywhere is 4.3× short of the loosest bar and 60× short of the primary.

## The replication itself is excellent

Joined against her delivered results on `MarkerID`:

| cohort | hers | ours | shared | r of −log10(p) |
|---|---:|---:|---:|---:|
| combined | 70,429 | 70,366 | 70,333 (99.9%) | **0.989** |
| EUR | 74,528 | 74,536 | 74,483 (99.9%) | **0.979** |
| AFR | 117,943 | 117,953 | 117,723 (99.8%) | **0.990** |

This is the closest agreement any arm of this project has produced. The variants that
appear in only one set number 33–230 per cohort and are explained by the cohort
difference: our corrected cohort has 40 more cases and 174 fewer controls, which moves
a handful of variants across the 1% threshold in each direction.

The largest per-variant disagreements all sit in the MHC (chr6:31–32 Mb), where LD is
extreme and our PC counts differ from hers (4/3 against her 9/10).

Her own minimum p-values — 2.57e-06 combined, 1.99e-05 EUR, 3.70e-06 AFR — sit within
a factor of 1.2 of ours, and her own `Significance` column reads "Not significant" for
every one of them.

---

## The poster reports six genome-wide significant results. They are not in the data.

`elena_publishes/Horizontal Draft LANS Poster - SUIP-PGSIP.pptx`, under "Genome Wide
Significant Results Total: 6", lists:

| gene | variant | cohort | p |
|---|---|---|---|
| KNL1 | chr15:40611486:A:G | Combined | 3.51e-12 |
| KNL1 | chr15:40611486:A:G | EUR | 7.21e-12 |
| ZFPM1 | chr16:88533289:A:C | EUR | 2.46e-10 |
| ZFPM1 | chr16:88533290:G:C | EUR | 2.52e-10 |

Three unique variants. **All three exist in the PMBB v4 release, and all three fail the
quality filter the poster itself describes.**

The poster's methods read: *"variants with a genotype call rate < 0.99 and a Minor Allele
Frequency (MAF) < 0.01 were removed."* Measured on the release:

| variant | missing | call rate | MAF | poster p |
|---|---:|---:|---:|---:|
| `chr15_40611486_A_G` | 56,646 / 70,925 | **0.2013** | 0.284 | 3.51e-12 |
| `chr16_88533289_A_C` | 2,717 / 70,925 | 0.9617 | 0.295 | 2.46e-10 |
| `chr16_88533290_G_C` | 3,736 / 70,925 | 0.9473 | 0.283 | 2.52e-10 |

**The KNL1 variant carrying the 3.51e-12 is genotyped in 20% of the cohort.** 79.9% of
participants have no call at that position. Differential missingness between cases and
controls is the standard way such a variant produces an enormous spurious p-value, and
excluding it is precisely what a call-rate filter is for.

### Where each was dropped

| stage | chr16 variants | poster variants present |
|---|---:|---|
| PMBB v4 release | 780,404 | yes |
| her `variant_qc/chr16_callrate.bim` | 670,730 | **yes** |
| her `plink_deduplicated/chr16_allele_ids.bim` | 662,670 | no |
| her `chr16_deduplicated.bim` | 662,667 | no |

The chr15 variant (call rate 0.20) was removed at the `variant_qc` stage. The two chr16
variants survived it and were dropped at the `allele_ids` stage.

A true 0.99 call-rate filter applied to the release leaves **475,376** variants on chr16.
Her `chr16_callrate.bim` holds **670,730**. So the stage named for call rate is not
filtering at 0.99 — it is considerably more permissive, which is how two variants at
0.96 and 0.95 passed through it.

The variants are absent from her raw SAIGE output, her annotated tables and our merged
results, which is consistent: by the time step 2 ran they were gone.

Where the genes *are* present, they carry no signal:

| | her delivered results | poster |
|---|---|---|
| KNL1 | 7 variants, 40,621,979–40,657,328, best p **0.33** | 3.51e-12 at 40,611,486 |
| ZFPM1 | 7 variants, 88,453,619–88,534,952, best p **0.26** | 2.46e-10 at 88,533,289 |

### What this is, and what it is not

**It is not a disagreement between our results and hers.** We reproduce her delivered
results at r ≈ 0.98–0.99, and her own `Significance` column reads "Not significant"
throughout. The poster disagrees with her own pipeline's output, not with ours.

The sequence that fits: the poster was prepared over the summer from a run that still
contained low-call-rate variants; the QC was later tightened and the variant set reworked;
the delivered files are dated 2 October, uploaded on the day of the meeting, and contain
nothing significant. Doug's reaction in that meeting — *"I wasn't excited about any of
those genes"* and *"if the variants that associated were not pathogenic, then I didn't
understand the value of that"* — turns out to have been right for a more concrete reason
than was available to him at the time.

**This is what a replication is for.** The arm Doug deprioritised as "just contributing
noise" is the one that produced the most consequential finding of the project so far, and
it did so by failing to reproduce a result rather than by finding one.

### Consequence for premise P7

P7 named KNL1 and ZFPM1 as "declared targets" whose p-values this arm would report
whatever they are. That is still done — they are reported above, at p 0.26 and 0.33 in
her data and with no variant of ours clearing any bar.

But the premise assumed a target that could be hit or missed. The three variants behind
it fail the poster's own stated QC, so there was never a result there to replicate. P7 is
amended to say that, rather than left standing as though the comparison were clean.

---

## Files

| | |
|---|---|
| `results/merged_{combined,EUR,AFR}.tsv` | one row per tested variant |
| `results/MERGE_REPORT.md` | the per-cohort table and the checks that passed |
| `scripts/01_saige_step2.bsub` | 66 tasks, step 1 inherited |
| `scripts/02_merge_results.py` | refuses on in-flight jobs, header drift, field counts, row counts above the plink2 ceiling or more than 0.5% below it, and p-values outside (0,1] |
