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

Three unique variants. Each was searched for in five places:

| searched | found |
|---|---|
| `plink_deduplicated/chr{15,16}_deduplicated.bim` | **no** |
| `plink_deduplicated/chr{15,16}_allele_ids.bim` (pre-deduplication) | **no** |
| her raw SAIGE output, `HL_only_ExWAS/saige_results/` | **no** |
| her annotated results, `HL_ExWAS_FullyAnnotated_AllVariants_*.csv` | **no** |
| our merged results | **no** |

Deduplication removed 2 variants on chr15 and 3 on chr16, so it is not the explanation.
The positions sit in dense regions — chr15 has tested variants at 40611381, 40611396,
40611397, 40611403 and 40611412, roughly 15 bp apart — so the region is covered while
the exact position is not.

Where the genes *are* present, they carry no signal:

| | her delivered results | poster |
|---|---|---|
| KNL1 | 7 variants, 40,621,979–40,657,328, best p **0.33** | 3.51e-12 at 40,611,486 |
| ZFPM1 | 7 variants, 88,453,619–88,534,952, best p **0.26** | 2.46e-10 at 88,533,289 |

The ZFPM1 poster position falls inside the range of her tested variants; the KNL1
position sits about 10 kb below hers.

### What this is, and what it is not

**It is not a failed replication.** We reproduce her delivered results at r ≈ 0.98–0.99.
What cannot be reproduced is a set of variants that is absent from the genotypes, from
her raw SAIGE output, and from her own annotated tables.

The most likely explanation is chronology rather than error: the poster was prepared
during the summer, and the delivered files are dated 2 October, uploaded by Nikki on the
day of the meeting — after Doug questioned the common-variant results and Nikki offered
to "adjust … and resend those values". A revised run on a reworked variant set would
produce exactly this pattern.

**The open question for the group is which variant set produced the poster's six hits**,
because it is not the one in `rarevariant_geneburden/plink_deduplicated/`. Until that is
answered, the poster's numbers and the delivered numbers cannot be placed on the same
axis, and neither this arm nor any future one can be said to have replicated or failed
to replicate them.

### Consequence for premise P7

P7 named KNL1 and ZFPM1 as "declared targets" whose p-values this arm would report
whatever they are. That is still done — they are reported above, at p 0.26 and 0.33 in
her data and with no variant of ours clearing any bar. But the premise assumed a target
that existed in the data to be hit or missed. It does not. P7 has been rewritten to say
so rather than left standing as though the comparison were clean.

---

## Files

| | |
|---|---|
| `results/merged_{combined,EUR,AFR}.tsv` | one row per tested variant |
| `results/MERGE_REPORT.md` | the per-cohort table and the checks that passed |
| `scripts/01_saige_step2.bsub` | 66 tasks, step 1 inherited |
| `scripts/02_merge_results.py` | refuses on in-flight jobs, header drift, field counts, row counts above the plink2 ceiling or more than 0.5% below it, and p-values outside (0,1] |
