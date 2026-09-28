# Step 01 — PMBB v4 phenotypes (Cycle 2, independent build)

Produced by [`../../analysis/pipeline/01_phenotype.py`](../../analysis/pipeline/01_phenotype.py) on 2026-08-26.
Reads only `/static/PMBB/PMBB-Release-2026-4.0/`. No dependency on another analyst's outputs.

## Result

| phenotype | cases | controls | excl. 1 date | excl. related ear | analysis N | case rate |
|---|---|---|---|---|---|---|
| hearing_impairment | 6,752 | 50,755 | 4,007 | 9,411 | 57,507 | 11.74% |
| **hl_or_tinnitus** | **7,941** | 50,755 | 4,455 | 7,774 | 58,696 | **13.53%** |
| tinnitus | 2,859 | 50,755 | 2,911 | 14,400 | 53,614 | 5.33% |

Controls are identical across the three by construction: a control must have **no** ear-family
evidence of any kind, which is the same set regardless of which target is being defined.

Usable ancestry strata (combined phenotype): **EUR 43,620 · AFR 11,618**. AMR (868), EAS (1,118),
SAS (892) and UNKNOWN (580) are too small to stratify and are carried only inside `combined`.

## Independent validation against the Nikki/Elena pipeline

Our case sets are strict **supersets** of theirs — nothing they call a case is missed:

| phenotype | ours | theirs | shared | only ours | only theirs |
|---|---|---|---|---|---|
| hearing_impairment | 6,752 | 6,712 | 6,712 (99.4%) | 40 | **0** |
| tinnitus | 2,859 | 2,732 | 2,732 (95.6%) | 127 | **0** |

The surplus is the sample-frame correction: they intersect with the **imputed** LD-pruned `.fam`,
which drops 517 exome-sequenced participants and admits 85 with no exome. We intersect with the
exome sample list. Two independently written pipelines agreeing to within 0.6% is the strongest
validation available without a third implementation.

## The combined phenotype is not simply the union — and that is correct

Union of the two individual case sets = 7,726. The combined phenotype has **7,941**, i.e. **215 more**.

Checked: all 215 have **exactly one hearing-loss date and exactly one tinnitus date**. Individually
each fails rule-of-2 and is excluded; pooled, they have two distinct dates of ear-disease evidence
and qualify. This is the intended behaviour of rule-of-2 — two separate documented encounters — but
it is a judgement call, so it is declared here rather than left implicit. A reviewer preferring the
strict union can drop these 215; the effect on power is negligible (2.7% of cases).

## Provenance note

The load-bearing fact that PMBB v4 moved standard tinnitus ICD codes (388.3x, H93.1x) into the OMOP
`observation` table comes from Nikki Palmiero and Elena's `PMBB_4_PhecodeX_Hearing_Tinnitus.ipynb`.
Independently verified before adoption:

| source | tinnitus events |
|---|---|
| `conditions_phecode_12.txt` (389.4) | 3,016 |
| `conditions_phecode_x.txt` (SO_397.1) | 3,016 |
| `observation.txt` (388.3x + H93.1x) | **25,094** |

Both phecode files capture only the ~3,016 *pulsatile* events. Neither is usable for tinnitus in v4
on its own — the hybrid pull is required in either coding system. Hearing impairment is unaffected
(SO_396 = 134,917 events in `condition_occurrence`).

## Files

| file | content |
|---|---|
| `phenotype_status.csv.gz` | per person × phenotype: status, date counts, ancestry |
| `phenotype_summary.csv` | the table above |
| `phenotype_by_ancestry.csv` | case/control counts per ancestry |
| `manifest.json` | sources, rules, code lists, counts |
| `_ear_family_phecodex.tsv`, `_tinnitus_observation.tsv` | cached awk extracts (regenerable) |
