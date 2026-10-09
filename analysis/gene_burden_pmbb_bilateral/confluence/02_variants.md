# Phase 2 — Which variants count

**Bilateral sensorineural hearing loss · PMBB Release 2026-4.0**
[← Phase 1 — cohort](01_cohort.md) · Part 2 of 5 · [Phase 3 — covariates →](03_covariates.md)

---

## In one line

Three masks decide which rare variants enter each gene's test — **462k, 890k and 1.34M entries**
across ~17,900 genes — built from the release VEP annotation and applied to protein-coding genes
only.

---

## What a mask is

A mask is a list: for each gene, which variants are summed into its test.

```
A1BG   var    19:58347031:T:C   19:58347349:T:C   19:58347352:C:A   ...
A1BG   anno   pLOF              pLOF              pLOF              ...
```

A participant carries ~20,000 rare exome variants, most of which do nothing. Summing all of them
buries any real signal, so the mask is where "which variants could plausibly break this gene" is
decided — and it is as consequential as the phenotype.

---

## The three masks

| mask | mechanism | genes | variants | per gene |
|---|---|---:|---:|---:|
| `pLOF` | the protein is **broken** | 17,841 | 462,144 | 26 |
| `pDM` | one amino acid **swapped**, predicted damaging | 17,623 | 890,332 | 51 |
| `pLOF_pDM` | the union | 17,945 | 1,340,936 | 75 |

A fourth mask, `ALL` — every rare variant regardless of consequence — was **not run**. At 884
variants per gene it pools synonymous and intronic variants with the ones that have a mechanism.
The clinical lead reached the same position independently:

> "that's the confusing part to me, why you would ever use the all category […] it's going to give
> you a lot of noise" — **D. Epstein**, 2026-10-02

---

## The rules, as implemented

**`pLOF`** — a variant enters if either holds:

1. the consequence is **exactly** one of `stop_gained`, `frameshift_variant`,
   `splice_acceptor_variant`, `splice_donor_variant`, `start_lost`, `stop_lost`,
   `transcript_ablation`
2. **or** the consequence mentions splice in any other way **and SpliceAI ≥ 0.2**

Condition 2 is the gate that separates splice variants that matter from ones that merely sit nearby.
Several low-impact splice terms — `splice_donor_region`, `splice_polypyrimidine_tract` — describe
variants *near* the splice site which usually do nothing; SpliceAI decides case by case.

Matching is **exact, not substring**: searching for `splice_donor_variant` inside the consequence
text would also catch `splice_donor_region_variant`, which is a different thing.

**`pDM`** — AlphaMissense classifies the variant `pathogenic` or `likely_pathogenic`, **or** the
**maximum REVEL across transcripts** is ≥ 0.5. The maximum matters: VEP returns REVEL as a list,
one value per transcript (`".,0.65,0.712,."`), because a variant can be missense in one transcript
and something else in another.

**`pLOF_pDM`** — the union, with `pLOF` taking precedence where both apply. 11,517 variants are
both, being loss-of-function in one transcript and missense in another.

**All three require `BIOTYPE == protein_coding`.** A loss-of-function test asks whether losing the
protein's function associates with disease; in a lncRNA or pseudogene there is no protein to lose,
and the test returns a number that cannot mean anything. This removes **1,101 genes of 19,038**.

---

## Provenance

These masks were **not rebuilt for this analysis**. They depend on the genome annotation and the
damage rules, neither of which changes when the phenotype changes, so they were carried over from
`analysis/elena_replication` and verified byte-for-byte. Rebuilding would reproduce the same files.

The rules above are the study's original analysis plan, applied as written. Where the earlier
implementation diverged from them is documented in that folder and is not repeated here.

---

## Known limitations

- **No X chromosome.** Masks cover chromosomes 1–22, which costs **6 of the 100** ClinGen
  Definitive/Strong deafness genes: `AIFM1`, `COL4A5`, `POU3F4`, `PRPS1`, `SMPX`, `TIMM8A`. Coverage
  is 94%.
- **Thresholds are choices.** REVEL at 0.5 and SpliceAI at 0.2 are the analysis plan's values and
  are defensible, but moving either would reshuffle the ranking. Neither has been tested.
- **AlphaMissense coverage is unmeasured.** Where it has no score, REVEL alone decides. How many
  genes that affects has not been quantified.

---

*Detail: [`phase_2/results/FINDINGS.md`](../phase_2/results/FINDINGS.md) ·
code: [`phase_2/scripts/07_rebuild_masks.py`](../phase_2/scripts/07_rebuild_masks.py) ·
provenance: [`PROVENANCE.md`](../PROVENANCE.md)*
