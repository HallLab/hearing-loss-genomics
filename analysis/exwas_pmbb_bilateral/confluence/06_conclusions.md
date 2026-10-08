# Conclusions

**Bilateral sensorineural hearing loss · PMBB Release 2026-4.0**
[← Phase 5 — results](05_results.md) · [Introduction & pipeline](00_introduction.md)

---

## The verdict

**No gene reaches exome-wide significance.** Not in any cohort, under Bonferroni or FDR, at any
level of aggregation. The closest is `FBXO25` in EUR at 1.37 × 10⁻⁵ against a bar of 2.79 × 10⁻⁶ —
an order of magnitude short once the nine-cell search is charged for.

**The absence is real, not a deflated test.** λ runs 0.68 to 1.00, the QQ points stay inside the
null envelope, and nothing replicates across ancestries beyond chance.

**The limit is the cohort, not the method.** 3,164 cases is small for rare-variant burden, and 6,752
was not large either. This is what the three arms below agree on.

---

## Three arms, one question

The analysis exists alongside two others, and each changes **one thing** from the one before it:

| arm | phenotype | masks | PCs |
|---|---|---|---|
| **Elena** | `SO_396` broad | original | 5 / 9 / 10 |
| **Replication** | `SO_396` broad | rebuilt | 5 / 4 / 3 |
| **Bilateral** *(this work)* | bilateral sensorineural | rebuilt | 5 / 4 / 3 |

So **Elena → Replication** isolates the mask and PC corrections, and **Replication → Bilateral**
isolates the phenotype decision.

To compare them, Elena's arm was given a per-gene omnibus computed the same way — her pipeline emits
none, since its `Cauchy` rows combine annotations inside a cell and leave about nine per gene.

### What the comparison shows

![three arms](../figures/three_way_clingen.png)

**No arm has a significant gene.** Three phenotype-and-method combinations, none crosses the bar.

**The rankings are not noise.** ClinGen Definitive/Strong deafness genes cluster near the top of two
arms more than chance allows — `Replication / combined` at 5 against 0.26 expected
(p = 5.9 × 10⁻⁶), `Elena / EUR` at 4 (p = 1.1 × 10⁻⁴). **This is the one positive result in the
work:** these orderings carry real biology, which a null headline might otherwise suggest they
do not.

### What the comparison does **not** show

The figure invites two stories — *our corrections found the biology*, and *the restriction destroyed
it*. **Neither is supported.** Fisher's exact between every pair of arms:

| cohort | comparison | counts | p |
|---|---|---|---|
| combined | Elena vs Replication | 2 vs 5 | 0.436 |
| combined | Replication vs Bilateral | 5 vs 0 | 0.056 |
| EUR | Elena vs Replication | 4 vs 2 | 0.678 |
| EUR | Replication vs Bilateral | 2 vs 0 | 0.495 |

**Not one pair is distinguishable**, and that test is already optimistic — the arms share people,
genes and controls, so they are not independent.

The counts run 0 to 5 against an expectation of 0.26. Differences among numbers that small are not
measurable here, and reading the figure instead of the table is how a claim gets made that the data
cannot carry.

---

## What the phenotype decision cost

Tested properly rather than argued about. Restricting to bilateral sensorineural halves the cases,
and the ClinGen enrichment disappears — which looks like evidence the restriction discards real
biology. Two controls say otherwise:

| | cases | ClinGen in top 50 |
|---|---:|---:|
| broad phenotype | 6,752 | 5 |
| broad, drawn down to 3,164, five times | 3,164 | 2 · 0 · 1 · 2 · 3 |
| bilateral sensorineural | 3,164 | **0** |
| the cases it discards | 3,588 | 1 |

The five size-matched draws are the yardstick: with the broad definition and 3,164 cases, between 0
and 3 is what can be found, mean 1.6. The restricted arm's 0 sits **inside** that range —
P(0 | Poisson 1.6) = 0.20, and one of the five draws also returned zero. The complement, run as its
own case set, returns 1.

**Neither half carries the enrichment, and they do not add up: 0 + 1 against 5 for the whole.**

> **The restriction's null is lost power, not lost signal.** The phenotype decision is not shown to
> discard biology. It is also not shown to be *better* — only that it is not worse for the reason
> suspected.

---

## What is worth carrying forward

**`FBXO25`** — the closest thing to a signal. In EUR one cell reaches 1.52 × 10⁻⁶, below the
per-gene bar, but the omnibus is 1.37 × 10⁻⁵ with the maximum search penalty, and it rests on 18
alleles, 7 of them in cases. Not a finding; the first gene to check if the cohort grows.

**`SIX1`** — the most consistent gene across everything run. Top 50 in the full broad arm, in two of
five random draws, and the only one in the complement's top 50. Still two orders of magnitude from
significance.

**`GJB3`** — connexin 31, DFNA2B. Shows the behaviour a real gene would: signal concentrated in the
rarest variants, diluting as commoner ones enter.

All three are flags, not results. With ~18,000 genes tested, a known deafness gene landing high by
chance is unremarkable; what makes these worth noting is that identity and behaviour point the same
way.

---

## Next steps

**1 · Audiograms — the one that changes the question.** Diagnosis codes are a proxy.
`H91.90 Unspecified hearing loss, unspecified ear` alone covers 2,698 people who have hearing loss
the record does not characterise, and no phenotyping rule can recover what was never written down.
Audiograms are quantitative and settle type and laterality without inference. The blocker is a
REDCap-to-PMBB ID bridge that does not exist yet.

> "we certainly think the audiograms are going to be a safer bet" — **D. Epstein**, 2026-10-02

When that bridge exists, **only Phase 1 changes**; Phases 2 through 5 run unaltered.

**2 · A larger cohort.** All of Us or UK Biobank. This is what the null asks for, and it is the
third sense of "replication" — reproducing a finding in independent data — that neither this work
nor the replication addresses.

**3 · The X chromosome.** Currently uncovered, costing 6 of the 100 ClinGen Definitive/Strong genes:
`AIFM1`, `COL4A5`, `POU3F4`, `PRPS1`, `SMPX`, `TIMM8A`. Coverage is 94%. Closing it needs X masks
and SAIGE handling hemizygosity in males.

**4 · Open clinical question.** Sudden idiopathic hearing loss (`SO_396.5`, 203 people in the broad
definition) is currently excluded. Whether it belongs in a genetic study of age-related loss is a
question for the clinical side, not a statistical one.

**Deliberately not next:** sensitivity analyses on REVEL 0.5 vs 0.6 or SpliceAI 0.2 vs 0.5. They
would move the ranking but not the conclusion, and running sensitivity analyses on a null result
produces more null results. They become worth doing if the cohort grows.

---

## In three sentences

This cohort, at this size, does not support a rare-variant gene burden finding for bilateral
sensorineural hearing loss, and the test is calibrated well enough that the absence is informative
rather than inconclusive. The rankings underneath that null do carry known deafness biology, which
means the pipeline is working and the signal is simply below the detection floor. What stands
between this result and a finding is sample size and phenotype quality — not analysis choices.

---

*Full detail: [`phase_5/results/FINDINGS.md`](../phase_5/results/FINDINGS.md) ·
three-way notebook: [`phase_5/notebooks/01_three_way.ipynb`](../phase_5/notebooks/01_three_way.ipynb) ·
premises and their sources: [`PREMISES.md`](../PREMISES.md) ·
open items, triaged: [`andre_notes/99_pontos_em_aberto.md`](../andre_notes/99_pontos_em_aberto.md)*
