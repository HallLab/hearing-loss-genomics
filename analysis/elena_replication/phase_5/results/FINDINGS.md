# Phase 5 — reading the result

**Status:** complete · **Arm:** corrected only
**Notebooks:** [`01_results.ipynb`](../notebooks/01_results.ipynb) (canonical) · [`01_results.pt.ipynb`](../notebooks/01_results.pt.ipynb)
**Scripts:** `01_gene_positions.py` · `02_tables.py` · `burden_plots.py`

Phase 5 covers the pipeline's Steps 9 to 11. The chromosome merge (her Step 9.1) was pulled forward
into Phase 4, because it was needed to compare against her results, and is not repeated.

Phase 4 already reported that nothing is significant in either arm and that the ranking moves. What
Phase 5 does is **show** both, quantify how much of the result is load-bearing, and close the
multiple-testing question — Finding 0 — which Phase 4 had left resting on Bonferroni alone.

---

## The figures

| file | what it answers |
|---|---|
| `manhattan_<cohort>.png` | does anything cross the bar, anywhere in the grid? |
| `qq_<cohort>.png` | is the test calibrated, or is the null manufactured? |
| `comparison_ours_vs_hers.png` | which genes did the corrections move, and which disappeared? |

The comparison figure is the one this replication exists to produce, and she has no equivalent. One
point per gene, her p on x and ours on y, with the genes our Phase 2 rebuild removed drawn as a
strip along the bottom — they have no y value because we never tested them. **The rightmost mark in
the combined panel is `TMC3-AS1`**, her rank-1 gene.

---

## Finding 0 — nothing passes, by Bonferroni or by FDR, at any level

The question this phase has to answer without ambiguity, since it is the first
thing anyone will press on.

### Which Bonferroni bar the red line is

The line drawn on each Manhattan panel is `0.05 / genes in that panel` — that panel treated as if it
were the only analysis. It is **nine separate corrections per cohort**, not one correction over the
whole grid, and so it is the most permissive bar in play: it does not charge for having looked at
nine panels.

| bar | denominator | combined | what it assumes |
|---|---|---|---|
| the red line | genes in that panel, ~17,800 | 2.81 × 10⁻⁶ | this panel is the whole analysis |
| `0.05/genes` | genes in the cohort, 17,943 | 2.79 × 10⁻⁶ | the nine cells per gene count as one test |
| `0.05/tests` | all 159,960 tests | 3.13 × 10⁻⁷ | every cell is its own test |

The first two nearly coincide, but only by arithmetic accident — a panel holds almost every gene in
the cohort. They are different claims, and the honest bar for "did anything in this cohort pass" is
one of the latter two, because scanning nine panels and reporting the best is nine times the search.

The permissive line is what gets drawn **because nothing crosses even it.** Where that is not true,
the distinction would have to be made explicitly.

### And FDR finds nothing either

Bonferroni is the strict correction; Benjamini-Hochberg is the lenient one, and it is where a weak
but real signal would show up first. It was run at both levels:

| cohort | min q, within a panel | min q, over the cohort | genes at q < 0.05 |
|---|---|---|---|
| combined | 0.075 | 0.374 | **0** |
| EUR | 0.271 | 0.785 | **0** |
| AFR | 0.117 | 0.260 | **0** |

Zero in every cohort, at either level. And the q-values themselves are the informative part, more
than the count: the best gene in each cohort carries between a 26% and a 79% chance of being a false
positive once the cohort is taken as a whole. That is not "narrowly missed".

This closes the question in a way Bonferroni alone does not. A study with a weak real effect usually
shows *something* under FDR. Nothing here does, and not by a small margin.

No FDR line is drawn on the Manhattan plots, because BH's cutoff is the largest p with q < 0.05 and
there is none — a line would have nowhere to sit. Each panel is annotated with its smallest q
instead, which carries the same information without implying a threshold that does not exist.

The q-values are in `summary.tsv` and per gene in `top_hits_<cohort>.tsv`.

### No Cauchy rows, and what that costs — raised by Nikki Palmiero

SAIGE emits a `Cauchy` row per gene when several annotation groups are requested in one run,
combining them into one omnibus p-value. **Our results contain none**, by construction: we pass one
annotation per mask, so there is nothing for it to combine. Hers has 489,886 of them, 27% of her
table.

Two things had to be checked rather than assumed.

**Did her Cauchy rows contaminate the comparison in Phase 4?** No. My filter kept them — it selects
on the mask folder, not on `Group` — so 30% of the rows I compared against were omnibus p-values
mixed in with per-annotation ones. Recomputed with them excluded, all three top-50 lists are
**identical, 50 of 50**, the best p is unchanged in every cohort and so is the rank-1 gene. Only
246–294 genes out of ~19,000 have their best p coming from a Cauchy row, and none sit near the top.
The reason is structural: an omnibus across annotations lands between its components, so it rarely
beats the best one.

**Is losing the omnibus a loss?** It is the principled answer to the multiple-testing question
above — one test per gene, with no denominator left to argue about. So it was computed post hoc over
the 9 cells per gene, using the same ACAT statistic SAIGE uses:

| cohort | best single cell | Cauchy omnibus | bar 0.05/genes | significant |
|---|---|---|---|---|
| combined | 4.21 × 10⁻⁶ | 3.02 × 10⁻⁵ | 2.79 × 10⁻⁶ | no |
| EUR | 1.54 × 10⁻⁵ | 6.91 × 10⁻⁵ | 2.79 × 10⁻⁶ | no |
| AFR | 7.23 × 10⁻⁶ | 2.17 × 10⁻⁵ | 2.81 × 10⁻⁶ | no |

Nothing passes, and the omnibus is **further** from passing than the best cell — by a factor of 3 to
7. That is the expected direction and it is informative: our top genes carry signal in one cell
only, so averaging in the eight where they show nothing dilutes it. A gene with a consistent effect
across masks and cutoffs would move the other way.

It also exposes something about the headline numbers. Taking the smallest of 9 correlated p-values
per gene and setting it against `0.05/genes` is mildly anti-conservative — it does not charge for
the nine looks. The omnibus is the version that does charge, and under it the null is firmer still.

Her `Cauchy` rows carry a separate problem for anyone reading her tables: `Newmasks_HL_TopGenes_*.csv`
has no `Group` column, so an omnibus p-value and a single-annotation p-value sit in the same ranking
with nothing to tell them apart.

---

## Finding 1 — the null is calibrated across all 27 panels

λ, the genomic inflation factor, over every cohort × mask × MAF cell:

```
                  maxMAF    1e-4    1e-3    1e-2
combined  pLOF              1.040   1.020   1.016
          pDM               0.962   0.947   0.953
          pLOF_pDM          0.972   0.927   0.938
EUR       pLOF              1.048   1.023   1.012
          pDM               1.002   0.969   0.959
          pLOF_pDM          1.006   0.950   0.950
AFR       pLOF              0.947   0.982   0.962
          pDM               1.033   1.006   1.004
          pLOF_pDM          1.108   1.056   1.060
```

Range 0.927 to 1.108, median 1.002, **none outside 0.9–1.15**. The QQ points stay inside the 95%
null envelope and sit slightly below the diagonal through the middle, which is the same mild
conservatism Phase 4 measured by counting tail fractions. Two methods, same answer.

This is what licenses the claim "there is nothing here" rather than only "nothing passed the bar".

### λ on a SKAT-O p-value, and her constant

Two notes a reviewer would raise.

λ is built for single-variant tests, where the null statistic is chi-square on 1 df. SKAT-O's is
not. It still works as used, because what is being tested is whether the p-values are **uniform**:
if they are, the median of `qchisq(1-p, 1)` is `qchisq(0.5, 1)` whatever produced them. It is a
uniformity check, not a claim about the statistic's distribution.

Her `step10_1_QQplot.R` divides by `0.456`. The exact value is `qchisq(0.5, 1)` = 0.4549364, so her
λ is inflated by 0.23%. It changes no conclusion anywhere, and is recorded only because the
correction is free.

---

## Finding 2 — most of the top of each list is fragile

A p-value alone ranks a five-allele gene alongside a four-hundred-allele one. `top_hits_<cohort>.tsv`
flags a test as fragile when any of three things holds:

- the burden component disagrees with the combined p by more than 100-fold — the signal is all SKAT,
  with no consistent direction
- `Number_rare = 0`, so SAIGE collapsed the whole gene into one unit and Burden and SKAT became the
  same test
- MAC below 20

| cohort | flagged, of the top 30 |
|---|---:|
| combined | 16 |
| EUR | 18 |
| **AFR** | **26** |

AFR is the striking one, and the explanation is its size: 11,334 people and 1,285 cases. At that
scale the smallest p-values are reached by genes with a handful of alleles, where removing one
person would move the result. `CFH` tops that cohort on **five alleles**, all in cases.

The practical consequence: the flagged fraction should be read before the gene names. A top-30 table
from AFR is not a list of 30 candidates.

---

## Finding 3 — gene coordinates, and a check worth keeping

The Manhattan plots need a coordinate per gene, which SAIGE's output does not carry. Two sources
were used, in order: the release's own Ensembl metadata
(`Exome/metadata/Homo_sapiens.GRCh38.113.ENSG_locations_symbols.tsv`), and our own masks, whose
variant IDs give a build-matched position for any symbol the metadata lacks.

All 17,943 tested genes were placed: 17,799 from the metadata, 144 from our masks.

The second source doubles as a check on the first, and it is worth keeping for that reason alone:

- **chromosome disagreements: 0** — no symbol collision put a gene on the wrong chromosome
- **variants outside the metadata interval: 429** — inspected, and benign. All are on the correct
  chromosome and overshoot by 3 kb to 73 kb, the ordinary difference between one Ensembl gene record
  and the transcripts VEP used to assign variants to that symbol. On a 250 Mb chromosome that is
  invisible.

Had either number come back large, the plots would have been wrong in a way nothing downstream would
have caught.
