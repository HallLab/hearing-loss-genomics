# Cycle 2 · Step 1 — Gene set selection and disease→phecode mapping

**Author:** Andre Rico · **Date:** 2026-08-26 · **Status:** for review
**Repo:** `cycle_2/analysis/clingen_fetcher/` · **Data:** `cycle_2/data/clingen/`

**In one line.** We adopted the ClinGen Hearing Loss expert-panel gene list (100 genes) as the
pre-specified set for Q1, then asked which phecodes those genes' diseases correspond to. Only
**61 of 100** genes route to the hearing-loss phecode. That gap is the finding.

---

## Decision requested from this meeting

1. **Ratify or replace the gene set** — ClinGen HL GCEP, Definitive + Strong, 100 genes.
   Adopted provisionally so Q1 is not blocked; swapping tiers costs one flag and a re-run.
2. **Syndromic in or out?** The 100 genes mix syndromic and non-syndromic hearing loss.
   The charter motivates the set as *non-syndromic congenital* HL genes. These are not the same list.
   → **Doug**
3. **Does the phecode gap below change the Tier A / Tier B balance?** → **Molly, Doug, Nikki**
4. **Keep tinnitus in the combined phenotype?** ClinGen curates **no** tinnitus gene–disease pair —
   in any expert panel. The 2026-07-01 decision to test *hearing loss and/or tinnitus* pairs a
   pre-specified HL gene set with a phenotype half that set has no prior for. → **Molly, Doug**

---

## 1. Gene set

**Source:** ClinGen Hearing Loss Gene Curation Expert Panel (affiliate 40007), downloaded
2026-08-26 from `search.clinicalgenome.org/kb/gene-validity/download`.

Why this source and not ClinVar, DVD, or the Cycle-1 173-gene list: it is the only one that is
itself adjudicated. ClinVar is submitter-level and carries conflicts; DVD is a variant database,
not a gene-validity authority; the Cycle-1 list has no documented provenance. Only the expert
panel publishes an evidence tier per gene–disease pair, with a named panel and a date.

| ClinGen classification | Gene–disease pairs | Genes (best) | In primary set? |
|---|---:|---:|---|
| Definitive | 102 | 94 | **yes** |
| Strong | 7 | 6 | **yes** |
| Moderate | 23 | 20 | sensitivity only |
| Limited | 43 | 35 | no |
| Disputed | 12 | 9 | no |
| Refuted | 4 | 2 | no |

**Primary set = Definitive + Strong = 100 genes.** Sensitivity = + Moderate = 120 genes.
Moderate is excluded from the primary because it inflates the FDR denominator ~20% to admit
gene–disease pairs the panel itself flags as possibly not holding up.

**The list moves.** 164 pairs (2019) → 174 (Aug 2024) → 191 (Aug 2026). Any result must
name its snapshot date. Citation: Tshering KC et al., *Genet Med* 2025;27(5):101397 (PMID 39987489) —
a paper whose finding *is* that these classifications drift.

## 2. Disease list

The 100 genes carry **121 curated gene–disease pairs** (75 AR, 39 AD, 7 XL) — full list in
`cycle_2/data/clingen/hl_diseases.tsv`, with MONDO, OMIM and Orphanet identifiers per row.

## 3. Mapping those diseases to phecodes

There is no disease→phecode table, so this is an inference chain, not a lookup:

```
MONDO term --xref--> ICD-9-CM / ICD-10-CM --Phecode 1.2 map--> phecode
```

Phecode map and labels are the ones already in use in this project — `phecode_map12.csv` from
Chapter 2, and phecode descriptions from the PheWAS R package `pheinfo`, i.e. the same table
`createPhenotypes` uses (verified identical to the copy installed on LPC).

**Three caveats, all load-bearing:**

- Most Mendelian syndromes have no billing code of their own — MONDO carries an ICD-10-CM xref for
  only ~2,100 of its ~59,000 terms. Where a term has none, we walk up the `is_a` hierarchy to the
  nearest ancestor that does, and record the depth.
- Not every phecode reached is worth quoting, so each is triaged by rule into **clinical**
  (14 non-hearing phecodes), **wastebasket** (6: the mapping is correct but the phecode is a
  residual "other/unspecified" bin nobody is coded into for this reason), or **suspect**
  (8: the ontology walk went too far, or the underlying xref is wrong). Only *clinical* rows
  appear in the tables below; the rest are in the appendix. The rules live in
  `map_clingen_diseases_to_phecodes.py` and are versioned with the data.
- This is a *disease*→phecode map, not a *patient*→phecode map. Patients are coded for what they
  present with. It tells us where to look for syndromic manifestations; it does not predict who
  will be a case.

## 4. Result — the finding

| Outcome | Genes |
|---|---:|
| Reach phecode **389 (Hearing loss)** | **61** |
| Reach a phecode, but **never** a hearing one | **30** |
| Reach no phecode at all | **9** |
| | **100** |

### 389 is the only hearing phecode reachable — for two different reasons

Never 389.1 (sensorineural), 389.2 (conductive), 389.4 (tinnitus). The two causes are not the same
and have different consequences:

**Tinnitus — nothing to map.** The string "tinnitus" appears **zero times in the entire ClinGen
gene-disease validity table**, across every expert panel, not just hearing loss. ClinGen does not
curate tinnitus as a Mendelian disease entity. This is not an ontology limitation; there is simply
no gene–disease pair to map.

**Sensorineural / conductive — the ontology stops one level short.** Exactly one MONDO term in this
set carries a hearing ICD cross-reference: `MONDO:0005365 "hearing loss disorder"` → `H90`,
`389`, `389.8`, `389.9`. All four land on phecode **389**. The codes that would yield 389.1
(`H90.3`–`H90.5`, `389.1x`) and 389.2 (`H90.0`–`H90.2`, `389.0x`) are cross-referenced by no term in
this set — even for diseases whose *name* says sensorineural. `MONDO:0009968` ("renal tubular
acidosis, distal, 2, with progressive **sensorineural** hearing loss") xrefs `ICD9:389.8`, which
maps to 389, not 389.1.

**Consequence for the phenotype decision.** The 2026-07-01 meeting settled on *hearing loss **and/or
tinnitus*** as the combined Q1 phenotype. This gene set has **no a-priori support for the tinnitus
half** — not weak support, none. Whether that argues for dropping tinnitus from the primary
phenotype, or for keeping it and treating any tinnitus signal as unexpected, is a call for the room.
It does not affect the hearing-loss half.

### Genes whose only phecode is a real, non-hearing manifestation

Every gene below causes deafness. None routes to a hearing-loss phecode, because ClinGen curates
the **syndrome** and the ontology routes the syndrome to its dominant manifestation. These are the
15 with clinically meaningful targets:

| Gene | ClinGen disease | Phecode |
|---|---|---|
| `AIFM1` | X-linked hereditary sensory and autonomic neuropath… | 356 Hereditary and idiopathic peripheral neuropathy |
| `BSND` | Bartter disease type 4A | 255.12 Hyperaldosteronism |
| `BTD` | biotinidase deficiency | 277.6 Other deficiencies of circulating enzymes |
| `CISD2` | Wolfram syndrome | 250.2 Type 2 diabetes |
| `CLPP` | Perrault syndrome 3 | 627.5 Premature menopause and other ovarian failure |
| `DNMT1` | autosomal dominant cerebellar ataxia, deafness and … | 327 Sleep disorders |
| `DSPP` | dentinogenesis imperfecta | 520.1 Hereditary disturbances in tooth structure |
| `HSD17B4` | Perrault syndrome | 627.5 Premature menopause and other ovarian failure |
| `KCNQ1` | Jervell and Lange-Nielsen syndrome | 426.8 Other cardiac conduction disorders |
| `LARS2` | Perrault syndrome | 627.5 Premature menopause and other ovarian failure |
| `MYH9` | macrothrombocytopenia and granulocyte inclusions wi… | 287.31 Primary thrombocytopenia; 580.14 Chronic glomeruloneph… |
| `PRPS1` | PRPS1 deficiency disorder; phosphoribosylpyrophosph… | 356 Hereditary and idiopathic peripheral neuropathy |
| `SLC26A4` | Pendred syndrome | 244.5 Congenital hypothyroidism |
| `TCOF1` | Treacher-Collins syndrome | 749.2 Congenital anomalies of skull and face bones |
| `TWNK` | Perrault syndrome 5 | 627.5 Premature menopause and other ovarian failure |

`KCNQ1` is worth pausing on: a carrier with undiagnosed long QT is a finding with immediate
clinical consequence, independent of hearing.

### Non-hearing phecodes implicated, clinical tier only

This view is phecode-centric, so it also lists genes that reach 389 *as well* — `WFS1` and
`ATP6V1B1` each have a second, non-syndromic curated pair — which is why the gene counts here do
not match the table above.

| Phecode | Description | Genes | Evidence |
|---|---|---|---|
| **627.5** | Premature menopause and other ovarian failure | CLPP, HSD17B4, LARS2, TWNK | ontology rollup |
| **250.2** | Type 2 diabetes | CISD2, WFS1 | direct xref |
| **356** | Hereditary and idiopathic peripheral neuropathy | AIFM1, PRPS1 | ontology rollup |
| **691** | Congenital anomalies of skin | GJB3, GJB6 | ontology rollup |
| **244.5** | Congenital hypothyroidism | SLC26A4 | ontology rollup |
| **255.12** | Hyperaldosteronism | BSND | ontology rollup |
| **277.6** | Other deficiencies of circulating enzymes | BTD | direct xref |
| **287.31** | Primary thrombocytopenia | MYH9 | direct xref |
| **327** | Sleep disorders | DNMT1 | ontology rollup |
| **426.8** | Other cardiac conduction disorders | KCNQ1 | ontology rollup |
| **520.1** | Hereditary disturbances in tooth structure | DSPP | direct xref |
| **580.14** | Chronic glomerulonephritis, NOS | MYH9 | direct xref |
| **588** | Disorders resulting from impaired renal function | ATP6V1B1 | direct xref |
| **749.2** | Congenital anomalies of skull and face bones | TCOF1 | ontology rollup |

### The 16 disease entries that reach no phecode

Two different counts. **16 gene–disease pairs** fail to map, but some of those genes are still
covered through a *different* curated disease. Only **9 genes** reach no phecode at all:
`BCS1L`, `CEP250`, `CEP78`, `CLRN1`, `GATA3`, `GPSM2`, `MITF`, `SLITRK6`, `USH2A`.

| Gene | ClinGen disease | Gene still reaches 389 via another curated disease? |
|---|---|---|
| `BCS1L` | Bjornstad syndrome | **no — no phecode at all** |
| `GPSM2` | Chudley-McCullough syndrome | **no — no phecode at all** |
| `TBC1D24` | DOORS syndrome | yes |
| `CIB2` | Usher syndrome type 1 | yes |
| `MYO7A` | Usher syndrome type 1 | yes |
| `USH1C` | Usher syndrome type 1 | yes |
| `ADGRV1` | Usher syndrome type 2 | yes |
| `USH2A` | Usher syndrome type 2 | **no — no phecode at all** |
| `WHRN` | Usher syndrome type 2D | yes |
| `CLRN1` | Usher syndrome type 3 | **no — no phecode at all** |
| `MITF` | Waardenburg syndrome type 2 | **no — no phecode at all** |
| `CEP78` | cone-rod dystrophy and hearing loss | **no — no phecode at all** |
| `CEP250` | cone-rod dystrophy and hearing loss 2 | **no — no phecode at all** |
| `CDC14A` | hearing impairment and infertile male syndrome | yes |
| `SLITRK6` | high myopia-sensorineural deafness syndrome | **no — no phecode at all** |
| `GATA3` | hypoparathyroidism-deafness-renal disease syndrome | **no — no phecode at all** |

7 of these entries are **Usher syndrome** — the most important syndromic hearing-loss entity,
and it has no ICD code of its own.

## 5. What this implies for the design

The ICD/phecode layer **systematically under-represents the hearing phenotype of syndromic
hearing-loss genes**. This was already the argument for the Tier B audiometric phenotype in the
Cycle 2 charter (§3); it is now quantified, with gene names, rather than asserted.

1. **A phecode-based secondary phenotype cannot be built from this mapping alone.** For
   39 of 100 genes it would point at the wrong organ system or at nothing.
2. **Tier A's `hearing impairment` phecode stays the discovery phenotype** — but its known
   coarseness now has a second, independent line of evidence behind it.
3. **The syndromic manifestations are a real, testable secondary hypothesis.** If carriers of
   `SLC26A4`, `MYH9`, `KCNQ1` or `TCOF1` variants are enriched for the phecodes above, that is a
   positive control for the carrier set that does not assume the hearing answer.

## 6. Open questions

- **Syndromic vs non-syndromic** (→ Doug). Restricting to non-syndromic means filtering on
  `DISEASE LABEL` / MONDO, not on the gene — a different operation, and a biology call.
- **16 genes carry gene–disease pairs at more than one evidence tier** (e.g. `CIB2` = Refuted +
  Definitive; `GJB6` = Definitive + Refuted). The current rule keeps each gene at its *best*
  classification. Defensible for a disease-agnostic burden test; revisit if the set is narrowed.
- **Do the syndromic phecodes above have usable N in PMBB v4?** Not yet counted.

## 7. Appendix — phecodes excluded from the tables above

Kept for completeness and for anyone auditing the mapping. Excluded because the phecode carries no
clinical information, or because the evidence does not survive review — not because the gene is
uninteresting.

| Gene | ClinGen disease | Phecode |
|---|---|---|
| `ABHD12` | PHARC syndrome | 272.9 Unspecified disorder of lipoid metabolism; 277.5 Other … |
| `ALMS1` | Alstrom syndrome | 759 Other and unspecified congenital anomalies |
| `CHD7` | CHARGE syndrome | 759 Other and unspecified congenital anomalies |
| `COL4A5` | Alport syndrome | 759 Other and unspecified congenital anomalies |
| `COL9A1` | Stickler syndrome, type 4 | 759 Other and unspecified congenital anomalies |
| `COL9A3` | Stickler syndrome | 759 Other and unspecified congenital anomalies |
| `DIAPH1` | DIAPH1-related sensorineural hearing loss-thrombocy… | 287.3 Thrombocytopenia |
| `EYA1` | branchio-oto-renal syndrome | 759 Other and unspecified congenital anomalies |
| `FGF3` | deafness with labyrinthine aplasia, microtia, and m… | 759 Other and unspecified congenital anomalies |
| `PAX3` | Waardenburg syndrome | 759 Other and unspecified congenital anomalies |
| `SIX1` | branchio-oto-renal syndrome | 759 Other and unspecified congenital anomalies |
| `SLC52A2` | Brown-Vialetto-van Laere syndrome 2 | 334.2 Anterior horn cell disease |
| `SLC52A3` | Brown-Vialetto-van Laere syndrome 1 | 334.2 Anterior horn cell disease |
| `SOX10` | Waardenburg syndrome type 4C | 564 Functional digestive disorders; 565 Anal and rectal condi… |
| `TIMM8A` | deafness dystonia syndrome | 759 Other and unspecified congenital anomalies |

**One known-bad mapping, hard-coded as such:** MONDO cross-references "primary ovarian failure" to
`ICD9:253.4`, which is *disorders of the anterior pituitary*. Primary ovarian failure is gonadal,
not pituitary. This affects the four Perrault-syndrome genes, whose correct phecode — 627.5 — the
same MONDO term also yields. The mapper flags `253.4` as suspect for this reason.

**A known false positive in the other direction:** the rule marks anything resolved three ontology
hops up as *suspect*, which demotes `SLC52A2` / `SLC52A3` → 334.2 (anterior horn cell disease).
Brown-Vialetto-van Laere syndrome genuinely is a motor neuron disorder. *Suspect* means "check
before quoting", not "wrong".

## 8. Reproducibility

| | |
|---|---|
| Gene set | `cycle_2/analysis/clingen_fetcher/fetch_clingen_hl_genes.py` |
| Phecode map | `cycle_2/analysis/clingen_fetcher/map_clingen_diseases_to_phecodes.py` |
| This page | `cycle_2/analysis/clingen_fetcher/render_confluence_page.py` |
| Outputs | `cycle_2/data/clingen/` (git-ignored; regenerate with the scripts above) |
| ClinGen snapshot | 2026-08-26, sha256 `753fb2c70af9bd69…` — pinned in `manifest.json` |

Every number and every table row on this page is read from those outputs at render time.
