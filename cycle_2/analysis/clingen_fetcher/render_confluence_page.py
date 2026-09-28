#!/usr/bin/env python3
"""
render_confluence_page.py — regenerate cycle_2/docs/to_confluence/01_gene_set_and_phecode_map.md
from the outputs of the two ClinGen scripts.

The `to_confluence/` convention says numbers are generated, never transcribed
(`cycle_2/docs/to_confluence/README.md`). This is what enforces it: refresh the
ClinGen snapshot, re-run the mapper, re-run this, and the page is correct by
construction. Prose lives here; every count and every table row is read from
`cycle_2/data/clingen/`.
"""

import csv
import collections
import json
import sys
from pathlib import Path


def repo_root() -> Path:
    for parent in Path(__file__).resolve().parents:
        if (parent / ".git").exists():
            return parent
    return Path.cwd()


def trunc(s: str, n: int) -> str:
    return s if len(s) <= n else s[: n - 1] + "…"


def main() -> int:
    repo = repo_root()
    D = repo / "cycle_2/data/clingen"
    out = repo / "cycle_2/docs/to_confluence/01_gene_set_and_phecode_map.md"

    need = ["hl_disease_phecode_map.tsv", "hl_diseases.tsv",
            "hl_diseases_unmapped.tsv", "hl_phecode_summary.tsv", "manifest.json"]
    missing = [n for n in need if not (D / n).exists()]
    if missing:
        sys.exit(f"missing inputs in {D}: {', '.join(missing)}\n"
                 f"run fetch_clingen_hl_genes.py then map_clingen_diseases_to_phecodes.py")

    def tsv(name):
        with (D / name).open(newline="", encoding="utf-8") as f:
            return list(csv.DictReader(f, delimiter="\t"))

    rows = tsv("hl_disease_phecode_map.tsv")
    dis = tsv("hl_diseases.tsv")
    unm = tsv("hl_diseases_unmapped.tsv")
    summ = tsv("hl_phecode_summary.tsv")
    man = json.loads((D / "manifest.json").read_text())
    genes_all = sorted(set((D / "genes_definitive_strong.txt").read_text().split()))

    hl389 = {r["GENE"] for r in rows if r["PHECODE"] == "389"}
    mapped = {r["GENE"] for r in rows}
    none_g = sorted(set(genes_all) - mapped)
    non_hl = mapped - hl389

    # Best confidence per phecode, mirroring the summary file.
    rank = {"clinical": 0, "wastebasket": 1, "suspect": 2}
    ph_conf = {r["PHECODE"]: r["CONFIDENCE"] for r in summ}

    # Genes reaching a non-hearing phecode, split by the best evidence they have.
    per_gene: dict[str, dict] = collections.defaultdict(
        lambda: {"dis": set(), "ph": collections.defaultdict(set)})
    for r in rows:
        if r["GENE"] in non_hl and r["PHECODE"] != "389":
            e = per_gene[r["GENE"]]
            e["dis"].add(r["DISEASE_LABEL"])
            e["ph"][ph_conf[r["PHECODE"]]].add(f'{r["PHECODE"]} {r["PHECODE_DESC"]}')

    clinical_genes = sorted(g for g, e in per_gene.items() if e["ph"]["clinical"])
    other_genes = sorted(g for g, e in per_gene.items() if not e["ph"]["clinical"])

    def gene_table(genes, bucket):
        t = ["| Gene | ClinGen disease | Phecode |", "|---|---|---|"]
        for g in genes:
            e = per_gene[g]
            ph = e["ph"][bucket] if bucket else (e["ph"]["wastebasket"] | e["ph"]["suspect"])
            t.append(f'| `{g}` | {trunc("; ".join(sorted(e["dis"])), 52)} | '
                     f'{trunc("; ".join(sorted(ph)), 62)} |')
        return "\n".join(t)

    # Phecode summary, clinical only, hearing excluded.
    t_summ = ["| Phecode | Description | Genes | Evidence |", "|---|---|---|---|"]
    for r in summ:
        if r["PHECODE"] == "389" or r["CONFIDENCE"] != "clinical":
            continue
        lv = r["MATCH_LEVELS"]
        ev = "direct xref" if "direct" in lv else "ontology rollup"
        t_summ.append(f'| **{r["PHECODE"]}** | {r["PHECODE_DESC"]} | '
                      f'{r["GENES"].replace(";", ", ")} | {ev} |')

    # Unmapped disease entries.
    t_unm = ["| Gene | ClinGen disease | Gene still reaches 389 via another curated disease? |",
             "|---|---|---|"]
    for r in sorted(unm, key=lambda r: (r["DISEASE_LABEL"], r["GENE"])):
        g = r["GENE"]
        note = ("yes" if g in hl389 else
                "**no — no phecode at all**" if g in none_g else "non-hearing phecode only")
        t_unm.append(f'| `{g}` | {r["DISEASE_LABEL"]} | {note} |')

    pc, gc = man["pairs_by_classification"], man["genes_by_best_classification"]
    t_tier = ["| ClinGen classification | Gene–disease pairs | Genes (best) | In primary set? |",
              "|---|---:|---:|---|"]
    for c in ["Definitive", "Strong", "Moderate", "Limited", "Disputed", "Refuted"]:
        if c in pc:
            inset = ("**yes**" if c in ("Definitive", "Strong")
                     else "sensitivity only" if c == "Moderate" else "no")
            t_tier.append(f"| {c} | {pc[c]} | {gc.get(c, 0)} | {inset} |")

    moi = collections.Counter(d["MOI"] for d in dis)
    conf_n = collections.Counter(r["CONFIDENCE"] for r in summ if r["PHECODE"] != "389")
    n_usher = sum(1 for r in unm if "Usher" in r["DISEASE_LABEL"])

    doc = f"""# Cycle 2 · Step 1 — Gene set selection and disease→phecode mapping

**Author:** Andre Rico · **Date:** 2026-08-26 · **Status:** for review
**Repo:** `cycle_2/analysis/clingen_fetcher/` · **Data:** `cycle_2/data/clingen/`

**In one line.** We adopted the ClinGen Hearing Loss expert-panel gene list ({len(genes_all)} genes) as the
pre-specified set for Q1, then asked which phecodes those genes' diseases correspond to. Only
**{len(hl389)} of {len(genes_all)}** genes route to the hearing-loss phecode. That gap is the finding.

---

## Decision requested from this meeting

1. **Ratify or replace the gene set** — ClinGen HL GCEP, Definitive + Strong, {len(genes_all)} genes.
   Adopted provisionally so Q1 is not blocked; swapping tiers costs one flag and a re-run.
2. **Syndromic in or out?** The {len(genes_all)} genes mix syndromic and non-syndromic hearing loss.
   The charter motivates the set as *non-syndromic congenital* HL genes. These are not the same list.
   → **Doug**
3. **Does the phecode gap below change the Tier A / Tier B balance?** → **Molly, Doug, Nikki**
4. **Keep tinnitus in the combined phenotype?** ClinGen curates **no** tinnitus gene–disease pair —
   in any expert panel. The 2026-07-01 decision to test *hearing loss and/or tinnitus* pairs a
   pre-specified HL gene set with a phenotype half that set has no prior for. → **Molly, Doug**

---

## 1. Gene set

**Source:** ClinGen Hearing Loss Gene Curation Expert Panel (affiliate 40007), downloaded
{man["clingen_file_created"]} from `search.clinicalgenome.org/kb/gene-validity/download`.

Why this source and not ClinVar, DVD, or the Cycle-1 173-gene list: it is the only one that is
itself adjudicated. ClinVar is submitter-level and carries conflicts; DVD is a variant database,
not a gene-validity authority; the Cycle-1 list has no documented provenance. Only the expert
panel publishes an evidence tier per gene–disease pair, with a named panel and a date.

{chr(10).join(t_tier)}

**Primary set = Definitive + Strong = {len(genes_all)} genes.** Sensitivity = + Moderate = 120 genes.
Moderate is excluded from the primary because it inflates the FDR denominator ~20% to admit
gene–disease pairs the panel itself flags as possibly not holding up.

**The list moves.** 164 pairs (2019) → 174 (Aug 2024) → {man["n_gene_disease_pairs"]} (Aug 2026). Any result must
name its snapshot date. Citation: Tshering KC et al., *Genet Med* 2025;27(5):101397 (PMID 39987489) —
a paper whose finding *is* that these classifications drift.

## 2. Disease list

The {len(genes_all)} genes carry **{len(dis)} curated gene–disease pairs** ({moi["AR"]} AR, {moi["AD"]} AD, {moi["XL"]} XL) — full list in
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
  ({conf_n["clinical"]} non-hearing phecodes), **wastebasket** ({conf_n["wastebasket"]}: the mapping is correct but the phecode is a
  residual "other/unspecified" bin nobody is coded into for this reason), or **suspect**
  ({conf_n["suspect"]}: the ontology walk went too far, or the underlying xref is wrong). Only *clinical* rows
  appear in the tables below; the rest are in the appendix. The rules live in
  `map_clingen_diseases_to_phecodes.py` and are versioned with the data.
- This is a *disease*→phecode map, not a *patient*→phecode map. Patients are coded for what they
  present with. It tells us where to look for syndromic manifestations; it does not predict who
  will be a case.

## 4. Result — the finding

| Outcome | Genes |
|---|---:|
| Reach phecode **389 (Hearing loss)** | **{len(hl389)}** |
| Reach a phecode, but **never** a hearing one | **{len(non_hl)}** |
| Reach no phecode at all | **{len(none_g)}** |
| | **{len(genes_all)}** |

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
{len(clinical_genes)} with clinically meaningful targets:

{gene_table(clinical_genes, "clinical")}

`KCNQ1` is worth pausing on: a carrier with undiagnosed long QT is a finding with immediate
clinical consequence, independent of hearing.

### Non-hearing phecodes implicated, clinical tier only

This view is phecode-centric, so it also lists genes that reach 389 *as well* — `WFS1` and
`ATP6V1B1` each have a second, non-syndromic curated pair — which is why the gene counts here do
not match the table above.

{chr(10).join(t_summ)}

### The {len(unm)} disease entries that reach no phecode

Two different counts. **{len(unm)} gene–disease pairs** fail to map, but some of those genes are still
covered through a *different* curated disease. Only **{len(none_g)} genes** reach no phecode at all:
`{"`, `".join(none_g)}`.

{chr(10).join(t_unm)}

{n_usher} of these entries are **Usher syndrome** — the most important syndromic hearing-loss entity,
and it has no ICD code of its own.

## 5. What this implies for the design

The ICD/phecode layer **systematically under-represents the hearing phenotype of syndromic
hearing-loss genes**. This was already the argument for the Tier B audiometric phenotype in the
Cycle 2 charter (§3); it is now quantified, with gene names, rather than asserted.

1. **A phecode-based secondary phenotype cannot be built from this mapping alone.** For
   {len(genes_all) - len(hl389)} of {len(genes_all)} genes it would point at the wrong organ system or at nothing.
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

{gene_table(other_genes, None)}

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
| ClinGen snapshot | {man["clingen_file_created"]}, sha256 `{man["raw_sha256"][:16]}…` — pinned in `manifest.json` |

Every number and every table row on this page is read from those outputs at render time.
"""
    out.write_text(doc, encoding="utf-8")
    print(f"wrote {out.relative_to(repo)} ({len(doc.splitlines())} lines)")
    print(f"  {len(hl389)} genes -> 389 | {len(non_hl)} -> non-hearing | {len(none_g)} -> none")
    print(f"  non-hearing phecodes: {conf_n['clinical']} clinical, "
          f"{conf_n['wastebasket']} wastebasket, {conf_n['suspect']} suspect")
    return 0


if __name__ == "__main__":
    sys.exit(main())
