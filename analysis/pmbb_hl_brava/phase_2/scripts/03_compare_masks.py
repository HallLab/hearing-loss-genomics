#!/usr/bin/env python3
"""
Phase 2 (BRaVa) -- how far apart are the two annotations, measured rather than asserted.

The point of this arm is that a result of ours and a BRaVa-annotated result from another
biobank cannot be placed side by side unless the distance between the annotations is
known. This measures it.

Compared against `pmbb_hl`, NOT `pmbb_hl_v1`. The v1 folder is the retired first pass --
its own PCA, 5/4/3 PCs, no tinnitus, no deduplication -- so a difference measured
against it would mix annotation with ancestry correction. See PREMISES.md P2.

WHAT IS REPORTED

  per mask   genes and entries in each, and the Jaccard of the VARIANT sets
  per gene   the distribution of per-gene Jaccard, which is what a burden test
             actually sees -- a gene can keep its variant count and swap its contents
  ClinGen    the same numbers restricted to the 93 Definitive/Strong genes, because
             those are the genes PREMISES.md P5 pre-declared the test on
  lost genes genes one annotation tests and the other does not

A GENE-LEVEL JACCARD IS THE HONEST FIGURE. An exome-wide Jaccard pools a few huge genes
with 17,000 small ones and reads better than the analysis experiences. The burden test
is per gene, so the per-gene distribution is the one that says what changed.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[1]
OURS = (HERE.parent.parent / "pmbb_hl/rare_variant/phase_2/results/masks")
BRAVA = HERE / "results/masks"
CLINGEN = (HERE.parent.parent / "pmbb_hl/rare_variant/phase_5/data/clingen_hl_genes.tsv")
OUT = HERE / "results"
SHARED = ["pLOF", "pDM", "pLOF_pDM"]


def read_mask(path):
    """{gene: set(variant ids)} from the alternating var/anno format."""
    genes = {}
    with path.open() as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 2 or f[1] != "var":
                continue
            genes[f[0]] = set(f[2:])
    return genes


def jaccard(a, b):
    u = a | b
    return len(a & b) / len(u) if u else np.nan


def main():
    for d in (OURS, BRAVA):
        if not d.exists():
            sys.exit(f"REFUSING: {d} not found")
    OUT.mkdir(parents=True, exist_ok=True)

    cg = pd.read_csv(CLINGEN, sep="\t")
    clingen = set(cg[cg.best_classification.isin(["Definitive", "Strong"])].gene)

    report = {"compared_against": "analysis/pmbb_hl (the current rebuild), "
                                  "not pmbb_hl_v1", "masks": {}}
    print(f"  {'mask':<10}{'genes ours':>12}{'genes BRaVa':>13}{'shared':>9}"
          f"{'Jaccard, all variants':>24}")
    per_gene = {}
    for m in SHARED:
        ours, brava = read_mask(OURS / f"{m}.txt"), read_mask(BRAVA / f"{m}.txt")
        av = set().union(*ours.values()) if ours else set()
        bv = set().union(*brava.values()) if brava else set()
        shared_genes = set(ours) & set(brava)
        j = jaccard(av, bv)
        per_gene[m] = pd.Series(
            {g: jaccard(ours[g], brava[g]) for g in shared_genes}, dtype=float)
        report["masks"][m] = {
            "genes_ours": len(ours), "genes_brava": len(brava),
            "genes_shared": len(shared_genes),
            "genes_only_ours": sorted(set(ours) - set(brava))[:50],
            "n_genes_only_ours": len(set(ours) - set(brava)),
            "n_genes_only_brava": len(set(brava) - set(ours)),
            "variants_ours": len(av), "variants_brava": len(bv),
            "jaccard_all_variants": j,
            "per_gene_jaccard": {
                "median": float(per_gene[m].median()),
                "q25": float(per_gene[m].quantile(0.25)),
                "q75": float(per_gene[m].quantile(0.75)),
                "share_identical": float((per_gene[m] == 1).mean()),
                "share_disjoint": float((per_gene[m] == 0).mean()),
            },
        }
        print(f"  {m:<10}{len(ours):>12,}{len(brava):>13,}{len(shared_genes):>9,}"
              f"{j:>24.3f}")

    print(f"\n  per-gene Jaccard (what a burden test actually sees)\n")
    print(f"  {'mask':<10}{'median':>9}{'q25':>8}{'q75':>8}"
          f"{'identical':>12}{'disjoint':>11}")
    for m in SHARED:
        e = report["masks"][m]["per_gene_jaccard"]
        print(f"  {m:<10}{e['median']:>9.3f}{e['q25']:>8.3f}{e['q75']:>8.3f}"
              f"{e['share_identical']:>11.1%}{e['share_disjoint']:>11.1%}")

    print(f"\n  the 93 ClinGen Definitive/Strong genes\n")
    print(f"  {'mask':<10}{'tested in both':>16}{'median Jaccard':>17}"
          f"{'lost by BRaVa':>15}")
    for m in SHARED:
        ours, brava = read_mask(OURS / f"{m}.txt"), read_mask(BRAVA / f"{m}.txt")
        both = clingen & set(ours) & set(brava)
        lost = sorted(clingen & set(ours) - set(brava))
        med = float(np.median([jaccard(ours[g], brava[g]) for g in both])) if both \
            else np.nan
        report["masks"][m]["clingen"] = {
            "in_both": len(both), "median_jaccard": med,
            "lost_by_brava": lost, "gained_by_brava": sorted(
                clingen & set(brava) - set(ours)),
        }
        print(f"  {m:<10}{len(both):>16}{med:>17.3f}{len(lost):>15}")
        if lost:
            print(f"             lost: {', '.join(lost)}")

    syn = read_mask(BRAVA / "synonymous.txt")
    report["synonymous_control"] = {
        "genes": len(syn), "variants": len(set().union(*syn.values())),
        "clingen_genes_covered": len(clingen & set(syn)),
    }
    print(f"\n  synonymous control set: {len(syn):,} genes, "
          f"{report['synonymous_control']['variants']:,} variants")

    (OUT / "03_mask_comparison.json").write_text(json.dumps(report, indent=2) + "\n")
    print(f"\nwritten: results/03_mask_comparison.json")


if __name__ == "__main__":
    main()
