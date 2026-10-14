#!/usr/bin/env python3
"""
Phase 5 -- do known deafness genes cluster at the top of the rankings?

THE ONE POSITIVE RESULT of the first pass was that they do, in the broad combined
cohort: 5 ClinGen Definitive/Strong genes in the top 50 against 0.26 expected,
p = 5.9e-06. This rebuild changed the principal components, so the question is whether
that survives better ancestry correction or was an artefact of the old ones.

THE TEST is hypergeometric. Of the N genes tested, K are on the ClinGen list. If the
ranking carried no information, the number of ClinGen genes in the top n would follow
a hypergeometric distribution. The p-value asks how surprising the observed count is.

WHAT IT DOES AND DOES NOT SAY. It is a statement about the RANKING, not about any gene.
A significant enrichment with no significant gene means the pipeline is finding real
biology that sits below the detection floor. It does not license naming any of the
genes as a finding.

REPORTED AT THREE DEPTHS because a result that appears only at one is a result about
that cutoff. The first pass's headline was the top 50.
"""
import json
from pathlib import Path

import pandas as pd
from scipy.stats import hypergeom

HERE = Path(__file__).resolve().parents[1]
PHASE4 = HERE.parent / "phase_4/results"
OUT = HERE / "results"
ARMS = [(p, c) for p in ("broad", "bilateral")
        for c in ("combined", "EUR", "AFR")]
DEPTHS = [50, 100, 250]
TIERS = {"definitive_strong": ["Definitive", "Strong"],
         "plus_moderate": ["Definitive", "Strong", "Moderate"]}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    cg = pd.read_csv(HERE / "data/clingen_hl_genes.tsv", sep="\t")
    manifest = json.loads((HERE / "data/clingen_manifest.json").read_text())

    report = {"clingen_sha256": manifest.get("sha256"),
              "clingen_genes": manifest.get("genes"), "arms": {}}

    for tier, classes in TIERS.items():
        gene_set = set(cg[cg.best_classification.isin(classes)].gene)
        print(f"=== {tier}: {len(gene_set)} genes on the list ===\n")
        for pheno, cohort in ARMS:
            tag = f"{pheno}_{cohort}"
            d = pd.read_csv(PHASE4 / f"omnibus_{tag}.tsv", sep="\t")
            d = d.sort_values("omnibus_p").reset_index(drop=True)
            N = len(d)
            tested = d.Region.isin(gene_set)
            K = int(tested.sum())

            entry = report["arms"].setdefault(tag, {})[tier] = {
                "genes_tested": N, "clingen_tested": K, "depths": {}}

            line = f"  {tag:<20} {K:>3} of {len(gene_set)} on the list are tested  "
            for n in DEPTHS:
                obs = int(tested.head(n).sum())
                exp = K * n / N
                p = hypergeom.sf(obs - 1, N, K, n) if obs else 1.0
                entry["depths"][str(n)] = {
                    "observed": obs, "expected": round(exp, 3), "p": p,
                    "genes": list(d.head(n).Region[tested.head(n)])}
                line += f"| top {n}: {obs} vs {exp:.2f}  p={p:.2e}  "
            print(line)
        print()

    dest = OUT / "02_enrichment.json"
    dest.write_text(json.dumps(report, indent=2) + "\n")

    print("=== which ClinGen genes, top 50, definitive/strong ===")
    for pheno, cohort in ARMS:
        tag = f"{pheno}_{cohort}"
        g = report["arms"][tag]["definitive_strong"]["depths"]["50"]["genes"]
        print(f"  {tag:<20} {', '.join(g) if g else '—'}")
    print(f"\nwritten: {dest.relative_to(HERE)}")


if __name__ == "__main__":
    main()
