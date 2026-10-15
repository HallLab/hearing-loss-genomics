#!/usr/bin/env python3
"""
Phase 6 -- are the known deafness genes in the part of the study that has power?

THE QUESTION THIS ANSWERS. Phase 5 found two things that sit oddly together: nothing is
significant anywhere, and known deafness genes cluster at the top of the broad combined
ranking at p = 5.9e-06. Phase 5's power curves explain the first. This asks whether the
second happened in spite of the first.

If the ClinGen genes carry as few variants as everything else, the enrichment is
appearing WITHOUT the study being powered to detect any of them individually. That
strengthens the reading considerably: it means the ranking is picking up real biology
from almost no evidence per gene, rather than being powered for those genes and not
for others.

If instead the ClinGen genes are unusually well covered, the enrichment could be partly
an artefact of where power happens to sit, and the claim needs softening.

Nothing in phase 5 distinguishes these.

HOW. For every gene, the carrier count is known from the test itself, so the detectable
odds ratio is computable per gene rather than assumed. Then compare the distribution
for ClinGen genes against everything else, and report how many of either group the
study could have found at all.

The power model is phase 5's, unchanged: the gene as one super-variant at its observed
cumulative frequency, every masked variant assumed causal. Optimistic by construction.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import brentq
from scipy.stats import chi2, mannwhitneyu, ncx2

HERE = Path(__file__).resolve().parents[1]
RV = HERE.parent
PHASE4 = RV / "phase_4/results"
PHASE5 = RV / "phase_5"
OUT = HERE / "results"

# N AND CASES ARE READ, NOT WRITTEN DOWN. They were hardcoded at the broad arm's
# values and reused for bilateral, which gave that arm 57,498 people when its model
# was fitted on 53,910 -- 3,588 people it does not have. The error was conservative,
# not optimistic: information goes as mac * phi * (1 - phi) with phi = cases/n, so an
# inflated n dilutes the case fraction and makes the arm look LESS powered than it is
# (bilateral combined: detectable OR 3.20 reported against 3.10 actual). Each
# phenotype excludes a different set of "other ear" patients, so every arm has its
# own N. Phase 3 already records what SAIGE was actually given; that is the source.
PHASE3 = Path(__file__).resolve().parents[2] / "phase_3/results/01_covariates.json"
POWER_TARGET = 0.80
PHENOTYPES = ("broad", "bilateral", "tinnitus")


def cohort_sizes():
    arms = json.loads(PHASE3.read_text())["arms"]
    out = {}
    for pheno, cohorts in arms.items():
        for cohort, e in cohorts.items():
            # 'written' is the covariate file SAIGE read; 'analysable' still counts
            # the people dropped for a missing AGE.
            out[(pheno, cohort)] = (e["written"], e["cases"])
    return out


def detectable_or(mac, n, cases, alpha):
    """Smallest OR this gene could have been found at, given its own carrier count."""
    f = mac / (2 * n)
    if f <= 0 or f >= 1:
        return np.inf
    phi = cases / n
    info = 2 * n * f * (1 - f) * phi * (1 - phi)
    crit = chi2.isf(alpha, 1)
    try:
        ncp = brentq(lambda l: ncx2.sf(crit, 1, l) - POWER_TARGET, 1e-6, 1e4)
    except ValueError:
        return np.inf
    return float(np.exp(np.sqrt(ncp / info)))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    cg = pd.read_csv(PHASE5 / "data/clingen_hl_genes.tsv", sep="\t")
    gene_set = set(cg[cg.best_classification.isin(["Definitive", "Strong"])].gene)

    report = {"power_target": POWER_TARGET, "arms": {}}
    rows = []

    sizes = cohort_sizes()
    for pheno in PHENOTYPES:
        for cohort in ("combined", "EUR", "AFR"):
            tag = f"{pheno}_{cohort}"
            cells = pd.read_csv(PHASE4 / f"cells_{tag}.tsv", sep="\t")
            om = pd.read_csv(PHASE4 / f"omnibus_{tag}.tsv", sep="\t")
            alpha = float(om.bonferroni_bar.iloc[0])

            # the widest mask at the loosest cutoff: the most evidence any gene gets
            g = cells[(cells.Group == "pLOF;pDM") & (cells.max_MAF == 0.01)].copy()
            if g.empty:
                sys.exit(f"REFUSING: no pLOF;pDM cells at MAF 0.01 for {tag}")
            g = g.groupby("Region", as_index=False).MAC.max()

            if (pheno, cohort) not in sizes:
                sys.exit(f"REFUSING: phase 3 records no cohort size for {tag}")
            n, ca = sizes[(pheno, cohort)]
            g["detectable_or"] = [detectable_or(m, n, ca, alpha) for m in g.MAC]
            g["clingen"] = g.Region.isin(gene_set)
            g["phenotype"], g["cohort"] = pheno, cohort
            rows.append(g)

            cgg, other = g[g.clingen], g[~g.clingen]
            u = mannwhitneyu(cgg.MAC, other.MAC, alternative="two-sided")
            ent = {
                "genes": len(g), "clingen_genes": int(len(cgg)),
                "median_mac_clingen": float(cgg.MAC.median()),
                "median_mac_other": float(other.MAC.median()),
                "mannwhitney_p": float(u.pvalue),
                "median_detectable_or_clingen": float(cgg.detectable_or.median()),
                "median_detectable_or_other": float(other.detectable_or.median()),
                "clingen_powered_for_or5": int((cgg.detectable_or <= 5).sum()),
                "clingen_powered_for_or3": int((cgg.detectable_or <= 3).sum()),
                "other_powered_for_or3_pct": float((other.detectable_or <= 3).mean()),
                "bonferroni_bar": alpha,
            }
            report["arms"][tag] = ent

            print(f"{tag}")
            print(f"  median MAC   ClinGen {ent['median_mac_clingen']:>8,.0f}   "
                  f"others {ent['median_mac_other']:>8,.0f}   "
                  f"Mann-Whitney p = {ent['mannwhitney_p']:.3f}")
            print(f"  median detectable OR   ClinGen "
                  f"{ent['median_detectable_or_clingen']:>6.2f}   "
                  f"others {ent['median_detectable_or_other']:>6.2f}")
            print(f"  ClinGen genes the study could find at OR 3: "
                  f"{ent['clingen_powered_for_or3']} of {ent['clingen_genes']}"
                  f"   (others: {100*ent['other_powered_for_or3_pct']:.1f}%)")
            print()

    pd.concat(rows, ignore_index=True).to_csv(
        OUT / "gene_power.tsv", sep="\t", index=False)
    (OUT / "01_power_vs_clingen.json").write_text(json.dumps(report, indent=2) + "\n")
    print("written: results/gene_power.tsv, 01_power_vs_clingen.json")


if __name__ == "__main__":
    main()
