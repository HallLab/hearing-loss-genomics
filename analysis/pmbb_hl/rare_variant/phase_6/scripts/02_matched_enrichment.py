#!/usr/bin/env python3
"""
Phase 6 -- is the ClinGen enrichment biology, or is it where the power happens to be?

WHY THIS EXISTS. 01_power_vs_clingen.py found that ClinGen genes are substantially
better covered than the average gene: median MAC 453 against 190 in the broad combined
cohort, Mann-Whitney p < 0.001 in all six arms. Well-studied genes tend to be larger and
more variable, which is partly why they were found in the first place.

That raises a specific doubt about phase 5's one positive result. If better-covered genes
rank higher for reasons unrelated to biology, the enrichment would be measuring coverage.

TWO CHECKS, AND ONLY THE SECOND IS DECISIVE.

First, calibration: under the null a p-value is uniform regardless of power, so coverage
should not shift rank at all. Measured across MAC sextiles in the broad combined arm,
the fraction below 0.05 is 0.041 to 0.054 everywhere -- flat. The MEDIAN p rises with MAC
(0.522 to 0.643), which is conservativeness at high coverage, the opposite of the bias
being looked for. So nothing in the tail, which is what ranks 1 to 50 depend on, favours
well-covered genes.

Second, and this is the test: hold coverage FIXED. Draw 93 non-ClinGen genes matched to
the ClinGen genes on log MAC, and count how many land in the top 50. Repeat. If the real
set beats matched draws, the enrichment is not about coverage.

SCOPE: the two arms with any ClinGen gene in their top 50. The other four have zero, so
there is nothing for a matched draw to beat.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[1]
PHASE4 = HERE.parent / "phase_4/results"
OUT = HERE / "results"
DRAWS = 2000
SEED = 20261014
WINDOW = 0.10          # log10 MAC, so a factor of about 1.26
TOP_N = 50
# ONLY THE ARMS THAT HAVE SOMETHING TO EXPLAIN AWAY. This test asks whether the
# ClinGen genes' presence in the top 50 survives matching on coverage. It is only
# meaningful where they are present: bilateral and tinnitus put 0 ClinGen genes in
# the top 50 of every cohort, and a matched null for an observed count of 0 answers
# nothing. Running it there would manufacture three more p-values of 1.0.
ARMS = [("broad", "combined"), ("broad", "EUR")]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    power = pd.read_csv(OUT / "gene_power.tsv", sep="\t")
    rng = np.random.default_rng(SEED)
    report = {"draws": DRAWS, "seed": SEED, "log10_mac_window": WINDOW,
              "top_n": TOP_N, "arms": {}}

    for pheno, cohort in ARMS:
        tag = f"{pheno}_{cohort}"
        om = pd.read_csv(PHASE4 / f"omnibus_{tag}.tsv", sep="\t")[
            ["Region", "omnibus_p"]]
        gp = power[(power.phenotype == pheno) & (power.cohort == cohort)][
            ["Region", "MAC", "clingen"]]
        d = om.merge(gp, on="Region").sort_values("omnibus_p").reset_index(drop=True)
        d["rank"] = np.arange(1, len(d) + 1)
        d["lmac"] = np.log10(d.MAC.clip(lower=1))

        cg = d[d.clingen]
        obs = int((cg["rank"] <= TOP_N).sum())
        pool = d[~d.clingen].reset_index(drop=True)
        targets = cg.lmac.to_numpy()

        counts = np.empty(DRAWS, dtype=int)
        for i in range(DRAWS):
            used = np.zeros(len(pool), dtype=bool)
            hits = 0
            for t in targets:
                near = np.where((np.abs(pool.lmac.to_numpy() - t) < WINDOW) & ~used)[0]
                if near.size == 0:
                    near = np.where(~used)[0]
                j = near[rng.integers(near.size)]
                used[j] = True
                hits += pool["rank"].iat[j] <= TOP_N
            counts[i] = hits

        p_emp = float((counts >= obs).mean())
        report["arms"][tag] = {
            "observed": obs, "draw_mean": float(counts.mean()),
            "draw_max": int(counts.max()),
            "draw_distribution": {int(k): int(v) for k, v in
                                  zip(*np.unique(counts, return_counts=True))},
            "empirical_p": p_emp,
            "clingen_median_mac": float(cg.MAC.median()),
            "others_median_mac": float(pool.MAC.median()),
        }
        r = report["arms"][tag]
        print(f"{tag}")
        print(f"  observed                 {obs} ClinGen genes in the top {TOP_N}")
        print(f"  MAC-matched draws        mean {r['draw_mean']:.2f}, "
              f"max {r['draw_max']}, n = {DRAWS:,}")
        print(f"  distribution             {r['draw_distribution']}")
        print(f"  empirical p              "
              f"{'< ' + format(1/DRAWS, '.4f') if p_emp == 0 else format(p_emp, '.4f')}")
        print(f"  (ClinGen median MAC {r['clingen_median_mac']:,.0f} vs "
              f"{r['others_median_mac']:,.0f} — matched away by construction)")
        print()

    (OUT / "02_matched_enrichment.json").write_text(
        json.dumps(report, indent=2) + "\n")
    print("written: results/02_matched_enrichment.json")


if __name__ == "__main__":
    main()
