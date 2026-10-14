#!/usr/bin/env python3
"""
ExWAS phase 5 -- calibration, and the comparison between the two phenotypes.

LAMBDA IS LEGITIMATE HERE in a way it is not in the burden arm. The genomic inflation
factor assumes a null statistic that is chi-square on 1 df, which is exactly what a
single-variant score test produces. The burden arm's omnibus is Cauchy-combined, which
deflates lambda to 0.68-0.84 throughout and makes it useless as a readout there.

The exact qchisq(0.5, 1) = 0.4549364 is used, not the 0.456 rounded in the original
run's step10_1_QQplot.R, which inflates lambda by 0.23%.

BY REGION, because a whole-genome lambda hides what matters. The first pass found
lambda = 1.34 inside the MHC in EUR with 4 PCs of its own PCA against the original's
0.97 with 9 release PCs, and could not tell whether the cause was the count or the
source. This rebuild uses 15 release PCs, so the question is answerable.

THE TWO ARMS ARE NESTED, and the comparison must be read that way: every bilateral case
is a broad case and the controls are identical, so agreement between them is the same
data twice and measures nothing about reproducibility. A variant getting WEAKER says
nothing either -- losing 53% of cases does that by itself. Only strengthening is
informative, and the first pass's power control showed that the COUNT of strengthening
variants is within what random size-matched case removal produces.
"""
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2, spearmanr

HERE = Path(__file__).resolve().parents[1]
PHASE4 = HERE.parent / "phase_4/results"
OUT = HERE / "results"
COHORTS = ["combined", "EUR", "AFR"]
CHI2_MEDIAN = chi2.ppf(0.5, 1)
MHC = (28_510_000, 33_400_000)
PRIMARY = 5e-8


def lam(p):
    p = np.asarray([v for v in p if 0 < v <= 1])
    return float(np.median(chi2.isf(p, 1)) / CHI2_MEDIAN)


def load(tag):
    d = pd.read_csv(PHASE4 / f"merged_{tag}.tsv", sep="\t",
                    usecols=["CHR", "POS", "MarkerID", "p.value", "BETA",
                             "AF_Allele2", "MissingRate"])
    return d[(d["p.value"] > 0) & (d["p.value"] <= 1)].copy()


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    report = {"arms": {}, "cross_phenotype": {}}

    print("calibration, by region\n")
    print(f"  {'arm':<22}{'λ all':>9}{'λ outside MHC':>16}{'λ in MHC':>11}"
          f"{'n MHC':>9}{'min p':>11}")
    data = {}
    for pheno in ("broad", "bilateral"):
        for cohort in COHORTS:
            tag = f"{pheno}_{cohort}"
            d = load(tag)
            data[tag] = d
            mhc = (d.CHR.astype(str) == "6") & d.POS.between(*MHC)
            e = {
                "tested": len(d), "lambda_all": lam(d["p.value"]),
                "lambda_outside_mhc": lam(d.loc[~mhc, "p.value"]),
                "lambda_inside_mhc": lam(d.loc[mhc, "p.value"]),
                "n_mhc": int(mhc.sum()),
                "min_p": float(d["p.value"].min()),
                "secondary_bar": 0.05 / len(d),
                "hits_primary": int((d["p.value"] < PRIMARY).sum()),
            }
            e["hits_secondary"] = int((d["p.value"] < e["secondary_bar"]).sum())
            report["arms"][tag] = e
            print(f"  {tag:<22}{e['lambda_all']:>9.4f}{e['lambda_outside_mhc']:>16.4f}"
                  f"{e['lambda_inside_mhc']:>11.4f}{e['n_mhc']:>9,}"
                  f"{e['min_p']:>11.2e}")

    print("\nbroad against bilateral, per variant\n")
    for cohort in COHORTS:
        b, l = data[f"broad_{cohort}"], data[f"bilateral_{cohort}"]
        m = b.merge(l[["MarkerID", "p.value", "BETA"]], on="MarkerID",
                    suffixes=("_broad", "_bil"))
        x = -np.log10(m["p.value_broad"])
        y = -np.log10(m["p.value_bil"])
        gained = (y > x + 1.0) & (m["p.value_bil"] < 1e-4)
        rho = float(spearmanr(x, y).statistic)
        report["cross_phenotype"][cohort] = {
            "shared": len(m), "spearman": rho,
            "stronger_on_half_the_cases": int(gained.sum()),
            "stronger_variants": list(m.loc[gained, "MarkerID"])[:20],
        }
        print(f"  {cohort:<10} shared {len(m):>8,}   ρ = {rho:.3f}   "
              f"stronger on half the cases: {int(gained.sum())}")
    print("\n  ρ is reported only to make explicit that it is NOT a replication\n"
          "  statistic: the arms share every control and the bilateral cases are a\n"
          "  strict subset, so agreement is the same data twice.")

    (OUT / "01_calibration.json").write_text(json.dumps(report, indent=2) + "\n")
    print(f"\nwritten: results/01_calibration.json")


if __name__ == "__main__":
    main()
