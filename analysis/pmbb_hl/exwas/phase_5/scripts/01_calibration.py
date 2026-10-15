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

THE ARMS ARE NOT INDEPENDENT, and each pair must be read for what it is. Broad and
bilateral are nested -- every bilateral case is a broad case, controls identical -- so
agreement between them is the same data twice. Broad and tinnitus overlap at 66% but
are not nested: 974 tinnitus cases are not HL cases, so a third of that case set is
new. Both pairs share all 50,755 controls, which induces agreement by itself. A variant getting WEAKER says
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
    for pheno in ("broad", "bilateral", "tinnitus"):
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

    # TWO PAIRS, AND THEY ARE NOT THE SAME KIND OF COMPARISON.
    #   broad vs bilateral: nested. Every bilateral case is a broad case, controls
    #     identical. Agreement is the same data twice and measures nothing.
    #   broad vs tinnitus: overlapping but NOT nested. 974 of 2,859 tinnitus cases
    #     are not broad cases (34%), so the case sets are partly independent even
    #     though the controls are shared. Agreement here is weak evidence; it is
    #     still not replication, because the shared controls alone induce it.
    PAIRS = [("broad", "bilateral", "nested: bilateral cases are a strict subset"),
             ("broad", "tinnitus", "overlapping: 974 of 2,859 tinnitus cases "
                                   "are not broad cases")]
    for a, b_ph, relation in PAIRS:
        print(f"\n{a} against {b_ph}, per variant  ({relation})\n")
        for cohort in COHORTS:
            left, right = data[f"{a}_{cohort}"], data[f"{b_ph}_{cohort}"]
            m = left.merge(right[["MarkerID", "p.value", "BETA"]], on="MarkerID",
                           suffixes=("_l", "_r"))
            x = -np.log10(m["p.value_l"])
            y = -np.log10(m["p.value_r"])
            gained = (y > x + 1.0) & (m["p.value_r"] < 1e-4)
            # Sign agreement is only meaningful among variants that reached
            # anything at all in the left arm; over 70,000 near-null variants the
            # betas are noise and the share sits at 0.5 by construction.
            lead = m[m["p.value_l"] < 1e-3]
            concord = (float((np.sign(lead.BETA_l) == np.sign(lead.BETA_r)).mean())
                       if len(lead) else None)
            rho = float(spearmanr(x, y).statistic)
            report["cross_phenotype"][f"{a}_vs_{b_ph}_{cohort}"] = {
                "relation": relation, "shared": len(m), "spearman": rho,
                "stronger_in_second": int(gained.sum()),
                "stronger_variants": list(m.loc[gained, "MarkerID"])[:20],
                "lead_variants_left": len(lead), "sign_concordance": concord,
            }
            cc = "    n/a" if concord is None else f"{concord:>7.2f}"
            print(f"  {cohort:<10} shared {len(m):>8,}   ρ = {rho:.3f}   "
                  f"stronger in {b_ph}: {int(gained.sum()):>2}   "
                  f"sign concordance (p_l < 1e-3, n={len(lead)}):{cc}")
    print("\n  ρ is NOT a replication statistic in either pair: both share every\n"
          "  control, which induces agreement on its own. For the nested pair it is\n"
          "  the same data twice. For the overlapping pair a third of the tinnitus\n"
          "  cases are new, so it is weak evidence and nothing more.")

    (OUT / "01_calibration.json").write_text(json.dumps(report, indent=2) + "\n")
    print(f"\nwritten: results/01_calibration.json")


if __name__ == "__main__":
    main()
