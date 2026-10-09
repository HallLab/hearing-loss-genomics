#!/usr/bin/env python3
"""
Phase 5, step 1 -- calibration, and comparison against the broad-phenotype ExWAS.

THE TWO ARMS ARE NOT INDEPENDENT, and the comparison has to be read that way.
Verified below rather than assumed:

    all 3,164 bilateral cases are also cases in the broad arm
    the 50,746 controls are the same individuals in both

So the bilateral analysis is the broad analysis with 3,588 cases removed. A p-value
that agrees between the two is not a replication -- it is the same data twice.
Correlation here measures nothing about reproducibility.

WHAT THE COMPARISON CAN SHOW is the opposite direction: a variant that gets STRONGER
on half the cases. That cannot come from shared data, because removing cases removes
information. It means the association concentrates in the bilateral sensorineural
subtype. Those variants, if any, are the only output of this comparison worth reading.

A variant that gets weaker tells us nothing on its own -- losing 53% of cases does
that by itself. The burden arm's power control
(gene_burden_pmbb_bilateral/phase_5/power_control/) already established that for this
phenotype, and the same caution applies here.
"""
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2, spearmanr

BASE = Path(__file__).resolve().parents[2]
BROAD = BASE.parent / "exwas_elena_replication"
BURDEN = BASE.parent / "gene_burden_pmbb_bilateral"
BROAD_BURDEN = BASE.parent / "gene_burden_elena_replication"
ELENA = BASE.parent / "elena/HL_only_ExWAS"
OUT = BASE / "phase_5/results"
COHORTS = ["combined", "EUR", "AFR"]
CHI2_MEDIAN = chi2.ppf(0.5, 1)


def lam(p):
    p = np.asarray([v for v in p if 0 < v <= 1])
    return float(np.median(chi2.isf(p, 1)) / CHI2_MEDIAN)


def results(root, cohort):
    d = pd.read_csv(root / f"phase_4/results/merged_{cohort}.tsv", sep="\t")
    return d[(d["p.value"] > 0) & (d["p.value"] <= 1)].copy()


def check_nesting():
    def read(p):
        cases, ctrls = set(), set()
        for line in Path(p).read_text().splitlines()[1:]:
            f = line.split("\t")
            (cases if f[1] == "1" else ctrls).add(f[0])
        return cases, ctrls

    bc, bk = read(BROAD_BURDEN / "phase_3/results/covariates/"
                  "corrected_combined_covariates.txt")
    lc, lk = read(BURDEN / "phase_3/results/covariates/combined.txt")
    out = {
        "broad_cases": len(bc), "broad_controls": len(bk),
        "bilateral_cases": len(lc), "bilateral_controls": len(lk),
        "bilateral_cases_also_broad": len(lc & bc),
        "bilateral_cases_not_broad": len(lc - bc),
        "controls_identical": lk == bk,
        "cases_dropped": len(bc - lc),
    }
    if out["bilateral_cases_not_broad"] or not out["controls_identical"]:
        raise SystemExit("REFUSING: the two arms are not nested as assumed. "
                         "The whole interpretation below depends on that. " + str(out))
    return out


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    nest = check_nesting()
    print("Cohort nesting, verified:")
    print(f"  broad       {nest['broad_cases']:,} cases, "
          f"{nest['broad_controls']:,} controls")
    print(f"  bilateral   {nest['bilateral_cases']:,} cases, "
          f"{nest['bilateral_controls']:,} controls")
    print(f"  every bilateral case is a broad case, controls identical")
    print(f"  cases removed by the restriction: {nest['cases_dropped']:,}\n")

    report = {"nesting": nest, "cohorts": {}}

    for c in COHORTS:
        ours = results(BASE, c)
        broad = results(BROAD, c)
        l_ours, l_broad = lam(ours["p.value"]), lam(broad["p.value"])

        m = ours.merge(broad[["MarkerID", "p.value", "BETA"]], on="MarkerID",
                       suffixes=("_bil", "_broad"))
        x = -np.log10(m["p.value_broad"])
        y = -np.log10(m["p.value_bil"])
        rho = float(spearmanr(x, y).statistic)

        # The only informative direction: stronger on half the cases.
        stronger = m[(y > x + 1.0) & (m["p.value_bil"] < 1e-4)].copy()
        stronger["gain"] = (-np.log10(stronger["p.value_bil"])
                            + np.log10(stronger["p.value_broad"]))
        stronger = stronger.nlargest(20, "gain")

        report["cohorts"][c] = {
            "n_tested": len(ours), "n_broad": len(broad), "n_shared": len(m),
            "lambda": l_ours, "lambda_broad": l_broad,
            "min_p": float(ours["p.value"].min()),
            "min_p_broad": float(broad["p.value"].min()),
            "spearman_vs_broad": rho,
            "n_stronger_by_1_log_and_p_lt_1e4": int(len(
                m[(y > x + 1.0) & (m["p.value_bil"] < 1e-4)])),
            "stronger": [
            {"variant": row["MarkerID"],
             "p_bilateral": float(row["p.value_bil"]),
             "p_broad": float(row["p.value_broad"]),
             "or_bilateral": float(np.exp(row["BETA_bil"])),
             "or_broad": float(np.exp(row["BETA_broad"])),
             "af": float(row["AF_Allele2"]),
             "log10_gain": float(row["gain"])}
            for _, row in stronger.iterrows()],
        }

        print(f"{c}")
        print(f"  tested          {len(ours):>8,}   (broad arm {len(broad):,})")
        print(f"  lambda          {l_ours:>8.4f}   (broad arm {l_broad:.4f})")
        print(f"  min p           {ours['p.value'].min():>8.2e}   "
              f"(broad arm {broad['p.value'].min():.2e})")
        print(f"  Spearman vs broad  {rho:>5.3f}   "
              f"[nested data -- not a replication statistic]")
        n_str = report["cohorts"][c]["n_stronger_by_1_log_and_p_lt_1e4"]
        print(f"  stronger on half the cases (>1 log and p<1e-4): {n_str}")
        for s in report["cohorts"][c]["stronger"][:5]:
            print(f"      {s['variant']:<22} bilateral {s['p_bilateral']:.2e} "
                  f"(OR {s['or_bilateral']:.2f})  broad {s['p_broad']:.2e} "
                  f"(OR {s['or_broad']:.2f})  AF {s['af']:.3f}")
        print()

    dest = OUT / "calibration_crossarm.json"
    dest.write_text(json.dumps(report, indent=2) + "\n")
    print(f"written: {dest.relative_to(BASE)}")


if __name__ == "__main__":
    main()
