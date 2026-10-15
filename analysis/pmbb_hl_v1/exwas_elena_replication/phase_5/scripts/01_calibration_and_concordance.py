#!/usr/bin/env python3
"""
Phase 5, step 1 -- calibration and concordance against Elena's delivered results.

LAMBDA IS LEGITIMATE HERE, unlike in the burden arm. The genomic inflation factor is
built for single-variant tests, where the null statistic is chi-square on 1 df, and
that is exactly what this arm produces. The burden arm needed a paragraph of defence
for applying lambda to a SKAT-O p-value; this one does not.

The exact qchisq(0.5, 1) = 0.4549364 is used, not the 0.456 rounded in Elena's
step10_1_QQplot.R, which inflates lambda by 0.23%.

CONCORDANCE is the headline of this arm. Nothing is significant, so what the phase has
to report is whether the pipeline reproduces hers -- joined per variant, not compared
as summary statistics.
"""
import csv
import json
import math
from pathlib import Path

import numpy as np
from scipy.stats import chi2

BASE = Path(__file__).resolve().parents[2]
OURS = BASE / "phase_4/results/merged_{}.tsv"
HERS = (BASE.parent / "elena/HL_only_ExWAS/"
        "HL_ExWAS_FullyAnnotated_AllVariants_{}.csv")
OUT = BASE / "phase_5/results"
COHORTS = ["combined", "EUR", "AFR"]
CHI2_MEDIAN = chi2.ppf(0.5, 1)  # 0.4549364


def lam(p):
    p = np.asarray([v for v in p if 0 < v <= 1])
    return float(np.median(chi2.isf(p, 1)) / CHI2_MEDIAN)


def read_ours(cohort):
    path = Path(str(OURS).format(cohort))
    out = {}
    with path.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        i_id, i_p = header.index("MarkerID"), header.index("p.value")
        i_af, i_miss = header.index("AF_Allele2"), header.index("MissingRate")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            try:
                out[f[i_id]] = (float(f[i_p]), float(f[i_af]), float(f[i_miss]))
            except ValueError:
                continue
    return out


def read_hers(cohort):
    path = Path(str(HERS).format(cohort))
    out = {}
    with path.open() as fh:
        for row in csv.DictReader(fh):
            try:
                out[row["MarkerID"]] = (float(row["p.value"]),
                                        row.get("SYMBOL", ""),
                                        row.get("Consequence", ""))
            except (ValueError, KeyError):
                continue
    return out


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    report = {}

    for cohort in COHORTS:
        ours, hers = read_ours(cohort), read_hers(cohort)
        shared = sorted(set(ours) & set(hers))

        lo = lam([v[0] for v in ours.values()])
        lh = lam([v[0] for v in hers.values()])

        x = np.array([-math.log10(hers[v][0]) for v in shared])
        y = np.array([-math.log10(ours[v][0]) for v in shared])
        r = float(np.corrcoef(x, y)[0, 1])

        # Rank agreement matters more than r when nothing is significant: r is
        # dominated by the bulk near zero, where both are noise.
        from scipy.stats import spearmanr
        rho = float(spearmanr(x, y).statistic)

        # Where do they disagree most, and is it structural?
        gap = np.abs(x - y)
        order = np.argsort(gap)[::-1]
        worst = [{"variant": shared[i],
                  "p_hers": hers[shared[i]][0], "p_ours": ours[shared[i]][0],
                  "gene": hers[shared[i]][1], "consequence": hers[shared[i]][2],
                  "missing_rate_ours": ours[shared[i]][2]}
                 for i in order[:15]]

        # MHC is chr6:28.5-33.4 Mb (GRCh38). Expect disagreement to pile up there.
        def in_mhc(v):
            try:
                c, pos = v.split(":")[0], int(v.split(":")[1])
            except (IndexError, ValueError):
                return False
            return c == "6" and 28_510_000 <= pos <= 33_400_000

        big = [shared[i] for i in order if gap[i] > 1.0]
        report[cohort] = {
            "n_ours": len(ours), "n_hers": len(hers), "n_shared": len(shared),
            "only_ours": len(set(ours) - set(hers)),
            "only_hers": len(set(hers) - set(ours)),
            "lambda_ours": lo, "lambda_hers": lh,
            "pearson_r_neglog10p": r, "spearman_rho": rho,
            "min_p_ours": min(v[0] for v in ours.values()),
            "min_p_hers": min(v[0] for v in hers.values()),
            "n_disagree_over_1_log": len(big),
            "share_of_those_in_mhc": (sum(in_mhc(v) for v in big) / len(big)
                                      if big else None),
            "worst_15": worst,
        }

        print(f"{cohort}")
        print(f"  tested      ours {len(ours):>8,}   hers {len(hers):>8,}   "
              f"shared {len(shared):>8,}")
        print(f"  lambda      ours {lo:>8.4f}   hers {lh:>8.4f}")
        print(f"  min p       ours {min(v[0] for v in ours.values()):>8.2e}   "
              f"hers {min(v[0] for v in hers.values()):>8.2e}")
        print(f"  r / rho     {r:>8.4f} / {rho:.4f}")
        if big:
            print(f"  |diff|>1 log: {len(big):,}  of which in MHC: "
                  f"{100*report[cohort]['share_of_those_in_mhc']:.0f}%")
        print()

    dest = OUT / "calibration_concordance.json"
    dest.write_text(json.dumps(report, indent=2) + "\n")
    print(f"written: {dest.relative_to(BASE)}")


if __name__ == "__main__":
    main()
