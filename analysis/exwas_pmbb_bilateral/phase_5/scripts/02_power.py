#!/usr/bin/env python3
"""
Phase 5, step 4 -- what effect size could this study have detected?

"Nothing was significant" is only informative next to what would have been findable.
Without this, a null reads as "no association exists"; with it, it reads as "no
association of a size this cohort could see", which is a different and honest claim.

THE CALCULATION
For a single-variant additive test, the variance of the log-odds estimate is

    var(log OR)  =  1 / ( 2 N f(1-f) phi(1-phi) )

with N the sample size, f the minor allele frequency, phi the case fraction. The
score statistic is then chi-square on 1 df with non-centrality

    ncp  =  (log OR)^2 / var(log OR)  =  (log OR)^2 * 2 N f(1-f) phi(1-phi)

Eighty-percent power at a threshold alpha means finding the ncp for which a
non-central chi-square exceeds the alpha critical value 80% of the time, then solving
back for log OR.

This is the standard large-sample approximation. It assumes an additive model, no
covariate adjustment loss and no relatedness, so it is mildly OPTIMISTIC -- the true
detectable effect is a little larger than reported here. Stated that way round on
purpose: an optimistic power curve makes the null harder to explain away, not easier.
"""
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.stats import chi2, ncx2

BASE = Path(__file__).resolve().parents[2]
OUT = BASE / "phase_5/results"

# From PREMISES P1, the cohorts this arm inherited.
COHORTS = {
    "combined": {"n": 53_910, "cases": 3_164},
    "EUR": {"n": 40_143, "cases": 2_547},
    "AFR": {"n": 10_539, "cases": 490},
}
BARS = {"primary 5e-8": 5e-8}
POWER = 0.80


def ncp_for_power(alpha, power=POWER):
    crit = chi2.isf(alpha, 1)
    return brentq(lambda l: ncx2.sf(crit, 1, l) - power, 1e-6, 1e4)


def detectable_or(maf, n, cases, alpha, power=POWER):
    phi = cases / n
    info = 2 * n * maf * (1 - maf) * phi * (1 - phi)
    return float(np.exp(np.sqrt(ncp_for_power(alpha, power) / info)))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    mafs = np.concatenate([np.linspace(0.01, 0.05, 40), np.linspace(0.05, 0.5, 90)])
    out = {"power": POWER, "assumptions":
           "additive, large-sample, no covariate or relatedness loss -- optimistic",
           "cohorts": {}}

    print(f"Detectable odds ratio at {POWER:.0%} power, alpha = 5e-8\n")
    print(f"  {'cohort':<10}{'N':>9}{'cases':>8}   " +
          "".join(f"MAF {m:>5.0%}" .replace("MAF ", "MAF=") + "  " for m in
                  (0.01, 0.05, 0.10, 0.25, 0.50)))
    for name, c in COHORTS.items():
        curve = [detectable_or(m, c["n"], c["cases"], 5e-8) for m in mafs]
        out["cohorts"][name] = {
            **c, "maf": mafs.tolist(), "detectable_or": curve,
            "at": {f"{m}": detectable_or(m, c["n"], c["cases"], 5e-8)
                   for m in (0.01, 0.05, 0.10, 0.25, 0.50)},
        }
        row = "".join(f"{detectable_or(m, c['n'], c['cases'], 5e-8):>9.2f}  "
                      for m in (0.01, 0.05, 0.10, 0.25, 0.50))
        print(f"  {name:<10}{c['n']:>9,}{c['cases']:>8,}   {row}")

    dest = OUT / "power.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print(f"\nwritten: {dest.relative_to(BASE)}")


if __name__ == "__main__":
    main()
