#!/usr/bin/env python3
"""
Power control, step 4 -- apply the decision rule in README.md.

The rule was written before any of this ran and is not restated here as a choice:

  k_obs <= max(draws)                     -> selection, not subtype. Not reported.
  k_obs > max(draws) and Poisson p < 0.05 -> suggestive subtype concentration.
  otherwise                               -> inconclusive, both counts shown.

CRITERION, identical for the observed arm and for every draw: a variant counts when
its p-value in the restricted run is more than one log unit smaller than in the broad
run AND below 1e-4. The broad run is the same file in all comparisons.
"""
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import poisson

PC = Path(__file__).resolve().parents[1]
ARM = PC.parents[1]
BROAD = ARM.parent / "exwas_elena_replication"
COHORTS = ["combined", "EUR", "AFR"]
CHROMS = range(1, 23)
N_DRAWS = 5
GAIN = 1.0
P_CEIL = 1e-4


def in_flight():
    live = []
    for n in ("pc_s1", "pc_s2", "pc_s2_retry"):
        try:
            out = subprocess.run(["bjobs", "-J", n], capture_output=True,
                                 text=True, timeout=60).stdout
        except (FileNotFoundError, subprocess.TimeoutExpired):
            return None
        live += [l for l in out.splitlines()[1:] if l.strip()]
    return live


def load_broad(cohort):
    d = pd.read_csv(BROAD / f"phase_4/results/merged_{cohort}.tsv", sep="\t",
                    usecols=["MarkerID", "p.value"])
    return d[(d["p.value"] > 0) & (d["p.value"] <= 1)].rename(
        columns={"p.value": "p_broad"})


def load_restricted(cohort):
    d = pd.read_csv(ARM / f"phase_4/results/merged_{cohort}.tsv", sep="\t",
                    usecols=["MarkerID", "p.value"])
    return d[(d["p.value"] > 0) & (d["p.value"] <= 1)].rename(
        columns={"p.value": "p_r"})


def load_draw(cohort, draw):
    parts = []
    for c in CHROMS:
        f = PC / f"results/step2/{cohort}_draw{draw}/{cohort}_draw{draw}_chr{c}.txt"
        if not f.exists():
            raise FileNotFoundError(f)
        d = pd.read_csv(f, sep="\t", usecols=["MarkerID", "p.value"])
        if d.empty:
            raise ValueError(f"{f} has no rows")
        parts.append(d)
    d = pd.concat(parts, ignore_index=True)
    return d[(d["p.value"] > 0) & (d["p.value"] <= 1)].rename(
        columns={"p.value": "p_r"})


def count(restricted, broad):
    m = restricted.merge(broad, on="MarkerID")
    gain = -np.log10(m["p_r"]) + np.log10(m["p_broad"])
    return int(((gain > GAIN) & (m["p_r"] < P_CEIL)).sum()), len(m)


def main():
    live = in_flight()
    if live:
        sys.exit(f"REFUSING: {len(live)} power-control tasks still in flight.")
    if live is None:
        print("WARNING: could not query LSF\n")

    out = {"criterion": {"log10_gain_over": GAIN, "p_below": P_CEIL},
           "n_draws": N_DRAWS, "cohorts": {}}
    missing = []

    for cohort in COHORTS:
        broad = load_broad(cohort)
        k_obs, n_obs = count(load_restricted(cohort), broad)

        draws = []
        for d in range(1, N_DRAWS + 1):
            try:
                k, n = count(load_draw(cohort, d), broad)
            except (FileNotFoundError, ValueError) as e:
                missing.append(f"{cohort} draw {d}: {e}")
                continue
            draws.append(k)

        if len(draws) < N_DRAWS:
            continue

        pooled_rate = float(np.mean(draws))
        # one-sided Poisson: how surprising is k_obs against the draws' rate?
        p_pois = float(poisson.sf(k_obs - 1, max(pooled_rate, 1e-9)))

        if k_obs <= max(draws):
            verdict = "SELECTION -- within what random case removal produces"
        elif p_pois < 0.05:
            verdict = "SUGGESTIVE subtype concentration (never significant: no 5e-8)"
        else:
            verdict = "INCONCLUSIVE"

        out["cohorts"][cohort] = {
            "k_observed": k_obs, "n_compared": n_obs,
            "k_draws": draws, "draw_mean": pooled_rate,
            "draw_max": max(draws), "poisson_p": p_pois, "verdict": verdict,
        }
        print(f"{cohort}")
        print(f"  bilateral restriction  k = {k_obs}")
        print(f"  five random draws      k = {draws}   mean {pooled_rate:.1f}  "
              f"max {max(draws)}")
        print(f"  Poisson p              {p_pois:.4f}")
        print(f"  -> {verdict}\n")

    if missing:
        print("REFUSING to write a verdict -- incomplete draws:", file=sys.stderr)
        for m in missing[:10]:
            print(f"  {m}", file=sys.stderr)
        sys.exit(1)

    dest = PC / "results/04_readout.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print(f"written: {dest.relative_to(ARM)}")


if __name__ == "__main__":
    main()
