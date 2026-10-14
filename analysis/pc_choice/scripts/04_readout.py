#!/usr/bin/env python3
"""
PC choice experiment, step 4 -- lambda by region, for the three arms.

THE READOUT. Lambda inside the MHC against lambda outside it. The current arm sits at
1.339 inside and 0.995 outside. Whichever change brings the inside figure to ~1.0
without moving the outside one has identified the problem.

A CAVEAT THE TABLE CANNOT CARRY. Arm B tests 196 fewer people, because they have no
imputed PCs. So B's p-values are not comparable to the others variant by variant, and
B's variant set differs slightly too (MAF is computed on a different sample). Lambda is
a distributional property and remains comparable; individual p-values are not, and the
minimum p column below is reported for orientation rather than for comparison.

Arms A and the current one share their sample set exactly, so those two ARE comparable
throughout -- which is the cleaner of the two contrasts.
"""
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parent
CURRENT = ROOT / "exwas_elena_replication/phase_4/results/merged_EUR.tsv"
CHI2_MEDIAN = chi2.ppf(0.5, 1)
MHC = (28_510_000, 33_400_000)
ARMS = [("A_exome10", "10 PCs, our exome PCA"),
        ("B_imputed9", "9 PCs, release imputed PCA")]


def lam(p):
    p = np.asarray([v for v in p if 0 < v <= 1])
    return float(np.median(chi2.isf(p, 1)) / CHI2_MEDIAN)


def in_flight():
    live = []
    for n in ("pcx_s1", "pcx_s2"):
        try:
            out = subprocess.run(["bjobs", "-J", n], capture_output=True,
                                 text=True, timeout=60).stdout
        except (FileNotFoundError, subprocess.TimeoutExpired):
            return None
        live += [l for l in out.splitlines()[1:] if l.strip()]
    return live


def load_arm(arm):
    parts = []
    for c in range(1, 23):
        f = HERE / f"results/step2/{arm}/{arm}_chr{c}.txt"
        if not f.exists():
            raise FileNotFoundError(f)
        d = pd.read_csv(f, sep="\t", usecols=["CHR", "POS", "MarkerID", "p.value"])
        if d.empty:
            raise ValueError(f"{f} is empty")
        parts.append(d)
    return pd.concat(parts, ignore_index=True)


def regions(d):
    d = d[(d["p.value"] > 0) & (d["p.value"] <= 1)]
    mhc = (d["CHR"].astype(str) == "6") & d["POS"].between(*MHC)
    return {
        "tested": int(len(d)),
        "min_p": float(d["p.value"].min()),
        "lambda_all": lam(d["p.value"]),
        "lambda_outside_mhc": lam(d.loc[~mhc, "p.value"]),
        "lambda_inside_mhc": lam(d.loc[mhc, "p.value"]),
        "n_mhc": int(mhc.sum()),
    }


def main():
    live = in_flight()
    if live:
        sys.exit(f"REFUSING: {len(live)} experiment tasks still in flight.")

    out = {}
    cur = pd.read_csv(CURRENT, sep="\t", usecols=["CHR", "POS", "MarkerID", "p.value"])
    out["current_exome4"] = {**regions(cur), "label": "4 PCs, our exome PCA (current)"}

    for arm, label in ARMS:
        try:
            out[arm] = {**regions(load_arm(arm)), "label": label}
        except (FileNotFoundError, ValueError) as e:
            sys.exit(f"REFUSING: {arm} incomplete -- {e}")

    print(f"  {'arm':<34}{'tested':>9}{'λ all':>9}{'λ outside':>11}"
          f"{'λ MHC':>9}{'min p':>11}")
    for k, v in out.items():
        print(f"  {v['label']:<34}{v['tested']:>9,}{v['lambda_all']:>9.4f}"
              f"{v['lambda_outside_mhc']:>11.4f}{v['lambda_inside_mhc']:>9.4f}"
              f"{v['min_p']:>11.2e}")

    base = out["current_exome4"]["lambda_inside_mhc"]
    print(f"\n  MHC inflation, current: {base:.4f}")
    for k, v in out.items():
        if k == "current_exome4":
            continue
        moved = base - v["lambda_inside_mhc"]
        verdict = ("resolves it" if abs(v["lambda_inside_mhc"] - 1) < 0.08
                   else "improves it" if moved > 0.1 else "does not fix it")
        print(f"    {v['label']:<34} {v['lambda_inside_mhc']:.4f}  ({verdict})")

    dest = HERE / "results/04_readout.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print(f"\n  written: {dest.relative_to(HERE)}")


if __name__ == "__main__":
    main()
