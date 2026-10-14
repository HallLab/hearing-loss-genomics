#!/usr/bin/env python3
"""
Phase 5 -- what would it take to find something? Power against cohort size.

WHY THIS BELONGS TO THE BURDEN ARM AND NOT THE SYNTHESIS. Power for a gene burden test
is not the single-variant calculation with a different label. A single-variant test is
parameterised by one allele's frequency; a burden test is parameterised by the gene's
CUMULATIVE carrier count, which varies over four orders of magnitude across genes. The
ExWAS arm's curve answers a different question and the two cannot share a figure.

THE MODEL. The gene is treated as one super-variant at its observed cumulative
frequency. Then the standard case-control approximation applies:

    ncp  =  (log OR)^2 * 2N * f(1-f) * phi(1-phi)

with f the cumulative minor allele frequency across the mask, phi the case fraction.
Power is the probability a non-central chi-square on 1 df exceeds the critical value.

THIS IS OPTIMISTIC, and deliberately so. It assumes every variant in the mask is causal
and carries the same effect. A real mask contains neutral variants, which dilute the
burden signal, so the true requirement is larger than the curve shows. An optimistic
curve makes a null harder to explain away, not easier.

FREQUENCIES ARE OURS, not assumed. From the pLOF mask at MAF <= 0.01 in the broad
combined cohort, cumulative MAC per gene:

    median gene      49 of 114,996 alleles    f = 4.3e-4
    90th percentile 373                       f = 3.2e-3
    top gene      5,602                       f = 4.9e-2

so the three scenarios are a typical gene, a well-carried one, and the best case.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import brentq
from scipy.stats import chi2, ncx2

HERE = Path(__file__).resolve().parents[1]
PHASE4 = HERE.parent / "phase_4/results"
OUT = HERE / "results"

CASE_FRACTION = 6752 / 57498          # the broad combined prevalence
BAR = 0.05 / 17943                    # the Bonferroni bar this arm uses
POWER_TARGET = 0.80


def power(n, f, log_or, alpha=BAR, phi=CASE_FRACTION):
    ncp = log_or ** 2 * 2 * n * f * (1 - f) * phi * (1 - phi)
    return float(ncx2.sf(chi2.isf(alpha, 1), 1, ncp))


def n_for_power(f, log_or, target=POWER_TARGET):
    lo, hi = 1e3, 1e9
    if power(hi, f, log_or) < target:
        return None
    return brentq(lambda n: power(n, f, log_or) - target, lo, hi)


def main():
    OUT.mkdir(parents=True, exist_ok=True)

    cells = pd.read_csv(PHASE4 / "cells_broad_combined.tsv", sep="\t")
    g = cells[(cells.Group == "pLOF") & (cells.max_MAF == 0.01)]
    two_n = 2 * 57498
    scen = {
        "median gene": g.MAC.median() / two_n,
        "90th percentile gene": g.MAC.quantile(0.90) / two_n,
        "best-carried gene": g.MAC.max() / two_n,
    }

    print("cumulative carrier frequency, from our own pLOF mask:\n")
    for k, f in scen.items():
        print(f"  {k:<24} MAC {f*two_n:>8,.0f}   f = {f:.2e}")

    ors = [1.5, 2.0, 3.0, 5.0]
    out = {"case_fraction": CASE_FRACTION, "bonferroni_bar": BAR,
           "scenarios": {}, "note": "gene treated as one super-variant; assumes every "
                                    "masked variant is causal, so optimistic"}

    print(f"\ncohort size needed for {POWER_TARGET:.0%} power at p < {BAR:.2e}\n")
    print(f"  {'scenario':<24}" + "".join(f"{'OR '+str(o):>14}" for o in ors))
    for k, f in scen.items():
        row = []
        out["scenarios"][k] = {"f": f, "mac": float(f * two_n), "needed": {}}
        for o in ors:
            n = n_for_power(f, np.log(o))
            out["scenarios"][k]["needed"][str(o)] = n
            row.append("—" if n is None else
                       (f"{n/1e6:.1f}M" if n > 1e6 else f"{n/1e3:,.0f}k"))
        print(f"  {k:<24}" + "".join(f"{v:>14}" for v in row))

    print(f"\n  for reference, this cohort is 57,498")

    (OUT / "04_power.json").write_text(json.dumps(out, indent=2) + "\n")
    print(f"\nwritten: results/04_power.json")




# --------------------------------------------------------------------------
def figure():
    """Power against cohort size, in the form Nikki shared, adapted three ways.

    Hers is a linear-regression PheWAS: beta coefficients, samples up to 1,000,
    alpha = 0.05. Ours is a binary burden test, so the effect is an odds ratio; this
    cohort is already 57,498, so the interesting range is 50k to 3M; and the bar is
    0.05/genes = 2.79e-06, not 0.05.

    The fourth adaptation is the one that matters: burden power is not a function of
    sample size alone. It depends on how many carriers the GENE accumulates, which
    spans four orders of magnitude across our genes, so one panel per carrier scenario
    rather than one figure.

    COLOUR. Odds ratio is ORDERED, not categorical, so it takes the documented blue
    ramp rather than four separate hues. These are discrete marks, so the ordinal rule
    applies and the lightest step used is 250, which clears 2:1 on the light surface.
    The validator could not be run here; the steps are the documented ones unchanged.

    BIOBANK MARKERS are approximate and labelled as such. They are drawn to show the
    scale of what a replication cohort would offer, not as exact figures -- confirm
    before quoting either in a slide.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    FIG = HERE / "figures"
    FIG.mkdir(parents=True, exist_ok=True)
    rep = json.loads((OUT / "04_power.json").read_text())

    RAMP = {1.5: "#86b6ef", 2.0: "#3987e5", 3.0: "#1c5cab", 5.0: "#0d366b"}
    INK, INK_2, INK_3, SURFACE = "#1a1a19", "#5c5b55", "#8a8980", "#ffffff"
    RED = "#c0392b"

    plt.rcParams.update({
        "figure.dpi": 110, "savefig.dpi": 190, "savefig.bbox": "tight",
        "font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.alpha": 0.22, "grid.linewidth": 0.5,
        "legend.frameon": False, "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE, "text.color": INK,
        "axes.edgecolor": INK_3, "xtick.color": INK_2, "ytick.color": INK_2,
    })

    ns = np.logspace(np.log10(2e4), np.log10(3e6), 240)
    cohorts = [("PMBB v4\n57,498", 57498, True),
               ("All of Us\n~250k", 250_000, False),
               ("UK Biobank\n~470k", 470_000, False)]

    fig, axes = plt.subplots(1, 3, figsize=(13, 4.3), sharey=True)
    for ax, (name, s) in zip(axes, rep["scenarios"].items()):
        f = s["f"]
        for o, col in RAMP.items():
            ax.plot(ns, [power(n, f, np.log(o)) for n in ns], color=col, lw=2,
                    label=f"OR {o:g}")
        ax.axhline(0.80, color=INK, lw=1.0, ls="--")
        for label, x, ours in cohorts:
            ax.axvline(x, color=RED if ours else INK_3, lw=1.2 if ours else 0.9,
                       ls="-" if ours else ":")
            ax.text(x * 1.08, 0.04, label, fontsize=6.5, rotation=90,
                    color=RED if ours else INK_3, va="bottom")
        ax.set_xscale("log")
        ax.set_xlim(2e4, 3e6)
        ax.set_ylim(0, 1.04)
        ax.set_xlabel("cohort size")
        ax.set_title(f"{name}\n{s['mac']:,.0f} rare alleles in 57,498 people",
                     fontsize=9)
        if ax is axes[0]:
            ax.set_ylabel("power")
            ax.text(2.2e4, 0.835, "80%", fontsize=7.5, color=INK)
    axes[-1].legend(loc="lower right", fontsize=8, title="effect size",
                    title_fontsize=8)

    fig.suptitle("What cohort size would make a gene findable, at p < 2.79e-06",
                 fontsize=11.5, x=0.055, ha="left", y=1.045)
    fig.text(0.055, 0.985,
             "burden power depends on how many carriers the GENE has, not on sample "
             "size alone — so one panel per carrier scenario, taken from our own pLOF "
             "mask. Optimistic: every masked variant is assumed causal.",
             fontsize=8, color=INK_2, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(FIG / "power.png")
    plt.close(fig)
    print("  power.png")


if __name__ == "__main__":
    main()
    figure()
