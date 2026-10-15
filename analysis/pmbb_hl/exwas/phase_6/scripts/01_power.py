#!/usr/bin/env python3
"""
ExWAS phase 6 -- what would it take to find something? Power against cohort size.

WHY THIS IS SEPARATE FROM THE BURDEN ARM'S. A single-variant test is parameterised by
ONE allele's frequency; a burden test by a gene's cumulative carrier count. The two
curves answer different questions and cannot share a figure. Here the x of the
parameterisation is simply MAF, and every tested variant is above 1% by construction --
so unlike the burden arm, there is no long tail of untestable units.

THE MODEL. The gene is treated as one super-variant at its observed cumulative
frequency. Then the standard case-control approximation applies:

    ncp  =  (log OR)^2 * 2N * f(1-f) * phi(1-phi)

with f the cumulative minor allele frequency across the mask, phi the case fraction.
Power is the probability a non-central chi-square on 1 df exceeds the critical value.

NO OPTIMISM CORRECTION IS NEEDED HERE, unlike the burden arm. There is no mask to
dilute: the test is on one variant at its own frequency, so the calculation is the
standard one with no causal-fraction assumption buried in it.

PREVALENCE DIFFERS BY ARM and is reported for all three, because phi(1-phi) enters
the non-centrality directly: the same variant at the same odds ratio needs roughly
twice the cohort for tinnitus (5.3% cases) as for broad hearing loss (11.7%).

FREQUENCIES ARE OURS, taken from what the ExWAS actually tested: the 1% floor itself,
the median tested variant, and a common one. Nothing below 1% is tested here, so the
scenarios span the real range rather than a hypothetical one.
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

# PREVALENCE IS NOT SHARED ACROSS THE ARMS, and it drives the answer. Information
# goes as phi(1-phi), so tinnitus at 5.3% carries 0.0505 against broad's 0.1035 --
# a factor of 2.05 in required cohort size, not a rounding difference. The curve is
# drawn at the broad prevalence and the table is reported at all three.
PREVALENCE = {"broad": 6752 / 57498,
              "bilateral": 3164 / 53910,
              "tinnitus": 2859 / 53605}
CASE_FRACTION = PREVALENCE["broad"]
BAR = 5e-8                            # the primary bar: the convention
POWER_TARGET = 0.80


def power(n, f, log_or, alpha=BAR, phi=CASE_FRACTION):
    ncp = log_or ** 2 * 2 * n * f * (1 - f) * phi * (1 - phi)
    return float(ncx2.sf(chi2.isf(alpha, 1), 1, ncp))


def n_for_power(f, log_or, target=POWER_TARGET, phi=CASE_FRACTION):
    lo, hi = 1e3, 1e9
    if power(hi, f, log_or, phi=phi) < target:
        return None
    return brentq(lambda n: power(n, f, log_or, phi=phi) - target, lo, hi)


def main():
    OUT.mkdir(parents=True, exist_ok=True)

    d = pd.read_csv(PHASE4 / "merged_broad_combined.tsv", sep="\t",
                    usecols=["AF_Allele2"])
    maf = np.minimum(d.AF_Allele2, 1 - d.AF_Allele2)
    scen = {
        "at the 1% floor": 0.01,
        "median tested variant": float(maf.median()),
        "common variant (MAF 25%)": 0.25,
    }

    print("minor allele frequency, from what the ExWAS actually tested:\n")
    for k, f in scen.items():
        print(f"  {k:<28} MAF = {f:.4f}")

    ors = [1.05, 1.1, 1.2, 1.5]
    out = {"case_fraction": CASE_FRACTION, "bonferroni_bar": BAR,
           "scenarios": {}, "note": "gene treated as one super-variant; assumes every "
                                    "masked variant is causal, so optimistic"}

    out["prevalence"] = PREVALENCE
    for pheno, phi in PREVALENCE.items():
        print(f"\ncohort size needed for {POWER_TARGET:.0%} power at p < {BAR:.2e}"
              f"   [{pheno}, {phi:.1%} cases]\n")
        print(f"  {'scenario':<24}" + "".join(f"{'OR '+str(o):>14}" for o in ors))
        for k, f in scen.items():
            row = []
            out["scenarios"].setdefault(k, {"f": f, "needed": {}})
            out["scenarios"][k]["needed"].setdefault(pheno, {})
            for o in ors:
                n = n_for_power(f, np.log(o), phi=phi)
                out["scenarios"][k]["needed"][pheno][str(o)] = n
                row.append("—" if n is None else
                           (f"{n/1e6:.1f}M" if n > 1e6 else f"{n/1e3:,.0f}k"))
            print(f"  {k:<24}" + "".join(f"{v:>14}" for v in row))

    print(f"\n  for reference, the arms as run are 57,498 / 53,910 / 53,605")

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

    RAMP = {1.05: "#86b6ef", 1.1: "#3987e5", 1.2: "#1c5cab", 1.5: "#0d366b"}
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
        ax.set_title(f"{name}\nMAF {s['f']:.3f}", fontsize=9)
        if ax is axes[0]:
            ax.set_ylabel("power")
            ax.text(2.2e4, 0.835, "80%", fontsize=7.5, color=INK)
    axes[-1].legend(loc="lower right", fontsize=8, title="effect size",
                    title_fontsize=8)

    fig.suptitle("What cohort size would make a common variant findable, at p < 5e-8",
                 fontsize=11.5, x=0.055, ha="left", y=1.045)
    fig.text(0.055, 0.985,
             "one panel per allele frequency, taken from what this ExWAS tested. "
             "Common-variant effects on complex age-related traits sit at OR 1.05 to "
             "1.15, which is the range these curves are drawn over.",
             fontsize=8, color=INK_2, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(FIG / "power.png")
    plt.close(fig)
    print("  power.png")


if __name__ == "__main__":
    main()
    figure()
