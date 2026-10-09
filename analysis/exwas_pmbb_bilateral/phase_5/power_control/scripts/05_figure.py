#!/usr/bin/env python3
"""
Power control, step 5 -- the figure for the readout.

One panel per cohort: the observed count against the five random draws. The decision
rule lives in README.md and was fixed before any of this ran; this only draws it.
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

PC = Path(__file__).resolve().parents[1]
FIG = PC.parents[1] / "figures"
BLUE, PALE, RED, GREY = "#3b6ea5", "#9ab4d1", "#c0392b", "#888888"


def main():
    d = json.loads((PC / "results/04_readout.json").read_text())
    cohorts = list(d["cohorts"])
    fig, axes = plt.subplots(1, len(cohorts), figsize=(10.5, 3.6), sharey=False)
    if len(cohorts) == 1:
        axes = [axes]

    plt.rcParams.update({"figure.dpi": 110, "savefig.dpi": 160,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.grid": True, "grid.alpha": 0.25,
                         "legend.frameon": False, "font.size": 9})

    for ax, c in zip(axes, cohorts):
        k = d["cohorts"][c]
        draws = k["k_draws"]
        x = np.arange(len(draws))
        ax.bar(x, draws, color=PALE, width=0.62,
               label="random size-matched draws")
        ax.axhline(k["k_observed"], color=RED, lw=1.6,
                   label=f"bilateral restriction = {k['k_observed']}")
        ax.axhline(k["draw_mean"], color=GREY, lw=1.0, ls="--",
                   label=f"draw mean = {k['draw_mean']:.1f}")
        ax.set_xticks(x)
        ax.set_xticklabels([f"draw {i+1}" for i in x], fontsize=7)
        ax.set_ylabel("variants stronger by >1 log, p < 1e-4"
                      if c == cohorts[0] else "")
        top = max(max(draws), k["k_observed"]) * 1.45 + 1
        ax.set_ylim(0, top)
        for xi, v in zip(x, draws):
            ax.text(xi, v + top * 0.02, str(v), ha="center", fontsize=7)
        verdict = k["verdict"].split(" -- ")[0].split(" (")[0]
        ax.set_title(f"{c}\nPoisson p = {k['poisson_p']:.3f}  ·  {verdict}",
                     fontsize=8.5)
        ax.legend(loc="upper right", fontsize=6.5)

    fig.suptitle("Does restricting to bilateral sensorineural sharpen signal, or is "
                 "it what removing cases does anyway?", fontsize=10)
    fig.tight_layout()
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / "power_control.png")
    print("  power_control.png")


if __name__ == "__main__":
    main()
