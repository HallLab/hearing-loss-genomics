#!/usr/bin/env python3
"""
Synthesis -- the one figure.

FORM. The question is "how close did anything come?", so the quantity is a ratio: each
cell's best p divided by the bar it is judged against. 1.0 is the threshold; everything
is above it, which is the headline. A shared log axis across both panels makes the
second headline visible -- the burden cells sit an order of magnitude closer than the
ExWAS ones.

A ratio rather than raw p, because the two arms are not judged by the same rule: the
burden bar is 0.05/genes and the ExWAS bar is the 5e-8 convention. Raw minima would
compare a minimum over 17,900 genes with one over 70,000-118,000 variants.

COLOUR. Two phenotypes, two documented categorical slots, values unchanged. The
validator could not be run -- node is not installed on this cluster -- and with two
series the all-pairs limit is not approached.
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parents[1]
COHORTS = ["combined", "EUR", "AFR"]
PHENOS = ["broad", "bilateral"]
BLUE, ORANGE = "#2a78d6", "#eb6834"
INK, INK_2, INK_3, SURFACE = "#1a1a19", "#5c5b55", "#8a8980", "#ffffff"
PCOL = {"broad": BLUE, "bilateral": ORANGE}


def main():
    d = json.loads((HERE / "results/consolidated.json").read_text())["cells"]
    plt.rcParams.update({
        "figure.dpi": 110, "savefig.dpi": 200, "savefig.bbox": "tight",
        "font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
        "axes.spines.left": False, "axes.edgecolor": INK_3, "axes.linewidth": 0.8,
        "xtick.color": INK_2, "ytick.color": INK, "text.color": INK,
        "legend.frameon": False, "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
    })

    fig, axes = plt.subplots(2, 1, figsize=(9.8, 5.2), sharex=True,
                             gridspec_kw={"hspace": 0.34})
    panels = [("Rare-variant gene burden", "burden",
               "SAIGE-GENE · MAF ≤ 0.01 / 0.001 / 0.0001 · "
               "bar 0.05/genes, plus FDR q < 0.1"),
              ("Common-variant ExWAS", "exwas",
               "SAIGE single variant · MAF ≥ 0.01 · bar 5e-8")]

    for ax, (title, arm, note) in zip(axes, panels):
        for ci, cohort in enumerate(COHORTS[::-1]):
            vals = [(p, d[f"{arm}_{p}_{cohort}"]["short_by"]) for p in PHENOS]
            ax.plot([min(v for _, v in vals), max(v for _, v in vals)],
                    [ci, ci], color=INK_3, lw=2, alpha=0.35,
                    solid_capstyle="round", zorder=1)
            for j, (p, v) in enumerate(sorted(vals, key=lambda t: t[1])):
                ax.plot(v, ci, "o", ms=9, color=PCOL[p], markeredgecolor=SURFACE,
                        markeredgewidth=2, zorder=3,
                        label=p if (ci == 2 and arm == "burden") else None)
                ax.annotate(f"{v:.0f}×" if v >= 10 else f"{v:.1f}×",
                            (v, ci), textcoords="offset points",
                            xytext=(0, 11 if j == 0 else -17),
                            ha="center", fontsize=7.5, color=INK_2)
        ax.axvline(1.0, color=INK, lw=1.2, zorder=2)
        ax.annotate("threshold", xy=(1.0, -0.52), ha="center", va="center",
                    fontsize=7, color=INK)
        ax.set_yticks(range(3))
        ax.set_yticklabels(COHORTS[::-1])
        ax.set_xscale("log")
        ax.set_xlim(0.72, 1400)
        ax.set_ylim(-0.62, 2.62)
        ax.grid(axis="x", alpha=0.18, lw=0.5)
        ax.set_title(title, loc="left", fontweight="bold", pad=26)
        ax.text(0, 1.035, note, transform=ax.transAxes, fontsize=7.5,
                color=INK_3, va="bottom")

    axes[-1].set_xlabel("how far the best result fell short of its own bar "
                        "(fold, log scale)")
    axes[0].legend(loc="lower right", bbox_to_anchor=(1.0, 1.06), ncol=2,
                   fontsize=8.5, handletextpad=0.4, columnspacing=1.4)
    fig.suptitle("Twelve cells. Nothing reaches its threshold in any of them.",
                 fontsize=11.5, x=0.062, ha="left", y=1.105)
    fig.text(0.062, 1.040, "and the rare-variant arms land an order of magnitude "
             "closer than the common-variant ones", fontsize=9, color=INK_2, ha="left")

    out = HERE / "figures/how_close.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out)
    plt.close(fig)
    print(f"  {out.relative_to(HERE)}")


if __name__ == "__main__":
    main()
