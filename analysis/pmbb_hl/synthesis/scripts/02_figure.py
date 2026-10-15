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

COLOUR. Four series, four documented categorical slots, values unchanged. The validator
could not be run -- node is not installed on this cluster -- and four series stays
inside the all-pairs limit.
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parents[1]
COHORTS = ["combined", "EUR", "AFR"]
# Four series, not four phenotypes. The original run and our broad arm share a
# phenotype and differ in pipeline; bilateral and tinnitus are phenotypes we added, and
# the original ran neither, so there is no second original series.
SERIES = ["original run", "rebuild, broad", "rebuild, bilateral", "rebuild, tinnitus"]
BLUE, ORANGE, AQUA, PLUM = "#2a78d6", "#eb6834", "#1baf7a", "#a855a0"
INK, INK_2, INK_3, SURFACE = "#1a1a19", "#5c5b55", "#8a8980", "#ffffff"
PCOL = {"original run": BLUE, "rebuild, broad": ORANGE,
        "rebuild, bilateral": AQUA, "rebuild, tinnitus": PLUM}


DY = (26, -27, 11, -12)
assert len(DY) == len(SERIES), \
    f"REFUSING: {len(DY)} label heights for {len(SERIES)} series -- labels would " \
    "share an offset and overlap wherever two arms agree closely"


def main():
    raw = json.loads((HERE / "results/consolidated.json").read_text())
    d, elena = raw["cells"], raw["elena"]

    def short(series, arm, cohort):
        if series == "original run":
            return elena[cohort][arm]["short_by"]
        ph = series.split(", ")[1]
        return d[f"{arm}_{ph}_{cohort}"]["short_by"]
    plt.rcParams.update({
        "figure.dpi": 110, "savefig.dpi": 200, "savefig.bbox": "tight",
        "font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
        "axes.spines.left": False, "axes.edgecolor": INK_3, "axes.linewidth": 0.8,
        "xtick.color": INK_2, "ytick.color": INK, "text.color": INK,
        "legend.frameon": False, "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
    })

    fig, axes = plt.subplots(2, 1, figsize=(10.4, 6.2), sharex=True,
                             gridspec_kw={"hspace": 0.42})
    panels = [("Rare-variant gene burden", "burden",
               "SAIGE-GENE · MAF ≤ 0.01 / 0.001 / 0.0001 · "
               "bar 0.05/genes, plus FDR q < 0.1"),
              ("Common-variant ExWAS", "exwas",
               "SAIGE single variant · MAF ≥ 0.01 · bar 5e-8")]

    for ax, (title, arm, note) in zip(axes, panels):
        for ci, cohort in enumerate(COHORTS[::-1]):
            vals = [(s, short(s, arm, cohort)) for s in SERIES]
            ax.plot([min(v for _, v in vals), max(v for _, v in vals)],
                    [ci, ci], color=INK_3, lw=2, alpha=0.35,
                    solid_capstyle="round", zorder=1)
            for j, (p, v) in enumerate(sorted(vals, key=lambda t: t[1])):
                ax.plot(v, ci, "o", ms=9, color=PCOL[p], markeredgecolor=SURFACE,
                        markeredgewidth=2, zorder=3,
                        label=None)
                # ONE DISTINCT HEIGHT PER SERIES. With three series sharing two
                # offsets, runs that agree closely -- 397x against 495x, 74x against
                # 82x -- overlapped into an unreadable smear. A fourth series needs a
                # fourth height for the same reason, so this tuple is indexed by rank
                # and must stay as long as SERIES.
                dy = DY[j]
                ax.annotate(f"{v:.0f}×" if v >= 10 else f"{v:.1f}×",
                            (v, ci), textcoords="offset points", xytext=(0, dy),
                            ha="center", fontsize=7.5, color=INK_2,
                            arrowprops=dict(arrowstyle="-", lw=0.4, color=INK_3,
                                            shrinkA=0, shrinkB=4)
                            if abs(dy) > 15 else None)
        ax.axvline(1.0, color=INK, lw=1.2, zorder=2)
        # The y floor has to clear the deepest label offset (-27pt) plus the
        # "threshold" tick below it, or both are cut off by the axes box.
        ax.annotate("threshold", xy=(1.0, -0.95), ha="center", va="center",
                    fontsize=7, color=INK)
        ax.set_yticks(range(3))
        ax.set_yticklabels(COHORTS[::-1])
        ax.set_xscale("log")
        ax.set_xlim(0.72, 1400)
        ax.set_ylim(-1.12, 2.80)
        ax.grid(axis="x", alpha=0.18, lw=0.5)
        ax.set_title(title, loc="left", fontweight="bold", pad=30)
        ax.text(0, 1.035, note, transform=ax.transAxes, fontsize=7.5,
                color=INK_3, va="bottom")

    axes[-1].set_xlabel("how far the best result fell short of its own bar "
                        "(fold, log scale)")
    # Built by hand so the key reads in SERIES order. Drawing order is by value
    # within each row, so labelling the markers gave whatever order the first row
    # happened to sort into.
    from matplotlib.lines import Line2D
    axes[0].legend(handles=[Line2D([], [], marker="o", ls="", ms=8,
                                   color=PCOL[s_], label=s_) for s_ in SERIES],
                   loc="lower right", bbox_to_anchor=(1.0, 1.07), ncol=4,
                   fontsize=8.5, handletextpad=0.4, columnspacing=1.3)
    fig.suptitle("Twenty-four cells across four arms. Nothing reaches its threshold "
                 "in any of them.", fontsize=11.5, x=0.062, ha="left", y=1.135)
    fig.text(0.062, 1.072, "the rare-variant arms land an order of magnitude closer "
             "than the common-variant ones, in every run",
             fontsize=9, color=INK_2, ha="left")
    fig.text(0.062, 1.020, "the original run\u2019s burden figure is its omnibus "
             "recomputed by our method \u2014 its delivered tables report per-cell "
             "p-values corrected within one MAF family, which is not the same quantity",
             fontsize=7.5, color=INK_3, ha="left")

    out = HERE / "figures/how_close.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out)
    plt.close(fig)
    print(f"  {out.relative_to(HERE)}")


if __name__ == "__main__":
    main()
