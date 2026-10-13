#!/usr/bin/env python3
"""
Synthesis, step 2 -- the one figure for the presentation.

FORM. The question the figure answers is "how close did anything come?", so the
quantity plotted is a ratio: each run's headline p-value divided by the bar that run
would have been judged against. 1.0 means exactly at the threshold. Everything is
above 1, which is the headline, and a shared log axis across both panels is what makes
the second headline visible -- the burden arms sit an order of magnitude closer to
their bar than the ExWAS arms do.

A ratio is used rather than raw p-values because the bars differ: the ExWAS bar is the
5e-8 convention while the burden bar is 0.05/genes. Plotting raw p would compare numbers
that were never judged by the same rule.

THE BURDEN PANEL USES LIKE-FOR-LIKE NUMBERS, not each run's own reporting. Elena's
delivered tables report per-cell p-values corrected within one MAF family
(0.05/53,583 for combined), while our arms report a Cauchy omnibus per gene corrected
across genes. Comparing her minimum cell p against our minimum omnibus p would be
comparing a minimum over nine times more tests. So the panel uses the omnibus
recomputed on her data by our method, from 09_three_way.json, where all three arms are
judged by the same rule. Under her own delivered bar her best cell is 3.6x short in
combined rather than 1.9x; that number belongs in the table, not in a figure that
places three runs on one axis.

COLOR. Three series, so the first three categorical slots of the documented palette,
in fixed order, values unchanged. This is a small-multiples form, which uses the
all-pairs pairlist; the palette reference states that the first three slots validate
all-pairs in both modes (worst CVD Delta E 9.2 light, normal-vision 24.0). The
validator could not be run here -- node is not installed on this cluster -- so no
re-stepping was attempted and no fourth slot is used; the documented three-slot
guarantee is the whole basis for the choice.

Slot 3 (aqua) sits below 3:1 on the light surface, which triggers the relief rule:
every mark carries a visible direct label, and the document this figure belongs to is
itself the table view.
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parents[1]
COHORTS = ["combined", "EUR", "AFR"]
RUNS = ["Elena", "Replication", "Bilateral"]

# Documented categorical slots 1-3, light surface. Not re-stepped.
SERIES = {"Elena": "#2a78d6", "Replication": "#eb6834", "Bilateral": "#1baf7a"}
SURFACE = "#ffffff"
INK, INK_2, INK_3 = "#1a1a19", "#5c5b55", "#8a8980"


def ratio(cell):
    p = cell.get("min_omnibus_p") or cell.get("min_cell_p") or cell.get("min_p")
    bar = cell.get("bar") or cell.get("bar_primary")
    return p / bar


def burden_like_for_like():
    """All three arms' omnibus under one bar rule, from the burden arm's own
    three-way comparison. Not each arm's own reporting -- see the module docstring."""
    f = (HERE.parent / "gene_burden_pmbb_bilateral/phase_5/results/09_three_way.json")
    d = json.loads(f.read_text())
    return {c: {a: k["min_omnibus_p"] / k["bar"] for a, k in v["arms"].items()}
            for c, v in d["cohorts"].items()}


def main():
    d = json.loads((HERE / "results/consolidated.json").read_text())
    burden = burden_like_for_like()

    plt.rcParams.update({
        "figure.dpi": 110, "savefig.dpi": 200, "savefig.bbox": "tight",
        "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.spines.left": False,
        "axes.edgecolor": INK_3, "axes.linewidth": 0.8,
        "xtick.color": INK_2, "ytick.color": INK, "text.color": INK,
        "legend.frameon": False, "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
    })

    fig, axes = plt.subplots(2, 1, figsize=(9.6, 5.4), sharex=True,
                             gridspec_kw={"hspace": 0.28})

    for ax, (title, key, note) in zip(axes, [
        ("Rare-variant gene burden", "burden",
         "SAIGE-GENE · MAF ≤ 0.01 / 0.001 / 0.0001 · "
         "all three judged by one bar rule, 0.05/genes"),
        ("Common-variant ExWAS", "exwas",
         "SAIGE single variant · MAF ≥ 0.01 · bar 5e-8"),
    ]):
        ypos = {c: i for i, c in enumerate(COHORTS[::-1])}
        for ci, c in enumerate(COHORTS[::-1]):
            vals = ([(r, burden[c][r]) for r in RUNS] if key == "burden"
                    else [(r, ratio(d[key][c][r])) for r in RUNS])
            # connector, so the three runs in a cohort read as one comparison
            ax.plot([min(v for _, v in vals), max(v for _, v in vals)],
                    [ci, ci], color=INK_3, lw=2, solid_capstyle="round",
                    zorder=1, alpha=0.35)
            for r, v in vals:
                ax.plot(v, ci, "o", ms=9, color=SERIES[r],
                        markeredgecolor=SURFACE, markeredgewidth=2, zorder=3)
            # direct labels: relief rule for slot 3, and the numbers are the point
            order = sorted(vals, key=lambda t: t[1])
            for j, (r, v) in enumerate(order):
                ax.annotate(f"{v:.0f}×" if v >= 10 else f"{v:.1f}×",
                            (v, ci), textcoords="offset points",
                            xytext=(0, 11 if j % 2 == 0 else -17),
                            ha="center", fontsize=7.5, color=INK_2)
        ax.axvline(1.0, color=INK, lw=1.2, zorder=2)
        ax.set_yticks(range(len(COHORTS)))
        ax.set_yticklabels(COHORTS[::-1])
        ax.set_xscale("log")
        ax.set_xlim(0.72, 1500)
        ax.set_ylim(-0.62, len(COHORTS) - 0.38)
        ax.grid(axis="x", alpha=0.18, lw=0.5)
        ax.set_title(f"{title}", loc="left", fontweight="bold", pad=26)
        ax.text(0, 1.035, note, transform=ax.transAxes, fontsize=7.5,
                color=INK_3, va="bottom")

    for ax in axes:
        ax.annotate("threshold", xy=(1.0, -0.52), xytext=(1.0, -0.52),
                    ha="center", va="center", fontsize=7, color=INK)
    axes[-1].set_xlabel("how far the best result fell short of its significance "
                        "bar  (fold, log scale)")

    handles = [plt.Line2D([], [], marker="o", ls="", ms=9, color=SERIES[r],
                          markeredgecolor=SURFACE, markeredgewidth=2, label=r)
               for r in RUNS]
    axes[0].legend(handles=handles, loc="lower right", bbox_to_anchor=(1.0, 1.06),
                   ncol=3, fontsize=8.5, handletextpad=0.4, columnspacing=1.4)

    fig.suptitle("No run, in either analysis, in any cohort, reaches its own "
                 "significance threshold", fontsize=11.5, x=0.062, ha="left", y=1.105)
    fig.text(0.062, 1.040, "and the rare-variant arms come roughly ten times "
             "closer than the common-variant ones", fontsize=9, color=INK_2, ha="left")

    out = HERE / "figures/how_close.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out)
    plt.close(fig)
    print(f"  {out.relative_to(HERE.parent)}")


if __name__ == "__main__":
    main()
