#!/usr/bin/env python3
"""
Synthesis -- the clean comparison: one phenotype, one method, three annotations.

WHY THIS EXISTS BESIDE how_close.png. That figure answers "did anything reach its
threshold" across every arm the project ran, which means it mixes phenotypes, methods
and annotations in one picture. It is the right figure for that question and the wrong
one for this: how much does the ANNOTATION change the answer? That needs everything
else held fixed.

Held fixed here: the broad hearing-loss phenotype, the rare-variant burden method, the
same genotypes, the same covariates, the same GRM, the same fitted null models. The
only thing that differs between the three series is which variants the mask puts in
each gene.

BURDEN ONLY, AND THAT IS NOT A GAP. The common-variant ExWAS has no mask at all -- it
tests every variant above 1% regardless of consequence -- so an annotation standard
cannot change it. An ExWAS panel here would show a comparison that does not exist.

AFR IS MISSING FROM THE BRaVa SERIES and is drawn as missing. pmbb_hl_brava/PREMISES.md
P3 scoped that arm to combined and EUR, the two cells where the ClinGen enrichment
lives, before the run. The hole is a scope decision, not a result, and the figure says
so rather than leaving a reader to guess.

THE BARS DIFFER AND THE RATIO IS WHY. Each arm's Bonferroni bar is 0.05 divided by the
genes it tested, and the annotations test different numbers of genes (17,943 ours,
19,042 the original run, 16,742 BRaVa). Plotting raw minimum p would compare minima
taken over different numbers of tests. The ratio to each arm's OWN bar is the quantity
that travels.

COLOUR. Slots 1-3 of the documented categorical order, in order. Validated with the
skill's Python validator (node is absent here): all checks pass, worst adjacent CVD
deltaE 9.2. The contrast WARN on slot 3 is relieved by the direct label every point
carries.
"""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parents[1]
COHORTS = ["combined", "EUR", "AFR"]
SERIES = ["original run", "rebuild", "BRaVa annotation"]
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK_2, INK_3, SURFACE = "#1a1a19", "#5c5b55", "#8a8980", "#ffffff"
PCOL = dict(zip(SERIES, (BLUE, ORANGE, AQUA)))
DY = (24, -25, 10)
assert len(DY) == len(SERIES), "one label height per series, or labels overlap"


def main():
    raw = json.loads((HERE / "results/consolidated.json").read_text())
    if "brava" not in raw:
        sys.exit("REFUSING: no BRaVa arm in consolidated.json -- run 01_collect.py "
                 "after pmbb_hl_brava/phase_4 has merged.")
    bv = raw["brava"]["cells"]

    def cell(series, cohort):
        """(short_by, units) or None where the series has no such cell."""
        if series == "original run":
            e = raw["elena"][cohort]["burden"]
            return e["short_by"], e["units"]
        if series == "rebuild":
            e = raw["cells"][f"burden_broad_{cohort}"]
            return e["short_by"], e["units"]
        if cohort not in bv:
            return None
        return bv[cohort]["short_by"], bv[cohort]["units"]

    plt.rcParams.update({
        "figure.dpi": 110, "savefig.dpi": 200, "savefig.bbox": "tight",
        "font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
        "axes.spines.left": False, "axes.edgecolor": INK_3, "axes.linewidth": 0.8,
        "xtick.color": INK_2, "ytick.color": INK, "text.color": INK,
        "legend.frameon": False, "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
    })
    fig, ax = plt.subplots(figsize=(9.2, 3.9))

    for ci, cohort in enumerate(COHORTS[::-1]):
        vals = [(s, cell(s, cohort)) for s in SERIES]
        present = [(s, v[0]) for s, v in vals if v is not None]
        missing = [s for s, v in vals if v is None]
        ax.plot([min(v for _, v in present), max(v for _, v in present)],
                [ci, ci], color=INK_3, lw=2, alpha=0.35,
                solid_capstyle="round", zorder=1)
        for j, (s, v) in enumerate(sorted(present, key=lambda t: t[1])):
            ax.plot(v, ci, "o", ms=10, color=PCOL[s], markeredgecolor=SURFACE,
                    markeredgewidth=2, zorder=3)
            dy = DY[j]
            ax.annotate(f"{v:.1f}×" if v < 10 else f"{v:.0f}×", (v, ci),
                        textcoords="offset points", xytext=(0, dy), ha="center",
                        fontsize=8, color=INK_2,
                        arrowprops=dict(arrowstyle="-", lw=0.4, color=INK_3,
                                        shrinkA=0, shrinkB=4)
                        if abs(dy) > 15 else None)
        # A missing series is stated, not left blank. "not run" reads as a decision;
        # an empty row reads as a result that failed to appear.
        if missing:
            ax.annotate(f"{' and '.join(missing)}: not run",
                        xy=(max(v for _, v in present) * 1.6, ci),
                        va="center", fontsize=7.5, color=INK_3, style="italic")

    ax.axvline(1.0, color=INK, lw=1.2, zorder=2)
    ax.annotate("threshold", xy=(1.0, -0.78), ha="center", va="center",
                fontsize=7, color=INK)
    ax.set_yticks(range(3))
    ax.set_yticklabels(COHORTS[::-1])
    ax.set_xscale("log")
    ax.set_xlim(0.8, 260)
    ax.set_ylim(-0.95, 2.72)
    ax.grid(axis="x", alpha=0.18, lw=0.5)
    ax.set_xlabel("how far the best gene fell short of that arm's own Bonferroni bar "
                  "(fold, log scale)")
    ax.legend(handles=[Line2D([], [], marker="o", ls="", ms=8, color=PCOL[s], label=s)
                       for s in SERIES],
              loc="lower right", bbox_to_anchor=(1.0, 1.02), ncol=3,
              fontsize=8.5, handletextpad=0.4, columnspacing=1.6)

    fig.suptitle("Three annotations of the same analysis. None reaches its threshold.",
                 fontsize=11.5, x=0.055, ha="left", y=1.16)
    fig.text(0.055, 1.085,
             "broad hearing loss · rare-variant gene burden · same cohort, "
             "covariates, GRM and null models — only the mask differs",
             fontsize=9, color=INK_2, ha="left")
    units = {s: [cell(s, c)[1] for c in COHORTS if cell(s, c)] for s in SERIES}
    fig.text(0.055, 1.030,
             f"each arm is judged against 0.05 divided by the genes IT tested "
             f"({units['original run'][0]:,} · {units['rebuild'][0]:,} · "
             f"{units['BRaVa annotation'][0]:,} in combined), so the ratio travels "
             f"where a raw p-value would not",
             fontsize=7.5, color=INK_3, ha="left")

    out = HERE / "figures/annotation_comparison.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out)
    plt.close(fig)
    print(f"  {out.relative_to(HERE)}")


if __name__ == "__main__":
    main()
