#!/usr/bin/env python3
"""
Phase 3 -- figures for the PCA and the GRM, to explain them rather than assert them.

FORM, picked per question:

  variance.png     "how many PCs is enough?" -- a cumulative curve per cohort, with the
                   80% line and the three runs' choices marked on it. Cumulative
                   magnitude over an ordered sequence is a line; the three runs are
                   points on it, not their own lines, because they differ only in where
                   they stop.

  what_pcs_see.png "what is a PC actually measuring?" -- the thing a scree plot cannot
                   show. Two panels: the combined PCA separating ancestries, and the
                   within-EUR PCA where that separation is gone and what remains is
                   finer structure. Identity in two dimensions is a scatter.

COLOR. Six ancestry classes exist, and a scatter uses the all-pairs pairlist where the
documented palette validates only its first three slots. So the four small classes fold
to one labelled group, as the palette reference requires, rather than taking a fourth
and fifth hue that would fail the floors. EUR and AFR are the analysis cohorts and keep
their own slots.

The validator could not be run -- node is not installed on this cluster -- so the
documented slots are used unchanged and no re-stepping was attempted. Slot 3 sits below
3:1 on the light surface, which triggers the relief rule: the legend is present and the
marks are large enough to read without relying on hue alone.
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[1]
REL = Path("/static/PMBB/PMBB-Release-2026-4.0")
PCA = REL / "Exome/PCA"
FIG = HERE / "figures"

BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK_2, INK_3, SURFACE = "#1a1a19", "#5c5b55", "#8a8980", "#ffffff"

# PCs each run used. The first pass is plotted only for `combined`: for EUR and AFR it
# used a PCA OF ITS OWN, so its count cannot be placed on the release's curve. Its 4 PCs
# reached 66.7% on its own EUR eigenvalues, not the 39% this curve would imply.
RUNS = {"original run": (5, 9, 10), "first pass": (5, None, None),
        "this rebuild": (11, 15, 16)}
RUN_COLOR = {"original run": BLUE, "first pass": ORANGE, "this rebuild": AQUA}
COHORTS = ["combined", "EUR", "AFR"]


def style():
    plt.rcParams.update({
        "figure.dpi": 110, "savefig.dpi": 190, "savefig.bbox": "tight",
        "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.alpha": 0.22, "grid.linewidth": 0.5,
        "legend.frameon": False, "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE, "text.color": INK,
        "axes.edgecolor": INK_3, "xtick.color": INK_2, "ytick.color": INK_2,
    })


def eigenvalues(cohort):
    f = (PCA / "combined/PMBB-Release-2026-4.0_genetic_exome.no1KG.eigenvalues.tsv"
         if cohort == "combined" else
         PCA / f"by_ancestry/PMBB-Release-2026-4.0_genetic_exome.no1KG"
               f".eigenvalues.{cohort}.tsv")
    with f.open() as fh:
        fh.readline()
        return np.array([float(x) for x in fh.readline().split()])


def fig_variance():
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.8), sharey=True)
    for ax, c in zip(axes, COHORTS):
        cum = np.cumsum(100 * eigenvalues(c) / eigenvalues(c).sum())
        x = np.arange(1, len(cum) + 1)
        ax.plot(x, cum, color=INK_3, lw=2, zorder=2)
        ax.axhline(80, color=INK, lw=1.0, ls="--", zorder=1)
        ax.text(1.2, 83, "80% — the rule", fontsize=7.5, color=INK)
        for run, counts in RUNS.items():
            n = counts[COHORTS.index(c)]
            if n is None:
                continue
            ax.plot(n, cum[n - 1], "o", ms=9, color=RUN_COLOR[run],
                    markeredgecolor=SURFACE, markeredgewidth=2, zorder=3,
                    label=run if c == "combined" else None)
            ax.annotate(f"{n} PCs\n{cum[n-1]:.0f}%", (n, cum[n - 1]),
                        textcoords="offset points", xytext=(0, -26),
                        ha="center", fontsize=7, color=INK_2)
        ax.set_xlim(0, 21)
        ax.set_ylim(0, 102)
        ax.set_xticks([1, 5, 10, 15, 20])
        ax.set_xlabel("principal components used")
        ax.set_title(c, fontsize=9.5)
        if c == "combined":
            ax.set_ylabel("variance explained, cumulative")
            ax.legend(loc="lower right", fontsize=8, markerscale=0.9)
        else:
            ax.text(0.97, 0.06, "the first pass is absent here:\nit used a PCA of its "
                    "own,\nso its count is not on this curve",
                    transform=ax.transAxes, ha="right", va="bottom",
                    fontsize=7, color=INK_3)
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    fig.suptitle("How many principal components is enough — and where each run stopped",
                 fontsize=11, x=0.055, ha="left", y=0.995)
    fig.text(0.055, 0.935,
             "on the release's exome eigenvalues. “Variance” is of what the top twenty "
             "components capture, since plink publishes twenty — not of total genotypic "
             "variance",
             fontsize=8, color=INK_2, ha="left")
    fig.savefig(FIG / "variance.png")
    plt.close(fig)
    return "variance.png"


def fig_what_pcs_see():
    # NOT the eigenvector file's own `Class` column: in the combined file it reads
    # UNK for all 70,925 people. It is a placeholder, and reading it silently folds
    # every ancestry into one group while the figure still renders.
    comb = pd.read_csv(
        PCA / "combined/PMBB-Release-2026-4.0_genetic_exome.no1KG.eigenvectors.tsv",
        sep="\t", usecols=["person_id", "exome_PC1", "exome_PC2"])
    anc = pd.read_csv(
        PCA / "combined/PMBB-Release-2026-4.0_genetic_exome.commonsnps"
              ".samples_ancestries.tsv",
        sep="\t", usecols=["IID", "Class"]).rename(columns={"IID": "person_id"})
    comb = comb.merge(anc, on="person_id", how="left")
    if comb.Class.nunique() < 2:
        raise SystemExit("REFUSING: ancestry column has one value; the join failed "
                         "and the figure would claim every person is one group.")
    eur = pd.read_csv(
        PCA / "by_ancestry/PMBB-Release-2026-4.0_genetic_exome.no1KG"
              ".eigenvectors.EUR.tsv",
        sep="\t", usecols=["person_id", "exome_EUR_PC1", "exome_EUR_PC2"])

    # four small classes fold to one, per the palette's three-slot all-pairs limit
    small = ["EAS", "SAS", "AMR", "UNK"]
    comb["grp"] = np.where(comb.Class.isin(small), "other ancestries", comb.Class)
    order = [("EUR", BLUE), ("AFR", ORANGE), ("other ancestries", AQUA)]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.4))
    for g, col in order:
        s = comb[comb.grp == g]
        ax1.scatter(s.exome_PC1, s.exome_PC2, s=2.5, linewidths=0, color=col,
                    alpha=0.5, rasterized=True, label=f"{g}  ({len(s):,})")
    ax1.set_xlabel("PC1")
    ax1.set_ylabel("PC2")
    ax1.set_title("Combined PCA — PC1 and PC2 separate ancestries", fontsize=9.5)
    ax1.legend(loc="best", fontsize=7.5, markerscale=4)

    ax2.scatter(eur.exome_EUR_PC1, eur.exome_EUR_PC2, s=2.5, linewidths=0,
                color=BLUE, alpha=0.45, rasterized=True)
    ax2.set_xlabel("PC1")
    ax2.set_ylabel("PC2")
    ax2.set_title(f"Within EUR — that separation is gone; what is left is finer "
                  f"structure\n({len(eur):,} people, one colour because they are "
                  f"one group)", fontsize=9.5)

    fig.tight_layout(rect=(0, 0, 1, 0.90))
    fig.suptitle("What a principal component is measuring", fontsize=11,
                 x=0.055, ha="left", y=0.995)
    fig.text(0.055, 0.945,
             "this is why the analysis uses a within-ancestry PCA for EUR and AFR, and "
             "the combined one only for the combined cohort",
             fontsize=8, color=INK_2, ha="left")
    fig.savefig(FIG / "what_pcs_see.png")
    plt.close(fig)
    return "what_pcs_see.png"


def main():
    style()
    FIG.mkdir(parents=True, exist_ok=True)
    for f in (fig_variance, fig_what_pcs_see):
        print(f"  {f()}")


if __name__ == "__main__":
    main()
