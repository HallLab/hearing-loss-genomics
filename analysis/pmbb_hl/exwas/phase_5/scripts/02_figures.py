#!/usr/bin/env python3
"""
ExWAS phase 5 -- figures.

  manhattan.png   the exome, both bars, 2 phenotypes x 3 cohorts
  lambda_mhc.png  the MHC question the first pass could not answer

COLOUR. Two phenotypes, so two of the documented categorical slots; cohorts are facets.
Values unchanged; the validator could not be run, node is not installed here.
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[1]
PHASE4 = HERE.parent / "phase_4/results"
FIG = HERE / "figures"
COHORTS = ["combined", "EUR", "AFR"]
PHENOS = ["broad", "bilateral"]
BLUE, ORANGE, PALE = "#2a78d6", "#eb6834", "#9ab4d1"
INK, INK_2, INK_3, SURFACE, RED = "#1a1a19", "#5c5b55", "#8a8980", "#ffffff", "#c0392b"
PCOL = {"broad": BLUE, "bilateral": ORANGE}
PRIMARY = 5e-8


def style():
    plt.rcParams.update({
        "figure.dpi": 110, "savefig.dpi": 190, "savefig.bbox": "tight",
        "font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.alpha": 0.22, "grid.linewidth": 0.5,
        "legend.frameon": False, "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE, "text.color": INK, "axes.edgecolor": INK_3,
        "xtick.color": INK_2, "ytick.color": INK_2,
    })


def load(tag):
    d = pd.read_csv(PHASE4 / f"merged_{tag}.tsv", sep="\t",
                    usecols=["CHR", "POS", "MarkerID", "p.value"])
    return d[(d["p.value"] > 0) & (d["p.value"] <= 1)]


def fig_manhattan():
    rep = json.loads((HERE / "results/01_calibration.json").read_text())["arms"]
    fig, axes = plt.subplots(3, 2, figsize=(13, 9), sharex=True, sharey=True)
    for r, cohort in enumerate(COHORTS):
        for c, pheno in enumerate(PHENOS):
            ax = axes[r, c]
            tag = f"{pheno}_{cohort}"
            d = load(tag).copy()
            d["CHR"] = pd.to_numeric(d.CHR, errors="coerce")
            d = d.dropna(subset=["CHR"]).sort_values(["CHR", "POS"])
            off, ticks, labels, x0 = {}, [], [], 0
            for ch, grp in d.groupby("CHR"):
                off[ch] = x0
                ticks.append(x0 + grp.POS.max() / 2)
                labels.append(int(ch))
                x0 += grp.POS.max() + 2e7
            x = d.POS + d.CHR.map(off)
            shade = np.where(d.CHR.astype(int) % 2 == 1, PCOL[pheno], PALE)
            ax.scatter(x, -np.log10(d["p.value"]), s=2, linewidths=0, c=shade,
                       rasterized=True)
            e = rep[tag]
            ax.axhline(-np.log10(PRIMARY), color=RED, lw=1.0, label="5e-8")
            ax.axhline(-np.log10(e["secondary_bar"]), color=RED, lw=1.0, ls="--",
                       label=f"0.05/{e['tested']:,}")
            if r == 0:
                ax.set_title(pheno, fontweight="bold", fontsize=10)
            if c == 0:
                ax.set_ylabel(f"{cohort}\n$-\\log_{{10}}p$")
            ax.set_ylim(0, 8.4)
            ax.legend(loc="upper right", fontsize=7, ncol=2)
    for ax in axes[-1]:
        ax.set_xticks(ticks)
        ax.set_xticklabels(labels, fontsize=7)
        ax.set_xlabel("chromosome")
    fig.suptitle("Common-variant ExWAS — no variant clears either bar in any of the "
                 "six cells", fontsize=11.5, x=0.06, ha="left", y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(FIG / "manhattan.png")
    plt.close(fig)
    return "manhattan.png"


def fig_lambda_mhc():
    rep = json.loads((HERE / "results/01_calibration.json").read_text())["arms"]
    # the first pass's EUR broad arm, 4 PCs of its own PCA
    FIRST_PASS = {"outside": 0.9949, "inside": 1.3386}

    labels, out_v, in_v, cols = [], [], [], []
    for pheno in PHENOS:
        for cohort in COHORTS:
            e = rep[f"{pheno}_{cohort}"]
            labels.append(f"{pheno}\n{cohort}")
            out_v.append(e["lambda_outside_mhc"])
            in_v.append(e["lambda_inside_mhc"])
            cols.append(PCOL[pheno])

    fig, ax = plt.subplots(figsize=(10.5, 4.4))
    x = np.arange(len(labels))
    w = 0.38
    ax.bar(x - w / 2, out_v, w, color=PALE, label="outside the MHC")
    ax.bar(x + w / 2, in_v, w, color=cols, label="inside the MHC")
    ax.axhline(1.0, color=INK, lw=1.0, ls="--")
    ax.axhline(FIRST_PASS["inside"], color=RED, lw=1.1, ls=":",
               label=f"first pass, EUR broad in MHC ({FIRST_PASS['inside']:.2f})")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("genomic inflation factor λ")
    ax.set_ylim(0.8, 1.48)
    for xi, v in zip(x - w / 2, out_v):
        ax.text(xi, v + 0.008, f"{v:.3f}", ha="center", fontsize=6.5)
    for xi, v in zip(x + w / 2, in_v):
        ax.text(xi, v + 0.008, f"{v:.3f}", ha="center", fontsize=7,
                fontweight="bold" if v > 1.10 else "normal")
    ax.legend(loc="upper left", fontsize=8)
    fig.suptitle("The MHC question, answered — and a new one opened",
                 fontsize=11.5, x=0.055, ha="left", y=1.03)
    fig.text(0.055, 0.975,
             "the broad arms' inflation is gone with the release's PCs. The bilateral "
             "arms show it instead, with the same PCs, the same GRM and the same "
             "region — only the phenotype differs.",
             fontsize=8, color=INK_2, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(FIG / "lambda_mhc.png")
    plt.close(fig)
    return "lambda_mhc.png"


def main():
    style()
    FIG.mkdir(parents=True, exist_ok=True)
    for f in (fig_manhattan, fig_lambda_mhc):
        print(f"  {f()}")


if __name__ == "__main__":
    main()
