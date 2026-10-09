#!/usr/bin/env python3
"""
Phase 5, step 4 -- figures for the bilateral arm.

Same palette and rcParams as the other arms so all three sets can sit in one document.

The comparison figure here is NOT a concordance plot. The two arms share every control
and the bilateral cases are a strict subset of the broad ones, so agreement is
guaranteed by construction and measures nothing. What the plot is for is the other
direction: points above the diagonal are variants that got stronger on half the cases,
which shared data cannot produce.
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import chi2

BASE = Path(__file__).resolve().parents[2]
BROAD = BASE.parent / "exwas_elena_replication"
FIG = BASE / "figures"
RES = BASE / "phase_5/results"
COHORTS = ["combined", "EUR", "AFR"]
BLUE, PALE, RED, GREY = "#3b6ea5", "#9ab4d1", "#c0392b", "#888888"
CHI2_MEDIAN = chi2.ppf(0.5, 1)
PRIMARY = 5e-8


def style():
    plt.rcParams.update({
        "figure.dpi": 110, "savefig.dpi": 160, "savefig.bbox": "tight",
        "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.alpha": 0.25, "grid.linewidth": 0.5,
        "legend.frameon": False,
    })


def load(root, c):
    d = pd.read_csv(root / f"phase_4/results/merged_{c}.tsv", sep="\t")
    return d[(d["p.value"] > 0) & (d["p.value"] <= 1)].copy()


def lam(p):
    return float(np.median(chi2.isf(np.asarray(p), 1)) / CHI2_MEDIAN)


def fig_qq():
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.6))
    for ax, c in zip(axes, COHORTS):
        for root, col, label in ((BROAD, PALE, "broad (6,752 cases)"),
                                 (BASE, BLUE, "bilateral (3,164 cases)")):
            p = np.sort(load(root, c)["p.value"].to_numpy())
            n = len(p)
            exp = -np.log10((np.arange(1, n + 1) - 0.5) / n)
            ax.scatter(exp, -np.log10(p), s=3, linewidths=0, color=col, alpha=0.6,
                       label=f"{label}  λ={lam(p):.3f}")
        hi = max(ax.get_xlim()[1], ax.get_ylim()[1])
        ax.plot([0, hi], [0, hi], color=GREY, lw=0.8, ls="--")
        bar = -np.log10(PRIMARY)
        ax.axhline(bar, color=RED, lw=0.9, ls=":")
        ax.set_ylim(top=bar + 1.1)
        ax.text(ax.get_xlim()[1], bar + 0.12, "5e-8", ha="right", va="bottom",
                fontsize=7, color=RED)
        ax.set_title(c)
        ax.set_xlabel("expected  $-\\log_{10}p$")
        if c == "combined":
            ax.set_ylabel("observed  $-\\log_{10}p$")
        ax.legend(loc="upper left", fontsize=7, markerscale=3,
                  bbox_to_anchor=(0.0, 0.88))
    fig.suptitle("Bilateral sensorineural against the broad phenotype — neither "
                 "approaches 5e-8", fontsize=10)
    fig.tight_layout()
    fig.savefig(FIG / "qq.png")
    plt.close(fig)
    return "qq.png"


def fig_manhattan():
    rep = json.loads((RES / "calibration_crossarm.json").read_text())["cohorts"]
    fig, axes = plt.subplots(3, 1, figsize=(11, 8), sharex=True)
    for ax, c in zip(axes, COHORTS):
        d = load(BASE, c)
        d["CHR"] = pd.to_numeric(d["CHR"], errors="coerce")
        d = d.dropna(subset=["CHR"]).sort_values(["CHR", "POS"])
        offset, ticks, labels, x0 = {}, [], [], 0
        for ch, grp in d.groupby("CHR"):
            offset[ch] = x0
            ticks.append(x0 + grp["POS"].max() / 2)
            labels.append(int(ch))
            x0 += grp["POS"].max() + 2e7
        x = d["POS"] + d["CHR"].map(offset)
        colors = np.where(d["CHR"].astype(int) % 2 == 1, BLUE, PALE)
        ax.scatter(x, -np.log10(d["p.value"]), s=2, linewidths=0, c=colors)
        sec = 0.05 / rep[c]["n_tested"]
        ax.axhline(-np.log10(PRIMARY), color=RED, lw=0.9, label="5e-8")
        ax.axhline(-np.log10(sec), color=RED, lw=0.9, ls="--",
                   label=f"0.05/{rep[c]['n_tested']:,} = {sec:.2e}")
        ax.set_ylabel(f"{c}\n$-\\log_{{10}}p$")
        ax.set_ylim(0, max(8.5, -np.log10(PRIMARY) + 1))
        ax.legend(loc="upper right", fontsize=7, ncol=2)
    axes[-1].set_xticks(ticks)
    axes[-1].set_xticklabels(labels, fontsize=7)
    axes[-1].set_xlabel("chromosome")
    fig.suptitle("Bilateral sensorineural ExWAS — no variant clears either bar",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(FIG / "manhattan.png")
    plt.close(fig)
    return "manhattan.png"


def fig_crossarm():
    rep = json.loads((RES / "calibration_crossarm.json").read_text())["cohorts"]
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 4.0))
    for ax, c in zip(axes, COHORTS):
        m = load(BASE, c).merge(load(BROAD, c)[["MarkerID", "p.value"]],
                                on="MarkerID", suffixes=("_bil", "_broad"))
        x = -np.log10(m["p.value_broad"])
        y = -np.log10(m["p.value_bil"])
        hi = float(max(x.max(), y.max())) * 1.05
        ax.plot([0, hi], [0, hi], color=GREY, lw=0.8, ls="--", zorder=1)
        gained = (y > x + 1.0) & (m["p.value_bil"] < 1e-4)
        ax.scatter(x[~gained], y[~gained], s=4, linewidths=0, color=PALE,
                   alpha=0.4, zorder=2, rasterized=True)
        ax.scatter(x[gained], y[gained], s=22, linewidths=0.7, facecolors="none",
                   edgecolors=RED, zorder=3,
                   label=f"stronger on half the cases ({int(gained.sum())})")
        ax.set_title(f"{c}\nρ = {rep[c]['spearman_vs_broad']:.3f}  "
                     f"(nested data — not a replication statistic)", fontsize=8.5)
        ax.set_xlabel("broad phenotype  $-\\log_{10}p$")
        if c == "combined":
            ax.set_ylabel("bilateral  $-\\log_{10}p$")
        ax.legend(loc="upper left", fontsize=7)
    fig.suptitle("Restricting the phenotype: the only informative direction is "
                 "upward, and the power control decides whether it means anything",
                 fontsize=9.5)
    fig.tight_layout()
    fig.savefig(FIG / "crossarm.png")
    plt.close(fig)
    return "crossarm.png"


def fig_power():
    pw = json.loads((RES / "power.json").read_text())
    bro = json.loads((BROAD / "phase_5/results/power.json").read_text())
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.9), sharey=True)
    for ax, c in zip(axes, COHORTS):
        d = load(BASE, c)
        maf = np.minimum(d["AF_Allele2"], 1 - d["AF_Allele2"])
        ax.scatter(maf, np.exp(np.abs(d["BETA"])), s=2, linewidths=0, color=PALE,
                   alpha=0.35, label="observed |effect|", rasterized=True)
        k, kb = pw["cohorts"][c], bro["cohorts"][c]
        ax.plot(kb["maf"], kb["detectable_or"], color=GREY, lw=1.2, ls="--",
                label=f"broad, {kb['cases']:,} cases")
        ax.plot(k["maf"], k["detectable_or"], color=RED, lw=1.7,
                label=f"bilateral, {k['cases']:,} cases")
        ax.set_xscale("log")
        ax.set_xlim(0.009, 0.55)
        ax.set_ylim(1.0, 5.0)
        ax.set_xlabel("minor allele frequency")
        ax.set_title(c, fontsize=9)
        if c == "combined":
            ax.set_ylabel("odds ratio (away from 1)")
            ax.legend(loc="upper right", fontsize=7, markerscale=4)
    fig.suptitle("What halving the cases costs: the detectable effect at 80% power, "
                 "5e-8", fontsize=10)
    fig.tight_layout()
    fig.savefig(FIG / "power.png")
    plt.close(fig)
    return "power.png"


def main():
    style()
    FIG.mkdir(parents=True, exist_ok=True)
    for f in (fig_qq, fig_manhattan, fig_crossarm, fig_power):
        print(f"  {f()}")


if __name__ == "__main__":
    main()
