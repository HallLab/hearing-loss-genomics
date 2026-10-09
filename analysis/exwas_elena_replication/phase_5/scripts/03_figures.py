#!/usr/bin/env python3
"""
Phase 5, step 3 -- figures.

Same palette and rcParams as the burden arm's burden_plots.py, so the two sets of
figures can sit in one document without looking like they came from two projects.

Four figures:
  qq.png             calibration, ours against hers, per cohort
  manhattan.png      the whole exome, with both bars drawn
  concordance.png    our -log10(p) against hers, per shared variant
  poster_variants.png the three poster variants' call rate against the exome

The last one is the figure that carries the finding. It is drawn from the release
measurements in poster_audit.json, not from our results.
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


def ours(cohort):
    d = pd.read_csv(BASE / f"phase_4/results/merged_{cohort}.tsv", sep="\t",
                    usecols=["CHR", "POS", "MarkerID", "p.value", "MissingRate"])
    return d[(d["p.value"] > 0) & (d["p.value"] <= 1)]


def hers(cohort):
    d = pd.read_csv(BASE.parent / "elena/HL_only_ExWAS/"
                    f"HL_ExWAS_FullyAnnotated_AllVariants_{cohort}.csv",
                    usecols=["MarkerID", "p.value", "SYMBOL"], low_memory=False)
    d["p.value"] = pd.to_numeric(d["p.value"], errors="coerce")
    return d[(d["p.value"] > 0) & (d["p.value"] <= 1)]


def lam(p):
    return float(np.median(chi2.isf(np.asarray(p), 1)) / CHI2_MEDIAN)


def fig_qq():
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.6))
    for ax, c in zip(axes, COHORTS):
        for d, col, label in ((hers(c), PALE, "Elena"), (ours(c), BLUE, "Replication")):
            p = np.sort(d["p.value"].to_numpy())
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
    fig.suptitle("Single-variant calibration: nothing approaches 5e-8 in either run",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(FIG / "qq.png")
    plt.close(fig)
    return "qq.png"


def fig_manhattan():
    rep = json.loads((RES / "calibration_concordance.json").read_text())
    fig, axes = plt.subplots(3, 1, figsize=(11, 8), sharex=True)
    for ax, c in zip(axes, COHORTS):
        d = ours(c).copy()
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
        secondary = 0.05 / rep[c]["n_ours"]
        ax.axhline(-np.log10(PRIMARY), color=RED, lw=0.9,
                   label=f"5e-8 (poster's bar)")
        ax.axhline(-np.log10(secondary), color=RED, lw=0.9, ls="--",
                   label=f"0.05/{rep[c]['n_ours']:,} = {secondary:.2e}")
        ax.set_ylabel(f"{c}\n$-\\log_{{10}}p$")
        ax.set_ylim(0, max(8.5, -np.log10(PRIMARY) + 1))
        ax.legend(loc="upper right", fontsize=7, ncol=2)
    axes[-1].set_xticks(ticks)
    axes[-1].set_xticklabels(labels, fontsize=7)
    axes[-1].set_xlabel("chromosome")
    fig.suptitle("Common-variant ExWAS, broad hearing loss — no variant clears either bar",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(FIG / "manhattan.png")
    plt.close(fig)
    return "manhattan.png"


def fig_concordance():
    rep = json.loads((RES / "calibration_concordance.json").read_text())
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.8))
    for ax, c in zip(axes, COHORTS):
        o, h = ours(c), hers(c)
        m = o.merge(h, on="MarkerID", suffixes=("_ours", "_hers"))
        x, y = -np.log10(m["p.value_hers"]), -np.log10(m["p.value_ours"])
        mhc = ((m["MarkerID"].str.split(":").str[0] == "6")
               & m["MarkerID"].str.split(":").str[1].astype(int).between(28_510_000,
                                                                        33_400_000))
        hi = float(max(x.max(), y.max())) * 1.05
        ax.plot([0, hi], [0, hi], color=GREY, lw=0.8, ls="--", zorder=1)
        ax.scatter(x[~mhc], y[~mhc], s=4, linewidths=0, color=BLUE, alpha=0.4, zorder=2)
        ax.scatter(x[mhc], y[mhc], s=9, linewidths=0, color=RED, alpha=0.8, zorder=3,
                   label="MHC, chr6:28.5–33.4 Mb")
        ax.set_title(f"{c}\nr={rep[c]['pearson_r_neglog10p']:.3f}  "
                     f"ρ={rep[c]['spearman_rho']:.3f}  "
                     f"n={rep[c]['n_shared']:,}", fontsize=9)
        ax.set_xlabel("Elena  $-\\log_{10}p$")
        if c == "combined":
            ax.set_ylabel("Replication  $-\\log_{10}p$")
        ax.legend(loc="upper left", fontsize=7, markerscale=1.5)
    fig.suptitle("Per-variant agreement with the delivered results; every large "
                 "disagreement is in the MHC", fontsize=10)
    fig.tight_layout()
    fig.savefig(FIG / "concordance.png")
    plt.close(fig)
    return "concordance.png"


def fig_poster_variants():
    """The finding: the poster's variants against the call rate of what was tested.

    Our MissingRate is computed WITHIN a cohort, while the QC that built these files
    was applied across all 70,925 samples. So a variant at 0.995 overall can read
    0.985 inside AFR, which is why a small tail falls under 0.99 (0.04% of combined,
    0.18% of EUR, 0.91% of AFR) without the filter having failed. Nothing tested
    anywhere falls below 0.95; the floor is 0.9895 / 0.9875 / 0.9534.

    That is the contrast the figure exists to make. The poster's worst variant is at
    0.2013 -- not a tail, a different regime.
    """
    audit = json.loads((RES / "poster_audit.json").read_text())

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 3.9),
                                   gridspec_kw={"width_ratios": [1.4, 1]})

    floors = {}
    bins = np.linspace(0.95, 1.0, 100)
    for c, col in zip(COHORTS, (BLUE, "#6f9bc6", PALE)):
        cr = 1 - ours(c)["MissingRate"].to_numpy()
        floors[c] = cr.min()
        ax1.hist(cr, bins=bins, histtype="step", linewidth=1.1, color=col,
                 label=f"{c}  (min {cr.min():.4f})")
    ax1.axvline(0.99, color=RED, lw=1.0, ls="--")
    ax1.set_yscale("log")
    ax1.set_xlim(0.95, 1.001)
    ax1.text(0.9905, ax1.get_ylim()[1] * 0.3, "0.99", fontsize=7, color=RED)
    ax1.set_xlabel("genotype call rate, within cohort")
    ax1.set_ylabel("variants tested (log scale)")
    ax1.set_title("What was actually tested: nothing below 0.95", fontsize=9)
    ax1.legend(loc="upper left", fontsize=7)

    rows = sorted(((v, m["call_rate"]) for v, m in audit["release_metrics"].items()
                   if "call_rate" in m), key=lambda r: r[1])
    ypos = np.arange(len(rows))
    ax2.barh(ypos, [r[1] for r in rows], color=RED, alpha=0.75, height=0.55)
    ax2.axvline(0.99, color=RED, lw=1.0, ls="--")
    ax2.axvspan(0.95, 1.02, color=PALE, alpha=0.45, zorder=0)
    ax2.annotate("range of what\nwas tested", xy=(0.975, len(rows) - 1),
                 xytext=(0.74, len(rows) - 0.75), fontsize=6.5, color="#1a1a1a",
                 ha="center", va="center",
                 arrowprops=dict(arrowstyle="->", lw=0.7, color="#1a1a1a"))
    ax2.set_ylim(-0.6, len(rows) - 0.1)
    ax2.set_yticks(ypos)
    ax2.set_yticklabels([r[0] for r in rows], fontsize=7)
    ax2.set_xlim(0, 1.02)
    ax2.set_xlabel("genotype call rate (measured on the PMBB v4 release)")
    ax2.set_title("The poster's three variants", fontsize=9)
    for y, (v, cr) in zip(ypos, rows):
        ax2.text(min(cr + 0.012, 0.86), y, f"{cr:.3f}", va="center", fontsize=7,
                 color="#1a1a1a")

    fig.suptitle("All three variants behind the poster's six results fall below the "
                 "QC floor the poster states", fontsize=10)
    fig.tight_layout()
    fig.savefig(FIG / "poster_variants.png")
    plt.close(fig)
    return "poster_variants.png"


def main():
    style()
    FIG.mkdir(parents=True, exist_ok=True)
    for f in (fig_qq, fig_manhattan, fig_concordance, fig_poster_variants):
        name = f()
        print(f"  {name}")
    print(f"\nwritten to {FIG.relative_to(BASE.parent.parent)}")


if __name__ == "__main__":
    main()
