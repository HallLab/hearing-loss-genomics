#!/usr/bin/env python3
"""
Phase 5 -- the data loading and drawing the notebooks use.

Kept out of the notebooks so that both the English and Portuguese versions call
the same code, and so that a figure can be regenerated without opening either.

The genomic inflation factor uses the exact qchisq(0.5, 1) = 0.4549364, not the
rounded 0.456 in Elena's step10_1_QQplot.R. Hers inflates lambda by 0.23%,
which changes no conclusion but is free to get right.

Lambda on a SKAT-O p-value needs one word of defence, since lambda is built for
single-variant tests where the null statistic is chi-square on 1 df and SKAT-O's
is not. It still works as written, because the quantity being tested is whether
the p-values are uniform: if they are, median(qchisq(1-p, 1)) is qchisq(0.5, 1)
whatever produced them. It is a uniformity check, not a claim about the
statistic's distribution.
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import chi2
from statsmodels.stats.multitest import multipletests

HERE = Path(__file__).resolve().parents[1]
REPL = HERE.parent

COHORTS = ["combined", "EUR", "AFR"]
MASKS = ["pLOF", "pDM", "pLOF_pDM"]
MAFS = [0.0001, 0.001, 0.01]
CHI2_MEDIAN = chi2.ppf(0.5, 1)          # 0.4549364, not 0.456

PCS = {"combined": 5, "EUR": 4, "AFR": 3}
N = {"combined": 57498, "EUR": 42779, "AFR": 11334}
CASES = {"combined": 6752, "EUR": 5183, "AFR": 1285}


def style():
    plt.rcParams.update({
        "figure.dpi": 110, "savefig.dpi": 160, "savefig.bbox": "tight",
        "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.alpha": 0.25, "grid.linewidth": 0.5,
        "legend.frameon": False,
    })


def burden():
    """Our results, with a genome coordinate per gene."""
    d = pd.read_csv(HERE / "results/burden_all_cohorts.tsv", sep="\t")
    d["Pvalue"] = pd.to_numeric(d["Pvalue"])
    pos = pd.read_csv(HERE / "results/gene_positions.tsv", sep="\t")
    return d.merge(pos[["Region", "POS", "coord_source"]], on="Region", how="left")


def hers(three_masks_only=True):
    """Her merged results, for the side-by-side.

    Her table carries two defects Phase 4 documented and this has to survive:
    5,391 rows whose columns are shifted (three files lost their header) and a
    '-' pseudo-gene. Both are dropped here rather than silently averaged in.
    """
    f = ("/project/hall/analysis/hearing-loss-genomics/analysis/elena/"
         "HL_only_rarevariant/newmasks_combined_gene_results.csv")
    d = pd.read_csv(f, low_memory=False, keep_default_na=False, dtype=str)
    if three_masks_only:
        d = d[d["Mask"].isin(MASKS)]
    d["Pvalue"] = pd.to_numeric(d["Pvalue"], errors="coerce")
    before = len(d)
    d = d.dropna(subset=["Pvalue"])
    d = d[~d["Region"].isin(["", "-"])]
    d.attrs["rows_dropped"] = before - len(d)
    return d.rename(columns={"Ancestry": "Cohort"})


def lam(p):
    p = np.asarray(p, dtype=float)
    p = p[(p > 0) & (p <= 1)]
    return float(np.median(chi2.ppf(1 - p, 1)) / CHI2_MEDIAN)


def fdr(p):
    """Benjamini-Hochberg q-values."""
    p = np.asarray(p, dtype=float)
    return multipletests(p, alpha=0.05, method="fdr_bh")[1]


def genome_axis(d):
    """Cumulative x coordinate and the per-chromosome tick positions."""
    sizes = d.groupby("CHR").POS.max()
    offset = sizes.reindex(range(1, 23)).fillna(0).cumsum().shift(1).fillna(0)
    x = d.POS + d.CHR.map(offset)
    ticks = [offset[c] + sizes.get(c, 0) / 2 for c in range(1, 23)]
    return x, ticks, offset


def manhattan(d, cohort, ax=None, mask=None, maf=None, label_top=5):
    s = d[d.Cohort == cohort]
    if mask:
        s = s[s.Mask == mask]
    if maf:
        s = s[s.max_MAF == maf]
    s = s.dropna(subset=["POS"]).copy()
    ax = ax or plt.gca()

    x, ticks, _ = genome_axis(s)
    s["x"] = x
    y = -np.log10(s.Pvalue)

    for c in range(1, 23):
        m = s.CHR == c
        ax.scatter(s.x[m], y[m], s=3, linewidths=0,
                   color="#3b6ea5" if c % 2 else "#9ab4d1")

    # Bonferroni over the genes IN THIS PANEL -- that is, this panel treated as
    # if it were the only analysis. It is the most permissive bar in play: it
    # does not charge for having looked at nine panels per cohort, which the
    # 0.05/tests bar in the findings does. Drawn because nothing crosses even
    # this one, so the conclusion does not turn on the choice.
    bar_genes = 0.05 / s.Region.nunique()
    ax.axhline(-np.log10(bar_genes), color="#c0392b", lw=0.9,
               label=f"Bonferroni, this panel alone: 0.05/{s.Region.nunique():,} genes "
                     f"= {bar_genes:.1e}")

    # No FDR line: Benjamini-Hochberg's cutoff is the largest p with q < 0.05,
    # and there is none, so a line would have nowhere to sit. The smallest q is
    # annotated instead, which carries the same information honestly.
    q = fdr(s.Pvalue.values)
    ax.annotate(f"min $q$ = {q.min():.2f}", (0.985, 0.93), xycoords="axes fraction",
                ha="right", fontsize=7, color="#c0392b")
    ax.set_xticks(ticks)
    ax.set_xticklabels([str(c) if c <= 12 or c % 2 else "" for c in range(1, 23)],
                       fontsize=6)
    ax.set_xlim(s.x.min() - 1e6, s.x.max() + 1e6)
    ax.set_ylim(0, max(6.2, y.max() * 1.18))
    ax.set_ylabel("$-\\log_{10}p$")
    ax.grid(axis="x", visible=False)

    # Labels near either edge get pushed inward; centring them there runs the
    # text under the y-axis, which a first render showed happening to THAP3.
    span = s.x.max() - s.x.min()
    for _, r in s.nsmallest(label_top, "Pvalue").iterrows():
        frac = (r.x - s.x.min()) / span
        ha, dx = "center", 0
        if frac < 0.06:
            ha, dx = "left", 2
        elif frac > 0.94:
            ha, dx = "right", -2
        ax.annotate(r.Region, (r.x, -np.log10(r.Pvalue)),
                    xytext=(dx, 4), textcoords="offset points",
                    ha=ha, fontsize=6, color="#1a1a1a")
    return ax


def qq(p, ax=None, title=None, color="#3b6ea5"):
    ax = ax or plt.gca()
    p = np.sort(np.asarray(p, dtype=float))
    p = p[(p > 0) & (p <= 1)]
    n = len(p)
    exp = -np.log10((np.arange(1, n + 1) - 0.5) / n)
    obs = -np.log10(p)

    # 95% band for the order statistics, from the beta distribution
    from scipy.stats import beta
    k = np.arange(1, n + 1)
    lo = -np.log10(beta.ppf(0.975, k, n - k + 1))
    hi = -np.log10(beta.ppf(0.025, k, n - k + 1))
    ax.fill_between(exp, lo, hi, color="#cccccc", alpha=0.45, linewidth=0)

    hi_lim = max(exp.max(), obs.max()) * 1.05
    ax.plot([0, hi_lim], [0, hi_lim], color="#888888", lw=0.8, ls="--")
    ax.scatter(exp, obs, s=3, linewidths=0, color=color)
    ax.set_xlabel("expected $-\\log_{10}p$")
    ax.set_ylabel("observed $-\\log_{10}p$")
    ax.set_xlim(0, hi_lim)
    ax.set_ylim(0, hi_lim)
    ax.annotate(f"$\\lambda$ = {lam(p):.3f}", (0.04, 0.92), xycoords="axes fraction",
                fontsize=9)
    if title:
        ax.set_title(title)
    return ax


def comparison(ours, theirs, cohort, ax=None, label_n=8):
    """Our -log10 p against hers, one point per gene, best p across the grid.

    Genes she tested and we did not cannot go in the scatter -- they have no y.
    They are the non-coding genes Phase 2 removed, and dropping them silently
    would hide the single most consequential difference between the two arms,
    so they are drawn as a strip along the bottom instead.
    """
    o = ours[ours.Cohort == cohort].groupby("Region").Pvalue.min()
    h = theirs[theirs.Cohort == cohort].groupby("Region").Pvalue.min()
    ax = ax or plt.gca()

    shared = pd.DataFrame({"hers": h, "ours": o}).dropna()
    only_hers = h[~h.index.isin(o.index)]

    hx, oy = -np.log10(shared.hers), -np.log10(shared.ours)
    hi = max(hx.max(), oy.max()) * 1.08

    ax.plot([0, hi], [0, hi], color="#999999", lw=0.8, ls="--", zorder=1)
    ax.scatter(hx, oy, s=5, linewidths=0, color="#3b6ea5", alpha=0.5, zorder=2)

    # the dropped genes, on a strip below the axis
    strip = -0.55
    dx = -np.log10(only_hers)
    ax.scatter(dx, np.full(len(dx), strip), s=7, marker="|",
               color="#c0392b", linewidths=0.8, zorder=3)
    ax.axhline(0, color="#dddddd", lw=0.6, zorder=1)

    # label the biggest movers among shared genes, and her top dropped gene
    moved = (oy - hx).abs().nlargest(label_n)
    for g in moved.index:
        ax.annotate(g, (hx[g], oy[g]), xytext=(3, 3), textcoords="offset points",
                    fontsize=6, color="#1a1a1a")
    if len(dx):
        top = dx.idxmax()
        ax.annotate(top, (dx[top], strip), xytext=(0, -11),
                    textcoords="offset points", ha="center", fontsize=7,
                    color="#c0392b", fontweight="bold")

    ax.set_xlim(0, hi)
    ax.set_ylim(strip - 0.9, hi)
    ax.set_xlabel("hers,  $-\\log_{10}p$")
    ax.set_ylabel("ours,  $-\\log_{10}p$")
    # The count goes in the title rather than a legend: a first render put the
    # legend box on top of a labelled gene in the EUR panel.
    ax.set_title(f"{cohort} — {len(shared):,} genes in both, "
                 f"{len(only_hers):,} only in hers")
    return ax
