#!/usr/bin/env python3
"""
Phase 5, step 6 -- the figures that explain the null rather than state it.

  power.png           observed effect sizes against what was detectable
  lambda_by_region.png  the MHC limitation, ours against hers
  qc_selfcheck.png    call rate against significance in our own tested set
  top_variants.png    the strongest signals, annotated, with their loci marked
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
ELENA = BASE.parent / "elena/HL_only_ExWAS"
COHORTS = ["combined", "EUR", "AFR"]
BLUE, PALE, RED, GREY = "#3b6ea5", "#9ab4d1", "#c0392b", "#888888"
CHI2_MEDIAN = chi2.ppf(0.5, 1)
MHC = (28_510_000, 33_400_000)


def style():
    plt.rcParams.update({
        "figure.dpi": 110, "savefig.dpi": 160, "savefig.bbox": "tight",
        "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.alpha": 0.25, "grid.linewidth": 0.5,
        "legend.frameon": False,
    })


def ours(c):
    d = pd.read_csv(BASE / f"phase_4/results/merged_{c}.tsv", sep="\t")
    return d[(d["p.value"] > 0) & (d["p.value"] <= 1)].copy()


def fig_power():
    """The figure that answers 'why nothing'."""
    pw = json.loads((RES / "power.json").read_text())
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.9), sharey=True)
    for ax, c in zip(axes, COHORTS):
        d = ours(c)
        maf = np.minimum(d["AF_Allele2"], 1 - d["AF_Allele2"])
        obs = np.exp(np.abs(d["BETA"]))
        ax.scatter(maf, obs, s=2, linewidths=0, color=PALE, alpha=0.35,
                   label="observed |effect|", rasterized=True)
        k = pw["cohorts"][c]
        ax.plot(k["maf"], k["detectable_or"], color=RED, lw=1.6,
                label="detectable at 80% power, 5e-8")
        ax.set_xscale("log")
        ax.set_xlim(0.009, 0.55)
        ax.set_ylim(1.0, 4.0)
        ax.set_xlabel("minor allele frequency")
        ax.set_title(f"{c}\n{k['cases']:,} cases of {k['n']:,}", fontsize=9)
        if c == "combined":
            ax.set_ylabel("odds ratio (away from 1)")
            ax.legend(loc="upper right", fontsize=7, markerscale=4)
    fig.suptitle("Nothing observed comes close to what this cohort could have "
                 "detected", fontsize=10)
    fig.tight_layout()
    fig.savefig(FIG / "power.png")
    plt.close(fig)
    return "power.png"


def fig_lambda_by_region():
    def lam(p):
        p = np.asarray([v for v in p if 0 < v <= 1])
        return float(np.median(chi2.isf(p, 1)) / CHI2_MEDIAN)

    rows = []
    for c in COHORTS:
        o = ours(c)[["MarkerID", "p.value"]]
        h = pd.read_csv(ELENA / f"HL_ExWAS_FullyAnnotated_AllVariants_{c}.csv",
                        usecols=["MarkerID", "p.value"], low_memory=False)
        h["p.value"] = pd.to_numeric(h["p.value"], errors="coerce")
        m = o.merge(h, on="MarkerID", suffixes=("_o", "_h")).dropna()
        pos = pd.to_numeric(m["MarkerID"].str.split(":").str[1], errors="coerce")
        in_mhc = (m["MarkerID"].str.split(":").str[0] == "6") & pos.between(*MHC)
        for region, sel in (("outside MHC", ~in_mhc), ("inside MHC", in_mhc)):
            rows.append({"cohort": c, "region": region, "n": int(sel.sum()),
                         "ours": lam(m.loc[sel, "p.value_o"]),
                         "hers": lam(m.loc[sel, "p.value_h"])})
    d = pd.DataFrame(rows)

    fig, ax = plt.subplots(figsize=(8.2, 3.9))
    x = np.arange(len(d))
    w = 0.38
    ax.bar(x - w / 2, d["hers"], w, color=PALE, label="Elena  (9 / 10 PCs)")
    ax.bar(x + w / 2, d["ours"], w, color=BLUE, label="Replication  (4 / 3 PCs)")
    ax.axhline(1.0, color=GREY, lw=1.0, ls="--")
    ax.set_xticks(x)
    ax.set_xticklabels([f"{r.cohort}\n{r.region}\nn={r.n:,}"
                        for r in d.itertuples()], fontsize=7.5)
    ax.set_ylabel("genomic inflation factor λ")
    ax.set_ylim(0.9, 1.42)
    ax.legend(loc="upper left", fontsize=8)
    for xi, v in zip(x - w / 2, d["hers"]):
        ax.text(xi, v + 0.008, f"{v:.3f}", ha="center", fontsize=6.5)
    for xi, v in zip(x + w / 2, d["ours"]):
        ax.text(xi, v + 0.008, f"{v:.3f}", ha="center", fontsize=6.5,
                fontweight="bold" if v > 1.05 else "normal")
    ax.set_title("Outside the MHC the two runs are indistinguishable. Inside it, "
                 "our fewer PCs cost us.", fontsize=9.5)
    fig.tight_layout()
    fig.savefig(FIG / "lambda_by_region.png")
    plt.close(fig)
    return "lambda_by_region.png"


def fig_qc_selfcheck():
    rep = json.loads((RES / "qc_and_top_hits.json").read_text())["self_check"]
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.7), sharey=True)
    for ax, c in zip(axes, COHORTS):
        d = ours(c)
        cr = 1 - d["MissingRate"]
        y = -np.log10(d["p.value"])
        ax.scatter(cr, y, s=2, linewidths=0, color=PALE, alpha=0.35, rasterized=True)
        # binned median, so a trend would be visible through the cloud
        bins = np.quantile(cr, np.linspace(0, 1, 21))
        bins = np.unique(bins)
        idx = np.digitize(cr, bins[1:-1])
        med = pd.DataFrame({"b": idx, "y": y}).groupby("b")["y"].median()
        ctr = [(bins[max(i, 0)] + bins[min(i + 1, len(bins) - 1)]) / 2
               for i in med.index]
        ax.plot(ctr, med.to_numpy(), color=RED, lw=1.4, marker="o", ms=2.5,
                label="median per call-rate bin")
        s = rep[c]
        ax.set_title(f"{c}\nρ = {s['spearman_callrate_vs_neglog10p']:+.4f}  "
                     f"(p = {s['p']:.2f})", fontsize=9)
        ax.set_xlabel("genotype call rate")
        ax.set_xlim(0.95, 1.001)
        if c == "combined":
            ax.set_ylabel("$-\\log_{10}p$")
            ax.legend(loc="upper left", fontsize=7)
    fig.suptitle("In our own tested set, a worse call rate buys no significance — "
                 "the pathology behind the poster's hits is absent here", fontsize=10)
    fig.tight_layout()
    fig.savefig(FIG / "qc_selfcheck.png")
    plt.close(fig)
    return "qc_selfcheck.png"


def fig_top_variants():
    rep = json.loads((RES / "qc_and_top_hits.json").read_text())["top"]
    sev = {"HIGH": "#7b241c", "MODERATE": "#c0392b", "LOW": "#3b6ea5",
           "MODIFIER": "#9ab4d1"}
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.6))
    for ax, c in zip(axes, COHORTS):
        rows = rep[c][:12][::-1]
        y = np.arange(len(rows))
        vals = [-np.log10(r["p"]) for r in rows]
        cols = [sev.get(r["impact"], GREY) for r in rows]
        ax.barh(y, vals, color=cols, height=0.6)
        ax.set_yticks(y)
        ax.set_yticklabels([f"{r['gene'] or '?'}  ·  {r['variant'].split(':')[0]}:"
                            f"{int(r['variant'].split(':')[1]):,}" for r in rows],
                           fontsize=6.8)
        ax.set_xlabel("$-\\log_{10}p$")
        ax.set_title(c, fontsize=9)
        ax.set_xlim(0, 6.2)
        for yi, r in zip(y, rows):
            ax.text(-np.log10(r["p"]) + 0.08, yi,
                    f"OR {np.exp(r['beta']):.2f}  AF {r['af']:.3f}",
                    va="center", fontsize=6, color="#1a1a1a")
    handles = [plt.Rectangle((0, 0), 1, 1, color=v) for v in sev.values()]
    fig.suptitle("The strongest signals in each cohort — all far short of either bar",
                 fontsize=10)
    fig.tight_layout(rect=(0, 0.10, 1, 1))
    fig.legend(handles, list(sev), loc="lower center", ncol=4, fontsize=7.5,
               title="VEP impact", title_fontsize=7.5, bbox_to_anchor=(0.5, 0.005))
    fig.savefig(FIG / "top_variants.png")
    plt.close(fig)
    return "top_variants.png"


def main():
    style()
    FIG.mkdir(parents=True, exist_ok=True)
    for f in (fig_power, fig_lambda_by_region, fig_qc_selfcheck, fig_top_variants):
        print(f"  {f()}")


if __name__ == "__main__":
    main()
