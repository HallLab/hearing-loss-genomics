#!/usr/bin/env python3
"""
Phase 5 -- figures.

  manhattan.png    the whole genome, both bars drawn, 2 phenotypes x 3 cohorts
  clingen.png      where the known deafness genes rank -- the one positive result
  vs_first_pass.png  did changing the PCs move anything?

COLOUR. Two phenotypes are the categorical dimension, cohorts are facets, so only two
of the documented palette's slots are needed and the all-pairs limit is not approached.
Values are the documented ones unchanged; the validator could not be run, as node is not
installed on this cluster.

FORM, per question. The Manhattan is position-vs-magnitude, so a scatter with reference
lines. The ClinGen figure asks "where in the ranking do these genes fall?", which is a
distribution over an ordered axis -- a strip plot on a log rank scale, because the whole
claim lives in the first few dozen ranks out of 17,900 and a linear axis would hide it.
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
PHENOS = ["broad", "bilateral"]
COHORTS = ["combined", "EUR", "AFR"]
BLUE, ORANGE, PALE = "#2a78d6", "#eb6834", "#9ab4d1"
INK, INK_2, INK_3, SURFACE = "#1a1a19", "#5c5b55", "#8a8980", "#ffffff"
PCOL = {"broad": BLUE, "bilateral": ORANGE}


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


def omnibus(pheno, cohort):
    return pd.read_csv(PHASE4 / f"omnibus_{pheno}_{cohort}.tsv", sep="\t")


def fig_manhattan():
    pos = pd.read_csv(HERE.parent.parent.parent
                      / "gene_burden_pmbb_bilateral/phase_5/results/gene_positions.tsv",
                      sep="\t")[["Region", "CHR", "POS"]]
    fig, axes = plt.subplots(3, 2, figsize=(13, 9), sharex=True, sharey=True)
    for r, cohort in enumerate(COHORTS):
        for c, pheno in enumerate(PHENOS):
            ax = axes[r, c]
            d = omnibus(pheno, cohort).merge(pos, on="Region", how="left",
                                             suffixes=("", "_p"))
            d = d.dropna(subset=["POS"]).copy()
            d["CHR"] = pd.to_numeric(d["CHR_p"] if "CHR_p" in d else d["CHR"],
                                     errors="coerce")
            d = d.dropna(subset=["CHR"]).sort_values(["CHR", "POS"])
            off, ticks, labels, x0 = {}, [], [], 0
            for ch, grp in d.groupby("CHR"):
                off[ch] = x0
                ticks.append(x0 + grp.POS.max() / 2)
                labels.append(int(ch))
                x0 += grp.POS.max() + 2e7
            x = d.POS + d.CHR.map(off)
            shade = np.where(d.CHR.astype(int) % 2 == 1, PCOL[pheno], PALE)
            ax.scatter(x, -np.log10(d.omnibus_p), s=3, linewidths=0, c=shade,
                       rasterized=True)
            bar = float(d.bonferroni_bar.iloc[0])
            ax.axhline(-np.log10(bar), color="#c0392b", lw=1.0,
                       label=f"Bonferroni {bar:.2e}")
            best = d.nsmallest(1, "omnibus_p").iloc[0]
            ax.annotate(best.Region, (x[best.name], -np.log10(best.omnibus_p)),
                        textcoords="offset points", xytext=(0, 7), ha="center",
                        fontsize=7, color=INK)
            if r == 0:
                ax.set_title(pheno, fontsize=10, fontweight="bold")
            if c == 0:
                ax.set_ylabel(f"{cohort}\n$-\\log_{{10}}$ omnibus $p$")
            ax.set_ylim(0, 6.6)
            ax.legend(loc="upper right", fontsize=7)
    for ax in axes[-1]:
        ax.set_xticks(ticks)
        ax.set_xticklabels(labels, fontsize=7)
        ax.set_xlabel("chromosome")
    fig.suptitle("Gene burden across the genome — no gene reaches its bar in any of "
                 "the six cells", fontsize=11.5, x=0.06, ha="left", y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(FIG / "manhattan.png")
    plt.close(fig)
    return "manhattan.png"


def fig_clingen():
    rep = json.loads((HERE / "results/02_enrichment.json").read_text())
    cg = pd.read_csv(HERE / "data/clingen_hl_genes.tsv", sep="\t")
    gene_set = set(cg[cg.best_classification.isin(["Definitive", "Strong"])].gene)

    fig, axes = plt.subplots(3, 1, figsize=(11, 6.4), sharex=True)
    for ax, cohort in zip(axes, COHORTS):
        for i, pheno in enumerate(PHENOS):
            d = omnibus(pheno, cohort).sort_values("omnibus_p").reset_index(drop=True)
            d["rank"] = np.arange(1, len(d) + 1)
            hit = d[d.Region.isin(gene_set)]
            y = 1 - i * 0.55
            ax.scatter(hit["rank"], np.full(len(hit), y), s=26, marker="|",
                       linewidths=1.3, color=PCOL[pheno])
            e = rep["arms"][f"{pheno}_{cohort}"]["definitive_strong"]["depths"]["50"]
            # label ABOVE the strip: inline, the ticks run through the text
            ax.text(1.05, y + 0.17, f"{pheno}: {e['observed']} in top 50 "
                                    f"(expected {e['expected']:.2f}, "
                                    f"p = {e['p']:.1e})",
                    fontsize=7.5, va="center", color=INK_2)
        ax.axvline(50, color="#c0392b", lw=1.0, ls="--")
        ax.set_xscale("log")
        ax.set_xlim(0.9, 2e4)
        ax.set_ylim(0.15, 1.42)
        ax.set_yticks([])
        ax.set_ylabel(cohort, rotation=0, ha="right", va="center", fontsize=9.5)
        ax.grid(axis="y", alpha=0)
    axes[0].text(55, 1.33, "top 50", fontsize=7.5, color="#c0392b")
    axes[-1].set_xlabel("rank of the gene among ~17,900, by omnibus p (log scale)")
    fig.suptitle("Where the 93 tested ClinGen deafness genes fall in each ranking",
                 fontsize=11.5, x=0.06, ha="left", y=1.03)
    fig.text(0.06, 0.985, "each tick is one gene. The claim is about the ranking, not "
             "about any gene: none of them is significant on its own.",
             fontsize=8, color=INK_2, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(FIG / "clingen.png")
    plt.close(fig)
    return "clingen.png"


def fig_vs_first_pass():
    old_base = {"broad": HERE.parent.parent.parent
                / "gene_burden_elena_replication/phase_4/results",
                "bilateral": HERE.parent.parent.parent
                / "gene_burden_pmbb_bilateral/phase_5/results"}
    fig, axes = plt.subplots(2, 3, figsize=(12, 7))
    for r, pheno in enumerate(PHENOS):
        for c, cohort in enumerate(COHORTS):
            ax = axes[r, c]
            new = omnibus(pheno, cohort)[["Region", "omnibus_p"]]
            old = pd.read_csv(old_base[pheno] / f"omnibus_{cohort}.tsv",
                              sep="\t")[["Region", "omnibus_p"]]
            m = new.merge(old, on="Region", suffixes=("_new", "_old"))
            x, y = -np.log10(m.omnibus_p_old), -np.log10(m.omnibus_p_new)
            hi = float(max(x.max(), y.max())) * 1.05
            ax.plot([0, hi], [0, hi], color=INK_3, lw=0.8, ls="--", zorder=1)
            ax.scatter(x, y, s=3, linewidths=0, color=PCOL[pheno], alpha=0.35,
                       rasterized=True, zorder=2)
            rho = float(pd.Series(x).corr(pd.Series(y), method="spearman"))
            ax.set_title(f"{pheno} · {cohort}\nρ = {rho:.3f}  n = {len(m):,}",
                         fontsize=8.5)
            if c == 0:
                ax.set_ylabel("rebuild  $-\\log_{10}p$")
            if r == 1:
                ax.set_xlabel("first pass  $-\\log_{10}p$")
    fig.suptitle("Changing the PCs: release exome PCA at 80% variance against the "
                 "first pass's own PCA", fontsize=11.5, x=0.06, ha="left", y=1.0)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(FIG / "vs_first_pass.png")
    plt.close(fig)
    return "vs_first_pass.png"


def main():
    style()
    FIG.mkdir(parents=True, exist_ok=True)
    for f in (fig_manhattan, fig_clingen, fig_vs_first_pass):
        print(f"  {f()}")


if __name__ == "__main__":
    main()
