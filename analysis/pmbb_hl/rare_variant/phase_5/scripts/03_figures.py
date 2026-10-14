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
    """Rank against p, which is the only axis pair on which BOTH bars are drawable.

    The earlier version used rank alone. On a rank axis Bonferroni has no position --
    it is a p-value threshold -- and FDR is not even a line: Benjamini-Hochberg requires
    p(k) <= (k/n)*q, so its boundary rises with rank. Here rank is x and p is y, which
    makes Bonferroni a horizontal line and the BH boundary the diagonal it actually is,
    and the distance from each ClinGen gene to both is readable directly.

    Only the ClinGen genes inside the top 50 are named: seven across the six arms, so
    the labels cost nothing, and naming all 93 would bury the point.
    """
    rep = json.loads((HERE / "results/02_enrichment.json").read_text())
    cg = pd.read_csv(HERE / "data/clingen_hl_genes.tsv", sep="\t")
    gene_set = set(cg[cg.best_classification.isin(["Definitive", "Strong"])].gene)
    Q = 0.10

    fig, axes = plt.subplots(2, 3, figsize=(13.5, 7.6), sharey=True)
    for r, pheno in enumerate(PHENOS):
        for c, cohort in enumerate(COHORTS):
            ax = axes[r, c]
            d = omnibus(pheno, cohort).sort_values("omnibus_p").reset_index(drop=True)
            d["rank"] = np.arange(1, len(d) + 1)
            n = len(d)
            bar = float(d.bonferroni_bar.iloc[0])

            ax.plot(d["rank"], -np.log10(d.omnibus_p), color=INK_3, lw=1.0,
                    zorder=1, label="all genes")
            hit = d[d.Region.isin(gene_set)]
            ax.scatter(hit["rank"], -np.log10(hit.omnibus_p), s=22, zorder=3,
                       color=PCOL[pheno], edgecolors=SURFACE, linewidths=0.6,
                       label="ClinGen Definitive/Strong")

            k = np.arange(1, n + 1)
            ax.plot(k, -np.log10(k / n * Q), color="#1baf7a", lw=1.3, ls="-",
                    zorder=2, label=f"FDR boundary, q = {Q:g}")
            ax.axhline(-np.log10(bar), color="#c0392b", lw=1.3, zorder=2,
                       label="Bonferroni, 0.05/genes")
            ax.axvline(50, color=INK_3, lw=0.9, ls=":", zorder=1)

            # staggered: in the combined arm five labels land within a few ranks
            # of each other and stack illegibly at a fixed offset
            top = hit[hit["rank"] <= 50].sort_values("rank")
            for j, (_, g) in enumerate(top.iterrows()):
                dy = 9 + (j % 3) * 11
                ax.annotate(g.Region, (g["rank"], -np.log10(g.omnibus_p)),
                            textcoords="offset points", xytext=(7, dy),
                            fontsize=7.5, color=INK, fontweight="bold",
                            arrowprops=dict(arrowstyle="-", lw=0.5, color=INK_3,
                                            shrinkA=0, shrinkB=2))

            e = rep["arms"][f"{pheno}_{cohort}"]["definitive_strong"]["depths"]["50"]
            ax.set_xscale("log")
            ax.set_xlim(0.8, 3e4)
            ax.set_ylim(0, 6.4)
            ax.set_title(f"{pheno} · {cohort}\n{e['observed']} ClinGen in top 50 "
                         f"(expected {e['expected']:.2f}, p = {e['p']:.1e})",
                         fontsize=8.5)
            if c == 0:
                ax.set_ylabel("$-\\log_{10}$ omnibus $p$")
            if r == 1:
                ax.set_xlabel("rank among ~17,900 genes (log)")
    axes[0, 0].text(52, 0.25, "top 50", fontsize=7, color=INK_3)
    axes[0, 2].legend(loc="upper right", fontsize=7)

    fig.suptitle("Known deafness genes against both bars, in all six cells",
                 fontsize=11.5, x=0.055, ha="left", y=1.0)
    fig.text(0.055, 0.962,
             "the FDR boundary rises with rank because Benjamini-Hochberg requires "
             "p(k) \u2264 (k/n)\u00b7q — it is not a single threshold. A gene is "
             "significant where its point sits ABOVE a line.",
             fontsize=8, color=INK_2, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.945))
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
