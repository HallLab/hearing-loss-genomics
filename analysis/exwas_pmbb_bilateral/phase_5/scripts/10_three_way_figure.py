#!/usr/bin/env python3
"""
The three-way figure -- how the known-gene enrichment moves between arms.

Small multiples by cohort, three arms each, same visual language as
clingen_enrichment.png. Categorical colour with three series, which is inside
the safe range, and every row carries its count and genes as a direct label so
identity never rests on colour.

Palette taken unchanged from the dataviz skill's validated reference instance
(slots 1-3 on the light surface); `node` is not installed here, so its validator
could not be run and no colours were invented to compensate.

Output: figures/three_way_clingen.png
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parents[2]
r = json.load(open(HERE / "phase_5/results/09_three_way.json"))

ARMS = ["Elena", "Replication", "Bilateral"]
COLOR = {"Elena": "#eb6834", "Replication": "#2a78d6", "Bilateral": "#1baf7a"}
SUB = {"Elena": "original masks · 5/9/10 PCs",
       "Replication": "rebuilt masks · 5/4/3 PCs",
       "Bilateral": "+ bilateral sensorineural"}
SURFACE, INK, INK2 = "#fcfcfb", "#0b0b0b", "#52514e"

plt.rcParams.update({
    "figure.dpi": 110, "savefig.dpi": 170, "savefig.bbox": "tight",
    "font.size": 10, "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "text.color": INK, "axes.labelcolor": INK2, "axes.edgecolor": "#dddddd"})

fig, axes = plt.subplots(3, 1, figsize=(10.5, 7.6), sharex=True)
for ax, cohort in zip(axes, ["combined", "EUR", "AFR"]):
    e = r["cohorts"][cohort]["arms"]
    for i, arm in enumerate(ARMS):
        y = len(ARMS) - i
        v = e[arm]["clingen_in_top50"]
        ax.plot([0, v], [y, y], color=COLOR[arm], lw=2.5, solid_capstyle="round")
        ax.scatter([v], [y], s=95, color=COLOR[arm], zorder=3,
                   edgecolors=SURFACE, linewidths=2)
        ax.text(-0.15, y, arm, ha="right", va="center", fontsize=10,
                fontweight="bold")
        ax.text(-0.15, y - 0.3, SUB[arm], ha="right", va="center",
                fontsize=7.5, color=INK2)
        ax.text(v + 0.1, y, str(v), ha="left", va="center",
                fontsize=10, fontweight="bold")
        g = e[arm]["clingen_genes"]
        if g:
            ax.text(v + 0.38, y, " · ".join(g), ha="left", va="center",
                    fontsize=8, color=INK2)
    exp = e["Elena"]["clingen_expected"]
    ax.axvline(exp, color="#c0392b", lw=1.3, zorder=1)
    ax.set_ylim(0.3, 3.9)
    ax.set_yticks([])
    ax.set_xlim(-0.05, 6.6)
    ax.set_title(f"{cohort}", loc="left", fontsize=11, x=0.0)
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    ax.grid(axis="x", color="#eeeeee", lw=0.8)
    ax.set_axisbelow(True)

axes[0].text(r["cohorts"]["combined"]["arms"]["Elena"]["clingen_expected"], 3.7,
             "  expected by chance", ha="left", va="center",
             fontsize=8.5, color="#c0392b")
axes[-1].set_xlabel("ClinGen Definitive/Strong genes in the top 50", color=INK2)
fig.suptitle("Known deafness genes at the top — the three arms",
             x=0.33, ha="left", fontsize=13, y=1.0)
fig.text(0.33, 0.955, "each arm changes one thing from the one above it",
         ha="left", fontsize=9, color=INK2)
fig.tight_layout()
fig.savefig(HERE / "figures/three_way_clingen.png")
print("figures/three_way_clingen.png")
