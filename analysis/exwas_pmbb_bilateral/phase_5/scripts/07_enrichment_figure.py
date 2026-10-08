#!/usr/bin/env python3
"""
The ClinGen figure -- the one piece of this analysis that is not a null.

Everything else in phase 5 says "nothing passed". This says something positive
and then takes it apart, and it was living only in tables.

FORM. Emphasis, not categorical: the five random draws are context (the null
distribution at this sample size) and the three real arms are the subject. One
accent hue plus a de-emphasis gray, which is the lowest-risk colour case there
is. Every row is directly labelled with its count and its genes, so identity
never rests on colour.

COLOUR. Taken unchanged from the dataviz skill's validated reference palette --
slot 1 blue #2a78d6 on surface #fcfcfb. The skill says to run its validator
rather than reason about it; `node` is not installed on this host, so that was
not possible. Using the documented values as-is, in the two-colour emphasis
form, is the mitigation.

Output: figures/clingen_enrichment.png
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parents[2]

ACCENT = "#2a78d6"
MUTED = "#9a9a95"
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"

pc = json.load(open(HERE / "phase_5/power_control/results/04_enrichment.json"))
co = json.load(open(HERE / "phase_5/complement/results/04_readout.json"))
by = {r["label"]: r for r in co["rows"]}

EXPECTED = 0.26
rows = [
    ("fenótipo amplo", "6.752 casos", 5,
     "SIX1 · GJB3 · COCH · MYO6 · TMPRSS3", ACCENT),
    (None, None, None, None, None),                      # separador
    ("amplo, sorteio 1", "3.164", 2, "COCH · SIX1", MUTED),
    ("amplo, sorteio 2", "3.164", 0, "", MUTED),
    ("amplo, sorteio 3", "3.164", 1, "ACTG1", MUTED),
    ("amplo, sorteio 4", "3.164", 2, "COCH · TMPRSS3", MUTED),
    ("amplo, sorteio 5", "3.164", 3, "GJB3 · MYO6 · SIX1", MUTED),
    (None, None, None, None, None),
    ("restrito  (bilateral neuro.)", "3.164", 0, "", ACCENT),
    ("complemento  (os descartados)", "3.588", 1, "SIX1", ACCENT),
]

plt.rcParams.update({
    "figure.dpi": 110, "savefig.dpi": 170, "savefig.bbox": "tight",
    "font.size": 10, "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "text.color": INK, "axes.labelcolor": INK2, "axes.edgecolor": "#dddddd",
})
fig, ax = plt.subplots(figsize=(11, 5.6))

ys = []
for i, (label, n, val, genes, color) in enumerate(rows):
    y = len(rows) - i
    if label is None:
        continue
    ys.append(y)
    ax.plot([0, val], [y, y], color=color, lw=2, solid_capstyle="round", zorder=2)
    ax.scatter([val], [y], s=90, color=color, zorder=3,
               edgecolors=SURFACE, linewidths=2)
    ax.text(-0.22, y, label, ha="right", va="center", fontsize=10,
            color=INK if color == ACCENT else INK2,
            fontweight="bold" if color == ACCENT else "normal")
    ax.text(-0.22, y - 0.33, n, ha="right", va="center", fontsize=8, color=INK2)
    ax.text(val + 0.13, y, str(val), ha="left", va="center", fontsize=10,
            fontweight="bold", color=INK)
    if genes:
        ax.text(val + 0.45, y, genes, ha="left", va="center", fontsize=8.5,
                color=INK2)

ax.axvline(EXPECTED, color="#c0392b", lw=1.4, zorder=1)
ax.text(EXPECTED, max(ys) + 0.75, f"  esperado por acaso: {EXPECTED}",
        ha="left", va="center", fontsize=9, color="#c0392b")

ax.set_xlim(-0.05, 6.9)
ax.set_ylim(min(ys) - 0.9, max(ys) + 1.2)
ax.set_yticks([])
ax.set_xticks(range(0, 7))
ax.set_xlabel("genes ClinGen Definitive/Strong entre os 50 primeiros", color=INK2)
for s in ("top", "right", "left"):
    ax.spines[s].set_visible(False)
ax.grid(axis="x", color="#eeeeee", lw=0.8, zorder=0)
ax.set_axisbelow(True)

fig.suptitle("Genes de surdez conhecidos no topo da lista — coorte combinada",
             x=0.125, ha="left", fontsize=13, y=1.02)
fig.text(0.125, 0.955,
         "cinza: metades aleatórias do fenótipo amplo, o que o acaso produz nesse tamanho",
         ha="left", fontsize=9, color=INK2)
fig.savefig(HERE / "figures/clingen_enrichment.png")
print("figures/clingen_enrichment.png")
