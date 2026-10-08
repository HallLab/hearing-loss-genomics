#!/usr/bin/env python3
"""
The phase 5 figures, organised the way the write-up reads.

One figure per question rather than per cohort: the Manhattan answers "does
anything cross the bar", the QQ answers "is the test calibrated". Each carries
all three cohorts, so a reader compares across them in one glance and a slide
can take either on its own.

An earlier version stacked Manhattan and QQ per cohort, which meant one image
served two sections of the document and neither could be used alone.

Output: figures/manhattan_omnibus.png, figures/qq_omnibus.png
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
import burden_plots as bp

FIG = Path(__file__).resolve().parents[2] / "figures"
FIG.mkdir(exist_ok=True)
bp.style()

LABEL = {"combined": "combinada", "EUR": "europeia", "AFR": "africana"}

# ---- Manhattan, one row per cohort ----
fig, axes = plt.subplots(3, 1, figsize=(13, 10))
for ax, cohort in zip(axes, bp.COHORTS):
    om = bp.omnibus(cohort).rename(columns={"omnibus_p": "Pvalue"})
    om["Cohort"] = cohort
    bp.manhattan(om, cohort, ax=ax, label_top=6)
    ax.set_title(f"{LABEL[cohort]} — {bp.N[cohort]:,} pessoas, "
                 f"{bp.CASES[cohort]:,} casos", fontsize=10)
    if ax is not axes[-1]:
        ax.set_xlabel("")
h, l = axes[0].get_legend_handles_labels()
fig.legend(h, l, loc="lower center", bbox_to_anchor=(0.5, -0.015))
fig.suptitle("Omnibus de Cauchy — um p-valor por gene", y=0.997)
fig.tight_layout()
fig.savefig(FIG / "manhattan_omnibus.png")
plt.close(fig)

# ---- QQ, three side by side ----
fig, axes = plt.subplots(1, 3, figsize=(13, 4.6))
for ax, cohort in zip(axes, bp.COHORTS):
    om = bp.omnibus(cohort)
    bp.qq(om.omnibus_p, ax=ax, title=f"{LABEL[cohort]} — {bp.CASES[cohort]:,} casos")
    if ax is not axes[0]:
        ax.set_ylabel("")
fig.suptitle("QQ do omnibus — a faixa cinza é o envelope nulo de 95%", y=1.0)
fig.tight_layout()
fig.savefig(FIG / "qq_omnibus.png")
plt.close(fig)

print("figures/manhattan_omnibus.png\nfigures/qq_omnibus.png")
