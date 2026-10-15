#!/usr/bin/env python3
"""
Synthesis -- the pipeline as one slide.

NOT A CHART. There is no data, so no axis, no legend, no palette validation. What
governs instead is that the reader must follow the flow without the document beside
them, and that the two arms read as one pipeline forking late.

LAYOUT IS COMPUTED, NOT GUESSED. The first version fixed box heights as fractions of
the axes while the text inside was sized in points, so nothing scaled together and
every box overflowed. Here each box's height is derived from its own line count in
points and converted once, and boxes are stacked top-down by a cursor. Adding a line to
any box now grows that box instead of spilling out of it.
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

HERE = Path(__file__).resolve().parents[1]
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK_2, INK_3 = "#1a1a19", "#5c5b55", "#8a8980"
SURFACE, SHADE = "#ffffff", "#f5f4f0"

FIG_W, FIG_H = 13.6, 10.4
BODY_PT, TITLE_PT = 7.6, 9.0
LINE_PT = BODY_PT * 1.62
PAD_PT = 9.0


def pt(v):
    """points -> axes fraction of the vertical span"""
    return v / (FIG_H * 72)


def box_height(n_lines):
    return pt(TITLE_PT * 1.9 + n_lines * LINE_PT + 2 * PAD_PT)


class Layout:
    def __init__(self, ax, top):
        self.ax, self.y = ax, top

    def place(self, x, w, title, lines, accent=INK_3, fill=SHADE, advance=True):
        h = box_height(len(lines))
        y = self.y - h
        self.ax.add_patch(FancyBboxPatch(
            (x, y), w, h, boxstyle="round,pad=0.004,rounding_size=0.012",
            linewidth=1.3, edgecolor=accent, facecolor=fill, zorder=2))
        self.ax.text(x + 0.014, y + h - pt(PAD_PT), title, fontsize=TITLE_PT,
                     fontweight="bold", color=INK, va="top", zorder=3)
        self.ax.text(x + 0.014, y + h - pt(PAD_PT + TITLE_PT * 1.75),
                     "\n".join(lines), fontsize=BODY_PT, color=INK_2, va="top",
                     linespacing=1.62, zorder=3, family="DejaVu Sans Mono")
        if advance:
            self.y = y
        return y, h

    def gap(self, points=16):
        self.y -= pt(points)

    def arrow(self, x, y0, y1, color=INK_3):
        self.ax.add_patch(FancyArrowPatch(
            (x, y0), (x, y1), arrowstyle="-|>", mutation_scale=11,
            linewidth=1.4, color=color, shrinkA=0, shrinkB=0, zorder=1))


def main():
    fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    fig.patch.set_facecolor(SURFACE)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.93, bottom=0.01)

    L = Layout(ax, 0.975)
    LEFT, RIGHT, FULL = 0.030, 0.515, 0.955
    HALF = 0.455

    y, h = L.place(LEFT, FULL, "PMBB-Release-2026-4.0", [
        "phecodeX conditions  655,946 SO_39x rows     OMOP observation  25,094 tinnitus rows",
        "covariates  age · sex · batch · exome_PC1-20  Exome/PCA  by-ancestry, 20 PCs each",
        "VEP annotations  106 GB                       genotypes  all_variants_NF"])
    L.gap(); L.arrow(0.5, y, L.y)

    y, h = L.place(LEFT, FULL, "PHASE 1   who is in the study", [
        "broad      SO_396 or any child, ≥ 2 distinct dates     6,752 cases · 57,507 analysable",
        "bilateral  SO_396.2 AND SO_396.8 on the SAME date      3,164 cases · 53,919 analysable",
        "tinnitus   SO_397 or ICD 388.3x/H93.1x, ≥ 2 dates      2,859 cases · 53,614 analysable",
        "controls   50,755 — no ear evidence from any source, SHARED by all three arms",
        "▷ refuses unless every bilateral case is also a broad case",
        "▷ tinnitus is NOT nested — 66% overlap, measured not asserted"], accent=INK)
    L.gap(); split_from = L.y + pt(16)
    L.arrow(0.5, y, L.y)

    top = L.y
    y1, _ = L.place(LEFT, HALF, "PHASE 2   burden masks", [
        "22 tasks over the release VEP",
        "pLOF 17,841 · pDM 17,623 · pLOF_pDM 17,945 genes",
        "one build serves all three phenotypes — a mask depends",
        "on annotation, not on who is a case"], accent=BLUE)
    L.y = top
    y2, _ = L.place(RIGHT, HALF, "PHASE 2   ExWAS MAF census", [
        "198 tasks — once per phenotype",
        "MAF ≥ 0.01 → 70,287 to 118,067 variants",
        "MAF is a property of the SAMPLE SET, so each",
        "cohort needs its own count"], accent=ORANGE)
    L.gap()
    L.arrow(0.25, y1, L.y); L.arrow(0.75, y2, L.y)

    y, h = L.place(LEFT, FULL, "PHASE 3   covariates and GRM        shared by all three arms", [
        "PCs from the RELEASE's exome PCA at 80% of the variance   combined 11 · EUR 15 · AFR 16",
        "GRM rebuilt for EUR and AFR; combined reused — it already covers all 70,925",
        "▷ refuses if anyone lacks a PC, which is the premise of using the release's"],
        accent=INK)
    L.gap(); L.arrow(0.5, y, L.y)

    y, h = L.place(LEFT, FULL, "PHASE 4   SAIGE  ·  step 1", [
        "9 null models, 3 phenotypes × 3 cohorts    logit P(case) = Xβ + b,   b ~ N(0, τ·GRM)",
        "fitted once and reused by every step 2 task — relatedness does not depend on the gene"],
        accent=INK)
    L.gap()
    L.arrow(0.25, y, L.y, color=BLUE); L.arrow(0.75, y, L.y, color=ORANGE)

    top = L.y
    y1, _ = L.place(LEFT, HALF, "PHASE 4   burden  ·  step 2", [
        "198 tasks, the whole grid per call",
        "3 masks × 3 MAF = 9 cells per gene",
        "deduplicate   31,765 cells removed",
        "Cauchy/ACAT   → one p per gene",
        "bars          0.05/genes  +  FDR q < 0.1",
        "17,943 genes tested"], accent=BLUE)
    L.y = top
    y2, _ = L.place(RIGHT, HALF, "PHASE 4   ExWAS  ·  step 2", [
        "198 tasks, one test per variant",
        "MAF ≥ 0.01, no mask, no annotation",
        "step 1 INHERITED — the null model does not",
        "depend on the MAF filter",
        "bars          5e-8  +  0.05/variants",
        "70,287 to 118,067 variants tested"], accent=ORANGE)
    L.gap()
    L.arrow(0.25, y1, L.y); L.arrow(0.75, y2, L.y)

    L.place(LEFT, FULL, "PHASES 5 and 6   results and power", [
        "calibration (λ by region) · ClinGen enrichment · figures",
        "detectable effect per gene and per variant · what cohort size would change it",
        "MAC-matched test of the enrichment — the check that tried to break the one result"],
        accent=AQUA)

    fig.text(0.030, 0.978, "The pipeline", fontsize=15.5, fontweight="bold", color=INK)
    fig.text(0.030, 0.952, "two analyses × three phenotypes × three cohorts = 18 cells. "
             "Everything above step 2 is shared; the arms fork there, and the ExWAS "
             "does not refit the null model.", fontsize=8.8, color=INK_2)

    out = HERE / "figures/pipeline.png"
    fig.savefig(out, dpi=190, facecolor=SURFACE)
    plt.close(fig)
    print(f"  {out.relative_to(HERE)}")


if __name__ == "__main__":
    main()
