#!/usr/bin/env python3
"""
Phase 4 -- merge, deduplicate, combine, and apply both bars.

THE THREE LEVELS OF MULTIPLICITY, and who pays for each:

  SKAT-O's rho sweep        SAIGE, internally
  the 9 cells of one gene   Cauchy/ACAT, here
  the ~17,900 genes         Bonferroni and FDR, here

The nine cells are not nine hypotheses. They are nine ways of asking one question about
one gene, and they are strongly correlated -- the MAF cutoffs are nested, so everything
in 0.0001 is also in 0.001 and 0.01. Bonferroni across them would charge nine times for
one question. ACAT handles dependent p-values without needing to know the dependence.

DEDUPLICATION, BEFORE THE COMBINATION. Raised by Nikki Palmiero on 2026-10-14:

  "Let's say these two have the exact same rows [...] all of the numbers are the same.
   What that means is there's no pDM. It's doing the test again. So we would drop this
   because it's a duplicate. SAIGE doesn't notice that it's a duplicate."

Two mechanisms produce them, and both are real:

  mask-level   a gene with no pDM variants has pLOF:pDM identical to pLOF. 445 of
               17,945 genes, measured on the phase 2 masks.
  MAF-level    after the cutoff, the surviving variants can coincide even when the full
               masks differ. Measured on the first pass: 1,606 duplicate cells against
               the 1,335 the mask-level count predicts, so 271 come from this.

WHY IT USUALLY DOES NOT MATTER, AND WHEN IT DOES. ACAT is a MEAN of transformed
p-values. If every cell is duplicated exactly once -- the common case, three masks
collapsing to two -- the mean is unchanged and the omnibus is identical. The effect
appears only when duplication is UNEVEN. Measured on the first pass: 550 genes affected,
median ratio 1.000, range 0.69 to 1.73, and no verdict changed either way.

It is corrected regardless, because the direction of the error is the dangerous one: if
the duplicated cell is a gene's most extreme, counting it twice pulls the omnibus toward
significance. Harmless while nothing is near a bar; not harmless when something is.

SAIGE'S OWN Cauchy ROW IS DISCARDED. It combines all nine cells including duplicates --
verified on the first pass's chr21, where it equals an ACAT over all cells in 36 of 36
genes that contain duplicates.
"""
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

HERE = Path(__file__).resolve().parents[1]
STEP2 = HERE / "results/step2"
OUT = HERE / "results"
PHENOS = ["broad", "bilateral", "tinnitus"]
COHORTS = ["combined", "EUR", "AFR"]
CHROMS = range(1, 23)
Q_BAR = 0.10
ALPHA = 0.05


def acat(p):
    p = np.clip(np.asarray(p, dtype=float), 1e-300, 1 - 1e-16)
    return float(0.5 - np.arctan(np.mean(np.tan((0.5 - p) * np.pi))) / np.pi)


def in_flight():
    try:
        out = subprocess.run(["bjobs"], capture_output=True,
                             text=True, timeout=60).stdout
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return None
    return [l for l in out.splitlines()[1:] if l.strip()]


def load(tag):
    """All 22 chromosomes for one phenotype/cohort, refusing on any gap."""
    parts = []
    for c in CHROMS:
        f = STEP2 / tag / f"{tag}_chr{c}.txt"
        if not f.exists():
            raise FileNotFoundError(f"{tag} chr{c}: absent")
        d = pd.read_csv(f, sep="\t")
        if d.empty:
            raise ValueError(f"{tag} chr{c}: header only")
        d["CHR"] = c
        parts.append(d)
    return pd.concat(parts, ignore_index=True)


def main():
    live = in_flight()
    if live:
        sys.exit(f"REFUSING: {len(live)} step 2 tasks still in flight.")
    if live is None:
        print("WARNING: could not query LSF; relying on the structural checks\n")

    report = {"q_bar": Q_BAR, "alpha": ALPHA, "arms": {}}

    for pheno in PHENOS:
        for cohort in COHORTS:
            tag = f"{pheno}_{cohort}"
            try:
                d = load(tag)
            except (FileNotFoundError, ValueError) as e:
                sys.exit(f"REFUSING: {e}")

            saige_cauchy = d[d.Group == "Cauchy"].set_index("Region").Pvalue
            cells = d[d.Group != "Cauchy"].copy()
            cells["Pvalue"] = pd.to_numeric(cells.Pvalue, errors="coerce")
            cells = cells[cells.Pvalue.notna() &
                          (cells.Pvalue > 0) & (cells.Pvalue <= 1)]
            n_cells = len(cells)

            # --- the deduplication, before any combination
            keyed = cells.assign(
                _k=list(zip(cells.max_MAF.astype(str),
                            cells.MAC.round(10), cells.Pvalue.round(12))))
            ded = keyed.drop_duplicates(subset=["Region", "_k"])
            n_dropped = n_cells - len(ded)

            # which mask pairs collided, for the record
            pairs = {}
            for (_, _), grp in keyed[keyed.duplicated(
                    subset=["Region", "_k"], keep=False)].groupby(["Region", "_k"]):
                key = " == ".join(sorted(grp.Group.unique()))
                pairs[key] = pairs.get(key, 0) + len(grp) - 1

            # --- omnibus per gene, over the deduplicated cells
            om = (ded.groupby("Region")
                     .agg(omnibus_p=("Pvalue", lambda s: acat(s.values)),
                          best_cell_p=("Pvalue", "min"),
                          cells_used=("Pvalue", "size"),
                          CHR=("CHR", "first"))
                     .reset_index())

            n_genes = len(om)
            bar = ALPHA / n_genes
            om["bonferroni_bar"] = bar
            om["bonferroni_significant"] = om.omnibus_p < bar
            rej, q, _, _ = multipletests(om.omnibus_p.values, alpha=Q_BAR,
                                         method="fdr_bh")
            om["fdr_q"] = q
            om["fdr_significant"] = rej
            om = om.sort_values("omnibus_p")

            dest = OUT / f"omnibus_{tag}.tsv"
            om.to_csv(dest, sep="\t", index=False)
            cells.drop(columns=[c for c in cells.columns if c.startswith("_")]) \
                 .to_csv(OUT / f"cells_{tag}.tsv", sep="\t", index=False)

            # --- does SAIGE's own Cauchy differ from ours, and where?
            cmp_n = cmp_diff = 0
            for g, s in saige_cauchy.items():
                try:
                    s = float(s)
                except (TypeError, ValueError):
                    continue
                mine = om.loc[om.Region == g, "omnibus_p"]
                if mine.empty or not (0 < s <= 1):
                    continue
                cmp_n += 1
                if abs(float(mine.iloc[0]) - s) / max(s, 1e-300) > 0.01:
                    cmp_diff += 1

            report["arms"][tag] = {
                "genes": n_genes, "cells": n_cells,
                "duplicate_cells_dropped": n_dropped,
                "duplicate_pairs": pairs,
                "bonferroni_bar": bar,
                "bonferroni_hits": int(om.bonferroni_significant.sum()),
                "fdr_hits": int(om.fdr_significant.sum()),
                "min_omnibus_p": float(om.omnibus_p.min()),
                "min_q": float(om.fdr_q.min()),
                "best_gene": om.Region.iloc[0],
                "differs_from_saige_cauchy": cmp_diff,
                "compared_to_saige_cauchy": cmp_n,
            }
            r = report["arms"][tag]
            print(f"{tag}")
            print(f"  genes {n_genes:>7,}   cells {n_cells:>8,}   "
                  f"duplicates dropped {n_dropped:>5,}")
            print(f"  bar {bar:.3e}   min p {r['min_omnibus_p']:.3e} ({r['best_gene']})"
                  f"   min q {r['min_q']:.4f}")
            print(f"  Bonferroni hits {r['bonferroni_hits']}   FDR hits "
                  f"{r['fdr_hits']}")
            print(f"  our omnibus differs from SAIGE's Cauchy in "
                  f"{cmp_diff:,} of {cmp_n:,} genes")
            if pairs:
                print(f"  collisions: {pairs}")
            print()

    (OUT / "04_merge.json").write_text(json.dumps(report, indent=2) + "\n")
    print(f"written: results/omnibus_*.tsv, cells_*.tsv, 04_merge.json")


if __name__ == "__main__":
    main()
