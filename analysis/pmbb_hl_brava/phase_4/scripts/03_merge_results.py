#!/usr/bin/env python3
"""
Phase 4 (BRaVa) -- merge, deduplicate, combine, apply the bars.

THE ARITHMETIC IS IMPORTED, NOT COPIED. `acat` comes from
`pmbb_hl/rare_variant/phase_4/scripts/04_merge_results.py` by import, and the
deduplication key and the Benjamini-Hochberg call are written here to match it exactly.
The whole point of this arm is that only the MASK differs; a merge that combined
p-values even slightly differently would contaminate the comparison with its own
arithmetic, and the difference would be invisible because both numbers would look
reasonable.

TWO SETS, ONE PASS

  burden      the pLOF/pDM grid -- 3 masks x 3 MAF = 9 cells per gene
  synonymous  the control -- 1 mask x 3 MAF = 3 cells per gene

The synonymous set is merged by the same code and judged against the same bars. It is
expected to produce nothing (PREMISES.md P4); running it through a separate path would
make a flat result less meaningful, because a different path could flatten it by itself.

DEDUPLICATION, as in pmbb_hl: a cell whose (max_MAF, MAC, Pvalue) matches another cell
of the same gene is the same test reported twice and enters the combination once. Note
this is the SAME-CUTOFF axis only; the open question about identical cells across
cutoffs (pmbb_hl phase 4 FINDINGS) is untouched here, deliberately, so that the two arms
are combined by identical rules.
"""
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import pandas as pd
from statsmodels.stats.multitest import multipletests

HERE = Path(__file__).resolve().parents[1]
STEP2 = HERE / "results/step2"
OUT = HERE / "results"
PMBB_HL_MERGE = (HERE.parent.parent / "pmbb_hl/rare_variant/phase_4/scripts"
                 / "04_merge_results.py")
COHORTS = ["combined", "EUR"]
SETS = {"burden": 9, "synonymous": 3}      # expected cells per gene, for the report
CHROMS = range(1, 23)
ALPHA = 0.05
Q_BAR = 0.10


def _borrow_acat():
    """Import pmbb_hl's acat rather than restating it. If that file moves or its
    function is renamed, this refuses instead of silently falling back to a local
    copy that could drift."""
    if not PMBB_HL_MERGE.exists():
        sys.exit(f"REFUSING: {PMBB_HL_MERGE} not found. The BRaVa arm combines "
                 f"p-values with pmbb_hl's function so the two are comparable; it "
                 f"will not substitute its own.")
    spec = importlib.util.spec_from_file_location("pmbb_hl_merge", PMBB_HL_MERGE)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    if not hasattr(mod, "acat"):
        sys.exit(f"REFUSING: {PMBB_HL_MERGE} no longer defines acat().")
    return mod.acat


acat = _borrow_acat()


def in_flight():
    try:
        out = subprocess.run(["bjobs", "-J", "bv_s2"], capture_output=True,
                             text=True, timeout=60).stdout
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return None
    return [l for l in out.splitlines() if " RUN " in l or " PEND " in l]


def load(sub, tag):
    parts = []
    for c in CHROMS:
        f = STEP2 / sub / tag / f"{tag}_chr{c}.txt"
        if not f.exists():
            raise FileNotFoundError(f"{sub}/{tag} chr{c}: absent")
        d = pd.read_csv(f, sep="\t")
        if d.empty:
            raise ValueError(f"{sub}/{tag} chr{c}: header only")
        d["CHR"] = c
        parts.append(d)
    return pd.concat(parts, ignore_index=True)


def main():
    live = in_flight()
    if live:
        sys.exit(f"REFUSING: {len(live)} step 2 tasks still in flight.")
    OUT.mkdir(parents=True, exist_ok=True)
    report = {"arms": {}, "acat_source": str(PMBB_HL_MERGE)}

    for sub, expect_cells in SETS.items():
        for cohort in COHORTS:
            tag = f"broad_{cohort}"
            key = f"{sub}_{tag}"
            d = load(sub, tag)

            # SAIGE's own Cauchy row combines the cells INCLUDING duplicates; the
            # omnibus is recomputed here after removing them, as in pmbb_hl.
            cells = d[d.Group != "Cauchy"].copy()
            n_cells = len(cells)
            keyed = cells.assign(
                _k=list(zip(cells.max_MAF.astype(str),
                            cells.MAC.round(10), cells.Pvalue.round(12))))
            ded = keyed.drop_duplicates(subset=["Region", "_k"])
            n_dropped = n_cells - len(ded)

            om = (ded.groupby("Region")
                     .agg(omnibus_p=("Pvalue", lambda s: acat(s.values)),
                          best_cell_p=("Pvalue", "min"),
                          cells_used=("Pvalue", "size"),
                          CHR=("CHR", "first"))
                     .reset_index())
            bar = ALPHA / len(om)
            om["bonferroni_bar"] = bar
            om["bonferroni_significant"] = om.omnibus_p < bar
            rej, q, _, _ = multipletests(om.omnibus_p.values, alpha=Q_BAR,
                                         method="fdr_bh")
            om["fdr_q"], om["fdr_significant"] = q, rej
            om = om.sort_values("omnibus_p")

            om.to_csv(OUT / f"omnibus_{key}.tsv", sep="\t", index=False)
            cells.drop(columns=[c for c in cells.columns if c.startswith("_")]) \
                 .to_csv(OUT / f"cells_{key}.tsv", sep="\t", index=False)

            best = om.iloc[0]
            report["arms"][key] = {
                "genes": len(om), "cells": n_cells,
                "duplicate_cells_dropped": n_dropped,
                "median_cells_per_gene": float(om.cells_used.median()),
                "expected_cells_per_gene": expect_cells,
                "min_omnibus_p": float(best.omnibus_p),
                "best_gene": str(best.Region),
                "bonferroni_bar": float(bar),
                "bonferroni_hits": int(om.bonferroni_significant.sum()),
                "min_q": float(om.fdr_q.min()),
                "fdr_hits": int(om.fdr_significant.sum()),
            }
            e = report["arms"][key]
            print(f"  {key:<26}{e['genes']:>7,} genes  min p {e['min_omnibus_p']:.3e}  "
                  f"{e['best_gene']:<10} q {e['min_q']:.4f}  "
                  f"hits {e['bonferroni_hits']}/{e['fdr_hits']}  "
                  f"dup {e['duplicate_cells_dropped']:,}")

    (OUT / "03_merge.json").write_text(json.dumps(report, indent=2) + "\n")
    print(f"\nwritten: results/omnibus_*.tsv, cells_*.tsv, 03_merge.json")

    syn = [report["arms"][k] for k in report["arms"] if k.startswith("synonymous")]
    if any(a["bonferroni_hits"] or a["fdr_hits"] for a in syn):
        print("\n  *** THE SYNONYMOUS CONTROL IS NOT FLAT. PREMISES.md P4 expected "
              "nothing here. Something structural separates cases from controls, and "
              "no other check in this project would have caught it.")
    else:
        print("\n  synonymous control: flat, as PREMISES.md P4 expected")


if __name__ == "__main__":
    main()
