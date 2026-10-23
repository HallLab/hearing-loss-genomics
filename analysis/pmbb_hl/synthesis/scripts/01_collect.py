#!/usr/bin/env python3
"""
Synthesis -- pull every reportable number out of the two arms' own results.

EIGHTEEN CELLS: {rare-variant burden, common-variant ExWAS} x {broad, bilateral,
tinnitus} x {combined, EUR, AFR}, plus the original run's six (broad only -- it never
ran the other two phenotypes), for twenty-four points in the figure. Nothing here is
retyped from a document; each number is read from the file the phase that produced it
wrote.

A COMPARABILITY WARNING THE TABLES CANNOT CARRY. The two arms do not report the same
quantity and are not judged by the same bar:

  burden   one omnibus p per GENE, Cauchy-combined over up to 9 mask x MAF cells,
           against 0.05/genes and an FDR boundary
  ExWAS    one p per VARIANT, no combination, against 5e-8 and 0.05/variants

Putting their minimum p-values side by side compares a minimum over 17,900 genes with
one over 70,000-118,000 variants. The ratio to each cell's OWN bar is the only figure
that travels between them, and that is what the synthesis figure plots.
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1].parent
OUT = Path(__file__).resolve().parents[1] / "results"
PHENOS = ["broad", "bilateral", "tinnitus"]
COHORTS = ["combined", "EUR", "AFR"]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    burden = json.loads((ROOT / "rare_variant/phase_4/results/04_merge.json").read_text())
    exwas = json.loads((ROOT / "exwas/phase_4/results/02_merge.json").read_text())
    enrich = json.loads((ROOT / "rare_variant/phase_5/results/02_enrichment.json").read_text())
    matched = json.loads((ROOT / "rare_variant/phase_6/results/02_matched_enrichment.json").read_text())
    calib = json.loads((ROOT / "exwas/phase_5/results/01_calibration.json").read_text())
    cohorts = json.loads((ROOT / "rare_variant/phase_1/results/01_cohorts.json").read_text())
    covar = json.loads((ROOT / "rare_variant/phase_3/results/01_covariates.json").read_text())

    out = {"cohorts": {}, "cells": {}}
    for ph in PHENOS:
        out["cohorts"][ph] = {
            "analysable": cohorts["arms"][ph]["analysable"],
            "cases": cohorts["arms"][ph]["cases"],
            "tested": {c: covar["arms"][ph][c]["written"] for c in COHORTS},
            "tested_cases": {c: covar["arms"][ph][c]["cases"] for c in COHORTS},
        }

    print(f"  {'cell':<26}{'units':>9}{'min p':>11}{'bar':>11}{'short by':>10}"
          f"{'hits':>6}")
    for ph in PHENOS:
        for c in COHORTS:
            tag = f"{ph}_{c}"
            b, e = burden["arms"][tag], exwas["arms"][tag]
            for arm, units, minp, bar, hits in (
                ("burden", b["genes"], b["min_omnibus_p"], b["bonferroni_bar"],
                 b["bonferroni_hits"] + b["fdr_hits"]),
                ("exwas", e["tested"], e["min_p"], 5e-8,
                 e["hits_primary"] + e["hits_secondary"]),
            ):
                key = f"{arm}_{tag}"
                out["cells"][key] = {
                    "arm": arm, "phenotype": ph, "cohort": c, "units": units,
                    "min_p": minp, "bar": bar, "short_by": minp / bar, "hits": hits,
                }
                if arm == "burden":
                    out["cells"][key].update({
                        "best_gene": b["best_gene"], "min_q": b["min_q"],
                        "duplicates_dropped": b["duplicate_cells_dropped"],
                        "clingen_top50": enrich["arms"][tag]["definitive_strong"]
                                               ["depths"]["50"]["observed"],
                        "clingen_p": enrich["arms"][tag]["definitive_strong"]
                                           ["depths"]["50"]["p"],
                    })
                else:
                    out["cells"][key].update({
                        "best_variant": e["best_variant"],
                        "secondary_bar": e["secondary_bar"],
                        "lambda_in_mhc": calib["arms"][tag]["lambda_inside_mhc"],
                        "lambda_outside_mhc": calib["arms"][tag]["lambda_outside_mhc"],
                    })
                print(f"  {key:<26}{units:>9,}{minp:>11.2e}{bar:>11.2e}"
                      f"{minp/bar:>9.1f}x{hits:>6}")

    # --- the original run, for comparison. Broad phenotype only: it never ran
    # bilateral, so that column is new work rather than replication.
    #
    # Her burden omnibus is NOT her delivered p-value. Her tables report per-CELL
    # p-values corrected within one MAF family (0.05/53,583 for combined, verified
    # from her Pvalue_Bonferroni column) across the 160,672 cell tests she ran. Placing
    # that beside our omnibus would compare a minimum over nine times more tests. So
    # the figure uses her omnibus recomputed with the same ACAT over the same nine
    # cells, from the first pass's three-way comparison, where all three are one
    # quantity judged by one rule.
    # The first pass now lives under pmbb_hl_v1/ -- it was moved there when this
    # rebuild took the pmbb_hl name. Checked explicitly so a future move says what
    # is missing instead of raising on read.
    tw = (ROOT.parent / "pmbb_hl_v1/gene_burden_pmbb_bilateral/phase_5/results"
          / "09_three_way.json")
    if not tw.exists():
        sys.exit(f"REFUSING: the first pass's three-way comparison is not at {tw}. "
                 "It holds the original run's omnibus recomputed by our method, "
                 "which is the only figure comparable to ours.")
    three_way = json.loads(tw.read_text())
    import csv
    out["elena"] = {}
    for c in COHORTS:
        k = three_way["cohorts"][c]["arms"]["Elena"]
        ps = []
        f = (ROOT.parent / "elena/HL_only_ExWAS"
             / f"HL_ExWAS_FullyAnnotated_AllVariants_{c}.csv")
        with f.open() as fh:
            for row in csv.DictReader(fh):
                try:
                    v = float(row["p.value"])
                except (TypeError, ValueError):
                    continue
                if 0 < v <= 1:
                    ps.append(v)
        out["elena"][c] = {
            "burden": {"units": k["genes"], "min_p": k["min_omnibus_p"],
                       "bar": k["bar"], "short_by": k["min_omnibus_p"] / k["bar"]},
            "exwas": {"units": len(ps), "min_p": min(ps), "bar": 5e-8,
                      "short_by": min(ps) / 5e-8},
        }
        e = out["elena"][c]
        print(f"  elena_{c:<20}burden {e['burden']['short_by']:>6.1f}x   "
              f"exwas {e['exwas']['short_by']:>7.0f}x")

    # --- the BRaVa annotation arm. Burden only, and only the two cohorts where the
    # ClinGen enrichment exists, so it is a PARTIAL series by construction -- see
    # pmbb_hl_brava/PREMISES.md P3. Carried with its coverage stated so the figure can
    # draw the gaps rather than quietly close them.
    bv = ROOT.parent / "pmbb_hl_brava/phase_4/results/03_merge.json"
    if bv.exists():
        b = json.loads(bv.read_text())["arms"]
        # WHICH COHORTS ARE READ FROM THE MERGE, not written down here. That arm was
        # scoped to combined and EUR, then expanded to AFR on 2026-10-23; a hardcoded
        # pair would have silently dropped the new cell and the figure would still
        # have rendered.
        have = [c for c in COHORTS if f"burden_broad_{c}" in b]
        out["brava"] = {"arm": "burden", "cohorts": have,
                        "note": "burden only -- the ExWAS has no mask, so an "
                                "annotation standard cannot change it",
                        "cells": {}}
        for c in have:
            a = b[f"burden_broad_{c}"]
            out["brava"]["cells"][c] = {
                "units": a["genes"], "min_p": a["min_omnibus_p"],
                "bar": a["bonferroni_bar"],
                "short_by": a["min_omnibus_p"] / a["bonferroni_bar"],
                "best_gene": a["best_gene"], "min_q": a["min_q"],
                "hits": a["bonferroni_hits"] + a["fdr_hits"]}
            print(f"  brava_{c:<20}burden {out['brava']['cells'][c]['short_by']:>6.1f}x")
    else:
        print(f"  (no BRaVa arm at {bv})")

    out["matched_enrichment"] = matched["arms"]
    (OUT / "consolidated.json").write_text(json.dumps(out, indent=2) + "\n")
    print(f"\n  written: results/consolidated.json")


if __name__ == "__main__":
    main()
