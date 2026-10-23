#!/usr/bin/env python3
"""
Phase 5 (BRaVa) -- the pre-declared test, and the control.

PREMISES.md P5 fixed this before anything ran: ClinGen Definitive/Strong genes in the
TOP 50 of the broad combined ranking, hypergeometric p, exactly as
`pmbb_hl/rare_variant/phase_5/scripts/02_enrichment.py` computes it.

  p < 0.05   the enrichment survives a mask set that shares a fifth of its variants
  p >= 0.05  it does not, and that must be reported as prominently as the original
             finding, with the original re-described as annotation-dependent

Depths of 100 and 250 and the looser Moderate tier are computed and printed, because
reporting only the pre-declared number would hide whether it sits on a cliff. They are
NOT the test. P5 says so, and this script prints which is which so a reader cannot
mistake a secondary number for the declared one.

THE SYNONYMOUS CONTROL runs through the same function. It is expected to show nothing
(P4). A synonymous mask that enriches would mean something structural separates cases
from controls, and no other check in this project would catch it.

NOTHING IS RE-DERIVED. The ClinGen list is the file pmbb_hl fetched and used, read by
path rather than re-downloaded, so the two arms are judged against the same 93 genes.
A list refreshed between the two runs would change the answer for reasons that have
nothing to do with annotation.
"""
import json
import sys
from pathlib import Path

import pandas as pd
from scipy.stats import hypergeom

HERE = Path(__file__).resolve().parents[1]
PHASE4 = HERE.parent / "phase_4/results"
PMBB_HL = HERE.parent.parent / "pmbb_hl/rare_variant"
CLINGEN = PMBB_HL / "phase_5/data/clingen_hl_genes.tsv"
OUT = HERE / "results"

DECLARED_DEPTH = 50                 # PREMISES.md P5, fixed before the run
DECLARED_TIER = "definitive_strong"
DECLARED_ARM = "burden_broad_combined"
DEPTHS = [50, 100, 250]
TIERS = {"definitive_strong": ["Definitive", "Strong"],
         "plus_moderate": ["Definitive", "Strong", "Moderate"]}
ARMS = ["burden_broad_combined", "burden_broad_EUR",
        "synonymous_broad_combined", "synonymous_broad_EUR"]


def enrich(ranked, gene_set, depth):
    """Hypergeometric, identical in form to pmbb_hl's: of N tested genes, K are on the
    list; how surprising is the observed count in the top `depth`?"""
    N = len(ranked)
    K = len(gene_set & set(ranked.Region))
    obs = int(ranked.head(depth).Region.isin(gene_set).sum())
    exp = depth * K / N
    p = float(hypergeom.sf(obs - 1, N, K, depth)) if obs else 1.0
    return {"tested": N, "on_list": K, "observed": obs,
            "expected": float(exp), "p": p}


def main():
    if not CLINGEN.exists():
        sys.exit(f"REFUSING: {CLINGEN} not found. This arm must be judged against the "
                 f"same ClinGen list pmbb_hl used, not a freshly fetched one.")
    cg = pd.read_csv(CLINGEN, sep="\t")
    OUT.mkdir(parents=True, exist_ok=True)

    report = {"declared_test": {"arm": DECLARED_ARM, "tier": DECLARED_TIER,
                                "depth": DECLARED_DEPTH,
                                "source": "PREMISES.md P5"},
              "clingen_list": str(CLINGEN), "arms": {}}

    for arm in ARMS:
        f = PHASE4 / f"omnibus_{arm}.tsv"
        if not f.exists():
            sys.exit(f"REFUSING: {f} not found -- run phase 4 first")
        ranked = pd.read_csv(f, sep="\t").sort_values("omnibus_p")
        report["arms"][arm] = {}
        for tier, classes in TIERS.items():
            gs = set(cg[cg.best_classification.isin(classes)].gene)
            report["arms"][arm][tier] = {
                str(d): enrich(ranked, gs, d) for d in DEPTHS}

    # --- the declared test, alone and labelled
    d = report["arms"][DECLARED_ARM][DECLARED_TIER][str(DECLARED_DEPTH)]
    verdict = "SURVIVES" if d["p"] < 0.05 else "DOES NOT SURVIVE"
    report["declared_test"].update(d)
    report["declared_test"]["verdict"] = verdict

    print("=" * 72)
    print("THE PRE-DECLARED TEST  (PREMISES.md P5, fixed before the run)")
    print(f"  {DECLARED_ARM}, {DECLARED_TIER}, top {DECLARED_DEPTH}")
    print(f"  observed {d['observed']} of {d['on_list']} listed genes, "
          f"{d['expected']:.2f} expected among {d['tested']:,} tested")
    print(f"  p = {d['p']:.3e}   ->   {verdict}")
    print("=" * 72)

    # --- everything else, explicitly secondary
    print("\nsecondary, NOT the declared test:\n")
    print(f"  {'arm':<28}{'tier':<20}{'depth':>7}{'obs':>6}{'exp':>8}{'p':>12}")
    for arm in ARMS:
        for tier in TIERS:
            for dep in DEPTHS:
                e = report["arms"][arm][tier][str(dep)]
                mark = " <-- declared" if (arm == DECLARED_ARM
                                           and tier == DECLARED_TIER
                                           and dep == DECLARED_DEPTH) else ""
                print(f"  {arm:<28}{tier:<20}{dep:>7}{e['observed']:>6}"
                      f"{e['expected']:>8.2f}{e['p']:>12.2e}{mark}")

    # --- the control
    syn = [report["arms"][a][DECLARED_TIER][str(DECLARED_DEPTH)]
           for a in ARMS if a.startswith("synonymous")]
    flat = all(s["p"] >= 0.05 for s in syn)
    report["synonymous_control_flat"] = flat
    print(f"\nsynonymous control: " + (
        "flat, as PREMISES.md P4 expected"
        if flat else
        "*** NOT FLAT. P4 expected nothing here. Something structural separates "
        "cases from controls."))

    (OUT / "01_enrichment.json").write_text(json.dumps(report, indent=2) + "\n")
    print(f"\nwritten: results/01_enrichment.json")


if __name__ == "__main__":
    main()
