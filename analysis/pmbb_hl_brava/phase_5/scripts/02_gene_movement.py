#!/usr/bin/env python3
"""
Phase 5 (BRaVa) -- why each ClinGen gene moved, gene by gene.

The enrichment survived the change of annotation but THREE OF FIVE genes are the same.
That is the finding people will ask about, and "the masks differ" is not an answer. This
writes the answer: for every ClinGen Definitive/Strong gene, its rank and p under both
annotations, the cell that drove each, and the case/control allele counts underneath.

WHY THE ALLELE COUNTS ARE THE POINT. A burden p is a function of how the gene's alleles
split between cases and controls. At these counts -- tens of case alleles for most of
these genes -- swapping a handful of variants swaps a handful of alleles, and a handful
of alleles moves the p by an order of magnitude. The gene movement is not evidence that
one annotation is right; it is the same fragility phase 7 of the main write-up argues
from, shown directly.

OUTPUT is a table of all 93, plus a narrative file naming the seven genes that appear in
either top 50 and what happened to each.
"""
import json
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parents[1]
BRAVA4 = HERE.parent / "phase_4/results"
PMBB_HL = HERE.parent.parent / "pmbb_hl/rare_variant"
CLINGEN = PMBB_HL / "phase_5/data/clingen_hl_genes.tsv"
OUT = HERE / "results"
DEPTH = 50


def ranked(path):
    d = pd.read_csv(path, sep="\t").sort_values("omnibus_p").reset_index(drop=True)
    d["rank"] = d.index + 1
    return d.set_index("Region")


def cell_detail(cells, gene):
    """The cell that drove the gene, and the widest cell for context. The widest is
    worth carrying because a gene can look strong on a narrow cell and have nothing in
    the broad one -- which is a different situation from being strong throughout."""
    s = cells[cells.Region == gene]
    if s.empty:
        return None
    best = s.loc[s.Pvalue.idxmin()]
    wide = s[(s.Group.isin(["pLOF;pDM", "pLOF_pDM"])) & (s.max_MAF == 0.01)]
    w = wide.iloc[0] if len(wide) else best
    return {"mask": str(best.Group), "max_maf": float(best.max_MAF),
            "p_cell": float(best.Pvalue), "mac": int(best.MAC),
            "mac_case": int(best.MAC_case), "mac_control": int(best.MAC_control),
            "widest_mac": int(w.MAC)}


def main():
    for f in (CLINGEN, BRAVA4 / "omnibus_burden_broad_combined.tsv"):
        if not f.exists():
            sys.exit(f"REFUSING: {f} not found")
    OUT.mkdir(parents=True, exist_ok=True)

    cg = pd.read_csv(CLINGEN, sep="\t")
    genes = sorted(set(cg[cg.best_classification.isin(["Definitive", "Strong"])].gene))

    ours = ranked(PMBB_HL / "phase_4/results/omnibus_broad_combined.tsv")
    brava = ranked(BRAVA4 / "omnibus_burden_broad_combined.tsv")
    oc = pd.read_csv(PMBB_HL / "phase_4/results/cells_broad_combined.tsv", sep="\t")
    bc = pd.read_csv(BRAVA4 / "cells_burden_broad_combined.tsv", sep="\t")

    rows = []
    for g in genes:
        o = ours.loc[g] if g in ours.index else None
        b = brava.loc[g] if g in brava.index else None
        od, bd = cell_detail(oc, g), cell_detail(bc, g)
        rows.append({
            "gene": g,
            "rank_ours": int(o["rank"]) if o is not None else None,
            "p_ours": float(o.omnibus_p) if o is not None else None,
            "rank_brava": int(b["rank"]) if b is not None else None,
            "p_brava": float(b.omnibus_p) if b is not None else None,
            "in_top50_ours": bool(o is not None and o["rank"] <= DEPTH),
            "in_top50_brava": bool(b is not None and b["rank"] <= DEPTH),
            **{f"ours_{k}": v for k, v in (od or {}).items()},
            **{f"brava_{k}": v for k, v in (bd or {}).items()},
        })
    t = pd.DataFrame(rows)
    t["rank_change"] = t.rank_brava - t.rank_ours
    t = t.sort_values("rank_ours", na_position="last")
    t.to_csv(OUT / "02_gene_movement.tsv", sep="\t", index=False)

    moved = t[t.in_top50_ours | t.in_top50_brava].copy()
    print(f"genes in either top {DEPTH}: {len(moved)}\n")
    print(f"  {'gene':<10}{'ours':>8}{'p ours':>11}{'BRaVa':>8}{'p BRaVa':>11}"
          f"{'moved':>10}  what happened")
    for _, r in moved.sort_values("rank_ours").iterrows():
        both = r.in_top50_ours and r.in_top50_brava
        tag = ("held" if both else
               "LOST by BRaVa" if r.in_top50_ours else "GAINED by BRaVa")
        orank = f"{int(r.rank_ours):,}" if pd.notna(r.rank_ours) else "absent"
        brank = f"{int(r.rank_brava):,}" if pd.notna(r.rank_brava) else "absent"
        mv = f"{int(r.rank_change):+,}" if pd.notna(r.rank_change) else "—"
        print(f"  {r.gene:<10}{orank:>8}{r.p_ours:>11.2e}{brank:>8}"
              f"{r.p_brava:>11.2e}{mv:>10}  {tag}")

    print(f"\n  the allele counts underneath -- best cell in each annotation\n")
    print(f"  {'gene':<10}{'ours mask':>12}{'MAF':>8}{'case/ctrl':>12}"
          f"{'  ':>4}{'BRaVa mask':>12}{'MAF':>8}{'case/ctrl':>12}")
    for _, r in moved.sort_values("rank_ours").iterrows():
        print(f"  {r.gene:<10}{str(r.ours_mask):>12}{r.ours_max_maf:>8}"
              f"{str(int(r.ours_mac_case)) + '/' + str(int(r.ours_mac_control)):>12}"
              f"{'  ':>4}{str(r.brava_mask):>12}{r.brava_max_maf:>8}"
              f"{str(int(r.brava_mac_case)) + '/' + str(int(r.brava_mac_control)):>12}")

    summary = {
        "depth": DEPTH,
        "clingen_genes": len(genes),
        "in_top50_ours": sorted(t[t.in_top50_ours].gene),
        "in_top50_brava": sorted(t[t.in_top50_brava].gene),
        "held": sorted(t[t.in_top50_ours & t.in_top50_brava].gene),
        "lost_by_brava": sorted(t[t.in_top50_ours & ~t.in_top50_brava].gene),
        "gained_by_brava": sorted(t[~t.in_top50_ours & t.in_top50_brava].gene),
        "median_abs_rank_change": float(t.rank_change.abs().median()),
        "clingen_list_size": len(genes),
    }
    (OUT / "02_gene_movement.json").write_text(json.dumps(summary, indent=2) + "\n")
    both = int(t.rank_change.notna().sum())
    summary["ranked_by_both"] = both
    print(f"\n  median |rank change| across the {both} ClinGen genes ranked by both: "
          f"{summary['median_abs_rank_change']:,.0f} positions")
    print(f"\nwritten: results/02_gene_movement.tsv, 02_gene_movement.json")


if __name__ == "__main__":
    main()
