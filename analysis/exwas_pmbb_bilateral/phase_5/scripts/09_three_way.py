#!/usr/bin/env python3
"""
The three-way comparison -- Elena, Replication, Bilateral.

Each arm changes ONE thing from the one before it, which is what makes the
comparison readable:

    Elena         SO_396 broad phenotype   ·  her masks        ·  5/9/10 PCs
    Replication   SO_396 broad phenotype   ·  rebuilt masks    ·  5/4/3 PCs
    Bilateral     bilateral sensorineural  ·  rebuilt masks    ·  5/4/3 PCs
                  └─ the only change ──┘      └─ the only change from Elena ─┘

So Elena -> Replication isolates the mask and PC corrections, and
Replication -> Bilateral isolates the phenotype decision.

COMPARABILITY. Elena's pipeline emits no per-gene omnibus -- its Cauchy rows
combine annotations inside one cell, leaving ~9 per gene. Hers is therefore
computed here with the same ACAT over the same nine cells, so all three arms are
the same quantity. Her ALL mask is excluded, since the other two never ran it.

Her table also carries two defects documented in the replication: 5,391 rows
whose columns are shifted, from three files that lost their header, and a '-'
pseudo-gene pooling every variant VEP could not name. Both are dropped loudly
rather than averaged in.

Output: phase_5/results/09_three_way.json
        figures/three_way_*.png
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hypergeom

HERE = Path(__file__).resolve().parents[2]
REPL = HERE.parent / "elena_replication"
ELENA = (HERE.parent / "elena" / "HL_only_rarevariant" /
         "newmasks_combined_gene_results.csv")
MASKS = ["pLOF", "pDM", "pLOF_pDM"]
COHORTS = ["combined", "EUR", "AFR"]

ARMS = ["Elena", "Replication", "Bilateral"]
DESIGN = {
    "Elena":       dict(phenotype="SO_396 (broad)", masks="original",
                        pcs="5 / 9 / 10", cases={"combined": 6712}),
    "Replication": dict(phenotype="SO_396 (broad)", masks="rebuilt",
                        pcs="5 / 4 / 3", cases={"combined": 6752}),
    "Bilateral":   dict(phenotype="bilateral sensorineural", masks="rebuilt",
                        pcs="5 / 4 / 3", cases={"combined": 3164}),
}


def acat(p):
    p = np.clip(np.asarray(p, dtype=float), 1e-300, 1 - 1e-16)
    return float(0.5 - np.arctan(np.mean(np.tan((0.5 - p) * np.pi))) / np.pi)


def elena_omnibus():
    """Her cells, collapsed to one p per gene with the same statistic."""
    h = pd.read_csv(ELENA, low_memory=False, keep_default_na=False, dtype=str)
    before = len(h)
    h = h[h["Mask"].isin(MASKS)]
    h = h[h["Group"] != "Cauchy"]          # her per-cell Cauchy, not the grid
    h["Pvalue"] = pd.to_numeric(h["Pvalue"], errors="coerce")
    h = h.dropna(subset=["Pvalue"])
    h = h[~h["Region"].isin(["", "-"])]
    out = {}
    for c in COHORTS:
        s = h[h.Ancestry == c]
        out[c] = (s.groupby("Region").Pvalue.apply(lambda x: acat(x.values))
                   .rename("omnibus_p").reset_index())
    return out, before - len(h)


def read_omnibus(root):
    return {c: pd.read_csv(root / f"omnibus_{c}.tsv", sep="\t")[
        ["Region", "omnibus_p"]] for c in COHORTS}


el, dropped = elena_omnibus()
arms = {"Elena": el,
        "Replication": read_omnibus(REPL / "phase_4/results"),
        "Bilateral": read_omnibus(HERE / "phase_5/results")}

cg = pd.read_csv(HERE / "phase_5/data/clingen_hl_genes.tsv", sep="\t")
GENES = set(cg[cg.best_classification.isin(["Definitive", "Strong"])].gene)

report = {"design": DESIGN, "her_rows_dropped": int(dropped), "cohorts": {}}

for c in COHORTS:
    entry = {"arms": {}, "top50_overlap": {}}
    tops = {}
    for a in ARMS:
        d = arms[a][c].dropna().sort_values("omnibus_p").reset_index(drop=True)
        d["rank"] = np.arange(1, len(d) + 1)
        inset = d.Region.isin(GENES)
        n_in = int(inset.sum())
        hit = int((d["rank"][inset] <= 50).sum())
        exp = 50 * n_in / len(d)
        entry["arms"][a] = {
            "genes": int(len(d)),
            "min_omnibus_p": float(d.omnibus_p.min()),
            "bar": 0.05 / len(d),
            "significant": bool(d.omnibus_p.min() < 0.05 / len(d)),
            "clingen_tested": n_in,
            "clingen_in_top50": hit,
            "clingen_expected": round(exp, 2),
            "clingen_p": float(hypergeom.sf(hit - 1, len(d), n_in, 50)),
            "clingen_genes": sorted(d.Region[inset & (d["rank"] <= 50)]),
            "top5": d.head(5).Region.tolist(),
        }
        tops[a] = set(d.head(50).Region)
    for i, a in enumerate(ARMS):
        for b in ARMS[i + 1:]:
            entry["top50_overlap"][f"{a} vs {b}"] = len(tops[a] & tops[b])
    report["cohorts"][c] = entry

(HERE / "phase_5/results/09_three_way.json").write_text(json.dumps(report, indent=2))

print(f"Elena's table: {dropped:,} rows dropped (shifted columns, '-' pseudo-gene)\n")
for c in COHORTS:
    e = report["cohorts"][c]
    print(f"=== {c} ===")
    print(f"{'arm':13s} {'genes':>7s} {'min omnibus':>12s} {'sig':>4s}  "
          f"{'ClinGen top50':>13s}  {'p':>9s}")
    for a in ARMS:
        v = e["arms"][a]
        print(f"{a:13s} {v['genes']:>7,} {v['min_omnibus_p']:>12.3g} "
              f"{str(v['significant']):>4}  {v['clingen_in_top50']:>6}/"
              f"{v['clingen_expected']:<6} {v['clingen_p']:>9.2g}")
    print("  top-50 shared: " + " · ".join(f"{k} {v}" for k, v in
                                            e["top50_overlap"].items()))
    print()
