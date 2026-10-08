#!/usr/bin/env python3
"""
Top-50 gene tables for each arm, with ClinGen membership marked.

The reference table the Confluence conclusions page ends on. Three arms side by
side so a reader can scan what each ranking puts at the top and see, without
cross-referencing anything, which of those are genes ClinGen has adjudicated as
causing hearing loss.

p-values are the per-gene omnibus in every arm, including Elena's, which is
recomputed with the same ACAT -- see 09_three_way.py.

Output: phase_5/results/top50_<cohort>.tsv   (all three cohorts)
        and the combined cohort printed as markdown for the page
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[2]
REPL = HERE.parent / "elena_replication"
ELENA = (HERE.parent / "elena" / "HL_only_rarevariant" /
         "newmasks_combined_gene_results.csv")
MASKS = ["pLOF", "pDM", "pLOF_pDM"]
COHORTS = ["combined", "EUR", "AFR"]
ARMS = ["Elena", "Replication", "Bilateral"]

cg = pd.read_csv(HERE / "phase_5/data/clingen_hl_genes.tsv", sep="\t")
CLINGEN = set(cg[cg.best_classification.isin(["Definitive", "Strong"])].gene)
CG_ANY = dict(zip(cg.gene, cg.best_classification))


def acat(p):
    p = np.clip(np.asarray(p, dtype=float), 1e-300, 1 - 1e-16)
    return float(0.5 - np.arctan(np.mean(np.tan((0.5 - p) * np.pi))) / np.pi)


h = pd.read_csv(ELENA, low_memory=False, keep_default_na=False, dtype=str)
h = h[h["Mask"].isin(MASKS) & (h["Group"] != "Cauchy")]
h["Pvalue"] = pd.to_numeric(h["Pvalue"], errors="coerce")
h = h.dropna(subset=["Pvalue"])
h = h[~h["Region"].isin(["", "-"])]


def arm_table(arm, cohort):
    if arm == "Elena":
        s = h[h.Ancestry == cohort]
        d = (s.groupby("Region").Pvalue.apply(lambda x: acat(x.values))
              .rename("p").reset_index())
    else:
        root = (REPL / "phase_4/results" if arm == "Replication"
                else HERE / "phase_5/results")
        d = pd.read_csv(root / f"omnibus_{cohort}.tsv", sep="\t")[
            ["Region", "omnibus_p"]].rename(columns={"omnibus_p": "p"})
    d = d.dropna().sort_values("p").head(50).reset_index(drop=True)
    d.insert(0, "rank", np.arange(1, len(d) + 1))
    d["clingen"] = d.Region.map(lambda g: CG_ANY.get(g, ""))
    d["arm"] = arm
    return d


for cohort in COHORTS:
    out = pd.concat([arm_table(a, cohort) for a in ARMS], ignore_index=True)
    out.to_csv(HERE / f"phase_5/results/top50_{cohort}.tsv", sep="\t", index=False)

# markdown for the page -- combined cohort, three arms side by side
t = {a: arm_table(a, "combined") for a in ARMS}
print("| # | Elena | p | Replication | p | Bilateral | p |")
print("|---:|---|---|---|---|---|---|")
for i in range(50):
    cells = []
    for a in ARMS:
        r = t[a].iloc[i]
        star = " **‡**" if r.Region in CLINGEN else (" †" if r.clingen else "")
        cells += [f"`{r.Region}`{star}", f"{r.p:.2g}"]
    print(f"| {i+1} | " + " | ".join(cells) + " |")

print("\n--- counts ---")
for a in ARMS:
    ds = sum(1 for g in t[a].Region if g in CLINGEN)
    any_ = sum(1 for g in t[a].Region if g in CG_ANY)
    print(f"  {a:12s} {ds} Definitive/Strong, {any_} at any ClinGen classification")
