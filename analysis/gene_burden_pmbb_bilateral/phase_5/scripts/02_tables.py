#!/usr/bin/env python3
"""
Phase 5, step 2 -- the tables (her Step 11).

Four of them, each answering a question someone will ask about Phase 4:

  top_hits_<cohort>.tsv   what came closest, and whether it holds up
  dropped_from_top.tsv    genes high in her list that we never tested, and why
  biggest_movers.tsv      genes the corrections moved most, both directions
  summary.tsv             one row per cohort: counts, bars, best p

top_hits carries a `fragile` column, because a p-value alone ranks a 5-allele
gene alongside a 400-allele one. A test is flagged when the burden component
disagrees with the combined p (signal is all SKAT, no consistent direction) or
when the whole gene was collapsed into a single unit (Number_rare = 0), which
is where Burden and SKAT become the same test.

Input : phase_4/results/burden_all_cohorts.tsv, phase_2/data/symbol_biotype.tsv
Output: phase_5/results/*.tsv
"""
from pathlib import Path

import pandas as pd
from statsmodels.stats.multitest import multipletests

HERE = Path(__file__).resolve().parents[1]
REPL = HERE.parent
OUT = HERE / "results"
COHORTS = ["combined", "EUR", "AFR"]

d = pd.read_csv(HERE / "results/burden_all_cohorts.tsv", sep="\t")
# In this analysis one SAIGE call covers the grid, so the mask is the Group
# column rather than a separate folder name (premise P10).
d["Mask"] = d["Group"].replace({"pLOF;pDM": "pLOF_pDM"})
d["Pvalue"] = pd.to_numeric(d["Pvalue"])

bt = pd.read_csv(HERE.parent / "phase_2/data/symbol_biotype.tsv", sep="\t", header=None,
                 names=["symbol", "biotype"], dtype=str)
coding = set(bt[bt.biotype.str.contains("protein_coding", na=False)].symbol)
biotype = bt.groupby("symbol").biotype.first()

# Benjamini-Hochberg q-values, at both levels, because they answer different
# questions. Per panel: would this gene survive if that panel were the only
# analysis. Per cohort: does it survive having looked at all nine.
d["q_panel"] = d.groupby(["Cohort", "Mask", "max_MAF"]).Pvalue.transform(
    lambda x: multipletests(x.values, method="fdr_bh")[1])
d["q_cohort"] = d.groupby("Cohort").Pvalue.transform(
    lambda x: multipletests(x.values, method="fdr_bh")[1])

# ---- 1. top hits ----
for c in COHORTS:
    s = d[d.Cohort == c].copy()
    best = s.loc[s.groupby("Region").Pvalue.idxmin()].nsmallest(30, "Pvalue")
    best["bar_per_gene"] = 0.05 / s.Region.nunique()
    best["all_skat"] = best.Pvalue_Burden > 100 * best.Pvalue
    best["single_unit"] = best.Number_rare == 0
    best["fragile"] = best.all_skat | best.single_unit | (best.MAC < 20)
    best["why_fragile"] = [
        ", ".join(filter(None, [
            "burden disagrees" if r.all_skat else "",
            "collapsed to one unit" if r.single_unit else "",
            f"MAC={int(r.MAC)}" if r.MAC < 20 else "",
        ])) or "-"
        for r in best.itertuples()]
    best.insert(1, "rank", range(1, len(best) + 1))
    cols = ["rank", "Region", "Mask", "max_MAF", "Pvalue", "q_panel", "q_cohort",
            "Pvalue_Burden", "Pvalue_SKAT", "BETA_Burden", "MAC", "MAC_case",
            "MAC_control", "Number_rare", "Number_ultra_rare", "bar_per_gene",
            "fragile", "why_fragile"]
    best[cols].to_csv(OUT / f"top_hits_{c}.tsv", sep="\t", index=False)
    print(f"top_hits_{c}.tsv  — {int(best.fragile.sum())} of 30 flagged fragile")

# Sections comparing against another arm were dropped: this analysis has one.

# ---- 4. summary ----
rows = []
for c in COHORTS:
    s = d[d.Cohort == c]
    rows.append({
        "cohort": c, "tests": len(s), "genes": s.Region.nunique(),
        "min_pvalue": s.Pvalue.min(),
        "bar_per_test": 0.05 / len(s), "bar_per_gene": 0.05 / s.Region.nunique(),
        "any_significant": bool(s.Pvalue.min() < 0.05 / s.Region.nunique()),
        "min_q_cohort": s.q_cohort.min(), "min_q_panel": s.q_panel.min(),
        "any_fdr_significant": bool(s.q_cohort.min() < 0.05),
    })
pd.DataFrame(rows).to_csv(OUT / "summary.tsv", sep="\t", index=False)
print("summary.tsv")
print(pd.DataFrame(rows)[["cohort", "genes", "min_pvalue", "bar_per_gene",
                          "any_significant"]].to_string(index=False))
