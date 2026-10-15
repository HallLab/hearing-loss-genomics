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

d = pd.read_csv(REPL / "phase_4/results/burden_all_cohorts.tsv", sep="\t")
d["Pvalue"] = pd.to_numeric(d["Pvalue"])

bt = pd.read_csv(REPL / "phase_2/data/symbol_biotype.tsv", sep="\t", header=None,
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

# ---- 2 and 3 need her results ----
f = ("/project/hall/analysis/hearing-loss-genomics/analysis/elena/"
     "HL_only_rarevariant/newmasks_combined_gene_results.csv")
h = pd.read_csv(f, low_memory=False, keep_default_na=False, dtype=str)
h = h[h["Mask"].isin(["pLOF", "pDM", "pLOF_pDM"])]
h["Pvalue"] = pd.to_numeric(h["Pvalue"], errors="coerce")
h = h.dropna(subset=["Pvalue"])
h = h[~h["Region"].isin(["", "-"])].rename(columns={"Ancestry": "Cohort"})

dropped, movers = [], []
for c in COHORTS:
    o = d[d.Cohort == c].groupby("Region").Pvalue.min()
    hh = h[h.Cohort == c].groupby("Region").Pvalue.min()

    gone = hh[~hh.index.isin(o.index)].nsmallest(25)
    dropped.append(pd.DataFrame({
        "Cohort": c, "Region": gone.index, "her_pvalue": gone.values,
        "her_rank": [int((hh < p).sum()) + 1 for p in gone.values],
        "biotype": [biotype.get(g, "not in release VEP") for g in gone.index],
        "reason": "not protein-coding; removed by the Phase 2 mask rebuild",
    }))

    both = pd.DataFrame({"hers": hh, "ours": o}).dropna()
    import numpy as np
    both["shift_log10"] = -np.log10(both.ours) + np.log10(both.hers)
    big = pd.concat([both.nlargest(15, "shift_log10"), both.nsmallest(15, "shift_log10")])
    movers.append(big.assign(Cohort=c, direction=np.where(big.shift_log10 > 0,
                                                          "stronger in ours",
                                                          "weaker in ours"))
                     .reset_index().rename(columns={"index": "Region"}))

pd.concat(dropped).to_csv(OUT / "dropped_from_top.tsv", sep="\t", index=False)
pd.concat(movers).to_csv(OUT / "biggest_movers.tsv", sep="\t", index=False)
print("dropped_from_top.tsv, biggest_movers.tsv")

# ---- 4. summary ----
rows = []
for c in COHORTS:
    s = d[d.Cohort == c]
    hh = h[h.Cohort == c]
    rows.append({
        "cohort": c, "tests": len(s), "genes": s.Region.nunique(),
        "min_pvalue": s.Pvalue.min(),
        "bar_per_test": 0.05 / len(s), "bar_per_gene": 0.05 / s.Region.nunique(),
        "any_significant": bool(s.Pvalue.min() < 0.05 / s.Region.nunique()),
        "min_q_cohort": s.q_cohort.min(), "min_q_panel": s.q_panel.min(),
        "any_fdr_significant": bool(s.q_cohort.min() < 0.05),
        "her_genes": hh.Region.nunique(), "her_min_pvalue": hh.Pvalue.min(),
        "her_bar_per_gene": 0.05 / hh.Region.nunique(),
        "her_any_significant": bool(hh.Pvalue.min() < 0.05 / hh.Region.nunique()),
    })
pd.DataFrame(rows).to_csv(OUT / "summary.tsv", sep="\t", index=False)
print("summary.tsv")
print(pd.DataFrame(rows)[["cohort", "genes", "min_pvalue", "bar_per_gene",
                          "any_significant", "her_any_significant"]].to_string(index=False))
