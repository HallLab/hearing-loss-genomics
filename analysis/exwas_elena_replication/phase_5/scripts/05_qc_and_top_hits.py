#!/usr/bin/env python3
"""
Phase 5, step 5 -- a QC self-check, and the top variants with annotation.

THE SELF-CHECK exists because of what the poster audit found. Having shown that the
poster's hits sit on variants with 20% call rate, this arm owes the same test of
itself: among the variants WE tested, is a lower call rate associated with a smaller
p-value? If it were, our own nulls would be suspect in the same way.

Measured three ways, because a single correlation can hide a tail effect:
  - Spearman rho between call rate and -log10(p) across everything tested
  - the same, restricted to the worst-call-rate decile
  - the call rate of the top 100 variants against the rest

THE TOP VARIANTS are annotated by joining her delivered tables on MarkerID, which
carry VEP SYMBOL, Consequence, IMPACT, REVEL and CADD. Premise P9 would have run VEP
on anything clearing a bar; nothing did, so this is the weaker version -- the best
signals described, not hits explained.
"""
import csv
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

BASE = Path(__file__).resolve().parents[2]
ELENA = BASE.parent / "elena/HL_only_ExWAS"
OUT = BASE / "phase_5/results"
COHORTS = ["combined", "EUR", "AFR"]
TOP_N = 100


def ours(cohort):
    d = pd.read_csv(BASE / f"phase_4/results/merged_{cohort}.tsv", sep="\t")
    d = d[(d["p.value"] > 0) & (d["p.value"] <= 1)].copy()
    d["call_rate"] = 1 - d["MissingRate"]
    d["neglog10p"] = -np.log10(d["p.value"])
    return d


def annotation(cohort):
    cols = ["MarkerID", "SYMBOL", "Consequence", "IMPACT", "REVEL_score", "CADD_PHRED"]
    d = pd.read_csv(ELENA / f"HL_ExWAS_FullyAnnotated_AllVariants_{cohort}.csv",
                    usecols=cols, low_memory=False)
    return d.drop_duplicates("MarkerID")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    report = {"self_check": {}, "top": {}}

    print("QC SELF-CHECK -- does low call rate buy significance in OUR tested set?\n")
    for c in COHORTS:
        d = ours(c)
        rho_all = spearmanr(d["call_rate"], d["neglog10p"])
        cut = d["call_rate"].quantile(0.10)
        worst = d[d["call_rate"] <= cut]
        rho_worst = spearmanr(worst["call_rate"], worst["neglog10p"])

        top = d.nlargest(TOP_N, "neglog10p")
        rest = d.drop(top.index)
        u = mannwhitneyu(top["call_rate"], rest["call_rate"], alternative="two-sided")

        report["self_check"][c] = {
            "spearman_callrate_vs_neglog10p": float(rho_all.statistic),
            "p": float(rho_all.pvalue),
            "worst_decile_cut": float(cut),
            "spearman_worst_decile": float(rho_worst.statistic),
            "p_worst_decile": float(rho_worst.pvalue),
            "median_callrate_top100": float(top["call_rate"].median()),
            "median_callrate_rest": float(rest["call_rate"].median()),
            "mannwhitney_p": float(u.pvalue),
        }
        print(f"  {c}")
        print(f"    rho(call rate, -log10 p), all    {rho_all.statistic:+.4f}  "
              f"p={rho_all.pvalue:.2g}")
        print(f"    rho, worst decile (<= {cut:.4f})   "
              f"{rho_worst.statistic:+.4f}  p={rho_worst.pvalue:.2g}")
        print(f"    median call rate  top {TOP_N}: {top['call_rate'].median():.5f}   "
              f"rest: {rest['call_rate'].median():.5f}   "
              f"Mann-Whitney p={u.pvalue:.2g}")
        print()

    print("TOP VARIANTS, annotated from her delivered tables\n")
    for c in COHORTS:
        d = ours(c).nlargest(15, "neglog10p")
        ann = annotation(c)
        m = d.merge(ann, on="MarkerID", how="left")
        rows = []
        for _, r in m.iterrows():
            rows.append({
                "variant": r["MarkerID"], "p": float(r["p.value"]),
                "beta": float(r["BETA"]), "af": float(r["AF_Allele2"]),
                "call_rate": float(r["call_rate"]),
                "gene": r.get("SYMBOL") if pd.notna(r.get("SYMBOL")) else "",
                "consequence": (r.get("Consequence")
                                if pd.notna(r.get("Consequence")) else ""),
                "impact": r.get("IMPACT") if pd.notna(r.get("IMPACT")) else "",
                "cadd": (float(r["CADD_PHRED"]) if pd.notna(r.get("CADD_PHRED"))
                         else None),
            })
        report["top"][c] = rows
        print(f"  === {c} ===")
        print(f"  {'variant':<22}{'p':>10}{'OR':>7}{'AF':>7}{'call':>8}  "
              f"{'gene':<12}{'consequence':<26}{'CADD':>6}")
        for r in rows[:8]:
            cadd = f"{r['cadd']:.1f}" if r["cadd"] is not None else "--"
            print(f"  {r['variant']:<22}{r['p']:>10.2e}"
                  f"{np.exp(r['beta']):>7.2f}{r['af']:>7.3f}{r['call_rate']:>8.4f}  "
                  f"{r['gene'][:11]:<12}{r['consequence'][:25]:<26}{cadd:>6}")
        print()

    dest = OUT / "qc_and_top_hits.json"
    dest.write_text(json.dumps(report, indent=2) + "\n")
    print(f"written: {dest.relative_to(BASE)}")


if __name__ == "__main__":
    main()
