#!/usr/bin/env python3
"""
Phase 4, check 05 -- is "nothing significant" a real absence or a broken analysis?

"No gene passes the threshold" and "there is nothing here" are different claims.
A deflated analysis produces the first without supporting the second. Two checks:

  A. Calibration. Under the null, the fraction of tests below p is p. A ratio
     far below 1 means the test is conservative and would miss real signal; far
     above 1 means inflation and the p-values cannot be trusted at all. Measured
     at one MAF cutoff so the same gene is not counted three times.

  B. Cross-ancestry consistency. EUR and AFR are independent samples. If a real
     effect existed, genes near the top in one should be near the top in the
     other more often than chance allows.

  C. How many tests are actually distinct. The three max-MAF cutoffs are nested,
     and in exome data a pLOF or damaging-missense variant is nearly always very
     rare, so the cutoffs often select the SAME variant set and repeat the same
     test. This decides whether 0.05/n_tests is the honest bar or a conservative
     one.

  D. The Cauchy (ACAT) omnibus, one p-value per gene. Raised by Nikki Palmiero:
     SAIGE emits a Cauchy row combining annotation groups when several are
     requested, and ours have none, because we pass one annotation per mask so
     there is nothing for it to combine. The omnibus is computed here instead,
     over the 9 cells per gene. It is the principled answer to the question C
     raises -- one test per gene, with no denominator left to argue about.

  E. Where the MAF cutoff does change the answer. A gene whose signal sits only
     in the very rarest variants, and dilutes as slightly commoner ones enter,
     is behaving the way a real gene would.

A WARNING ABOUT CHECK B, because the first version of it was wrong.
Per gene we take the minimum p over 9 correlated tests (3 masks x 3 MAF). The
marginal rate of "min-p < t" per gene is therefore NOT t -- it is roughly 3x
larger. Using t for the expectation produced an apparent 12-fold excess that
does not exist. The expectation must use the MEASURED marginal rates.

Input : phase_4/results/burden_all_cohorts.tsv
Output: phase_4/results/05_calibration_and_consistency.json
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest

OUT = Path(__file__).resolve().parents[1] / "results"
d = pd.read_csv(OUT / "burden_all_cohorts.tsv", sep="\t")
d["Pvalue"] = pd.to_numeric(d["Pvalue"])

report = {"calibration": {}, "cross_ancestry": {}}

# ---- A. calibration ----
for c in ["combined", "EUR", "AFR"]:
    s = d[(d.Cohort == c) & (d.max_MAF == 0.01)]
    report["calibration"][c] = {
        "n_tests": int(len(s)),
        "ratios": {str(t): round(float((s.Pvalue < t).mean() / t), 3)
                   for t in (0.05, 0.01, 0.001)},
    }

# ---- B. cross-ancestry consistency ----
eur = d[d.Cohort == "EUR"].groupby("Region").Pvalue.min()
afr = d[d.Cohort == "AFR"].groupby("Region").Pvalue.min()
both = pd.DataFrame({"EUR": eur, "AFR": afr}).dropna()
n = len(both)
report["cross_ancestry"]["genes_in_both"] = int(n)

for t in (0.01, 0.005, 0.001):
    rate_e = float((both.EUR < t).mean())
    rate_a = float((both.AFR < t).mean())
    obs = int(((both.EUR < t) & (both.AFR < t)).sum())
    exp = n * rate_e * rate_a
    report["cross_ancestry"][f"min_p_lt_{t}"] = {
        "marginal_rate_EUR": round(rate_e, 4),
        "marginal_rate_AFR": round(rate_a, 4),
        "note": "marginal rate, not t -- min over 9 correlated tests inflates it",
        "observed": obs,
        "expected": round(exp, 1),
        "ratio": round(obs / exp, 2) if exp else None,
        "binomial_p_greater": float(binomtest(obs, n, rate_e * rate_a,
                                              alternative="greater").pvalue),
    }

# ---- C. how many tests are actually distinct ----
piv = d.pivot_table(index=["Cohort", "Mask", "Region"], columns="max_MAF",
                    values="Pvalue").dropna()
same_all = (piv[0.0001] == piv[0.001]) & (piv[0.001] == piv[0.01])
same_two = (piv[0.001] == piv[0.01]) & ~same_all
report["maf_nesting"] = {
    "triples_with_all_three_cutoffs": int(len(piv)),
    "pct_all_three_identical": round(100 * float(same_all.mean()), 1),
    "pct_only_two_identical": round(100 * float(same_two.mean()), 1),
    "pct_all_three_distinct": round(100 * float((~same_all & ~same_two).mean()), 1),
    "mean_variants_in_set": {
        str(m): round(float((d[d.max_MAF == m].Number_rare
                             + d[d.max_MAF == m].Number_ultra_rare).mean()), 2)
        for m in (0.0001, 0.001, 0.01)
    },
    "bonferroni": {},
}
for c in ["combined", "EUR", "AFR"]:
    s_ = d[d.Cohort == c]
    nominal = len(s_)
    distinct = int(s_.groupby(["Mask", "Region"]).Pvalue.nunique().sum())
    genes = int(s_.Region.nunique())
    report["maf_nesting"]["bonferroni"][c] = {
        "n_nominal": int(nominal), "bar_nominal": 0.05 / nominal,
        "n_distinct": distinct, "bar_distinct": 0.05 / distinct,
        "n_genes": genes, "bar_per_gene": 0.05 / genes,
        "min_pvalue": float(s_.Pvalue.min()),
        "passes_most_permissive_bar": bool(s_.Pvalue.min() < 0.05 / genes),
    }

# ---- D. Cauchy / ACAT omnibus, one p per gene ----
def acat(p):
    p = np.clip(np.asarray(p, dtype=float), 1e-300, 1 - 1e-16)
    return 0.5 - np.arctan(np.mean(np.tan((0.5 - p) * np.pi))) / np.pi


report["cauchy_omnibus"] = {
    "why": "we emit no Cauchy rows (one annotation per mask), so this is computed "
           "post hoc over the 9 cells per gene",
    "by_cohort": {},
}
for c in ["combined", "EUR", "AFR"]:
    s_ = d[d.Cohort == c]
    om = s_.groupby("Region").Pvalue.apply(lambda x: acat(x.values))
    bar = 0.05 / len(om)
    report["cauchy_omnibus"]["by_cohort"][c] = {
        "genes": int(len(om)),
        "min_omnibus_p": float(om.min()),
        "min_single_cell_p": float(s_.Pvalue.min()),
        "bar_per_gene": bar,
        "any_significant": bool(om.min() < bar),
        "top5": {g: float(v) for g, v in om.nsmallest(5).items()},
    }

# ---- E. genes where the cutoff changes the answer ----
r = piv[piv[0.01] > 0]
ratio = r[0.01] / r[0.0001]
sens = r[ratio > 100].assign(fold=ratio[ratio > 100]).sort_values(0.0001)
report["maf_sensitive_genes"] = {
    "n_with_100x_improvement_at_strictest": int(len(sens)),
    "top": [
        {"cohort": i[0], "mask": i[1], "gene": i[2],
         "p_at_1e-4": float(row[0.0001]), "p_at_1e-2": float(row[0.01]),
         "fold": round(float(row["fold"]), 1)}
        for i, row in sens.head(10).iterrows()
    ],
}

(OUT / "05_calibration_and_consistency.json").write_text(json.dumps(report, indent=2))

print("A. calibration (observed/expected; 1.0 = perfectly calibrated)")
for c, v in report["calibration"].items():
    print(f"   {c:9s} " + "  ".join(f"p<{t} {r}" for t, r in v["ratios"].items())
          + f"   (n={v['n_tests']:,})")
print("\nB. cross-ancestry, EUR and AFR independent")
for t in (0.01, 0.005, 0.001):
    v = report["cross_ancestry"][f"min_p_lt_{t}"]
    print(f"   min-p<{t:<6} observed {v['observed']:>3}  expected {v['expected']:>5}  "
          f"ratio {v['ratio']}  p={v['binomial_p_greater']:.3g}")

m = report["maf_nesting"]
print(f"\nC. the three MAF cutoffs, over {m['triples_with_all_three_cutoffs']:,} "
      f"(cohort, mask, gene) triples")
print(f"   all three identical {m['pct_all_three_identical']}%  "
      f"two identical {m['pct_only_two_identical']}%  "
      f"all distinct {m['pct_all_three_distinct']}%")
for c, v in m["bonferroni"].items():
    print(f"   {c:9s} min p {v['min_pvalue']:.3g} | bar nominal {v['bar_nominal']:.2e} "
          f"| distinct {v['bar_distinct']:.2e} | per gene {v['bar_per_gene']:.2e} "
          f"| passes most permissive: {v['passes_most_permissive_bar']}")

print("\nD. Cauchy omnibus, one p per gene over the 9 cells")
for c, v in report["cauchy_omnibus"]["by_cohort"].items():
    print(f"   {c:9s} omnibus {v['min_omnibus_p']:.3g}  vs best single cell "
          f"{v['min_single_cell_p']:.3g}  bar {v['bar_per_gene']:.2e}  "
          f"significant: {v['any_significant']}")

print(f"\nE. {report['maf_sensitive_genes']['n_with_100x_improvement_at_strictest']} "
      f"tests improve >100x at the strictest cutoff; top 5:")
for e in report["maf_sensitive_genes"]["top"][:5]:
    print(f"   {e['cohort']:9s} {e['mask']:9s} {e['gene']:10s} "
          f"{e['p_at_1e-4']:.3g} -> {e['p_at_1e-2']:.3g}  ({e['fold']:.0f}x)")
