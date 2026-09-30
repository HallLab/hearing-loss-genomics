#!/usr/bin/env python3
"""
Phase 1 / check 04 — close the phase.

Three things checks 01-03 left open.

  A. Step 5, Ancestry Stratification, had never been checked at all. Phase 1 is
     Step 1 + Step 5 by our own mapping, and only Step 1 had been done.

  B. Whether the 517 excluded participants are a random 0.7% of the cohort. The
     review page states the power loss is small; that holds only if they do not
     differ systematically from everyone else. Recorded as open since Phase 1
     began and never tested.

  C. Tinnitus. Step 1 produces two phenotypes; checks 01-03 covered only hearing
     impairment. This script does not close that gap -- it declares it, so the
     scope is explicit rather than silently narrow.

On (B): with 517 against ~70,000, a significance test flags differences far too
small to matter. Proportions and standardised differences are reported alongside
the p-values, and the reading should follow the effect size.

Writes only under phase_1/results/.
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent.parent
REL  = Path("/static/PMBB/PMBB-Release-2026-4.0")
E    = Path("/project/hall/analysis/hearing-loss-genomics/analysis/elena")
OUT  = HERE / "results"; OUT.mkdir(parents=True, exist_ok=True)

COV = REL / "Phenotype/4.0/PMBB-Release-2026-4.0_phenotype_covariates.txt"
ANC = REL / "Exome/PCA/combined/PMBB-Release-2026-4.0_genetic_exome.commonsnps.samples_ancestries.tsv"
SAIGE = E / "HL_only_rarevariant/new_SAIGE_covariates"

report = {}

# ---------------------------------------------------------------- A. Step 5
anc = pd.read_csv(ANC, sep="\t", usecols=["IID", "Class"])
amap = dict(zip(anc.IID, anc.Class))

cohorts = {}
for k in ["combined", "EUR", "AFR"]:
    d = pd.read_csv(SAIGE / f"hearing_impairment_{k}_SAIGE_Step1_covariates.txt",
                    sep="\t", usecols=["IID", "PHENO"])
    cohorts[k] = d

step5 = {"cohorts": {}, "purity": {}, "combined_composition": {}}
for k, d in cohorts.items():
    step5["cohorts"][k] = {"N": len(d), "cases": int(d.PHENO.sum()),
                           "controls": int((d.PHENO == 0).sum()),
                           "case_rate_pct": round(100 * d.PHENO.mean(), 2)}
for k in ["EUR", "AFR"]:
    cl = pd.Series([amap.get(i, "?") for i in cohorts[k].IID]).value_counts()
    step5["purity"][k] = {"expected_class": k, "all_match": bool(set(cl.index) == {k}),
                          "classes_present": cl.to_dict()}
step5["combined_composition"] = pd.Series(
    [amap.get(i, "?") for i in cohorts["combined"].IID]).value_counts().to_dict()
step5["EUR_plus_AFR"] = len(cohorts["EUR"]) + len(cohorts["AFR"])
step5["combined"] = len(cohorts["combined"])
step5["carried_only_in_combined"] = step5["combined"] - step5["EUR_plus_AFR"]
# --- the same stratification on the corrected arm (57,507), which nobody has run.
# The consumed cohort above is the reproduction arm; this is what the strata would be
# if the 431 came back and the 556 were excluded as the rules say.
stat = pd.read_csv(E / "rarevariantExWAS/PMBBv4_phecodex/PMBBv4_hearing_tinnitus_all_statuses.csv.gz")
stat = stat[stat.phenotype == "hearing_impairment"].set_index("person_id")["status"]
corrected = stat[stat.isin(["case", "control"])]

cdf = pd.DataFrame({"IID": corrected.index,
                    "PHENO": (corrected.values == "case").astype(int)})
cdf["ancestry"] = cdf.IID.map(amap)

arm = {"combined": {"N": len(cdf), "cases": int(cdf.PHENO.sum()),
                    "controls": int((cdf.PHENO == 0).sum()),
                    "case_rate_pct": round(100 * cdf.PHENO.mean(), 2)}}
for k in ["EUR", "AFR"]:
    d = cdf[cdf.ancestry == k]
    arm[k] = {"N": len(d), "cases": int(d.PHENO.sum()),
              "controls": int((d.PHENO == 0).sum()),
              "case_rate_pct": round(100 * d.PHENO.mean(), 2)}
arm["composition"] = cdf.ancestry.value_counts().to_dict()
step5["corrected_arm_57507"] = arm

# what changes between the arms, per stratum
delta = {}
for k in ["combined", "EUR", "AFR"]:
    r, c = step5["cohorts"][k], arm[k]
    delta[k] = {"N": c["N"] - r["N"], "cases": c["cases"] - r["cases"],
                "controls": c["controls"] - r["controls"],
                "case_rate_pct_points": round(c["case_rate_pct"] - r["case_rate_pct"], 2)}
step5["corrected_minus_reproduction"] = delta
report["A_step5_ancestry_stratification"] = step5

# ---------------------------------------------------------------- B. the 517
cov = pd.read_csv(COV, sep="\t", low_memory=False)
imp = [c for c in cov.columns if c.startswith("imputed_PC")]
cov["dropped"] = ~cov[imp].notna().all(axis=1)
cov["ancestry"] = cov.person_id.map(amap)
cov["year"] = pd.to_datetime(cov.sample_date, errors="coerce").dt.year

D, K = cov[cov.dropped], cov[~cov.dropped]
prof = {"n_dropped": int(len(D)), "n_kept": int(len(K))}

def smd(a, b):
    """standardised mean difference — effect size, not significance"""
    a, b = a.dropna(), b.dropna()
    if len(a) < 2 or len(b) < 2: return None
    sd = np.sqrt((a.var() + b.var()) / 2)
    return None if sd == 0 else round(float((a.mean() - b.mean()) / sd), 3)

for col, label in [("sample_age", "age"), ("year", "enrolment_year")]:
    a, b = D[col], K[col]
    u = stats.mannwhitneyu(a.dropna(), b.dropna(), alternative="two-sided")
    prof[label] = {"dropped_median": float(a.median()), "kept_median": float(b.median()),
                   "dropped_mean": round(float(a.mean()), 2), "kept_mean": round(float(b.mean()), 2),
                   "standardised_difference": smd(a, b), "mannwhitney_p": float(u.pvalue)}

for col in ["ancestry", "sequenced_gender", "batch"]:
    dv = D[col].value_counts(normalize=True).mul(100).round(2)
    kv = K[col].value_counts(normalize=True).mul(100).round(2)
    tab = pd.crosstab(cov[col], cov.dropped)
    chi = stats.chi2_contingency(tab)
    keys = sorted(set(dv.index) | set(kv.index), key=str)
    prof[col] = {"pct_dropped": {str(k): float(dv.get(k, 0.0)) for k in keys},
                 "pct_kept":    {str(k): float(kv.get(k, 0.0)) for k in keys},
                 "chi2_p": float(chi.pvalue),
                 "largest_pct_point_gap": round(float(max(abs(dv.get(k, 0.0) - kv.get(k, 0.0)) for k in keys)), 2)}
# exclusion rate WITHIN each group -- the interpretable framing
rates = {}
for col in ["ancestry", "sequenced_gender", "batch"]:
    g = cov.groupby(col)["dropped"].agg(["sum", "count"])
    g["pct_excluded"] = (100 * g["sum"] / g["count"]).round(2)
    rates[col] = {str(k): {"excluded": int(r["sum"]), "total": int(r["count"]),
                           "pct_excluded": float(r["pct_excluded"])}
                  for k, r in g.iterrows()}
prof["exclusion_rate_within_group"] = rates
prof["cohort_wide_exclusion_rate_pct"] = round(100 * len(D) / len(cov), 2)
report["B_profile_of_the_excluded"] = prof

# ---------------------------------------------------------------- C. scope
report["C_scope_not_covered"] = {
    "tinnitus": "Step 1 produces hearing_impairment and tinnitus. Checks 01-04 cover "
                "hearing impairment only. The tinnitus phenotype, and the combined "
                "HL-and/or-tinnitus phenotype decided on 2026-07-01, are unverified.",
}

(OUT / "04_step5_and_dropped_profile.json").write_text(json.dumps(report, indent=2, default=str))
print(json.dumps(report, indent=2, default=str))
