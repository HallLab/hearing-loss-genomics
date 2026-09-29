#!/usr/bin/env python3
"""
Phase 1 / check 02 — rebuild cases and controls from the release, independently.

Re-derives the hearing-impairment case/control assignment from
/static/PMBB/PMBB-Release-2026-4.0/ and compares person-by-person against the
assignment the pipeline consumed.

Two target definitions are computed, because the choice is not self-evident:

  exact   SO_396 only
  rollup  SO_396 and its children (SO_396.1, .2, .3, .8, .9)

PheWAS convention rolls a child phecode up into its parent. Whether this pipeline
did so is exactly what this check establishes.

Rules, as documented by the pipeline's own status categories:
  case                            evidence on >= 2 distinct dates
  excluded_one_date               evidence on exactly 1 date
  excluded_related_ear_phenotype  no target evidence, but other SO_39x evidence
  control                         no ear-family evidence at all

Writes only under phase_1/results/.
"""
import json
from pathlib import Path
import pandas as pd

HERE = Path(__file__).resolve().parent.parent
REL  = Path("/static/PMBB/PMBB-Release-2026-4.0")
PIPE = Path("/project/hall/analysis/hearing-loss-genomics/analysis/elena/rarevariantExWAS")
OUT  = HERE / "results"; OUT.mkdir(parents=True, exist_ok=True)

EAR = HERE / "data/_ear_family.tsv"          # cached extract, regenerable
ANC = REL / "Exome/PCA/combined/PMBB-Release-2026-4.0_genetic_exome.commonsnps.samples_ancestries.tsv"

OBS = HERE / "data/_tinnitus_obs.tsv"        # cached extract, regenerable

ear = pd.read_csv(EAR, sep="\t", names=["person_id", "date", "code"], dtype=str)

# PMBB v4 relocated the standard tinnitus ICD codes (388.3x, H93.1x) into the OMOP
# `observation` table. The phecode files keep only the ~3,016 pulsatile events, so
# ear-family evidence built from conditions_phecode_x alone is incomplete. Established
# here by check 02: 558 people the pipeline flags as having related ear evidence look
# like controls without this source, and 100% of them carry tinnitus in `observation`.
obs = pd.read_csv(OBS, sep="\t", names=["person_id", "date", "code"], dtype=str)
obs["code"] = "OBS_TINNITUS"
ear = pd.concat([ear, obs], ignore_index=True)

ear = ear[ear.date.notna() & (ear.date != "")]
cohort = set(pd.read_csv(ANC, sep="\t", usecols=["IID"])["IID"])
ear = ear[ear.person_id.isin(cohort)]

def classify(is_target):
    """is_target: boolean Series over `ear` rows."""
    tgt   = ear[is_target].groupby("person_id")["date"].nunique()
    other = ear[~is_target].groupby("person_id")["date"].nunique()
    st = pd.Series("control", index=sorted(cohort), dtype=object)
    st.loc[st.index.intersection(other.index)] = "excluded_related_ear_phenotype"
    one = tgt[tgt == 1].index
    st.loc[st.index.intersection(one)] = "excluded_one_date"
    two = tgt[tgt >= 2].index
    st.loc[st.index.intersection(two)] = "case"
    return st

defs = {
    "exact":  ear.code == "SO_396",
    "rollup": ear.code.str.match(r"^SO_396($|\.)"),
}
built = {k: classify(v) for k, v in defs.items()}

# ------------------------------------------------- the pipeline's own assignment
ref = pd.read_csv(PIPE / "PMBBv4_phecodex/PMBBv4_hearing_tinnitus_all_statuses.csv.gz")
ref = ref[ref.phenotype == "hearing_impairment"].set_index("person_id")["status"]

report = {"cohort": len(cohort), "reference_counts": ref.value_counts().to_dict()}
for name, st in built.items():
    aligned = st.reindex(ref.index)
    agree = (aligned == ref)
    report[name] = {
        "counts":      st.value_counts().to_dict(),
        "agree":       int(agree.sum()),
        "disagree":    int((~agree).sum()),
        "pct_agree":   round(100 * agree.mean(), 4),
        "confusion":   pd.crosstab(ref, aligned).to_dict(),
    }
    pd.DataFrame({"person_id": ref.index[~agree],
                  "pipeline":  ref[~agree].values,
                  "rebuilt":   aligned[~agree].values}) \
      .to_csv(OUT / f"02_disagreements_{name}.csv", index=False)

(OUT / "02_rebuild_cases_controls.json").write_text(json.dumps(report, indent=2, default=str))

print(f"cohort {len(cohort):,}")
print("pipeline :", {k: int(v) for k, v in report['reference_counts'].items()})
for name in defs:
    r = report[name]
    print(f"\n{name:7s} :", {k: int(v) for k, v in r['counts'].items()})
    print(f"        agree {r['agree']:,} / {len(ref):,}  ({r['pct_agree']}%)   disagree {r['disagree']:,}")
