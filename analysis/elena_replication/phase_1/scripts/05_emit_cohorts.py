#!/usr/bin/env python3
"""
Phase 1 / check 05 — emit the two cohorts, and close the phase.

Phase 1's job is to establish who is in the study. Checks 01-04 established it;
this writes it down, so later phases consume a file instead of re-deriving a set
and risking a different answer.

Two arms, per pipeline_plan.md section 6:

  reproduction  57,632  the set the SAIGE run actually consumed.
                        Taken verbatim from the covariate file that survives --
                        NOT rebuilt, because the phenotype file it came from was
                        overwritten and a rebuild would not reproduce it.

  corrected     57,507  the analysable cohort the phenotype supports.
                        Derived from the release and the pipeline's own documented
                        rules, which check 02 reproduced person-for-person.

Scope: this emits the COHORT -- who is in, their phenotype, their ancestry. It does
not build covariate files; that is Phase 3.

Writes only under phase_1/results/.
"""
import json
from pathlib import Path
import pandas as pd

HERE = Path(__file__).resolve().parent.parent
REL  = Path("/static/PMBB/PMBB-Release-2026-4.0")
E    = Path("/project/hall/analysis/hearing-loss-genomics/analysis/elena")
OUT  = HERE / "results" / "cohorts"; OUT.mkdir(parents=True, exist_ok=True)

ANC  = REL / "Exome/PCA/combined/PMBB-Release-2026-4.0_genetic_exome.commonsnps.samples_ancestries.tsv"
USED = E / "HL_only_rarevariant/new_SAIGE_covariates/hearing_impairment_combined_SAIGE_Step1_covariates.txt"
STAT = E / "rarevariantExWAS/PMBBv4_phecodex/PMBBv4_hearing_tinnitus_all_statuses.csv.gz"

amap = dict(pd.read_csv(ANC, sep="\t", usecols=["IID", "Class"]).values)

# ---- reproduction: verbatim from what ran
rep = pd.read_csv(USED, sep="\t", usecols=["IID", "PHENO"])
rep["ancestry"] = rep.IID.map(amap)

# ---- corrected: re-derived from the release-verified phenotype
st = pd.read_csv(STAT)
st = st[st.phenotype == "hearing_impairment"].set_index("person_id")["status"]
cor = st[st.isin(["case", "control"])]
cor = pd.DataFrame({"IID": cor.index, "PHENO": (cor.values == "case").astype(int)})
cor["ancestry"] = cor.IID.map(amap)

manifest = {
    "purpose": "Phase 1 output: the two cohorts later phases consume.",
    "arms": {}, "relationship": {}, "provenance": {},
}

for name, df, src in [("reproduction", rep, str(USED)), ("corrected", cor, str(STAT))]:
    strata = {}
    for k in ["combined", "EUR", "AFR"]:
        d = df if k == "combined" else df[df.ancestry == k]
        path = OUT / f"{name}_{k}.tsv"
        d[["IID", "PHENO"]].to_csv(path, sep="\t", index=False)
        strata[k] = {"N": len(d), "cases": int(d.PHENO.sum()),
                     "controls": int((d.PHENO == 0).sum()),
                     "case_rate_pct": round(100 * d.PHENO.mean(), 2),
                     "file": path.name}
    manifest["arms"][name] = {"strata": strata,
                              "composition": df.ancestry.value_counts().to_dict()}
    manifest["provenance"][name] = src

R, C = set(rep.IID), set(cor.IID)
manifest["relationship"] = {
    "corrected_only": len(C - R), "reproduction_only": len(R - C), "shared": len(R & C),
    "note": "Neither is a subset of the other. 57,507 - 431 + 556 = 57,632.",
}
manifest["not_included"] = ("Covariates are not built here -- that is Phase 3. "
                            "The reproduction arm is copied verbatim rather than rebuilt, "
                            "because the phenotype file it derives from was overwritten.")

(HERE / "results" / "05_cohorts_manifest.json").write_text(json.dumps(manifest, indent=2, default=str))
print(json.dumps({k: manifest[k] for k in ["arms", "relationship"]}, indent=2, default=str))
