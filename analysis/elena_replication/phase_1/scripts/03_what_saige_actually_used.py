#!/usr/bin/env python3
"""
Phase 1 / check 03 — what the SAIGE run actually consumed.

Checks 01 and 02 compared the release against the phenotype artifacts in
`PMBBv4_phecodex/`. This check asks a different question: is the phenotype file on
disk today the one the analysis was actually run on?

It is not. The covariate files SAIGE consumed were built at 2026-07-31 20:21 from a
phenotype file whose row count the build log records as 57,636. The phenotype file on
disk carries 57,080 rows and a timestamp of 20:57 -- 36 minutes later. It was
overwritten after being read.

Writes only under phase_1/results/.
"""
import json
from pathlib import Path
import pandas as pd

HERE = Path(__file__).resolve().parent.parent
E    = Path("/project/hall/analysis/hearing-loss-genomics/analysis/elena")
REL  = Path("/static/PMBB/PMBB-Release-2026-4.0")
OUT  = HERE / "results"; OUT.mkdir(parents=True, exist_ok=True)

USED  = E / "HL_only_rarevariant/new_SAIGE_covariates/hearing_impairment_combined_SAIGE_Step1_covariates.txt"
DISK  = E / "rarevariantExWAS/PMBBv4_phecodex/hearing_impairment_PMBBv4_SAIGE.txt"
STAT  = E / "rarevariantExWAS/PMBBv4_phecodex/PMBBv4_hearing_tinnitus_all_statuses.csv.gz"
ANC   = REL / "Exome/PCA/combined/PMBB-Release-2026-4.0_genetic_exome.commonsnps.samples_ancestries.tsv"

used = pd.read_csv(USED, sep="\t", usecols=["IID", "PHENO"])
disk = pd.read_csv(DISK, sep="\t", usecols=["IID", "PHENO"])
exome = set(pd.read_csv(ANC, sep="\t", usecols=["IID"])["IID"])

st = pd.read_csv(STAT)
st = st[st.phenotype == "hearing_impairment"].set_index("person_id")["status"]

U, D = set(used.IID), set(disk.IID)
uc = set(used.loc[used.PHENO == 1, "IID"])
dc = set(disk.loc[disk.PHENO == 1, "IID"])

only_used = U - D
only_disk = D - U

report = {
    "consumed_by_saige": {"n": len(U), "cases": len(uc), "controls": len(U) - len(uc)},
    "phenotype_file_on_disk": {"n": len(D), "cases": len(dc), "controls": len(D) - len(dc)},
    "case_sets_identical": uc == dc,
    "consumed_has_no_exome": len(U - exome),
    "only_in_consumed": {
        "n": len(only_used),
        "by_phenotype_status": st.reindex(sorted(only_used)).value_counts(dropna=False).to_dict(),
    },
    "only_on_disk": {
        "n": len(only_disk),
        "by_phenotype_status": st.reindex(sorted(only_disk)).value_counts(dropna=False).to_dict(),
    },
}

# improperly included: people the phenotype rules exclude, present as controls in the run
improper = {p for p in only_used if str(st.get(p, "")).startswith("excluded")}
report["improperly_included_as_controls"] = len(improper)

(OUT / "03_what_saige_actually_used.json").write_text(json.dumps(report, indent=2, default=str))
pd.DataFrame({"person_id": sorted(improper),
              "phenotype_status": [st.get(p) for p in sorted(improper)]}) \
  .to_csv(OUT / "03_improperly_included_controls.csv", index=False)

print(json.dumps(report, indent=2, default=str))
