#!/usr/bin/env python3
"""
Phase 1 / check 01 — who is in the study.

Traces the cohort from the institutional release through to the sample list the
SAIGE run actually used, and attributes every drop along the way.

Reads:
  - the release: covariates, exome ancestry list
  - the pipeline (read-only, as the reference being checked):
      PMBBv4_phecodex/PMBBv4_hearing_tinnitus_all_statuses.csv.gz
      PMBBv4_phecodex/hearing_impairment_PMBBv4_SAIGE.txt
      HL_TIN_PMBBv4_keep.txt
      HL_TIN_PMBBv4_SAIGE_samples.txt

Writes everything under phase_1/results/. Writes nothing outside this directory.
"""
import gzip, json
from pathlib import Path
import pandas as pd

REL  = Path("/static/PMBB/PMBB-Release-2026-4.0")
PIPE = Path("/project/hall/analysis/hearing-loss-genomics/analysis/elena/rarevariantExWAS")
OUT  = Path(__file__).resolve().parent.parent / "results"
OUT.mkdir(parents=True, exist_ok=True)

COV   = REL / "Phenotype/4.0/PMBB-Release-2026-4.0_phenotype_covariates.txt"
ANC   = REL / "Exome/PCA/combined/PMBB-Release-2026-4.0_genetic_exome.commonsnps.samples_ancestries.tsv"
PHECO = PIPE / "PMBBv4_phecodex"

# ---------------------------------------------------------------- release side
cov = pd.read_csv(COV, sep="\t", low_memory=False)
exome_ids = set(pd.read_csv(ANC, sep="\t", usecols=["IID"])["IID"])

imp_pcs   = [c for c in cov.columns if c.startswith("imputed_PC")]
exome_pcs = [c for c in cov.columns if c.startswith("exome_PC")]

cov["has_imputed_pcs"] = cov[imp_pcs].notna().all(axis=1)
cov["has_exome_pcs"]   = cov[exome_pcs].notna().all(axis=1)

# ---------------------------------------------------------------- pipeline side
status = pd.read_csv(PHECO / "PMBBv4_hearing_tinnitus_all_statuses.csv.gz")
hl = status[status.phenotype == "hearing_impairment"].set_index("person_id")["status"]

saige_hl = pd.read_csv(PHECO / "hearing_impairment_PMBBv4_SAIGE.txt", sep="\t", usecols=["IID", "PHENO"])
keep     = set(pd.read_csv(PIPE / "HL_TIN_PMBBv4_keep.txt", sep="\t", header=None)[0])
samples  = set(pd.read_csv(PIPE / "HL_TIN_PMBBv4_SAIGE_samples.txt", sep="\t", usecols=["IID"])["IID"])

# ---------------------------------------------------------------- the chain
analysable = set(hl[hl.isin(["case", "control"])].index)
delivered  = set(saige_hl.IID)
dropped    = analysable - delivered
dropped_by_status = hl.loc[sorted(dropped)].value_counts().to_dict()

lost_at_samples = keep - samples          # where the cohort actually shrinks

no_imp = set(cov.loc[~cov.has_imputed_pcs, "person_id"])
no_exo = set(cov.loc[~cov.has_exome_pcs,   "person_id"])

chain = {
    "exome_cohort_release":        len(exome_ids),
    "phenotyped_all_statuses":     int(hl.shape[0]),
    "analysable_case_or_control":  len(analysable),
    "delivered_to_saige":          len(delivered),
    "dropped":                     len(dropped),
    "dropped_by_status":           {k: int(v) for k, v in dropped_by_status.items()},
    "lost_between_keep_and_samples": len(lost_at_samples),
    "lost_all_have_exome":         len(lost_at_samples & exome_ids) == len(lost_at_samples),
    "lost_all_lack_imputed_pcs":   lost_at_samples.issubset(no_imp),
    "lost_all_have_exome_pcs":     len(lost_at_samples & no_exo) == 0,
    "cohort_lacking_imputed_pcs":  len(no_imp),
    "cohort_lacking_exome_pcs":    len(no_exo),
    "saige_file_uses_pcs":         "exome_PC1-6",
}

(OUT / "01_sample_frame.json").write_text(json.dumps(chain, indent=2))

pd.DataFrame({"person_id": sorted(dropped),
              "hl_status": hl.loc[sorted(dropped)].values}) \
  .to_csv(OUT / "01_dropped_participants.csv", index=False)

for k, v in chain.items():
    print(f"{k:32s} {v}")
