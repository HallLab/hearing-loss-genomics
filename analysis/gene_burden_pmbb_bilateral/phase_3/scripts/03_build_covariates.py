#!/usr/bin/env python3
"""
Phase 3, step 3.2 -- SAIGE covariate files, one per cohort.

One arm. This analysis has a position; it is not comparing two.

Columns follow the pipeline's convention exactly, so that a file from here can be
dropped into any of the existing scripts: IID PHENO AGE AGE2 SEX Batch PC1..PCn

PC counts, premise P7:

    combined   5    the release's exome PCA, which covers everyone
    EUR        4    our within-ancestry PCA   (the pipeline used 9)
    AFR        3    our within-ancestry PCA   (the pipeline used 10)

WHERE THE COMPONENTS COME FROM, and why the source differs by stratum. For the
combined cohort the release's own exome PCs are right: the structure being
adjusted for is between ancestries, which is what a global PCA measures. Inside
EUR or AFR a global PC1 is busy separating continents and says almost nothing
about structure within the group, so those use the within-ancestry PCA copied
into phase_3/data/pca.

Note what is NOT re-derived: that PCA ran over everyone with exome data per
ancestry, not over any analysis cohort, so changing the phenotype cannot
invalidate it. See PROVENANCE.md.

PEOPLE WITHOUT A RECORDED AGE DROP OUT, and that is correct. Age is a covariate
the model consumes; a row without it cannot be fitted. This is the opposite of
the cut premise P3 reverses, where people were dropped for missing an imputed
PC the model never touches. Same mechanism, opposite justification.

Input : phase_1/results/cohort.tsv, phase_3/data/pca/{EUR,AFR}.eigenvec, release covariates
Output: phase_3/results/covariates/{combined,EUR,AFR}.txt
        phase_3/results/03_covariates.json
"""
import json
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parents[1]
REL = Path("/static/PMBB/PMBB-Release-2026-4.0")
OUT = HERE / "results" / "covariates"
OUT.mkdir(parents=True, exist_ok=True)

PCA = HERE / "data/pca"
N_PCS = {"combined": 5, "EUR": 4, "AFR": 3}
COV = REL / "Phenotype/4.0/PMBB-Release-2026-4.0_phenotype_covariates.txt"

cohort = pd.read_csv(HERE.parent / "phase_1/results/cohort.tsv", sep="\t",
                     dtype={"IID": str})

base = pd.read_csv(COV, sep="\t", low_memory=False, dtype={"person_id": str},
                   usecols=["person_id", "batch", "sequenced_gender", "sample_age"]
                           + [f"exome_PC{i}" for i in range(1, 6)])
base = base.rename(columns={"person_id": "IID", "batch": "Batch", "sample_age": "AGE"})
base["SEX"] = base.sequenced_gender.map({"Male": 1, "Female": 2})
base["AGE2"] = base.AGE ** 2


def pcs_for(stratum):
    k = N_PCS[stratum]
    if stratum == "combined":
        d = base[["IID"] + [f"exome_PC{i}" for i in range(1, k + 1)]].copy()
        return d.rename(columns={f"exome_PC{i}": f"PC{i}" for i in range(1, k + 1)})
    ev = pd.read_csv(PCA / f"{stratum}.eigenvec", sep=r"\s+", header=None,
                     dtype={1: str})
    ev.columns = ["FID", "IID"] + [f"PC{i}" for i in range(1, ev.shape[1] - 1)]
    return ev[["IID"] + [f"PC{i}" for i in range(1, k + 1)]]


report = {"phenotype": "bilateral sensorineural, premise P1",
          "pc_counts": N_PCS, "cohorts": {}}

for stratum in ["combined", "EUR", "AFR"]:
    c = cohort if stratum == "combined" else cohort[cohort.ancestry == stratum]
    d = c[["IID", "PHENO"]].merge(base[["IID", "AGE", "AGE2", "SEX", "Batch"]],
                                  on="IID", how="left")
    d = d.merge(pcs_for(stratum), on="IID", how="left")
    cols = ["IID", "PHENO", "AGE", "AGE2", "SEX", "Batch"] + \
           [f"PC{i}" for i in range(1, N_PCS[stratum] + 1)]
    d = d[cols]
    missing = {col: int(d[col].isna().sum()) for col in cols if d[col].isna().any()}
    clean = d.dropna()
    clean = clean.astype({"PHENO": int, "SEX": int, "Batch": int})
    path = OUT / f"{stratum}.txt"
    clean.to_csv(path, sep="\t", index=False)
    report["cohorts"][stratum] = {
        "requested": int(len(c)), "written": int(len(clean)),
        "dropped_for_missing_covariate": int(len(d) - len(clean)),
        "missing_by_column": missing,
        "cases": int((clean.PHENO == 1).sum()),
        "controls": int((clean.PHENO == 0).sum()),
        "n_pcs": N_PCS[stratum], "file": path.name,
    }

(HERE / "results" / "03_covariates.json").write_text(json.dumps(report, indent=2))
for s, v in report["cohorts"].items():
    print(f"{s:9s} {v['written']:>7,} written  "
          f"({v['cases']:>5,} cases, {v['controls']:>7,} controls)  "
          f"{v['n_pcs']} PCs   dropped {v['dropped_for_missing_covariate']}")
