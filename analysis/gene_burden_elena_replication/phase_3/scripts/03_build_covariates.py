#!/usr/bin/env python3
"""
Phase 3 / Step 6 — build the SAIGE covariate files for both arms.

Mirrors the pipeline's Step 6 (6.2 combined, 6.4 by-ancestry, 6.5 batch) and closes
Phase 3 the way Phases 1 and 2 closed: with files, not numbers inside a script.

    reproduction   the pipeline's own covariates, copied verbatim
    corrected      rebuilt on the corrected cohort with the PC counts the scree supports

PC counts, from phase_3/scripts/01_pc_selection.ipynb:

    combined   5    the release's exome PCA -- what the pipeline used, and correct
    EUR        4    our re-run PCA   (the pipeline used 9)
    AFR        3    our re-run PCA   (the pipeline used 10)

WHERE THE PCs COME FROM, and why it differs by stratum:

  combined  the release's exome_PC1-5. The pipeline used these and they cover everyone,
            so there is nothing to re-derive.
  EUR, AFR  our own within-ancestry PCA (02_pca_by_ancestry.bsub). The pipeline's ran
            after the cohort was cut against the imputed .fam, so 203 EUR and 34 AFR
            participants restored by Phase 1 have no components in it.

The reproduction arm keeps the pipeline's PCs untouched. Ours are not identical even
for shared participants -- a PCA over a different sample set yields different
eigenvectors -- so mixing them would make a later difference un-attributable.

Covariate columns follow the pipeline exactly: IID PHENO AGE AGE2 SEX Batch PC1..PCn.

Writes only under phase_3/results/.
"""
import json
import shutil
from pathlib import Path
import pandas as pd

HERE = Path(__file__).resolve().parent.parent
REL  = Path("/static/PMBB/PMBB-Release-2026-4.0")
PIPE = Path("/project/hall/analysis/hearing-loss-genomics/analysis/elena")
OUT  = HERE / "results" / "covariates"; OUT.mkdir(parents=True, exist_ok=True)

COHORTS = HERE.parent / "phase_1/results/cohorts"
PCA     = HERE / "data/pca"
N_PCS   = {"combined": 5, "EUR": 4, "AFR": 3}
COV     = REL / "Phenotype/4.0/PMBB-Release-2026-4.0_phenotype_covariates.txt"

# ---------------------------------------------------------------- base covariates
base = pd.read_csv(COV, sep="\t", low_memory=False,
                   usecols=["person_id", "batch", "sequenced_gender", "sample_age"]
                          + [f"exome_PC{i}" for i in range(1, 6)])
base = base.rename(columns={"person_id": "IID", "batch": "Batch", "sample_age": "AGE"})
base["SEX"] = base.sequenced_gender.map({"Male": 1, "Female": 2})
base["AGE2"] = base.AGE ** 2

def pcs_for(stratum):
    """Per-person components, from the source appropriate to this stratum."""
    k = N_PCS[stratum]
    if stratum == "combined":
        d = base[["IID"] + [f"exome_PC{i}" for i in range(1, k + 1)]].copy()
        return d.rename(columns={f"exome_PC{i}": f"PC{i}" for i in range(1, k + 1)})
    ev = pd.read_csv(PCA / f"{stratum}.eigenvec", sep=r"\s+", header=None)
    ev.columns = ["FID", "IID"] + [f"PC{i}" for i in range(1, ev.shape[1] - 1)]
    return ev[["IID"] + [f"PC{i}" for i in range(1, k + 1)]]

manifest = {"pc_counts": N_PCS, "arms": {}, "pc_source": {
    "combined": "release exome_PC1-5",
    "EUR": "phase_3/data/pca/EUR.eigenvec (our re-run)",
    "AFR": "phase_3/data/pca/AFR.eigenvec (our re-run)"}}

# ---------------------------------------------------------------- corrected arm
built = {}
for stratum in ["combined", "EUR", "AFR"]:
    coh = pd.read_csv(COHORTS / f"corrected_{stratum}.tsv", sep="\t")
    d = coh.merge(base[["IID", "AGE", "AGE2", "SEX", "Batch"]], on="IID", how="left")
    d = d.merge(pcs_for(stratum), on="IID", how="left")
    cols = ["IID", "PHENO", "AGE", "AGE2", "SEX", "Batch"] + \
           [f"PC{i}" for i in range(1, N_PCS[stratum] + 1)]
    d = d[cols]
    missing = {c: int(d[c].isna().sum()) for c in cols if d[c].isna().any()}
    path = OUT / f"corrected_{stratum}_covariates.txt"
    d.dropna().to_csv(path, sep="\t", index=False)
    built[stratum] = {"requested": len(coh), "written": int(len(d.dropna())),
                      "dropped_for_missing": int(len(d) - len(d.dropna())),
                      "missing_by_column": missing, "n_pcs": N_PCS[stratum],
                      "file": path.name}
manifest["arms"]["corrected"] = built

# ---------------------------------------------------------------- reproduction arm
src = PIPE / "HL_only_rarevariant/new_SAIGE_covariates"
repro = {}
for stratum in ["combined", "EUR", "AFR"]:
    s = src / f"hearing_impairment_{stratum}_SAIGE_Step1_covariates.txt"
    dst = OUT / f"reproduction_{stratum}_covariates.txt"
    shutil.copyfile(s, dst)
    d = pd.read_csv(dst, sep="\t", nrows=1)
    repro[stratum] = {"rows": sum(1 for _ in dst.open()) - 1,
                      "n_pcs": sum(1 for c in d.columns if c.startswith("PC")),
                      "copied_from": str(s), "file": dst.name}
manifest["arms"]["reproduction"] = repro
manifest["note"] = ("The reproduction arm is copied verbatim, not rebuilt: the phenotype file it "
                    "derives from was overwritten, so a rebuild would not reproduce it.")

(HERE / "results" / "03_covariates_manifest.json").write_text(json.dumps(manifest, indent=2))
print(json.dumps(manifest, indent=2))
