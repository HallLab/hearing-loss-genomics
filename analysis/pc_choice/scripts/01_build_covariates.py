#!/usr/bin/env python3
"""
PC choice experiment, step 1 -- build the three covariate sets to compare.

THE QUESTION. Our EUR ExWAS leaves lambda = 1.34 inside the MHC where the original
run's leaves 0.97. Two explanations have been offered and nothing in the numbers
separates them:

  count    4 PCs is too few, and more of the SAME family fixes it
  family   the exome PCs are the wrong substrate past PC1, and the imputed ones fix it

The second explanation exists because the two PC families turn out to be unrelated
beyond the first component. Measured per person on the release's own eigenvectors:
EUR exome_PC1 correlates 0.979 with imputed_PC1, and every later exome PC correlates
below 0.27 with ANY of the twenty imputed PCs. So "use more PCs" is not a question of
how many -- it is a question of which variables.

THREE ARMS, EUR only, since that is where the problem is:

  current   4 PCs, our within-ancestry exome PCA        (already run)
  A        10 PCs, the same PCA, its own 80% count      no cost
  B         9 PCs, the release's imputed PCA, 80% count costs 427 people

B's cost is not incidental. Those 427 have exome data and complete exome_PC1-20 and no
imputed PCs at all, and dropping them is one of the two cohort defects phase 1
corrected. Choosing B means reinstating that.

READOUT: lambda inside and outside the MHC. Whichever brings 1.34 down to ~1.0 wins the
argument, and it does so on a measurement rather than on a convention.
"""
import sys
from pathlib import Path

import pandas as pd

ROOT = Path("/project/hall/analysis/hearing-loss-genomics/analysis")
HERE = Path(__file__).resolve().parents[1]
BASE_COV = (ROOT / "gene_burden_elena_replication/phase_3/results/covariates"
            / "corrected_EUR_covariates.txt")
OUR_EIGENVEC = ROOT / "gene_burden_elena_replication/phase_3/data/pca/EUR.eigenvec"
RELEASE_IMPUTED = Path(
    "/static/PMBB/PMBB-Release-2026-4.0/Imputed/PCA/by_ancestry"
    "/PMBB-Release-2026-4.0_genetic_imputed.no1KG.eigenvectors.EUR.tsv")
OUT = HERE / "data"

KEEP = ["IID", "PHENO", "AGE", "AGE2", "SEX", "Batch"]


def main():
    base = pd.read_csv(BASE_COV, sep="\t")
    n0 = len(base)
    print(f"base cohort (the corrected EUR arm): {n0:,} samples, "
          f"{int((base.PHENO == 1).sum()):,} cases\n")

    # ---- arm A: 10 PCs from our own within-ancestry exome PCA
    ev = pd.read_csv(OUR_EIGENVEC, sep=r"\s+", header=None)
    ev = ev.rename(columns={1: "IID"})
    pcs = [c for c in ev.columns if isinstance(c, int) and c >= 2]
    ev = ev[["IID"] + pcs[:10]]
    ev.columns = ["IID"] + [f"PC{i+1}" for i in range(10)]
    a = base[KEEP].merge(ev, on="IID", how="left")
    missing = a[[f"PC{i+1}" for i in range(10)]].isna().any(axis=1).sum()
    if missing:
        sys.exit(f"REFUSING: arm A loses {missing} samples to missing PCs; our PCA was "
                 f"re-run precisely so that nobody is missing one.")
    dest = OUT / "EUR_A_exome10.txt"
    a.to_csv(dest, sep="\t", index=False)
    print(f"A  10 PCs, our exome PCA      {len(a):>7,} samples  (lost {n0-len(a)})"
          f"  -> {dest.name}")

    # ---- arm B: 9 PCs from the release's imputed PCA
    imp = pd.read_csv(RELEASE_IMPUTED, sep="\t")
    cols = [f"imputed_PC_EUR{i+1}" for i in range(9)]
    missing_cols = [c for c in cols if c not in imp.columns]
    if missing_cols:
        sys.exit(f"REFUSING: release file lacks {missing_cols}")
    imp = imp[["person_id"] + cols].rename(columns={"person_id": "IID"})
    imp.columns = ["IID"] + [f"PC{i+1}" for i in range(9)]
    b = base[KEEP].merge(imp, on="IID", how="inner")
    lost = n0 - len(b)
    dest = OUT / "EUR_B_imputed9.txt"
    b.to_csv(dest, sep="\t", index=False)
    print(f"B   9 PCs, release imputed    {len(b):>7,} samples  (lost {lost})"
          f"  -> {dest.name}")

    if lost == 0:
        print("\nNOTE: arm B lost nobody, which contradicts phase 1's finding that "
              "some cohort members have no imputed PCs. Check before trusting B.")
    else:
        cases_lost = int((base.PHENO == 1).sum()) - int((b.PHENO == 1).sum())
        print(f"\nB drops {lost} people, of whom {cases_lost} are cases. Those are the "
              f"participants\nwith exome data and no imputed PCs -- the group phase 1 "
              f"restored.")


if __name__ == "__main__":
    main()
