#!/usr/bin/env python3
"""
Phase 3 -- covariate files, from the release's own exome PCA.

THE CHANGE THIS REBUILD MAKES. The first pass re-ran a within-ancestry PCA of its own.
The stated reason was that the original pipeline's PCA had a coverage gap. That was
wrong in a checkable way: the original computed no PCA at all -- its delivered PCs are
identical to the release's `exome_PC*` columns at max|difference| = 0.00e+00 across all
24 components -- and what its `step7_ldprune.bsub` builds is the GRM, with no `--pca`
anywhere in its scripts. A GRM step was read as a PCA step.

The release's by-ancestry exome PCA covers 100% of the analysable cohort. So there was
never a gap to close here, and this phase uses the release's PCs.

The GRM gap is a different thing and it is real -- 196 EUR and 34 AFR cohort members
are absent from the original's pruned genotypes. That is handled in 02_build_grm.bsub.

HOW MANY PCs: the count that reaches 80% of the variance in each cohort's own
eigenvalues, which is the lab convention Nikki stated on 2026-10-14.

    combined   11 PCs   80.3%
    EUR        15 PCs   83.3%
    AFR        16 PCs   84.1%

A caveat worth carrying rather than hiding: "80% of the variance" here means 80% of what
the top twenty eigenvalues capture, because plink publishes twenty. It is not 80% of
total genotypic variance, which twenty components come nowhere near. The convention is
still the convention -- it is just a proportion of a stated denominator rather than of
everything.

Note the direction of travel. The original used 5 / 9 / 10 and reached 64.2% / 60.8% /
58.6% on these same eigenvalues, so it did not meet the rule either. Neither did the
first pass, at 64.2% / 66.7% / 56.5% on its own PCA. This is the first run that does.

SOURCES
    combined PCs    Phenotype/4.0/..._phenotype_covariates.txt, exome_PC1..11
    EUR / AFR PCs   Exome/PCA/by_ancestry/..._eigenvectors.{EUR,AFR}.tsv
    AGE, SEX, Batch the same covariate table: sample_age, sequenced_gender, batch

Column layout follows the original exactly: IID PHENO AGE AGE2 SEX Batch PC1..PCn.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[1]
PHASE1 = HERE.parent / "phase_1/results"
REL = Path("/static/PMBB/PMBB-Release-2026-4.0")
COV = REL / "Phenotype/4.0/PMBB-Release-2026-4.0_phenotype_covariates.txt"
PCA_BY_ANC = REL / "Exome/PCA/by_ancestry"
OUT = HERE / "results/covariates"

PHENOTYPES = ["broad", "bilateral", "tinnitus"]
# cohort -> (how many PCs, where they come from)
PLAN = {
    "combined": (11, None),
    "EUR": (15, "EUR"),
    "AFR": (16, "AFR"),
}


def variance_reached(cohort, n):
    if cohort == "combined":
        f = (REL / "Exome/PCA/combined"
             / "PMBB-Release-2026-4.0_genetic_exome.no1KG.eigenvalues.tsv")
    else:
        f = (PCA_BY_ANC
             / f"PMBB-Release-2026-4.0_genetic_exome.no1KG.eigenvalues.{cohort}.tsv")
    with f.open() as fh:
        fh.readline()
        ev = np.array([float(x) for x in fh.readline().split()])
    return float(np.cumsum(100 * ev / ev.sum())[n - 1])


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    base = pd.read_csv(COV, sep="\t", low_memory=False)
    base = base.rename(columns={"person_id": "IID", "batch": "Batch",
                                "sample_age": "AGE"})
    base["SEX"] = base.sequenced_gender.map({"Male": 1, "Female": 2})
    base["AGE2"] = base.AGE ** 2

    report = {"pc_source": "release exome PCA", "rule": "80% of the top-20 variance",
              "arms": {}}
    problems = []

    for pheno in PHENOTYPES:
        coh = pd.read_csv(PHASE1 / f"cohort_{pheno}.tsv", sep="\t")
        report["arms"][pheno] = {}
        print(f"=== {pheno} ===")

        for cohort, (npcs, anc) in PLAN.items():
            sub = coh if cohort == "combined" else coh[coh.ancestry == anc]
            analysable = len(sub)

            if anc is None:
                pcs = base[["IID"] + [f"exome_PC{i+1}" for i in range(npcs)]]
            else:
                f = (PCA_BY_ANC / f"PMBB-Release-2026-4.0_genetic_exome.no1KG"
                     f".eigenvectors.{anc}.tsv")
                pcs = pd.read_csv(f, sep="\t").rename(columns={"person_id": "IID"})
                pcs = pcs[["IID"] + [f"exome_{anc}_PC{i+1}" for i in range(npcs)]]
            pcs.columns = ["IID"] + [f"PC{i+1}" for i in range(npcs)]

            # The whole point of using the release's PCA: nobody should be missing one.
            no_pc = analysable - len(sub[["IID"]].merge(pcs, on="IID"))
            if no_pc:
                problems.append(f"{pheno}/{cohort}: {no_pc} cohort members have no PC "
                                f"in the release PCA, which contradicts the premise "
                                f"of this phase")

            d = sub.merge(base[["IID", "AGE", "AGE2", "SEX", "Batch"]], on="IID",
                          how="left").merge(pcs, on="IID", how="left")
            cols = ["IID", "PHENO", "AGE", "AGE2", "SEX", "Batch"] + \
                   [f"PC{i+1}" for i in range(npcs)]
            before = len(d)
            missing_by_col = {c: int(d[c].isna().sum()) for c in cols
                              if d[c].isna().any()}
            d = d[cols].dropna()
            d["PHENO"] = d.PHENO.astype(int)

            dest = OUT / f"{pheno}_{cohort}.txt"
            d.to_csv(dest, sep="\t", index=False)

            var = variance_reached(cohort, npcs)
            report["arms"][pheno][cohort] = {
                "analysable": analysable, "written": len(d),
                "dropped_for_missing": before - len(d),
                "missing_by_column": missing_by_col,
                "cases": int((d.PHENO == 1).sum()),
                "controls": int((d.PHENO == 0).sum()),
                "n_pcs": npcs, "variance_pct": round(var, 1),
                "pc_source": "release exome combined" if anc is None
                             else f"release exome {anc}",
                "file": dest.name,
            }
            r = report["arms"][pheno][cohort]
            print(f"  {cohort:<9} analysable {analysable:>7,}  tested {len(d):>7,}  "
                  f"(dropped {r['dropped_for_missing']})  "
                  f"{r['cases']:>6,} cases  {npcs} PCs = {var:.1f}%")
            if missing_by_col:
                print(f"            missing: {missing_by_col}")
        print()

    if problems:
        print("REFUSING:", file=sys.stderr)
        for p in problems:
            print(f"  {p}", file=sys.stderr)
        sys.exit(1)

    (HERE / "results/01_covariates.json").write_text(
        json.dumps(report, indent=2) + "\n")
    print("every cohort member has a PC from the release — no one dropped for that")
    print(f"written: results/covariates/, 01_covariates.json")


if __name__ == "__main__":
    main()
