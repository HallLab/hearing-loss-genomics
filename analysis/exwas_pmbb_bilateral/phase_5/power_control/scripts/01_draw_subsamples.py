#!/usr/bin/env python3
"""
Power control -- is the restricted arm's null a lost signal or a lost half of the cases?

Finding 3 showed ClinGen deafness genes clustering at the top of the BROAD
analysis (5 of the top 50, p = 5.9e-06) and nowhere near the top of the
RESTRICTED one. Two explanations survive that:

  signal  the restriction removed the cases that carried the biology
  power   the restriction removed half the cases, and 3,164 is simply too few

They are separated by holding N fixed. This draws the broad phenotype down to
the restricted arm's exact case count and runs it again. If the enrichment
survives at matched N, the restriction removed signal. If it vanishes, the
restriction removed only cases and the phenotype decision is exonerated.

FIVE REPLICATES, because one draw can be lucky in either direction and a single
run would let me report whichever answer I drew. Each replicate is an
independent random subsample of cases; controls are held fixed, since they are
identical between the two arms anyway.

Only the combined cohort. That is where the enrichment is, and running EUR and
AFR would add 150 jobs to answer a question neither is powered to settle.

Input : elena_replication/phase_3/results/covariates/corrected_combined_covariates.txt
Output: phase_5/power_control/data/rep<N>.txt  (5 covariate files)
"""
import json
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parents[1]
BROAD = Path("/project/hall/analysis/hearing-loss-genomics/analysis/elena_replication/"
             "phase_3/results/covariates/corrected_combined_covariates.txt")
TARGET_CASES = 3164          # the restricted arm, exactly
N_REP = 5
SEED = 20261008

d = pd.read_csv(BROAD, sep="\t", dtype={"IID": str})
cases = d[d.PHENO == 1]
controls = d[d.PHENO == 0]
assert len(cases) > TARGET_CASES, "nothing to subsample"

report = {"source": str(BROAD), "broad_cases": int(len(cases)),
          "controls": int(len(controls)), "target_cases": TARGET_CASES,
          "replicates": N_REP, "seed": SEED, "draws": {}}

(HERE / "data").mkdir(parents=True, exist_ok=True)
for r in range(1, N_REP + 1):
    drawn = cases.sample(n=TARGET_CASES, random_state=SEED + r)
    sub = pd.concat([drawn, controls]).sort_values("IID")
    path = HERE / "data" / f"rep{r}.txt"
    sub.to_csv(path, sep="\t", index=False)
    report["draws"][f"rep{r}"] = {
        "cases": int((sub.PHENO == 1).sum()),
        "controls": int((sub.PHENO == 0).sum()),
        "n": int(len(sub)),
        "file": path.name,
    }
    print(f"rep{r}: {len(sub):,} samples  {(sub.PHENO==1).sum():,} cases")

(HERE / "results" / "01_subsamples.json").parent.mkdir(parents=True, exist_ok=True)
(HERE / "results" / "01_subsamples.json").write_text(json.dumps(report, indent=2))

# how much do the draws overlap each other, and the restricted arm?
restricted = set(pd.read_csv(HERE.parent.parent / "phase_1/results/cohort.tsv",
                             sep="\t", dtype={"IID": str}).query("PHENO==1").IID)
drawn_sets = {r: set(pd.read_csv(HERE / "data" / f"rep{r}.txt", sep="\t",
                                 dtype={"IID": str}).query("PHENO==1").IID)
              for r in range(1, N_REP + 1)}
print(f"\noverlap of each draw with the restricted arm's cases "
      f"(expected ~{TARGET_CASES * len(restricted) / len(cases):.0f} by chance):")
for r, s in drawn_sets.items():
    print(f"  rep{r}: {len(s & restricted):,}")
