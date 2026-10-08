#!/usr/bin/env python3
"""
Complement run -- do the cases the restriction EXCLUDES carry the signal?

The restricted cases are a strict subset of the broad ones: 3,164 of 6,752,
checked. The complement is therefore exactly the 3,588 the restriction throws
away -- hearing loss that is unilateral, conductive, mixed, unspecified, or
bilateral sensorineural on only one date.

WHAT THE TWO OUTCOMES WOULD MEAN, written before the run.

  enrichment in the complement
      Known deafness genes associating with hearing loss that is NOT bilateral
      sensorineural. Biologically odd, and the likely explanation is not
      biological: ICD coding imprecision. H91.90 'Unspecified hearing loss,
      unspecified ear' alone covers 2,698 people, and someone with genuine
      bilateral sensorineural loss is often coded exactly that way. This would
      say the restriction discards real cases that are merely badly labelled.

  nothing in the complement
      The signal sits in the cases the restriction keeps, and the restricted
      arm's null is lost power, full stop. The phenotype decision comes out
      clean.

Either way the question is whether ICD codes separate what we assume they
separate -- which is the case for audiograms, measured instead of asserted.

ON SIZE. The complement is 3,588 against the restricted arm's 3,164 and the
power replicates' 3,164 -- 13% more cases. It is run at its natural size rather
than cut, because discarding 424 real cases to make a round number is worse
than carrying the caveat. The direction matters and is favourable: more cases
biases TOWARD finding enrichment, so a null here is a stronger null than a
matched one would be.

Output: phase_5/complement/data/complement.txt
"""
import json
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parents[1]
BROAD = Path("/project/hall/analysis/hearing-loss-genomics/analysis/elena_replication/"
             "phase_3/results/covariates/corrected_combined_covariates.txt")
RESTRICTED = HERE.parent.parent / "phase_3/results/covariates/combined.txt"

b = pd.read_csv(BROAD, sep="\t", dtype={"IID": str})
r = pd.read_csv(RESTRICTED, sep="\t", dtype={"IID": str})
keep = set(r[r.PHENO == 1].IID)
broad_cases = set(b[b.PHENO == 1].IID)
assert keep <= broad_cases, "restricted cases are not a subset of broad cases"

comp = b[(b.PHENO == 0) | ((b.PHENO == 1) & ~b.IID.isin(keep))].copy()
(HERE / "data").mkdir(parents=True, exist_ok=True)
comp.sort_values("IID").to_csv(HERE / "data" / "complement.txt", sep="\t", index=False)

report = {
    "broad_cases": len(broad_cases), "restricted_cases": len(keep),
    "complement_cases": int((comp.PHENO == 1).sum()),
    "controls": int((comp.PHENO == 0).sum()),
    "n": int(len(comp)),
    "note": "complement is 13% larger than the 3,164 comparison points; "
            "that biases toward finding enrichment, so a null is a strong null",
}
(HERE / "results").mkdir(parents=True, exist_ok=True)
(HERE / "results" / "01_cohort.json").write_text(json.dumps(report, indent=2))
print(json.dumps(report, indent=2))
