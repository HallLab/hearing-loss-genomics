#!/usr/bin/env python3
"""
Phase 1 -- cases and controls, bilateral sensorineural only.

This is the one place this analysis differs from the replication it inherits
from. Everything downstream is the replication's corrected pipeline unchanged.

WHY. Douglas Epstein, 2026-10-02:

    "we typically exclude one-sided hearing loss, and we just focus on the
     bilateral sensorineural hearing loss"
    "unilateral hearing loss we exclude, because it's less likely genetic,
     more likely environment related"

The replication's case definition is `SO_396`, the phecodeX parent, which
bundles every kind of hearing loss: conductive, mixed, unilateral, unspecified.
That faithfully reproduces what the pipeline did, which is what a replication
owes. It is not what the clinical lead considers the phenotype.

HOW THE CODES MAKE THIS EXPRESSIBLE. The children of SO_396 are not mutually
exclusive subtypes -- they encode two orthogonal axes at once, and a single ICD
code maps to one of each:

    type        .1 conductive   .2 sensorineural   .3 mixed   .5 sudden idiopathic
    laterality  .8 bilateral    .9 unilateral

So `H90.3  Sensorineural hearing loss, bilateral` carries SO_396, SO_396.2 and
SO_396.8 together. "Bilateral sensorineural" is therefore the intersection
.2 AND .8, and it can be built from the release without touching ICD strings.

THE RULE OF 2 IS APPLIED TO THE RESTRICTED SET, NOT AFTERWARDS, AND THE TWO
AXES MUST MEET ON THE SAME DATE. Three readings of "bilateral sensorineural on
two dates" are possible, and they do not agree:

    A  >= 2 dates carrying .2 OR .8, among people who have both   4,029
    B  >= 2 dates of .2 AND >= 2 dates of .8, counted separately  3,259
    C  >= 2 dates on which the SAME date carries .2 and .8        3,164   <-- used

A is too loose: a date with only `H90.0 Conductive hearing loss, bilateral`
contributes .8 and counts toward the total, so conductive loss can carry a case
over the line. B is closer but still lets a sensorineural-unilateral visit and a
conductive-bilateral visit combine into a case that was never diagnosed with
bilateral sensorineural loss at all. C is what the phrase actually means: on at
least two distinct dates, the diagnosis recorded was bilateral sensorineural.

C is a subset of both. The 865 people A admits and C does not are the measure of
how much the loose reading would have cost.

CONTROLS ARE UNCHANGED, AND DELIBERATELY SO. A control is still someone with no
ear-family evidence at all, from either source. Narrowing the cases must not
widen the controls: the people now excluded as cases -- unilateral, conductive,
mixed -- still carry hearing loss, and putting them in the control group would
be worse than leaving them in the case group.

Input : phase_1/data/_ear_family.tsv, _tinnitus_obs.tsv  (script 00)
Output: phase_1/results/cohort.tsv, 01_phenotype.json
"""
import json
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parents[1]
REL = Path("/static/PMBB/PMBB-Release-2026-4.0")
ANC = REL / ("Exome/PCA/combined/PMBB-Release-2026-4.0_genetic_exome"
             ".commonsnps.samples_ancestries.tsv")
OUT = HERE / "results"
OUT.mkdir(parents=True, exist_ok=True)

TYPE_SENSORINEURAL = "SO_396.2"
LATERALITY_BILATERAL = "SO_396.8"

ear = pd.read_csv(HERE / "data/_ear_family.tsv", sep="\t",
                  names=["pid", "date", "code"], dtype=str)
obs = pd.read_csv(HERE / "data/_tinnitus_obs.tsv", sep="\t",
                  names=["pid", "date", "code"], dtype=str)
obs["code"] = "OBS_TINNITUS"
ear = pd.concat([ear, obs], ignore_index=True)
ear = ear[ear.date.notna() & (ear.date != "")]

anc = pd.read_csv(ANC, sep="\t", usecols=["IID", "Class"])
cohort = set(anc.IID)
ear = ear[ear.pid.isin(cohort)]

# ---- who is a case: reading C, both axes present on the same date ----
pair = ear[ear.code.isin([TYPE_SENSORINEURAL, LATERALITY_BILATERAL])]
per_date = pair.groupby(["pid", "date"]).code.nunique()
qualifying = per_date[per_date == 2].reset_index().groupby("pid").date.nunique()
cases = set(qualifying[qualifying >= 2].index)
one_date = set(qualifying[qualifying == 1].index)

# ---- who is a control: no ear-family evidence from any source ----
any_ear = set(ear.pid)
controls = cohort - any_ear

# ---- everyone else is excluded, and the reasons are counted ----
excluded = cohort - cases - controls
hearing_loss = set(ear[ear.code.str.startswith("SO_396")].pid)
# Split rather than lumped: 'one date' people DO have the target phenotype and
# fail only the rule of 2, which is a different reason from not having it.
bilat_sn_one_date = one_date - cases
hl_never_bilat_sn = hearing_loss - cases - bilat_sn_one_date
other_ear_only = excluded - hearing_loss

status = pd.Series("excluded", index=sorted(cohort), dtype=object)
status.loc[sorted(controls)] = "control"
status.loc[sorted(cases)] = "case"

frame = pd.DataFrame({"IID": status.index, "status": status.values})
frame["PHENO"] = frame.status.map({"case": 1, "control": 0})
frame = frame.merge(anc.rename(columns={"Class": "ancestry"}), on="IID", how="left")
frame[frame.PHENO.notna()].to_csv(OUT / "cohort.tsv", sep="\t", index=False)

report = {
    "definition": {
        "cases": f"{TYPE_SENSORINEURAL} and {LATERALITY_BILATERAL} present on "
                 f"the same date, on >= 2 distinct dates",
        "controls": "no ear-family evidence from conditions_phecode_x or observation",
        "source": "Douglas Epstein, meeting 2026-10-02",
    },
    "counts": {
        "cohort_with_exome": len(cohort),
        "cases": len(cases),
        "controls": len(controls),
        "excluded_total": len(excluded),
        "excluded_bilateral_sensorineural_but_one_date_only": len(bilat_sn_one_date),
        "excluded_hearing_loss_never_bilateral_sensorineural": len(hl_never_bilat_sn),
        "excluded_other_ear_evidence_only": len(other_ear_only),
    },
    "by_ancestry": {
        a: {"cases": int(((frame.ancestry == a) & (frame.PHENO == 1)).sum()),
            "controls": int(((frame.ancestry == a) & (frame.PHENO == 0)).sum())}
        for a in sorted(frame.ancestry.dropna().unique())
    },
}
(OUT / "01_phenotype.json").write_text(json.dumps(report, indent=2))

print(f"cohort with exome   {len(cohort):>7,}")
print(f"  cases             {len(cases):>7,}   bilateral sensorineural, >= 2 dates")
print(f"  controls          {len(controls):>7,}   no ear evidence at all")
print(f"  excluded          {len(excluded):>7,}")
print(f"    bilateral sensorineural, but on one date only        : "
      f"{len(bilat_sn_one_date):,}")
print(f"    hearing loss, never bilateral sensorineural          : "
      f"{len(hl_never_bilat_sn):,}")
print(f"    other ear evidence only                              : "
      f"{len(other_ear_only):,}")
print("\nby ancestry:")
for a, v in report["by_ancestry"].items():
    print(f"  {a:9s} {v['cases']:>6,} cases  {v['controls']:>7,} controls")
