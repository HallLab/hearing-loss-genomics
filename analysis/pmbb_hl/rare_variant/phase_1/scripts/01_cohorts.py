#!/usr/bin/env python3
"""
Phase 1 -- the two cohorts, built in one pass.

THREE PHENOTYPES, ONE SCRIPT, and that is the point. The broad and the bilateral arms
share a control definition and an ancestry source. Building them separately means the
controls agree only if nobody edits one copy, and the two arms then differ for reasons
nobody intended. Here they are the same objects by construction, and the nesting
between the arms is asserted rather than hoped for.

THE CASE RULES

  broad       SO_396 or any of its children, on >= 2 distinct dates.
              This is the phenotype the original PMBB v4 run used.

  bilateral   SO_396.2 (sensorineural) AND SO_396.8 (bilateral) recorded on the SAME
              date, on >= 2 distinct dates.

The bilateral rule needs its awkward shape because phecodeX splits one diagnosis across
two orthogonal axes -- type and laterality -- as separate rows. A person coded
conductive-bilateral on one visit and sensorineural-unilateral on another carries both
.2 and .8 without ever having had bilateral sensorineural hearing loss. Requiring the
two codes on the same date is what prevents that, and it is not a detail: the looser
reading gives 4,029 cases where this one gives 3,164.

THE CONTROL RULE, identical for both arms: no ear-family evidence from any source. Not
"no hearing loss" -- anyone with any SO_39x diagnosis, or tinnitus from the observation
table, is excluded rather than used as a control.

WHY THE OBSERVATION TABLE. PMBB v4 moved the standard tinnitus codes (388.3x, H93.1x)
out of the phecode files and into OMOP `observation`. The phecode files keep only the
~3,016 pulsatile events; the other 25,094 are only in `observation`. A control set built
without that source silently contains tinnitus patients.

WHAT THIS PHASE DOES NOT DO. It defines who is analysable. The count that reaches SAIGE
is smaller, because phase 3 drops anyone missing a covariate -- 9 people with no AGE in
the broad combined cohort, for instance. Both numbers are emitted here so the later
phases cannot quietly conflate them.
"""
import json
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parents[1]
REL = Path("/static/PMBB/PMBB-Release-2026-4.0")
ANC = REL / ("Exome/PCA/combined/PMBB-Release-2026-4.0_genetic_exome"
             ".commonsnps.samples_ancestries.tsv")
OUT = HERE / "results"

PHENOTYPES = ["broad", "bilateral", "tinnitus"]
SENSORINEURAL = "SO_396.2"
BILATERAL = "SO_396.8"
TINNITUS = "SO_397"
MIN_DATES = 2


def load_evidence(cohort):
    """Every ear-family row for people with exome data, from both sources."""
    ear = pd.read_csv(HERE / "data/_ear_family.tsv", sep="\t",
                      names=["pid", "date", "code"], dtype=str)
    obs = pd.read_csv(HERE / "data/_tinnitus_obs.tsv", sep="\t",
                      names=["pid", "date", "code"], dtype=str)
    obs["code"] = "OBS_TINNITUS"
    ev = pd.concat([ear, obs], ignore_index=True)
    ev = ev[ev.date.notna() & (ev.date != "")]
    return ev[ev.pid.isin(cohort)]


def cases_broad(ev):
    hl = ev[ev.code.str.match(r"^SO_396($|\.)")]
    n = hl.groupby("pid").date.nunique()
    return set(n[n >= MIN_DATES].index), set(n[n == 1].index)


def cases_tinnitus(ev):
    """SO_397 or any child, OR the observation-table codes, on >= 2 distinct dates.

    The second source is not optional. PMBB v4 moved the standard tinnitus ICD codes
    (388.3x, H93.1x) out of the phecode files, which keep 3,016 of 28,110 events; the
    other 25,094 exist only in OMOP `observation`. An arm built on the phecode files
    alone would miss 89% of its own phenotype.
    """
    tin = ev[ev.code.str.match(rf"^{TINNITUS}($|\.)") | (ev.code == "OBS_TINNITUS")]
    n = tin.groupby("pid").date.nunique()
    return set(n[n >= MIN_DATES].index), set(n[n == 1].index)


def cases_bilateral(ev):
    pair = ev[ev.code.isin([SENSORINEURAL, BILATERAL])]
    per_date = pair.groupby(["pid", "date"]).code.nunique()
    qualifying = per_date[per_date == 2].reset_index().groupby("pid").date.nunique()
    return (set(qualifying[qualifying >= MIN_DATES].index),
            set(qualifying[qualifying == 1].index))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    anc = pd.read_csv(ANC, sep="\t", usecols=["IID", "Class"])
    cohort = set(anc.IID)
    ev = load_evidence(cohort)

    # One control set, used by both arms.
    controls = cohort - set(ev.pid)
    print(f"exome cohort            {len(cohort):>8,}")
    print(f"with any ear evidence   {len(set(ev.pid)):>8,}")
    print(f"controls (both arms)    {len(controls):>8,}\n")

    report = {"release": REL.name, "exome_cohort": len(cohort),
              "controls": len(controls), "min_dates": MIN_DATES, "arms": {}}

    built = {}
    builders = {"broad": cases_broad, "bilateral": cases_bilateral,
                "tinnitus": cases_tinnitus}
    for name in PHENOTYPES:
        fn = builders[name]
        cases, one_date = fn(ev)
        frame = pd.DataFrame(
            [(i, 1) for i in sorted(cases)] + [(i, 0) for i in sorted(controls)],
            columns=["IID", "PHENO"])
        frame = frame.merge(anc.rename(columns={"Class": "ancestry"}),
                            on="IID", how="left")
        frame.to_csv(OUT / f"cohort_{name}.tsv", sep="\t", index=False)
        built[name] = cases

        by_anc = frame.groupby(["ancestry", "PHENO"]).size().unstack(fill_value=0)
        report["arms"][name] = {
            "cases": len(cases), "controls": len(controls),
            "analysable": len(frame),
            "excluded_one_date_only": len(one_date),
            "excluded_other_ear": len(cohort) - len(cases) - len(controls)
                                  - len(one_date),
            "by_ancestry": {a: {"cases": int(r.get(1, 0)),
                                "controls": int(r.get(0, 0))}
                            for a, r in by_anc.iterrows()},
        }
        print(f"{name}")
        print(f"  cases                 {len(cases):>8,}")
        print(f"  analysable            {len(frame):>8,}")
        print(f"  excluded, one date    {len(one_date):>8,}  (they DO have it)")
        for a, r in by_anc.iterrows():
            print(f"    {a:<8} {int(r.get(1,0)):>6,} cases  {int(r.get(0,0)):>7,} controls")
        print()

    # The nesting the later phases depend on. Asserted, not assumed.
    #
    # SCOPED TO ONE PAIR ON PURPOSE. Bilateral sensorineural is a kind of hearing loss,
    # so every bilateral case must satisfy the broad rule. Tinnitus is NOT a kind of
    # hearing loss -- it is a different symptom that often accompanies it -- so the
    # same check applied to it would fail correctly and for the wrong reason.
    nested = built["bilateral"] <= built["broad"]
    report["bilateral_nested_in_broad"] = bool(nested)
    if not nested:
        stray = len(built["bilateral"] - built["broad"])
        raise SystemExit(
            f"REFUSING: {stray} bilateral cases are not broad cases. Bilateral "
            f"sensorineural is a kind of hearing loss, so every bilateral case must "
            f"satisfy the broad rule. If it does not, one of the two rules is wrong.")
    print(f"nesting verified: all {len(built['bilateral']):,} bilateral cases are "
          f"also broad cases")

    # The tinnitus arm OVERLAPS rather than nests, and the amount matters for every
    # later comparison, so it is measured and recorded rather than left implicit.
    overlap = built["tinnitus"] & built["broad"]
    report["tinnitus_overlap_with_broad"] = {
        "tinnitus_cases": len(built["tinnitus"]),
        "also_broad_cases": len(overlap),
        "share": round(len(overlap) / len(built["tinnitus"]), 3),
        "tinnitus_only": len(built["tinnitus"] - built["broad"]),
    }
    o = report["tinnitus_overlap_with_broad"]
    print(f"tinnitus overlaps rather than nests: {o['also_broad_cases']:,} of "
          f"{o['tinnitus_cases']:,} ({100*o['share']:.0f}%) are also broad cases")

    (OUT / "01_cohorts.json").write_text(json.dumps(report, indent=2) + "\n")
    print(f"\nwritten: " + ", ".join(f"cohort_{p}.tsv" for p in PHENOTYPES) + ", 01_cohorts.json")


if __name__ == "__main__":
    main()
