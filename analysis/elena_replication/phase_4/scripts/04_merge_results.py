#!/usr/bin/env python3
"""
Phase 4, collation -- merge the 594 SAIGE step-2 outputs into one table.

This script refuses rather than globs.  Elena's equivalent
(HL_only_rarevariant/step3_2_merge.bsub) reads whatever `rglob("*.txt")`
returns, prints value_counts, and writes the CSV.  Two defects in her run got
through it unseen, and both are things a merge can check:

  * 11 of her 792 cells never produced output -- all combined/ALL at
    maf 0.001 and 0.01, on chromosomes 1, 2, 11, 12, 16, 17, 19.
  * 3 of her files lost their header line (AFR/pDM/chr8, all three MAF
    cutoffs), so pandas read the first data row as column names.  Concatenating
    by name parked 5,391 rows in phantom columns with Region, Group and Pvalue
    blank, and consumed gene ABRA's own result as the header.

So this script asserts three things before it writes anything: every expected
cell exists, every file carries the expected header, and no row has a blank or
unnamed gene.  Any failure stops the merge and names the files.

Input : phase_4/results/step2/<cohort>/<mask>/<cohort>_chr<N>_<mask>_maf<M>.txt
Output: phase_4/results/burden_<cohort>.tsv   (one per cohort)
        phase_4/results/burden_all_cohorts.tsv
        phase_4/results/04_merge_manifest.json
"""

import json
import re
import sys
from pathlib import Path

import pandas as pd

PHASE4 = Path(__file__).resolve().parents[1]
STEP2 = PHASE4 / "results" / "step2"
OUT = PHASE4 / "results"

COHORTS = ["combined", "EUR", "AFR"]
MASKS = ["pLOF", "pDM", "pLOF_pDM"]
MAFS = ["0.0001", "0.001", "0.01"]
CHROMS = list(range(1, 23))

EXPECTED_HEADER = [
    "Region", "Group", "max_MAF", "Pvalue", "Pvalue_Burden", "Pvalue_SKAT",
    "BETA_Burden", "SE_Burden", "MAC", "MAC_case", "MAC_control",
    "Number_rare", "Number_ultra_rare",
]

FNAME = re.compile(r"^(?P<cohort>\w+)_chr(?P<chr>\d+)_(?P<mask>\w+)_maf(?P<maf>[\d.]+)\.txt$")


def expected_cells():
    return {
        (c, m, f, k)
        for c in COHORTS for m in MASKS for f in MAFS for k in CHROMS
    }


def main():
    problems = []

    # ---- 1. every expected cell present, and nothing unexpected on disk ----
    found = {}
    for cohort in COHORTS:
        for mask in MASKS:
            d = STEP2 / cohort / mask
            if not d.is_dir():
                continue
            for f in d.glob("*.txt"):
                if f.name.endswith((".index", ".singleAssoc.txt")):
                    continue
                m = FNAME.match(f.name)
                if not m:
                    problems.append(f"unparseable filename: {f}")
                    continue
                key = (m["cohort"], mask, m["maf"], int(m["chr"]))
                if key in found:
                    problems.append(f"duplicate cell {key}: {f} and {found[key]}")
                found[key] = f

    missing = sorted(expected_cells() - set(found))
    if missing:
        problems.append(f"{len(missing)} of {len(expected_cells())} cells missing")
        for key in missing[:25]:
            problems.append(f"    missing: {key[0]} / {key[1]} / maf{key[2]} / chr{key[3]}")
        if len(missing) > 25:
            problems.append(f"    ... and {len(missing) - 25} more")

    extra = sorted(set(found) - expected_cells())
    for key in extra:
        problems.append(f"unexpected cell on disk: {key}")

    # ---- 2. every file carries the expected header ----
    for key, f in sorted(found.items()):
        with f.open() as fh:
            header = fh.readline().rstrip("\n").split("\t")
        if header != EXPECTED_HEADER:
            problems.append(
                f"header mismatch in {f.name}: got {header[:4]}... "
                f"({len(header)} cols, expected {len(EXPECTED_HEADER)})"
            )

    if problems:
        print("MERGE REFUSED\n", file=sys.stderr)
        for p in problems:
            print(f"  {p}", file=sys.stderr)
        sys.exit(1)

    # ---- 3. read, annotate, and check gene names ----
    frames = []
    for (cohort, mask, maf, chrom), f in sorted(found.items()):
        df = pd.read_csv(f, sep="\t", dtype={"Region": str, "Group": str})
        df["Cohort"] = cohort
        df["Mask"] = mask
        df["CHR"] = chrom
        # max_MAF is already a SAIGE column and equals the cutoff; assert rather
        # than overwrite it, which is what Elena's merge does
        if not df.empty:
            got = set(df["max_MAF"].astype(float))
            if got != {float(maf)}:
                print(f"MERGE REFUSED: {f.name} reports max_MAF {got}, "
                      f"filename says {maf}", file=sys.stderr)
                sys.exit(1)
        frames.append(df)

    merged = pd.concat(frames, ignore_index=True)

    blank = merged["Region"].isna() | merged["Region"].isin(["", "-"])
    if blank.any():
        print(f"MERGE REFUSED: {int(blank.sum())} rows with a blank or unnamed "
              f"gene. Elena's masks carried a '-' pseudo-gene pooling every "
              f"variant VEP could not assign to a symbol; ours should not.",
              file=sys.stderr)
        sys.exit(1)

    nonnum = pd.to_numeric(merged["Pvalue"], errors="coerce").isna()
    if nonnum.any():
        print(f"MERGE REFUSED: {int(nonnum.sum())} rows with a non-numeric "
              f"Pvalue", file=sys.stderr)
        sys.exit(1)

    # ---- 4. write ----
    OUT.mkdir(parents=True, exist_ok=True)
    manifest = {"cells": len(found), "rows": int(len(merged)), "by_cohort": {}}

    for cohort in COHORTS:
        sub = merged[merged["Cohort"] == cohort]
        path = OUT / f"burden_{cohort}.tsv"
        sub.to_csv(path, sep="\t", index=False)
        n_tests = len(sub)
        manifest["by_cohort"][cohort] = {
            "rows": int(n_tests),
            "genes": int(sub["Region"].nunique()),
            "bonferroni_0.05": 0.05 / n_tests if n_tests else None,
            "min_pvalue": float(pd.to_numeric(sub["Pvalue"]).min()) if n_tests else None,
            "file": path.name,
        }

    merged.to_csv(OUT / "burden_all_cohorts.tsv", sep="\t", index=False)
    (OUT / "04_merge_manifest.json").write_text(json.dumps(manifest, indent=2))

    print(f"{len(found)} cells, {len(merged):,} rows")
    for cohort, v in manifest["by_cohort"].items():
        print(f"  {cohort:9s} {v['rows']:>9,} tests  {v['genes']:>6,} genes  "
              f"min p {v['min_pvalue']:.3g}  bonferroni {v['bonferroni_0.05']:.3g}")


if __name__ == "__main__":
    main()
