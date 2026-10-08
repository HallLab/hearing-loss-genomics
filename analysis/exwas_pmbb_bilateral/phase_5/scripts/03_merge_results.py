#!/usr/bin/env python3
"""
Phase 5, step 5.1 -- merge the 66 step-2 outputs.

This script refuses rather than globs. The pipeline this analysis descends from
read whatever rglob returned, printed value_counts and wrote the CSV, and two
defects got through it unseen: 11 of 792 cells never produced output, and three
files lost their header line so 5,391 rows landed in phantom columns. Both are
things a merge can check.

WHAT THE FILES LOOK LIKE HERE, which differs from the replication. Premise P10
runs one SAIGE call per cohort and chromosome over the whole grid, so each file
holds every cell for that chromosome plus the omnibus:

    Group       max_MAF          what it is
    pLOF        1e-4/1e-3/1e-2   one annotation, one cutoff -- a cell
    pDM         1e-4/1e-3/1e-2   "
    pLOF;pDM    1e-4/1e-3/1e-2   "
    Cauchy      NA               the omnibus over all nine, per gene

THE CAUCHY ROW IS THE HEADLINE. It is one p-value per gene with the multiple
testing across the grid already paid, so the bar is 0.05/genes with no
denominator left to argue about. The replication had to compute this afterwards
because it sliced the grid into separate calls; here SAIGE produces it.

A gene's `search_penalty` -- omnibus over its best cell -- says how much of the
nine it took to find the signal. 1x means the effect shows everywhere; 9x means
it lives in one corner.

Input : phase_4/results/step2/<cohort>/<cohort>_chr<N>.txt
Output: phase_5/results/burden_all_cohorts.tsv   the cells
        phase_5/results/omnibus_<cohort>.tsv     one row per gene
        phase_5/results/03_merge_manifest.json
"""
import json
import re
import subprocess
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parents[1]
STEP2 = HERE.parent / "phase_4" / "results" / "step2"
OUT = HERE / "results"
OUT.mkdir(parents=True, exist_ok=True)

COHORTS = ["combined", "EUR", "AFR"]
CHROMS = list(range(1, 23))
ALLOW_RUNNING = "--allow-running" in sys.argv

EXPECTED_HEADER = [
    "Region", "Group", "max_MAF", "Pvalue", "Pvalue_Burden", "Pvalue_SKAT",
    "BETA_Burden", "SE_Burden", "MAC", "MAC_case", "MAC_control",
    "Number_rare", "Number_ultra_rare",
]
CELL_GROUPS = {"pLOF", "pDM", "pLOF;pDM"}
FNAME = re.compile(r"^(?P<cohort>\w+)_chr(?P<chr>\d+)\.txt$")


def jobs_in_flight():
    try:
        out = subprocess.run(["bjobs", "-a", "-J", "s2_bilat*"],
                             capture_output=True, text=True, timeout=60).stdout
    except (OSError, subprocess.SubprocessError):
        return None
    return [l for l in out.splitlines()[1:]
            if len(l.split()) > 2 and l.split()[2] in
            {"PEND", "RUN", "SSUSP", "USUSP", "PSUSP"}]


def main():
    problems = []

    flight = jobs_in_flight()
    if flight and not ALLOW_RUNNING:
        problems.append(f"{len(flight)} step-2 elements still in flight. SAIGE creates "
                        f"its output file when a task starts, so merging now would read "
                        f"files still being written. Pass --allow-running to override.")

    # ---- every expected cell, and nothing unexpected ----
    found = {}
    for cohort in COHORTS:
        d = STEP2 / cohort
        if not d.is_dir():
            continue
        for f in d.glob("*.txt"):
            if f.name.endswith((".index", ".singleAssoc.txt")):
                continue
            m = FNAME.match(f.name)
            if not m:
                problems.append(f"unparseable filename: {f}")
                continue
            found[(cohort, int(m["chr"]))] = f

    expected = {(c, k) for c in COHORTS for k in CHROMS}
    for key in sorted(expected - set(found)):
        problems.append(f"missing: {key[0]} / chr{key[1]}")
    for key in sorted(set(found) - expected):
        problems.append(f"unexpected on disk: {key}")

    # ---- structural integrity ----
    for key, f in sorted(found.items()):
        raw = f.read_bytes()
        if not raw:
            problems.append(f"empty: {f.name}")
            continue
        if not raw.endswith(b"\n"):
            problems.append(f"truncated, no trailing newline: {f.name}")
        lines = raw.decode().splitlines()
        if lines[0].split("\t") != EXPECTED_HEADER:
            problems.append(f"header mismatch: {f.name}")
            continue
        for i, line in enumerate(lines[1:], start=2):
            if len(line.split("\t")) != len(EXPECTED_HEADER):
                problems.append(f"{f.name} line {i}: wrong field count")
                break

    if problems:
        print("MERGE REFUSED\n", file=sys.stderr)
        for p in problems[:30]:
            print(f"  {p}", file=sys.stderr)
        if len(problems) > 30:
            print(f"  ... and {len(problems) - 30} more", file=sys.stderr)
        sys.exit(1)

    # ---- read ----
    frames = []
    for (cohort, chrom), f in sorted(found.items()):
        df = pd.read_csv(f, sep="\t", dtype={"Region": str, "Group": str})
        df["Cohort"] = cohort
        df["CHR"] = chrom
        frames.append(df)
    allrows = pd.concat(frames, ignore_index=True)

    blank = allrows.Region.isna() | allrows.Region.isin(["", "-"])
    if blank.any():
        print(f"MERGE REFUSED: {int(blank.sum())} rows with a blank or unnamed gene",
              file=sys.stderr)
        sys.exit(1)

    cells = allrows[allrows.Group.isin(CELL_GROUPS)].copy()
    omni = allrows[allrows.Group == "Cauchy"].copy()

    cells["Pvalue"] = pd.to_numeric(cells.Pvalue)
    omni["Pvalue"] = pd.to_numeric(omni.Pvalue)

    # SAIGE should give exactly one Cauchy row per gene per cohort. If it does
    # not, the headline is not what it claims to be, so this is checked.
    dup = omni.duplicated(["Cohort", "Region"]).sum()
    if dup:
        print(f"MERGE REFUSED: {dup} duplicate Cauchy rows -- the omnibus is "
              f"supposed to be one per gene", file=sys.stderr)
        sys.exit(1)

    cells.to_csv(OUT / "burden_all_cohorts.tsv", sep="\t", index=False)

    manifest = {"cells_files": len(found), "cell_rows": int(len(cells)),
                "by_cohort": {}}
    for cohort in COHORTS:
        o = omni[omni.Cohort == cohort].copy()
        best = (cells[cells.Cohort == cohort].groupby("Region")
                .Pvalue.min().rename("best_cell_p"))
        o = o.merge(best, on="Region", how="left")
        o["search_penalty"] = o.Pvalue / o.best_cell_p
        bar = 0.05 / len(o)
        o["bonferroni_bar"] = bar
        o["significant"] = o.Pvalue < bar
        o = o.rename(columns={"Pvalue": "omnibus_p"})
        o.sort_values("omnibus_p").to_csv(OUT / f"omnibus_{cohort}.tsv",
                                          sep="\t", index=False)
        manifest["by_cohort"][cohort] = {
            "genes": int(len(o)),
            "cell_tests": int((cells.Cohort == cohort).sum()),
            "min_omnibus_p": float(o.omnibus_p.min()),
            "min_cell_p": float(cells[cells.Cohort == cohort].Pvalue.min()),
            "bar": bar,
            "n_significant": int(o.significant.sum()),
            "median_search_penalty": float(o.search_penalty.median()),
        }

    (OUT / "03_merge_manifest.json").write_text(json.dumps(manifest, indent=2))

    print(f"{len(found)} files, {len(cells):,} cell tests, "
          f"{len(omni):,} omnibus rows\n")
    print(f"{'cohort':9s} {'genes':>7s} {'omnibus p':>11s} {'bar':>10s} "
          f"{'sig':>4s}  {'best cell':>10s}")
    for c, v in manifest["by_cohort"].items():
        print(f"{c:9s} {v['genes']:>7,} {v['min_omnibus_p']:>11.3g} "
              f"{v['bar']:>10.2e} {v['n_significant']:>4} "
              f"{v['min_cell_p']:>10.3g}")


if __name__ == "__main__":
    main()
