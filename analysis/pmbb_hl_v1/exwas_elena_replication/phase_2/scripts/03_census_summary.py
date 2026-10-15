#!/usr/bin/env python3
"""
Phase 2, step 3 -- aggregate the 66 histograms and fix the secondary bar.

This is the script that closes premise P6. It must run to completion BEFORE any
step 2 p-value exists, because its output is a significance threshold and a
threshold chosen after seeing the data is not a threshold.

It refuses rather than reporting a partial census. Three ways a census can be
wrong while looking fine:

  - a task is still running, so its histogram is on disk but short. Checked by
    asking LSF, not by trusting the file.
  - a task failed after writing a header, so the file parses and sums to zero.
    Checked against the .bim count for that chromosome.
  - the high-frequency bins are non-empty, in which case the [0.99, 0.999) bin
    straddles the MAF 0.01 boundary and the count is off by however many variants
    sit exactly at ALT1 frequency 0.99. Reported, not silently absorbed.
"""
import json
import subprocess
import sys
from pathlib import Path

BASE = Path("/project/hall/analysis/hearing-loss-genomics/analysis/exwas_elena_replication")
CENSUS = BASE / "phase_2/results/census"
GENO = Path("/project/hall/analysis/hearing-loss-genomics/analysis/elena/"
            "rarevariant_geneburden/plink_deduplicated")
COHORTS = ["combined", "EUR", "AFR"]
CHROMS = list(range(1, 23))

# A variant is tested when MAF >= 0.01, i.e. ALT1 frequency in [0.01, 0.99].
TESTED_BINS = ["0.01", "0.05", "0.95"]
# Bins that must be empty for TESTED_BINS to be an exact count rather than a
# near-count. See the docstring.
BOUNDARY_BINS = ["0.99", "0.999", "0.9999"]


def lsf_in_flight():
    """Any of our census tasks still queued or running?"""
    try:
        out = subprocess.run(["bjobs", "-J", "maf_census"], capture_output=True,
                             text=True, timeout=60).stdout
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return None  # cannot tell; caller decides
    return [l for l in out.splitlines()[1:] if l.strip()]


def read_hist(path):
    """{bin_start_as_written: count}. Raises on anything unexpected."""
    lines = [l for l in path.read_text().splitlines() if l.strip()]
    if not lines:
        raise ValueError("file is empty")
    if lines[0].split() != ["#BIN_START", "OBS_CT"]:
        raise ValueError(f"unexpected header: {lines[0]!r}")
    bins = {}
    for ln in lines[1:]:
        f = ln.split()
        if len(f) != 2:
            raise ValueError(f"expected 2 fields, got {len(f)}: {ln!r}")
        bins[f[0]] = int(f[1])
    return bins


def main():
    flight = lsf_in_flight()
    if flight:
        sys.exit(f"REFUSING: {len(flight)} census tasks still in flight.\n"
                 + "\n".join("  " + l for l in flight[:5]))
    if flight is None:
        print("WARNING: could not query LSF; relying on the per-file checks below\n")

    bim_counts = {}
    for c in CHROMS:
        bim = GENO / f"chr{c}_deduplicated.bim"
        if not bim.exists():
            sys.exit(f"REFUSING: {bim} not found")
        bim_counts[c] = sum(1 for _ in bim.open())

    problems, cells = [], {}
    for cohort in COHORTS:
        for c in CHROMS:
            path = CENSUS / f"{cohort}.chr{c}.afreq.alt1.bins"
            if not path.exists():
                problems.append(f"{cohort} chr{c}: histogram absent")
                continue
            try:
                bins = read_hist(path)
            except ValueError as e:
                problems.append(f"{cohort} chr{c}: {e}")
                continue
            total = sum(bins.values())
            if total != bim_counts[c]:
                problems.append(f"{cohort} chr{c}: bins sum to {total:,}, "
                                f".bim has {bim_counts[c]:,}")
                continue
            cells[(cohort, c)] = bins

    if problems:
        print(f"REFUSING -- {len(problems)} of 66 cells are not usable:", file=sys.stderr)
        for p in problems[:20]:
            print(f"  {p}", file=sys.stderr)
        if len(problems) > 20:
            print(f"  ... and {len(problems) - 20} more", file=sys.stderr)
        sys.exit(1)

    print(f"all 66 cells present, each summing to its .bim count\n")

    out = {"cohorts": {}, "premise": "P6", "tested_definition": "MAF >= 0.01"}
    straddle_total = 0

    for cohort in COHORTS:
        spectrum, tested, on_file = {}, 0, 0
        for c in CHROMS:
            for k, v in cells[(cohort, c)].items():
                spectrum[k] = spectrum.get(k, 0) + v
        on_file = sum(spectrum.values())
        tested = sum(spectrum.get(b, 0) for b in TESTED_BINS)
        straddle = sum(spectrum.get(b, 0) for b in BOUNDARY_BINS)
        straddle_total += straddle
        bar = 0.05 / tested if tested else None

        out["cohorts"][cohort] = {
            "variants_on_file": on_file,
            "variants_tested": tested,
            "secondary_bar": bar,
            "primary_bar": 5e-8,
            "alt1_spectrum": spectrum,
            "high_frequency_bins": straddle,
        }

        print(f"{cohort}")
        print(f"  on file           {on_file:>12,}")
        print(f"  tested (MAF>=1%)  {tested:>12,}   {100*tested/on_file:.2f}%")
        print(f"  secondary bar     {bar:>12.3e}   = 0.05 / {tested:,}")
        print(f"  primary bar       {5e-8:>12.3e}   "
              f"({'stricter' if 5e-8 < bar else 'LOOSER'} of the two)")
        print()

    if straddle_total:
        print(f"ATTENTION: {straddle_total:,} variants have ALT1 frequency >= 0.99.")
        print("The [0.99, 0.999) bin straddles MAF 0.01, so 'tested' is a near-count")
        print("for those. Resolve before quoting the bar as exact.\n")
    else:
        print("No variant has ALT1 frequency >= 0.99, so REF is the major allele")
        print("throughout and the tested count is exact rather than approximate.\n")

    dest = BASE / "phase_2/results/census_summary.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print(f"written: {dest.relative_to(BASE)}")


if __name__ == "__main__":
    main()
