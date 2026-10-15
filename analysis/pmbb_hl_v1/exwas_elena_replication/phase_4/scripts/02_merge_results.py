#!/usr/bin/env python3
"""
Phase 4, step 2 -- merge the 66 step 2 outputs, one table per cohort.

This script refuses rather than reporting a partial result. SAIGE creates its
output file when a task STARTS, so file existence proves nothing -- in the burden
arm 15 of 547 cells were empty files that a glob happily merged. The checks here
are the ones that lesson paid for:

  - LSF is asked whether anything is still in flight, rather than trusting mtimes
  - every cell must have the exact same header
  - every data row must have the same field count as the header
  - the row count must be at or below the phase 2 plink2 census for that cell,
    and within 0.5% of it (see phase_2/results/FINDINGS.md for why the census is
    an upper bound and not an equality: plink2 divides by observed genotypes,
    SAIGE by the whole cohort)
  - p-values must parse and sit in (0, 1]

Output: results/merged_<cohort>.tsv plus results/MERGE_REPORT.md
"""
import math
import subprocess
import sys
from pathlib import Path

BASE = Path("/project/hall/analysis/hearing-loss-genomics/analysis/exwas_elena_replication")
STEP2 = BASE / "phase_4/results/step2"
CENSUS = BASE / "phase_2/results/census"
OUT = BASE / "phase_4/results"
COHORTS = ["combined", "EUR", "AFR"]
CHROMS = list(range(1, 23))
PRIMARY_BAR = 5e-8
SHORTFALL_CAP = 0.005  # 0.5%; largest observed was 0.05%


def in_flight():
    names = ["exwas_s2", "exwas_s2_retry"]
    live = []
    for n in names:
        try:
            out = subprocess.run(["bjobs", "-J", n], capture_output=True,
                                 text=True, timeout=60).stdout
        except (FileNotFoundError, subprocess.TimeoutExpired):
            return None
        live += [l for l in out.splitlines()[1:] if l.strip()]
    return live


def ceiling(cohort, chrom):
    h = CENSUS / f"{cohort}.chr{chrom}.afreq.alt1.bins"
    if not h.exists():
        return None
    return sum(int(l.split()[1]) for l in h.read_text().splitlines()[1:]
               if l.split()[0] in ("0.01", "0.05", "0.95"))


def main():
    live = in_flight()
    if live:
        sys.exit(f"REFUSING: {len(live)} step 2 tasks still in flight.\n"
                 + "\n".join("  " + l for l in live[:5]))
    if live is None:
        print("WARNING: could not query LSF; relying on the structural checks\n")

    problems = []
    header_ref = None
    cells = {}

    for cohort in COHORTS:
        for chrom in CHROMS:
            path = STEP2 / cohort / f"{cohort}_chr{chrom}.txt"
            tag = f"{cohort} chr{chrom}"

            if not path.exists():
                problems.append(f"{tag}: output absent")
                continue
            lines = path.read_text().splitlines()
            if not lines:
                problems.append(f"{tag}: file is empty (task started, never finished)")
                continue

            header = lines[0].split("\t")
            if header_ref is None:
                header_ref = header
            elif header != header_ref:
                problems.append(f"{tag}: header differs from the other cells")
                continue

            body = [l for l in lines[1:] if l.strip()]
            if not body:
                problems.append(f"{tag}: header only, no tested variants")
                continue

            nwrong = sum(1 for l in body if len(l.split("\t")) != len(header))
            if nwrong:
                problems.append(f"{tag}: {nwrong} rows with the wrong field count")
                continue

            cap = ceiling(cohort, chrom)
            if cap is None:
                problems.append(f"{tag}: no phase 2 census histogram to check against")
                continue
            short = cap - len(body)
            if short < 0:
                problems.append(f"{tag}: {len(body)} rows exceeds the plink2 ceiling "
                                f"{cap}; SAIGE cannot test more than plink2 counts")
                continue
            if short > cap * SHORTFALL_CAP:
                problems.append(f"{tag}: {len(body)} of {cap} rows, shortfall {short} "
                                f"(>{SHORTFALL_CAP:.1%}) is too large for the "
                                f"missingness mechanism")
                continue

            pcol = header.index("p.value")
            bad_p = 0
            for l in body:
                try:
                    v = float(l.split("\t")[pcol])
                except ValueError:
                    bad_p += 1
                    continue
                if not (0 < v <= 1) or math.isnan(v):
                    bad_p += 1
            if bad_p:
                problems.append(f"{tag}: {bad_p} p-values outside (0, 1] or unparseable")
                continue

            cells[(cohort, chrom)] = (body, cap, short)

    if problems:
        print(f"REFUSING -- {len(problems)} of 66 cells unusable:", file=sys.stderr)
        for p in problems[:25]:
            print(f"  {p}", file=sys.stderr)
        sys.exit(1)

    print(f"all 66 cells pass: identical headers, consistent field counts,")
    print(f"row counts at or under the plink2 ceiling, p-values in range\n")

    report = ["# Phase 4 — merge report", "",
              "66 of 66 cells merged. Checks in `02_merge_results.py`.", ""]
    report += ["| cohort | tested | plink2 ceiling | shortfall | secondary bar | "
               "min p | hits at 5e-8 | hits at secondary |",
               "|---|---:|---:|---:|---:|---:|---:|---:|"]

    for cohort in COHORTS:
        rows, cap_tot, short_tot = [], 0, 0
        for chrom in CHROMS:
            body, cap, short = cells[(cohort, chrom)]
            rows += body
            cap_tot += cap
            short_tot += short

        pcol = header_ref.index("p.value")
        ps = [float(r.split("\t")[pcol]) for r in rows]
        n = len(rows)
        secondary = 0.05 / n

        dest = OUT / f"merged_{cohort}.tsv"
        with dest.open("w") as fh:
            fh.write("\t".join(header_ref) + "\n")
            for r in rows:
                fh.write(r + "\n")

        hits_p = sum(1 for v in ps if v < PRIMARY_BAR)
        hits_s = sum(1 for v in ps if v < secondary)
        print(f"{cohort}")
        print(f"  tested          {n:>12,}   (plink2 ceiling {cap_tot:,}, "
              f"shortfall {short_tot})")
        print(f"  secondary bar   {secondary:>12.3e}")
        print(f"  min p           {min(ps):>12.3e}")
        print(f"  hits  5e-8      {hits_p:>12,}")
        print(f"  hits  secondary {hits_s:>12,}")
        print(f"  written         {dest.name}")
        print()

        report.append(f"| {cohort} | {n:,} | {cap_tot:,} | {short_tot} | "
                      f"{secondary:.3e} | {min(ps):.3e} | {hits_p} | {hits_s} |")

    report += ["", "Shortfall is the plink2 census minus what SAIGE tested. plink2 "
               "divides by observed", "genotypes, SAIGE by the whole cohort, so the "
               "census is an upper bound. See", "`phase_2/results/FINDINGS.md`.", ""]
    (OUT / "MERGE_REPORT.md").write_text("\n".join(report) + "\n")
    print(f"written: {(OUT / 'MERGE_REPORT.md').relative_to(BASE)}")


if __name__ == "__main__":
    main()
