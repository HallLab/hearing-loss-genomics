#!/usr/bin/env python3
"""
ExWAS phase 4 -- merge the 132 cells, with the checks the first pass paid for.

Refuses rather than reporting a partial result. SAIGE creates its output file when a
task STARTS, so existence proves nothing -- in the burden arm's first pass 15 of 547
cells were empty files that a glob happily merged.

THE ROW-COUNT CHECK IS AN UPPER BOUND, NOT AN EQUALITY. plink2 computes MAF over
observed genotypes, SAIGE over the whole cohort, so a variant carrying missingness can
clear 1% for plink2 and not for SAIGE. The reverse is impossible. Verified in the first
pass on AFR chr1: all six differing variants straddled the threshold exactly as that
predicts, and all 41 discrepancies across 66 cells ran the same way. So SAIGE's count
must be at or below the census, and within 0.5% of it.

TWO BARS, both fixed before step 2 ran. 5e-8 is the convention and the primary; the
secondary is 0.05 over the variants actually tested, which for exome-only testing is
about ten times looser. No FDR: that was explicit on 2026-10-14, and the split is
principled -- FDR is a threshold over a count of tests we performed, while 5e-8 counts
nothing we did.
"""
import json
import math
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
STEP2 = HERE / "results/step2"
CENSUS = HERE.parent / "phase_2/results/census"
OUT = HERE / "results"
ARMS = [f"{p}_{c}" for p in ("broad", "bilateral", "tinnitus")
        for c in ("combined", "EUR", "AFR")]
CHROMS = range(1, 23)
PRIMARY = 5e-8
SHORTFALL_CAP = 0.005


def in_flight():
    try:
        out = subprocess.run(["bjobs"], capture_output=True, text=True,
                             timeout=60).stdout
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return None
    return [l for l in out.splitlines()[1:] if l.strip()]


def ceiling(tag, chrom):
    h = CENSUS / f"{tag}.chr{chrom}.afreq.alt1.bins"
    if not h.exists():
        return None
    return sum(int(l.split()[1]) for l in h.read_text().splitlines()[1:]
               if l.split()[0] in ("0.01", "0.05", "0.95"))


def main():
    live = in_flight()
    if live:
        sys.exit(f"REFUSING: {len(live)} jobs still in flight.")
    if live is None:
        print("WARNING: could not query LSF; relying on the structural checks\n")

    problems, header_ref, cells = [], None, {}
    for tag in ARMS:
        for chrom in CHROMS:
            path = STEP2 / tag / f"{tag}_chr{chrom}.txt"
            label = f"{tag} chr{chrom}"
            if not path.exists():
                problems.append(f"{label}: absent")
                continue
            lines = path.read_text().splitlines()
            if not lines:
                problems.append(f"{label}: empty (task started, never finished)")
                continue
            header = lines[0].split("\t")
            if header_ref is None:
                header_ref = header
            elif header != header_ref:
                problems.append(f"{label}: header differs from the other cells")
                continue
            body = [l for l in lines[1:] if l.strip()]
            if not body:
                problems.append(f"{label}: header only")
                continue
            if any(len(l.split("\t")) != len(header) for l in body):
                problems.append(f"{label}: rows with the wrong field count")
                continue
            cap = ceiling(tag, chrom)
            if cap is None:
                problems.append(f"{label}: no census histogram to check against")
                continue
            short = cap - len(body)
            if short < 0:
                problems.append(f"{label}: {len(body)} rows exceeds the plink2 "
                                f"ceiling {cap}")
                continue
            if short > cap * SHORTFALL_CAP:
                problems.append(f"{label}: shortfall {short} of {cap} exceeds "
                                f"{SHORTFALL_CAP:.1%}")
                continue
            pcol = header.index("p.value")
            bad = 0
            for l in body:
                try:
                    v = float(l.split("\t")[pcol])
                except ValueError:
                    bad += 1
                    continue
                if not (0 < v <= 1) or math.isnan(v):
                    bad += 1
            if bad:
                problems.append(f"{label}: {bad} p-values outside (0,1]")
                continue
            cells[(tag, chrom)] = (body, cap, short)

    if problems:
        print(f"REFUSING -- {len(problems)} of 198 cells unusable:", file=sys.stderr)
        for p in problems[:20]:
            print(f"  {p}", file=sys.stderr)
        sys.exit(1)

    print("all 198 cells pass: identical headers, consistent field counts,\n"
          "row counts at or under the plink2 ceiling, p-values in range\n")
    report = {"primary_bar": PRIMARY, "arms": {}}

    for tag in ARMS:
        rows, cap_tot, short_tot = [], 0, 0
        for chrom in CHROMS:
            body, cap, short = cells[(tag, chrom)]
            rows += body
            cap_tot += cap
            short_tot += short
        pcol = header_ref.index("p.value")
        ps = [float(r.split("\t")[pcol]) for r in rows]
        n = len(rows)
        secondary = 0.05 / n
        best_i = min(range(n), key=lambda i: ps[i])
        best_id = rows[best_i].split("\t")[header_ref.index("MarkerID")]

        with (OUT / f"merged_{tag}.tsv").open("w") as fh:
            fh.write("\t".join(header_ref) + "\n")
            for r in rows:
                fh.write(r + "\n")

        report["arms"][tag] = {
            "tested": n, "plink2_ceiling": cap_tot, "shortfall": short_tot,
            "secondary_bar": secondary, "min_p": min(ps), "best_variant": best_id,
            "hits_primary": sum(1 for v in ps if v < PRIMARY),
            "hits_secondary": sum(1 for v in ps if v < secondary),
        }
        r = report["arms"][tag]
        print(f"{tag}")
        print(f"  tested {n:>8,}  (ceiling {cap_tot:,}, shortfall {short_tot})")
        print(f"  secondary bar {secondary:.3e}   min p {min(ps):.3e} ({best_id})")
        print(f"  hits: 5e-8 {r['hits_primary']}   secondary {r['hits_secondary']}")
        print()

    (OUT / "02_merge.json").write_text(json.dumps(report, indent=2) + "\n")
    print("written: results/merged_*.tsv, 02_merge.json")


if __name__ == "__main__":
    main()
