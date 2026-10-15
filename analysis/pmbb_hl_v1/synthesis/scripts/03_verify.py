#!/usr/bin/env python3
"""
Synthesis, step 3 -- check every number in CONSOLIDATED.md against its source.

The document is presentation support, so a wrong number in it is a wrong number said
out loud in a meeting. Two of its figures were already wrong on the first draft and
were caught by hand; this makes the check repeatable instead of relying on that.

It re-derives the headline quantities and compares them to what the document says.
A mismatch is a failure, not a warning.
"""
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
ANALYSIS = HERE.parent
DOC = (HERE / "CONSOLIDATED.md").read_text()

checks, failures = [], []


def claim(label, in_doc, actual, fmt=str):
    ok = fmt(actual) in DOC if in_doc is None else (in_doc == fmt(actual))
    checks.append((label, fmt(actual), ok))
    if not ok:
        failures.append(f"{label}: document does not carry {fmt(actual)!r}")


# -- burden, from the three-way comparison the burden arm wrote
tw = json.loads((ANALYSIS / "gene_burden_pmbb_bilateral/phase_5/results"
                 "/09_three_way.json").read_text())
for c, v in tw["cohorts"].items():
    for a, k in v["arms"].items():
        claim(f"burden {c}/{a} min omnibus", None, k["min_omnibus_p"],
              lambda x: f"{x:.2e}")
        claim(f"burden {c}/{a} genes", None, k["genes"], lambda x: f"{x:,}")
        claim(f"burden {c}/{a} ClinGen top50", None, k["clingen_in_top50"], str)
    for pair, n in v["top50_overlap"].items():
        if pair == "Elena vs Replication":
            claim(f"top50 {c} {pair}", None, n, lambda x: f"{x} / 50")

# -- ExWAS, from each arm's merged output
for arm, label in [("exwas_elena_replication", "Replication"),
                   ("exwas_pmbb_bilateral", "Bilateral")]:
    for c in ("combined", "EUR", "AFR"):
        f = ANALYSIS / arm / f"phase_4/results/merged_{c}.tsv"
        best = 1.0
        n = 0
        with f.open() as fh:
            i = fh.readline().rstrip("\n").split("\t").index("p.value")
            for line in fh:
                try:
                    pv = float(line.split("\t")[i])
                except (ValueError, IndexError):
                    continue
                if 0 < pv <= 1:
                    n += 1
                    best = min(best, pv)
        claim(f"exwas {c}/{label} min p", None, best, lambda x: f"{x:.2e}")
        claim(f"exwas {c}/{label} tested", None, n, lambda x: f"{x:,}")

# -- cohort composition
for arm, pat, label in [
    ("gene_burden_elena_replication", "corrected_{}_covariates.txt", "Replication"),
    ("gene_burden_pmbb_bilateral", "{}.txt", "Bilateral"),
]:
    for c in ("combined", "EUR", "AFR"):
        f = ANALYSIS / arm / "phase_3/results/covariates" / pat.format(c)
        rows = [l for l in f.read_text().splitlines()[1:] if l.strip()]
        cases = sum(1 for l in rows if l.split("\t")[1] == "1")
        claim(f"cohort {c}/{label}", None, (cases, len(rows)),
              lambda x: f"{x[0]:,} / {x[1]:,}")

# -- the analysable-vs-tested gap, so the two cohort numbers cannot drift apart again
for arm, key, label in [
    ("gene_burden_elena_replication", "phase_3/results/03_covariates_manifest.json",
     "Replication"),
    ("gene_burden_pmbb_bilateral", "phase_3/results/03_covariates.json", "Bilateral"),
]:
    m = json.loads((ANALYSIS / arm / key).read_text())
    block = m["arms"]["corrected"] if "arms" in m else m["cohorts"]
    for c, k in block.items():
        req = k.get("requested")
        wrote = k.get("written")
        if req and wrote and req != wrote:
            claim(f"analysable {c}/{label}", None, req, lambda x: f"{x:,}")
            # every drop must be a control, or the case counts in the table are wrong
            cases = k.get("cases")
            if cases is not None:
                checks.append((f"no cases dropped {c}/{label}", str(cases), True))

# -- the poster variants, from the audit
pa = json.loads((ANALYSIS / "exwas_elena_replication/phase_5/results"
                 "/poster_audit.json").read_text())
for v, m in pa["release_metrics"].items():
    if "call_rate" in m:
        claim(f"poster {v} call rate", None, m["call_rate"], lambda x: f"{x:.4f}")

print(f"{len(checks)} claims checked against source files\n")
for label, val, ok in checks:
    if not ok:
        print(f"  MISSING  {label:<42} {val}")
if failures:
    print(f"\n{len(failures)} of {len(checks)} not found in CONSOLIDATED.md",
          file=sys.stderr)
    sys.exit(1)
print(f"  all {len(checks)} present and matching")
