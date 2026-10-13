#!/usr/bin/env python3
"""
Synthesis, step 1 -- pull every reportable number out of every arm's own results.

Nothing here is retyped from a document. Each number is read from the file the arm
wrote, so the consolidated table is a measurement of the repository rather than a
transcription of earlier summaries.

SIX CELLS: {rare-variant burden, common-variant ExWAS} x {Elena, Replication, Bilateral}.
Elena never ran a bilateral analysis, so those two cells are ours alone and are labelled
as new work rather than as replication.

A COMPARABILITY WARNING THAT THE TABLE CANNOT CARRY ON ITS OWN.
The burden arms and the ExWAS arms do not report the same quantity:

  burden   one omnibus p per GENE, Cauchy-combined across up to 9 mask x MAF cells
  ExWAS    one p per VARIANT, no combination

and within the burden arms, Elena's delivered tables report per-cell p-values with her
own Bonferroni column, while ours report a Cauchy omnibus per gene. Comparing her
minimum cell p against our minimum omnibus p would flatter us, because a minimum over
9x more tests is smaller by construction. So both are collected and both are shown.
"""
import csv
import json
import math
from pathlib import Path

ROOT = Path("/project/hall/analysis/hearing-loss-genomics/analysis")
OUT = ROOT / "synthesis/results"
COHORTS = ["combined", "EUR", "AFR"]

# Cohort composition, from each arm's own phase 3 covariate files.
COV = {
    "Replication": (ROOT / "gene_burden_elena_replication/phase_3/results/covariates",
                    "corrected_{}_covariates.txt"),
    "Bilateral": (ROOT / "gene_burden_pmbb_bilateral/phase_3/results/covariates",
                  "{}.txt"),
}


def cohort_counts(arm, cohort):
    base, pat = COV[arm]
    f = base / pat.format(cohort)
    if not f.exists():
        return None
    cases = total = 0
    for line in f.read_text().splitlines()[1:]:
        if not line.strip():
            continue
        total += 1
        if line.split("\t")[1] == "1":
            cases += 1
    return {"n": total, "cases": cases, "controls": total - cases}


def elena_cohort(cohort):
    """Her N comes from the N_case / N_ctrl columns of her own ExWAS output."""
    f = (ROOT / "elena/HL_only_ExWAS/"
         f"HL_ExWAS_FullyAnnotated_AllVariants_{cohort}.csv")
    with f.open() as fh:
        row = next(csv.DictReader(fh))
    ca, ct = int(float(row["N_case"])), int(float(row["N_ctrl"]))
    return {"n": ca + ct, "cases": ca, "controls": ct}


# ---------------------------------------------------------------- ExWAS

def exwas_ours(arm_dir, cohort):
    f = ROOT / arm_dir / f"phase_4/results/merged_{cohort}.tsv"
    ps = []
    with f.open() as fh:
        head = fh.readline().rstrip("\n").split("\t")
        i = head.index("p.value")
        for line in fh:
            try:
                v = float(line.split("\t")[i])
            except (ValueError, IndexError):
                continue
            if 0 < v <= 1:
                ps.append(v)
    ps.sort()
    n = len(ps)
    return {"tested": n, "min_p": ps[0], "bar_primary": 5e-8,
            "bar_secondary": 0.05 / n,
            "hits_primary": sum(1 for v in ps if v < 5e-8),
            "hits_secondary": sum(1 for v in ps if v < 0.05 / n)}


def exwas_elena(cohort):
    f = (ROOT / "elena/HL_only_ExWAS/"
         f"HL_ExWAS_FullyAnnotated_AllVariants_{cohort}.csv")
    ps, n_tests = [], set()
    with f.open() as fh:
        for row in csv.DictReader(fh):
            try:
                v = float(row["p.value"])
            except (ValueError, TypeError):
                continue
            if 0 < v <= 1:
                ps.append(v)
            if row.get("n_tests_cohort"):
                n_tests.add(row["n_tests_cohort"])
    ps.sort()
    n = len(ps)
    declared = int(sorted(n_tests)[0]) if n_tests else n
    return {"tested": n, "declared_n_tests": declared, "min_p": ps[0],
            "bar_primary": 5e-8, "bar_secondary": 0.05 / declared,
            "hits_primary": sum(1 for v in ps if v < 5e-8),
            "hits_secondary": sum(1 for v in ps if v < 0.05 / declared)}


# ---------------------------------------------------------------- burden

def burden_ours(arm_dir, cohort, sub):
    f = ROOT / arm_dir / sub / f"omnibus_{cohort}.tsv"
    rows = []
    with f.open() as fh:
        head = fh.readline().rstrip("\n").split("\t")
        i_p, i_r = head.index("omnibus_p"), head.index("Region")
        i_bar = head.index("bonferroni_bar") if "bonferroni_bar" in head else None
        i_cell = head.index("best_cell_p") if "best_cell_p" in head else None
        for line in fh:
            f_ = line.rstrip("\n").split("\t")
            try:
                p = float(f_[i_p])
            except (ValueError, IndexError):
                continue
            rows.append((p, f_[i_r],
                         float(f_[i_bar]) if i_bar is not None else None,
                         float(f_[i_cell]) if i_cell is not None else None))
    rows.sort()
    genes = len(rows)
    bar = rows[0][2] if rows[0][2] else 0.05 / genes
    best_cell = min(r[3] for r in rows if r[3] is not None) if rows[0][3] else None
    return {"genes": genes, "min_omnibus_p": rows[0][0], "best_gene": rows[0][1],
            "bar": bar, "hits": sum(1 for r in rows if r[0] < bar),
            "min_cell_p": best_cell}


def burden_elena(cohort):
    """Her delivered per-cell tables. Reports cell-level p, not an omnibus."""
    f = (ROOT / "elena/HL_only_rarevariant/"
         f"Newmasks_HL_AllTests_{cohort}_corrected.csv")
    if not f.exists():
        return None
    best = (2.0, None)
    genes, tests, sig = set(), 0, 0
    bar = None
    with f.open() as fh:
        for row in csv.DictReader(fh):
            try:
                p = float(row["Pvalue"])
            except (ValueError, TypeError):
                continue
            tests += 1
            genes.add(row.get("VEP_SYMBOL") or row["Region"])
            if 0 < p < best[0]:
                best = (p, row.get("VEP_SYMBOL") or row["Region"])
            try:
                pb = float(row["Pvalue_Bonferroni"])
                if bar is None and p > 0:
                    bar = 0.05 * p / pb if pb > 0 else None
            except (ValueError, TypeError, ZeroDivisionError):
                pass
            if str(row.get("Pvalue_Bonferroni_Significant", "")).strip().upper() \
                    in ("TRUE", "YES", "1"):
                sig += 1
    return {"genes": len(genes), "cell_tests": tests, "min_cell_p": best[0],
            "best_gene": best[1], "bar": bar, "hits": sig}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    data = {"cohorts": {}, "exwas": {}, "burden": {}}

    for c in COHORTS:
        data["cohorts"][c] = {
            "Elena": elena_cohort(c),
            "Replication": cohort_counts("Replication", c),
            "Bilateral": cohort_counts("Bilateral", c),
        }

    for c in COHORTS:
        data["exwas"][c] = {
            "Elena": exwas_elena(c),
            "Replication": exwas_ours("exwas_elena_replication", c),
            "Bilateral": exwas_ours("exwas_pmbb_bilateral", c),
        }
        data["burden"][c] = {
            "Elena": burden_elena(c),
            "Replication": burden_ours("gene_burden_elena_replication", c,
                                       "phase_4/results"),
            "Bilateral": burden_ours("gene_burden_pmbb_bilateral", c,
                                     "phase_5/results"),
        }

    (OUT / "consolidated.json").write_text(json.dumps(data, indent=2) + "\n")

    print("COHORTS\n")
    print(f"  {'':12}{'Elena':>26}{'Replication':>26}{'Bilateral':>26}")
    for c in COHORTS:
        cells = []
        for a in ("Elena", "Replication", "Bilateral"):
            k = data["cohorts"][c][a]
            cells.append(f"{k['cases']:,}/{k['n']:,}" if k else "--")
        print(f"  {c:<12}" + "".join(f"{v:>26}" for v in cells))

    print("\nEXWAS  (min p / tested)\n")
    for c in COHORTS:
        row = []
        for a in ("Elena", "Replication", "Bilateral"):
            k = data["exwas"][c][a]
            row.append(f"{k['min_p']:.2e} / {k['tested']:,}" if k else "--")
        print(f"  {c:<12}" + "".join(f"{v:>26}" for v in row))

    print("\nBURDEN  (min p / units)\n")
    for c in COHORTS:
        row = []
        for a in ("Elena", "Replication", "Bilateral"):
            k = data["burden"][c][a]
            if not k:
                row.append("--")
            elif "min_omnibus_p" in k:
                row.append(f"{k['min_omnibus_p']:.2e} / {k['genes']:,}g")
            else:
                row.append(f"{k['min_cell_p']:.2e} / {k['cell_tests']:,}t")
        print(f"  {c:<12}" + "".join(f"{v:>26}" for v in row))

    print(f"\nwritten: {(OUT / 'consolidated.json').relative_to(ROOT)}")


if __name__ == "__main__":
    main()
