#!/usr/bin/env python3
"""
build_all_genes_pvalues_csv.py — comprehensive CSV with all WES genes and
their p-values across Phases 5, 6, 7 (ours vs Daniel) + Paper published.

Master gene universe = Phase 7 STable (18,547 genes — Daniel's preserved
linear-regression output covers every gene tested exome-wide).

Inputs:
- data/PMBB_Exome/all_genes_including_ShadisList.txt.gz       (179 HL genes)
- analysis/daniel/outputs/phase5/biobin/...-bins.csv          (our Phase 5)
- data/PMBB_Exome/biobin/...-bins.csv.gz                      (Daniel Phase 5)
- analysis/daniel/outputs/phase6/results/all_chrom_meta_HL_needAud.txt (our Phase 6)
- data/PMBB_Exome/allGenes/HL_needAud/results_allChr_needAud_withBH.txt.gz (Daniel Phase 6)
- data/PMBB_Exome/allGenes/20PCs/degreeHL/results/allChr_STable_degHL.txt.gz (Daniel Phase 7)

Output:
- analysis/daniel/outputs/phase8/comparison/all_genes_pvalues.csv
"""

import argparse
import csv
import gzip
import sys
from pathlib import Path


PAPER = {
    "TCOF1":   {"p": 5.20e-6, "beta": 0.798,  "fdr": 7.2e-4,  "table": 3},
    "ESRRB":   {"p": 7.00e-4, "beta": 0.148,  "fdr": 0.06,    "table": 3},
    "COL5A1":  {"p": 3.31e-5, "beta": 0.0586, "fdr": 0.0123,  "table": 4},
    "HMMR":    {"p": 6.89e-5, "beta": 0.0974, "fdr": 0.0128,  "table": 4},
    "RAPGEF3": {"p": 1.88e-4, "beta": 0.0731, "fdr": 0.0189,  "table": 4},
    "NNT":     {"p": 2.03e-4, "beta": 0.0507, "fdr": 0.0189,  "table": 4},
}


def parse_bins_csv(path: Path) -> dict[str, float]:
    """Phase 5 bins.csv: wide format. Row 0 = gene names, row 8 = logistic p."""
    open_fn = gzip.open if str(path).endswith(".gz") else open
    with open_fn(path, "rt") as f:
        reader = csv.reader(f)
        rows = []
        for i, row in enumerate(reader):
            rows.append(row)
            if i > 10:
                break
    header_row = rows[0]
    pval_row = rows[8]
    out = {}
    for col_idx in range(2, len(header_row)):
        gene = header_row[col_idx]
        if not gene or gene.lower().startswith("loc"):
            continue
        try:
            p = float(pval_row[col_idx])
        except (ValueError, IndexError):
            continue
        if gene not in out or p < out[gene]:
            out[gene] = p
    return out


def parse_long_tsv(path: Path, gene_col: int, pval_col: int,
                   beta_col: int = -1, carriers_col: int = -1,
                   case_carriers_col: int = -1, has_header: bool = True) -> dict[str, dict]:
    open_fn = gzip.open if str(path).endswith(".gz") else open
    out = {}
    with open_fn(path, "rt") as f:
        if has_header:
            next(f)
        for line in f:
            parts = line.rstrip("\n").split("\t")
            try:
                gene = parts[gene_col].strip()
                p = float(parts[pval_col])
            except (ValueError, IndexError):
                continue
            entry = {"p": p}
            for name, idx in [("beta", beta_col), ("carriers", carriers_col),
                              ("case_carriers", case_carriers_col)]:
                if idx >= 0:
                    try:
                        entry[name] = float(parts[idx])
                    except (ValueError, IndexError):
                        pass
            if gene not in out or p < out[gene]["p"]:
                out[gene] = entry
    return out


def read_gene_set(path: Path) -> set[str]:
    open_fn = gzip.open if str(path).endswith(".gz") else open
    with open_fn(path, "rt") as f:
        return {line.strip() for line in f if line.strip()}


def fmt_p(x):
    if x is None or x == "":
        return ""
    if isinstance(x, dict):
        x = x.get("p")
    if x is None:
        return ""
    if x == 0:
        return "0"
    if abs(x) < 1e-3 or abs(x) > 1e3:
        return f"{x:.4e}"
    return f"{x:.6f}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--project-root", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()
    args.out.parent.mkdir(parents=True, exist_ok=True)

    root = args.project_root
    print("Loading all data sources...")

    hl_genes = read_gene_set(root / "data/PMBB_Exome/all_genes_including_ShadisList.txt.gz")
    print(f"  HL gene set: {len(hl_genes)} genes")

    p5_ours = parse_bins_csv(root / "analysis/daniel/outputs/phase5/biobin/merged_maf.001_noRels_keepHLcases-bins.csv")
    p5_daniel = parse_bins_csv(root / "data/PMBB_Exome/biobin/merged_maf.001_noRels_keepHLcases-bins.csv.gz")
    print(f"  Phase 5 ours: {len(p5_ours)} genes; Daniel: {len(p5_daniel)} genes")

    p6_ours = parse_long_tsv(
        root / "analysis/daniel/outputs/phase6/results/all_chrom_meta_HL_needAud.txt",
        gene_col=0, pval_col=3, has_header=True,
    )
    p6_daniel = parse_long_tsv(
        root / "data/PMBB_Exome/allGenes/HL_needAud/results_allChr_needAud_withBH.txt.gz",
        gene_col=1, pval_col=2, has_header=True,
    )
    print(f"  Phase 6 ours: {len(p6_ours)} genes; Daniel: {len(p6_daniel)} genes")

    p7_daniel = parse_long_tsv(
        root / "data/PMBB_Exome/allGenes/20PCs/degreeHL/results/allChr_STable_degHL.txt.gz",
        gene_col=0, pval_col=3, beta_col=1, carriers_col=4, case_carriers_col=5,
        has_header=True,
    )
    print(f"  Phase 7 Daniel STable: {len(p7_daniel)} genes")

    # Master universe = union of all genes seen across all phases
    all_genes = set(p7_daniel.keys()) | set(p6_ours.keys()) | set(p6_daniel.keys()) | set(p5_ours.keys()) | set(p5_daniel.keys())
    print(f"\n  Total unique genes in master CSV: {len(all_genes)}")

    # Sort: HL genes first (alphabetical), then non-HL (alphabetical)
    sorted_genes = sorted(all_genes, key=lambda g: (g not in hl_genes, g))

    # Build rows
    rows = []
    for gene in sorted_genes:
        p7 = p7_daniel.get(gene, {})
        is_loc_artifact = gene.startswith("LOC") or gene.startswith("LINC") or gene.startswith("MIR")
        row = {
            "Gene": gene,
            "Is_real_gene": "no" if is_loc_artifact else "yes",
            "In_HL_gene_set": "yes" if gene in hl_genes else "no",
            "In_paper_table_3_or_4": "yes" if gene in PAPER else "no",
            "Phase5_HL_4PC_binary_p_ours": fmt_p(p5_ours.get(gene)),
            "Phase5_HL_4PC_binary_p_daniel": fmt_p(p5_daniel.get(gene)),
            "Phase6_WES_4PC_binary_p_ours": fmt_p(p6_ours.get(gene, {}).get("p") if isinstance(p6_ours.get(gene), dict) else None),
            "Phase6_WES_4PC_binary_p_daniel": fmt_p(p6_daniel.get(gene, {}).get("p")),
            "Phase7_WES_20PC_linear_p_daniel_stable": fmt_p(p7.get("p")),
            "Phase7_WES_20PC_linear_beta_daniel_stable": fmt_p(p7.get("beta")),
            "Phase7_carriers_daniel_stable": int(p7["carriers"]) if "carriers" in p7 else "",
            "Phase7_case_carriers_daniel_stable": int(p7["case_carriers"]) if "case_carriers" in p7 else "",
            "Paper_p_published": fmt_p(PAPER.get(gene, {}).get("p")),
            "Paper_beta_published": fmt_p(PAPER.get(gene, {}).get("beta")),
            "Paper_fdr_published": fmt_p(PAPER.get(gene, {}).get("fdr")),
            "Joe_PMBB_v1_p": "",
            "Joe_PMBB_v1_notes": "",
        }
        # Handle p6_ours possibly being a dict or a flat float
        v = p6_ours.get(gene)
        if isinstance(v, dict):
            row["Phase6_WES_4PC_binary_p_ours"] = fmt_p(v.get("p"))
        elif isinstance(v, (int, float)):
            row["Phase6_WES_4PC_binary_p_ours"] = fmt_p(v)

        rows.append(row)

    cols = list(rows[0].keys())

    print(f"\nWriting {args.out}...")
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)

    # Summary
    hl_with_p5 = sum(1 for r in rows if r["In_HL_gene_set"] == "yes" and r["Phase5_HL_4PC_binary_p_ours"])
    wes_with_p6 = sum(1 for r in rows if r["Phase6_WES_4PC_binary_p_ours"])
    wes_with_p7 = sum(1 for r in rows if r["Phase7_WES_20PC_linear_p_daniel_stable"])
    paper_genes = sum(1 for r in rows if r["In_paper_table_3_or_4"] == "yes")
    real_genes = sum(1 for r in rows if r["Is_real_gene"] == "yes")
    loc_artifacts = sum(1 for r in rows if r["Is_real_gene"] == "no")

    print(f"\n=== Summary ===")
    print(f"  Total rows: {len(rows)}")
    print(f"    Real genes (Is_real_gene=yes): {real_genes}")
    print(f"    LOC/LINC/MIR artifacts (newer LOKI): {loc_artifacts}")
    print(f"  HL gene set genes: {sum(1 for r in rows if r['In_HL_gene_set'] == 'yes')}")
    print(f"  Paper-cited genes: {paper_genes}")
    print(f"  HL genes with Phase 5 data (ours): {hl_with_p5}")
    print(f"  Genes with Phase 6 data (ours): {wes_with_p6}")
    print(f"  Genes with Phase 7 data (Daniel STable): {wes_with_p7}")

    # Also save a "real-genes-only" subset
    real_out = args.out.parent / "all_real_genes_pvalues.csv"
    print(f"\nAlso writing real-genes-only subset → {real_out}")
    with open(real_out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows([r for r in rows if r["Is_real_gene"] == "yes"])

    print(f"\n✓ Done. Two files:")
    print(f"  - {args.out} (all {len(rows)} entries including LOC/LINC artifacts)")
    print(f"  - {real_out} (only {real_genes} real genes — use this for normal analysis)")


if __name__ == "__main__":
    sys.exit(main())
