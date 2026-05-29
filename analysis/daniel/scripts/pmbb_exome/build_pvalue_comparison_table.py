#!/usr/bin/env python3
"""
build_pvalue_comparison_table.py — consolidate per-gene p-values across
all our runs (Phases 5, 6) + Daniel's preserved outputs (Phase 5 eq,
Phase 6 eq, binary 20PC, degHL linear 20PC) + paper published values.

Outputs:
  - pvalue_comparison_hl_genes.tsv  (rows = the 173 known HL gene set)
  - pvalue_comparison_wes_top.tsv   (rows = union of top 30 by p in each source)
  - markdown summary tables for the 6 paper-cited genes
"""

import argparse
import csv
import gzip
import sys
from pathlib import Path


PAPER_VALUES = {
    "TCOF1":   {"p": 5.20e-6, "beta": 0.798,  "fdr": 7.2e-4,  "table": 3},
    "ESRRB":   {"p": 7.00e-4, "beta": 0.148,  "fdr": 0.06,    "table": 3},
    "COL5A1":  {"p": 3.31e-5, "beta": 0.0586, "fdr": 0.0123,  "table": 4},
    "HMMR":    {"p": 6.89e-5, "beta": 0.0974, "fdr": 0.0128,  "table": 4},
    "RAPGEF3": {"p": 1.88e-4, "beta": 0.0731, "fdr": 0.0189,  "table": 4},
    "NNT":     {"p": 2.03e-4, "beta": 0.0507, "fdr": 0.0189,  "table": 4},
}


def parse_bins_csv(path: Path) -> dict[str, float]:
    """Phase 5 / Daniel Phase 5 bins.csv — wide format.

    Row 0 (after ID,) = gene names (one per column, may repeat across bins).
    Row 8 = 'logistic p-value' row.

    Returns dict {gene: best (lowest) p across columns sharing that gene name}.
    """
    open_fn = gzip.open if str(path).endswith(".gz") else open
    with open_fn(path, "rt") as f:
        reader = csv.reader(f)
        rows = []
        for i, row in enumerate(reader):
            rows.append(row)
            if i > 10:
                break
    header_row = rows[0]                # 'ID', '', 'ESPN', 'ESPN', 'LOC...', ...
    pval_row = rows[8]                  # 'logistic p-value', 'nan', '0.94...', ...
    out = {}
    for col_idx in range(2, len(header_row)):  # skip ID and the blank col
        gene = header_row[col_idx]
        if not gene or gene.lower().startswith("loc"):
            continue                    # skip LOCs (loki-introduced; not real genes)
        try:
            p = float(pval_row[col_idx])
        except (ValueError, IndexError):
            continue
        if gene not in out or p < out[gene]:
            out[gene] = p
    # also keep LOCs separately if any caller wants them
    return out


def parse_long_tsv(path: Path, gene_col: int, pval_col: int, beta_col: int = -1,
                   has_header: bool = True, gz: bool = None) -> dict[str, dict]:
    """Generic long-format TSV parser. Returns {gene: {p, beta?}}."""
    if gz is None:
        gz = str(path).endswith(".gz")
    open_fn = gzip.open if gz else open
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
            if beta_col >= 0:
                try:
                    entry["beta"] = float(parts[beta_col])
                except (ValueError, IndexError):
                    pass
            # keep the lowest p if duplicates
            if gene not in out or p < out[gene]["p"]:
                out[gene] = entry
    return out


def read_gene_set(path: Path) -> set[str]:
    open_fn = gzip.open if str(path).endswith(".gz") else open
    with open_fn(path, "rt") as f:
        return {line.strip() for line in f if line.strip()}


def fmt_p(p):
    if p is None:
        return "—"
    if isinstance(p, str):
        return p
    if p < 1e-3 or p > 1e3:
        return f"{p:.2e}"
    return f"{p:.4f}"


def write_tsv(rows: list[dict], path: Path, cols: list[str]):
    with open(path, "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--project-root", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    root = args.project_root

    print("Loading sources...")

    # 1. Our Phase 5 — HL genes, 4PC, binary logistic (biobin)
    ours_p5 = parse_bins_csv(root / "analysis/daniel/outputs/phase5/biobin/merged_maf.001_noRels_keepHLcases-bins.csv")
    print(f"  Our Phase 5 (HL, 4PC, binary): {len(ours_p5)} genes")

    # 2. Daniel Phase 5 equivalent — same params
    daniel_p5 = parse_bins_csv(root / "data/PMBB_Exome/biobin/merged_maf.001_noRels_keepHLcases-bins.csv.gz")
    print(f"  Daniel Phase 5 (HL, 4PC, binary): {len(daniel_p5)} genes")

    # 3. Our Phase 6 — WES, 4PC, binary logistic
    ours_p6 = parse_long_tsv(
        root / "analysis/daniel/outputs/phase6/results/all_chrom_meta_HL_needAud.txt",
        gene_col=0, pval_col=3, has_header=True, gz=False,
    )
    print(f"  Our Phase 6 (WES, 4PC, binary): {len(ours_p6)} genes")

    # 4. Daniel Phase 6 equivalent — Chromosome_hg38 | Gene | Logistic_regression_beta_p | BH-corrected_p
    daniel_p6 = parse_long_tsv(
        root / "data/PMBB_Exome/allGenes/HL_needAud/results_allChr_needAud_withBH.txt.gz",
        gene_col=1, pval_col=2, has_header=True,
    )
    print(f"  Daniel Phase 6 (WES, 4PC, binary): {len(daniel_p6)} genes")

    # 5. Daniel binary 20PC STable — no header. Gene | Beta | SE | P | Carriers | Case_carriers
    daniel_binary_20pc = parse_long_tsv(
        root / "data/PMBB_Exome/allGenes/20PCs/binary/results/allChr_STable_binaryHL.txt.gz",
        gene_col=0, pval_col=3, beta_col=1, has_header=False,
    )
    print(f"  Daniel binary 20PC (WES, logistic): {len(daniel_binary_20pc)} genes")

    # 6. Daniel degHL linear 20PC STable — WITH header. Gene | Beta | SE | P | Carriers | Case_carriers
    daniel_linear_20pc = parse_long_tsv(
        root / "data/PMBB_Exome/allGenes/20PCs/degreeHL/results/allChr_STable_degHL.txt.gz",
        gene_col=0, pval_col=3, beta_col=1, has_header=True,
    )
    print(f"  Daniel degHL linear 20PC (WES): {len(daniel_linear_20pc)} genes")

    # HL gene set
    hl_genes = read_gene_set(root / "data/PMBB_Exome/all_genes_including_ShadisList.txt.gz")
    print(f"  HL gene set: {len(hl_genes)} genes")

    # ------------------------------------------------------------------
    # Build HL-gene table (rows = 173 known HL genes)
    # ------------------------------------------------------------------
    print("\nBuilding HL genes table...")
    hl_rows = []
    for gene in sorted(hl_genes):
        row = {
            "Gene": gene,
            "Ours_HL_4PC_binary_p": fmt_p(ours_p5.get(gene)),
            "Daniel_HL_4PC_binary_p": fmt_p(daniel_p5.get(gene)),
            "Ours_WES_4PC_binary_p": fmt_p(ours_p6.get(gene, {}).get("p") if isinstance(ours_p6.get(gene), dict) else ours_p6.get(gene)),
            "Daniel_WES_4PC_binary_p": fmt_p(daniel_p6.get(gene, {}).get("p")),
            "Daniel_WES_20PC_binary_p": fmt_p(daniel_binary_20pc.get(gene, {}).get("p")),
            "Daniel_WES_20PC_linear_p": fmt_p(daniel_linear_20pc.get(gene, {}).get("p")),
            "Daniel_WES_20PC_linear_beta": fmt_p(daniel_linear_20pc.get(gene, {}).get("beta")),
            "Paper_p":   fmt_p(PAPER_VALUES.get(gene, {}).get("p")),
            "Paper_beta": fmt_p(PAPER_VALUES.get(gene, {}).get("beta")),
        }
        hl_rows.append(row)

    cols = ["Gene", "Ours_HL_4PC_binary_p", "Daniel_HL_4PC_binary_p",
            "Ours_WES_4PC_binary_p", "Daniel_WES_4PC_binary_p",
            "Daniel_WES_20PC_binary_p", "Daniel_WES_20PC_linear_p",
            "Daniel_WES_20PC_linear_beta", "Paper_p", "Paper_beta"]
    write_tsv(hl_rows, args.out_dir / "pvalue_comparison_hl_genes.tsv", cols)
    print(f"  Wrote pvalue_comparison_hl_genes.tsv: {len(hl_rows)} rows")

    # ------------------------------------------------------------------
    # WES top genes — union of top 30 by p in each Daniel source + 6 paper
    # ------------------------------------------------------------------
    print("\nBuilding WES top table...")
    def top_genes(d: dict, n: int = 30) -> list[str]:
        return [g for g, _ in sorted(d.items(), key=lambda kv: kv[1]["p"] if isinstance(kv[1], dict) else kv[1])[:n]]
    top_set = set(PAPER_VALUES.keys())
    top_set.update(top_genes(daniel_p6, 20))
    top_set.update(top_genes(daniel_binary_20pc, 20))
    top_set.update(top_genes(daniel_linear_20pc, 30))
    wes_rows = []
    for gene in sorted(top_set):
        row = {
            "Gene": gene,
            "In_HL_set": "✓" if gene in hl_genes else "",
            "Ours_WES_4PC_binary_p": fmt_p(ours_p6.get(gene, {}).get("p") if isinstance(ours_p6.get(gene), dict) else ours_p6.get(gene)),
            "Daniel_WES_4PC_binary_p": fmt_p(daniel_p6.get(gene, {}).get("p")),
            "Daniel_WES_20PC_binary_p": fmt_p(daniel_binary_20pc.get(gene, {}).get("p")),
            "Daniel_WES_20PC_linear_p": fmt_p(daniel_linear_20pc.get(gene, {}).get("p")),
            "Daniel_WES_20PC_linear_beta": fmt_p(daniel_linear_20pc.get(gene, {}).get("beta")),
            "Paper_p":   fmt_p(PAPER_VALUES.get(gene, {}).get("p")),
            "Paper_beta": fmt_p(PAPER_VALUES.get(gene, {}).get("beta")),
        }
        wes_rows.append(row)

    cols_wes = ["Gene", "In_HL_set", "Ours_WES_4PC_binary_p", "Daniel_WES_4PC_binary_p",
                "Daniel_WES_20PC_binary_p", "Daniel_WES_20PC_linear_p",
                "Daniel_WES_20PC_linear_beta", "Paper_p", "Paper_beta"]
    write_tsv(wes_rows, args.out_dir / "pvalue_comparison_wes_top.tsv", cols_wes)
    print(f"  Wrote pvalue_comparison_wes_top.tsv: {len(wes_rows)} rows")

    # ------------------------------------------------------------------
    # Console summary for the 6 paper genes
    # ------------------------------------------------------------------
    print("\n" + "=" * 80)
    print("Paper-cited genes — full comparison")
    print("=" * 80)
    for gene in ["TCOF1", "ESRRB", "COL5A1", "HMMR", "RAPGEF3", "NNT"]:
        print(f"\n{gene}:")
        print(f"  Ours   HL,   4PC, binary  = {fmt_p(ours_p5.get(gene))}")
        print(f"  Daniel HL,   4PC, binary  = {fmt_p(daniel_p5.get(gene))}")
        ours_wes = ours_p6.get(gene)
        print(f"  Ours   WES,  4PC, binary  = {fmt_p(ours_wes['p'] if isinstance(ours_wes, dict) else ours_wes)}")
        print(f"  Daniel WES,  4PC, binary  = {fmt_p(daniel_p6.get(gene, {}).get('p'))}")
        print(f"  Daniel WES, 20PC, binary  = {fmt_p(daniel_binary_20pc.get(gene, {}).get('p'))}")
        print(f"  Daniel WES, 20PC, linear  = {fmt_p(daniel_linear_20pc.get(gene, {}).get('p'))} (β={fmt_p(daniel_linear_20pc.get(gene, {}).get('beta'))})")
        print(f"  Paper  WES, 20PC, linear  = {fmt_p(PAPER_VALUES[gene]['p'])} (β={fmt_p(PAPER_VALUES[gene]['beta'])})")

    print("\nDone.")


if __name__ == "__main__":
    sys.exit(main())
