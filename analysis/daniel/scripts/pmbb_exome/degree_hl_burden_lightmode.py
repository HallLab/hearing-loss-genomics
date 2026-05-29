#!/usr/bin/env python3
"""
degree_hl_burden_lightmode.py — Chapter 1 · Phase 7 light-mode analysis.

Replicates Tables 3 & 4 of Hui et al. 2023 (PLOS Genetics) using Daniel's
preserved per-gene supplementary table (linear regression of degreeHL on
burden count with 20 PCs).

Paper filters:
- Table 3 (known HL genes): restrict to the 173 HL gene set, case_carriers > 0.
- Table 4 (novel genes):    restrict to non-HL genes, case_carriers > 25
                             (paper reports 373 such genes).

BH-FDR computed independently within each subset.
"""

import argparse
import gzip
import sys
from pathlib import Path


PAPER_TABLE_3 = {
    "TCOF1": {"beta": 0.798,  "se": 0.175,   "p": 5.20e-6, "carriers": 8,    "fdr": 7.2e-4},
    "ESRRB": {"beta": 0.148,  "se": 0.0434,  "p": 7.00e-4, "carriers": 129,  "fdr": 0.06},
}
PAPER_TABLE_4 = {
    "COL5A1":  {"beta": 0.0586, "se": 0.0141, "p": 3.31e-5, "carriers": 1255, "fdr": 0.0123},
    "HMMR":    {"beta": 0.0974, "se": 0.0245, "p": 6.89e-5, "carriers": 407,  "fdr": 0.0128},
    "RAPGEF3": {"beta": 0.0731, "se": 0.0196, "p": 1.88e-4, "carriers": 647,  "fdr": 0.0189},
    "NNT":     {"beta": 0.0507, "se": 0.0136, "p": 2.03e-4, "carriers": 1013, "fdr": 0.0189},
}


def read_hl_genes(path: Path) -> set[str]:
    with gzip.open(path, "rt") as f:
        return {line.strip() for line in f if line.strip()}


def read_stable(path: Path):
    rows = []
    with gzip.open(path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        for line in f:
            parts = line.rstrip("\n").split("\t")
            try:
                rows.append({
                    "gene": parts[0],
                    "beta": float(parts[1]),
                    "se": float(parts[2]),
                    "p": float(parts[3]),
                    "carriers": int(parts[4]),
                    "case_carriers": int(parts[5]),
                })
            except (ValueError, IndexError):
                continue
    return header, rows


def bh_fdr(pvals: list[float]) -> list[float]:
    """Benjamini-Hochberg step-up FDR (monotonic non-decreasing in sorted p)."""
    n = len(pvals)
    indexed = sorted(enumerate(pvals), key=lambda x: x[1])
    raw = [p * n / (rank + 1) for rank, (_, p) in enumerate(indexed)]
    monotone = [0.0] * n
    cur_min = 1.0
    for i in range(n - 1, -1, -1):
        cur_min = min(cur_min, raw[i])
        monotone[i] = cur_min
    out = [0.0] * n
    for rank, (orig_idx, _) in enumerate(indexed):
        out[orig_idx] = monotone[rank]
    return out


def annotate_fdr(rows: list[dict]) -> list[dict]:
    if not rows:
        return rows
    sorted_rows = sorted(rows, key=lambda r: r["p"])
    fdrs = bh_fdr([r["p"] for r in sorted_rows])
    for r, fdr in zip(sorted_rows, fdrs):
        r["fdr"] = fdr
    return sorted_rows


def write_tsv(rows: list[dict], path: Path, header: list[str]):
    with open(path, "w") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join(str(r[k]) for k in header) + "\n")


def fmt(x: float) -> str:
    if x is None:
        return "NA"
    if abs(x) < 1e-3 or abs(x) > 1e3:
        return f"{x:.3e}"
    return f"{x:.4f}"


def compare(label: str, our: list[dict], paper: dict):
    print(f"\n=== {label} ===")
    by_gene = {r["gene"]: r for r in our}
    cols = ["Gene", "Paper β/p/carr", "Ours β/p/carr/FDR", "Δβ rel", "FDR<0.05?"]
    widths = [10, 28, 36, 10, 10]
    print("  ".join(c.ljust(w) for c, w in zip(cols, widths)))
    print("  ".join("-" * w for w in widths))
    for gene, p in paper.items():
        if gene not in by_gene:
            row = ["MISSING from STable subset", "—", "—", "—"]
            print(f"{gene:<10}  {row[0]:<28}  —")
            continue
        r = by_gene[gene]
        paper_str = f"{p['beta']:.4f} / {fmt(p['p'])} / {p['carriers']}"
        ours_str = f"{r['beta']:.4f} / {fmt(r['p'])} / {r['carriers']} / {fmt(r['fdr'])}"
        delta = 100.0 * (r["beta"] - p["beta"]) / p["beta"] if p["beta"] else 0.0
        ok = "✓" if r["fdr"] < 0.05 else "—"
        print("  ".join([
            gene.ljust(widths[0]),
            paper_str.ljust(widths[1]),
            ours_str.ljust(widths[2]),
            f"{delta:+.1f}%".ljust(widths[3]),
            ok.ljust(widths[4]),
        ]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stable", required=True, type=Path)
    ap.add_argument("--hl-genes", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    hl = read_hl_genes(args.hl_genes)
    header, rows = read_stable(args.stable)
    print(f"Loaded STable: {len(rows)} genes")
    print(f"Loaded HL gene list: {len(hl)} genes")

    table3 = [r for r in rows if r["gene"] in hl and r["case_carriers"] > 0]
    table4 = [r for r in rows if r["gene"] not in hl and r["case_carriers"] > 25]
    table3 = annotate_fdr(table3)
    table4 = annotate_fdr(table4)
    print(f"Table 3 (known HL ∩ case_carriers>0): {len(table3)} genes (paper: 173)")
    print(f"Table 4 (novel ∩ case_carriers>25):  {len(table4)} genes (paper: 373)")

    out_cols = ["gene", "beta", "se", "p", "carriers", "case_carriers", "fdr"]
    write_tsv(table3, args.out_dir / "table3_known_hl_genes.tsv", out_cols)
    write_tsv(table4, args.out_dir / "table4_novel_genes.tsv", out_cols)

    print("\n=== Table 3 — top 10 by p (known HL genes) ===")
    for r in table3[:10]:
        print(f"  {r['gene']:<12} β={r['beta']:+.4f}  SE={r['se']:.4f}  "
              f"p={fmt(r['p'])}  FDR={fmt(r['fdr'])}  carriers={r['carriers']}  case_carriers={r['case_carriers']}")

    print("\n=== Table 4 — top 10 by p (novel genes, case_carriers>25) ===")
    for r in table4[:10]:
        print(f"  {r['gene']:<12} β={r['beta']:+.4f}  SE={r['se']:.4f}  "
              f"p={fmt(r['p'])}  FDR={fmt(r['fdr'])}  carriers={r['carriers']}  case_carriers={r['case_carriers']}")

    compare("Paper Table 3 vs our Table 3", table3, PAPER_TABLE_3)
    compare("Paper Table 4 vs our Table 4", table4, PAPER_TABLE_4)

    sig3 = [r for r in table3 if r["fdr"] < 0.05]
    sig4 = [r for r in table4 if r["fdr"] < 0.05]
    print(f"\n=== Significance summary (FDR < 0.05) ===")
    print(f"  Known HL genes (Table 3): {len(sig3)} significant — {[r['gene'] for r in sig3]}")
    print(f"  Novel genes    (Table 4): {len(sig4)} significant — {[r['gene'] for r in sig4]}")

    paper_t3 = set(PAPER_TABLE_3.keys())
    paper_t4 = set(PAPER_TABLE_4.keys())
    our_t3 = {r["gene"] for r in sig3}
    our_t4 = {r["gene"] for r in sig4}
    print(f"\n  Paper T3 ∩ Ours T3: {sorted(paper_t3 & our_t3)} ({len(paper_t3 & our_t3)}/{len(paper_t3)} expected)")
    print(f"  Paper T4 ∩ Ours T4: {sorted(paper_t4 & our_t4)} ({len(paper_t4 & our_t4)}/{len(paper_t4)} expected)")


if __name__ == "__main__":
    sys.exit(main())
