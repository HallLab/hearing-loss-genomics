#!/usr/bin/env python3
"""
fetch_clingen_hl_genes.py — download the ClinGen gene-disease validity table
and extract the Hearing Loss Gene Curation Expert Panel (HL GCEP, affiliate
40007) curations.

Serves Cycle 2 open dependency #1 (`cycle_2/README.md` §5): fixing an
authoritative adjudicated hearing-loss gene list for the Q1 primary test.

The ClinGen table is a *living* resource, not a paper snapshot — Tshering KC
et al. 2025 (Genet Med 27:101397, PMID 39987489) exists precisely to show that
gene-disease validity classifications move over time. So every run stamps its
provenance (ClinGen's own FILE CREATED date + sha256 of the raw download) into
manifest.json. Cite the paper; pin the manifest.

Source: https://search.clinicalgenome.org/kb/gene-validity/download
        (single CSV covering every GCEP; we filter on the GCEP column)

Outputs, under --outdir (default cycle_2/data/clingen/):
  - clingen_gene_validity_all.csv          raw download, unmodified
  - clingen_hl_gcep.tsv                    HL GCEP rows only, tidy header
  - genes_definitive.txt                   one gene symbol per line
  - genes_definitive_strong.txt              "
  - genes_definitive_strong_moderate.txt     "
  - manifest.json                          provenance + counts
"""

import argparse
import csv
import hashlib
import json
import sys
import urllib.request
from collections import Counter
from pathlib import Path

CLINGEN_URL = "https://search.clinicalgenome.org/kb/gene-validity/download"
HL_GCEP = "Hearing Loss Gene Curation Expert Panel"

# ClinGen orders these by strength of evidence; anything below Moderate is not
# a candidate for a pre-specified gene set.
TIERS = {
    "definitive": ("Definitive",),
    "definitive_strong": ("Definitive", "Strong"),
    "definitive_strong_moderate": ("Definitive", "Strong", "Moderate"),
}
CLASS_ORDER = [
    "Definitive", "Strong", "Moderate",
    "Limited", "Disputed", "Refuted",
    "No Known Disease Relationship",
]


def download(url: str, dest: Path) -> None:
    """Fetch the ClinGen CSV to dest."""
    print(f"[fetch] {url}", file=sys.stderr)
    with urllib.request.urlopen(url, timeout=120) as resp:
        if resp.status != 200:
            sys.exit(f"ClinGen returned HTTP {resp.status}")
        dest.write_bytes(resp.read())
    print(f"[fetch] wrote {dest} ({dest.stat().st_size:,} bytes)", file=sys.stderr)


def parse_clingen_csv(path: Path) -> tuple[list[str], list[list[str]], str]:
    """Parse ClinGen's CSV, which is not a plain table.

    Layout: 3 banner lines, a '++++' separator, the real header, another
    '++++' separator, then data. Returns (header, rows, file_created).
    """
    with path.open(newline="", encoding="utf-8") as f:
        raw = list(csv.reader(f))

    file_created = ""
    header: list[str] | None = None
    rows: list[list[str]] = []

    for row in raw:
        if not row or not row[0]:
            continue
        first = row[0].strip()
        if first.startswith("+++"):
            continue
        if first.startswith("FILE CREATED:"):
            file_created = first.split(":", 1)[1].strip()
            continue
        if header is None:
            if first == "GENE SYMBOL":
                header = [c.strip() for c in row]
            continue  # still in the banner
        if len(row) == len(header):
            rows.append(row)

    if header is None:
        sys.exit(f"{path}: no 'GENE SYMBOL' header row — ClinGen format changed?")
    return header, rows, file_created


def repo_root() -> Path:
    """Walk up from this file to the repo root (marked by .git).

    Not a fixed number of parents[] — the script has already been moved once,
    and a hard-coded depth silently retargets --outdir when it moves again.
    """
    here = Path(__file__).resolve()
    for parent in here.parents:
        if (parent / ".git").exists():
            return parent
    return here.parent


def main() -> int:
    repo = repo_root()
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", type=Path, default=repo / "cycle_2/data/clingen",
                    help="output directory (default: cycle_2/data/clingen)")
    ap.add_argument("--url", default=CLINGEN_URL, help="ClinGen download URL")
    ap.add_argument("--gcep", default=HL_GCEP, help="GCEP name to filter on")
    ap.add_argument("--no-download", action="store_true",
                    help="reuse the raw CSV already in --outdir instead of refetching")
    args = ap.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    raw_csv = args.outdir / "clingen_gene_validity_all.csv"

    if args.no_download:
        if not raw_csv.exists():
            sys.exit(f"--no-download given but {raw_csv} does not exist")
        print(f"[fetch] reusing {raw_csv}", file=sys.stderr)
    else:
        download(args.url, raw_csv)

    sha256 = hashlib.sha256(raw_csv.read_bytes()).hexdigest()
    header, rows, file_created = parse_clingen_csv(raw_csv)

    col = {name: i for i, name in enumerate(header)}
    for required in ("GENE SYMBOL", "CLASSIFICATION", "GCEP"):
        if required not in col:
            sys.exit(f"missing expected column {required!r}; got {header}")

    hl = [r for r in rows if r[col["GCEP"]].strip() == args.gcep]
    if not hl:
        gceps = sorted({r[col["GCEP"]].strip() for r in rows})
        sys.exit(f"no rows for GCEP {args.gcep!r}. Available:\n  " + "\n  ".join(gceps))

    # HL GCEP rows, tidy — one gene-disease pair per line.
    pairs_tsv = args.outdir / "clingen_hl_gcep.tsv"
    with pairs_tsv.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(sorted(hl, key=lambda r: (r[col["GENE SYMBOL"]],
                                              r[col["DISEASE LABEL"]])))

    # A gene can hold several gene-disease pairs at different strengths
    # (e.g. syndromic Definitive + non-syndromic Limited). Take its *best*
    # classification, so a gene is never demoted by a weaker second pair.
    best: dict[str, int] = {}
    for r in hl:
        gene, cls = r[col["GENE SYMBOL"]], r[col["CLASSIFICATION"]].strip()
        rank = CLASS_ORDER.index(cls) if cls in CLASS_ORDER else len(CLASS_ORDER)
        best[gene] = min(best.get(gene, len(CLASS_ORDER)), rank)

    tier_files: dict[str, dict] = {}
    for tier, keep in TIERS.items():
        allowed = {CLASS_ORDER.index(c) for c in keep}
        genes = sorted(g for g, rank in best.items() if rank in allowed)
        out = args.outdir / f"genes_{tier}.txt"
        out.write_text("\n".join(genes) + "\n", encoding="utf-8")
        tier_files[tier] = {"file": out.name, "n_genes": len(genes)}

    pair_counts = Counter(r[col["CLASSIFICATION"]].strip() for r in hl)
    gene_counts = Counter(CLASS_ORDER[rank] if rank < len(CLASS_ORDER) else "Other"
                          for rank in best.values())
    moi_counts = Counter(r[col["MOI"]].strip() for r in hl) if "MOI" in col else {}

    manifest = {
        "source_url": args.url,
        "gcep": args.gcep,
        "clingen_file_created": file_created,
        "raw_csv": raw_csv.name,
        "raw_sha256": sha256,
        "n_curations_all_gceps": len(rows),
        "n_gene_disease_pairs": len(hl),
        "n_unique_genes": len(best),
        "pairs_by_classification": {c: pair_counts[c] for c in CLASS_ORDER if pair_counts[c]},
        "genes_by_best_classification": {c: gene_counts[c] for c in CLASS_ORDER if gene_counts[c]},
        "pairs_by_moi": dict(sorted(moi_counts.items())),
        "gene_lists": tier_files,
        "citation": ("Tshering KC et al. ClinGen recuration of hearing loss "
                     "associated-genes demonstrates significant changes in "
                     "gene-disease validity over time. Genet Med. 2025;27(5):101397. "
                     "PMID 39987489"),
    }
    (args.outdir / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    print(f"\nClinGen file created : {file_created}")
    print(f"GCEP                 : {args.gcep}")
    print(f"gene-disease pairs   : {len(hl)}   unique genes: {len(best)}")
    print(f"\n{'classification':<32}{'pairs':>7}{'genes':>7}")
    for c in CLASS_ORDER:
        if pair_counts[c]:
            print(f"  {c:<30}{pair_counts[c]:>7}{gene_counts[c]:>7}")
    print("\ngene lists:")
    for tier, meta in tier_files.items():
        print(f"  {meta['file']:<40}{meta['n_genes']:>4} genes")
    print(f"\nwrote {args.outdir}/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
