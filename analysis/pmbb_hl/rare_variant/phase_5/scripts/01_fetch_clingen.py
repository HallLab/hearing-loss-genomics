#!/usr/bin/env python3
"""
Phase 5, step 5.5a -- fetch the ClinGen hearing-loss gene list.

Fetched here rather than borrowed from cycle_2, so this folder stays
self-contained. It is a public endpoint and the list is a LIVING resource --
191 gene-disease pairs in Aug 2026 against 164 in 2019 -- so every run stamps
ClinGen's own file date and a sha256 of what it got. A result quoted against
"the ClinGen list" without saying which day's list is not reproducible.

Source: https://search.clinicalgenome.org/kb/gene-validity/download
Panel : Hearing Loss Gene Curation Expert Panel (affiliate 40007).
        NOT 50007, which is the Variant CEP -- a different body, often confused.

The CSV is not a plain table. Three banner lines, then the real header fenced
above and below by rows of `++++`. Those fences are QUOTED CSV FIELDS, not lines
beginning with `+`, so they have to be found after parsing, not before -- a
first version filtered on `line.startswith("++++")`, matched nothing, and
silently produced an empty gene list.

A gene can hold several pairs at different strengths -- CIB2 is both Refuted and
Definitive, GJB6 both Definitive and Refuted. Each gene is taken at its BEST
classification, so a weaker second pair never demotes it.

Output: phase_5/data/clingen_hl_genes.tsv, clingen_manifest.json
"""
import csv
import hashlib
import io
import json
import re
import urllib.request
from collections import defaultdict
from pathlib import Path

URL = "https://search.clinicalgenome.org/kb/gene-validity/download"
PANEL = "Hearing Loss Gene Curation Expert Panel"
ORDER = ["Definitive", "Strong", "Moderate", "Limited", "Disputed",
         "Refuted", "No Known Disease Relationship"]
OUT = Path(__file__).resolve().parents[1] / "data"
OUT.mkdir(parents=True, exist_ok=True)

raw = urllib.request.urlopen(URL, timeout=120).read()
sha = hashlib.sha256(raw).hexdigest()
text = raw.decode("utf-8-sig", errors="replace")

file_created = ""
for line in text.splitlines()[:6]:
    m = re.search(r"(\d{4}-\d{2}-\d{2})", line)
    if m:
        file_created = m.group(1)
        break

# Find the header by its content, and drop the fence rows, after parsing.
raw_rows = list(csv.reader(io.StringIO(text)))
hdr_i = next(i for i, r in enumerate(raw_rows) if r and r[0].strip() == "GENE SYMBOL")
header = [c.strip() for c in raw_rows[hdr_i]]
data = [r for r in raw_rows[hdr_i + 1:]
        if r and not r[0].startswith("+") and len(r) == len(header)]
rows = [dict(zip(header, r)) for r in data]
if not rows:
    raise SystemExit("no data rows parsed -- the CSV layout changed")

idx = {c: i for i, c in enumerate(header)}
pairs = [r for r in rows if PANEL.lower() in r.get("GCEP", "").lower()]
if not pairs:
    raise SystemExit(f"no rows for panel {PANEL!r} -- check the GCEP label")
best = {}
pair_counts = defaultdict(int)
for r in pairs:
    g = r.get("GENE SYMBOL", "").strip()
    c = r.get("CLASSIFICATION", "").strip()
    if not g:
        continue
    pair_counts[c] += 1
    rank = ORDER.index(c) if c in ORDER else len(ORDER)
    if g not in best or rank < best[g][1]:
        best[g] = (c, rank)

with (OUT / "clingen_hl_genes.tsv").open("w") as fh:
    fh.write("gene\tbest_classification\n")
    for g, (c, _) in sorted(best.items()):
        fh.write(f"{g}\t{c}\n")

gene_counts = defaultdict(int)
for _, (c, _) in best.items():
    gene_counts[c] += 1

manifest = {
    "source": URL, "panel": PANEL, "clingen_file_created": file_created,
    "sha256": sha, "pairs": len(pairs), "genes": len(best),
    "pairs_by_classification": {c: pair_counts[c] for c in ORDER if pair_counts[c]},
    "genes_by_best_classification": {c: gene_counts[c] for c in ORDER if gene_counts[c]},
    "tiers": {
        "definitive_strong": sum(gene_counts[c] for c in ["Definitive", "Strong"]),
        "plus_moderate": sum(gene_counts[c] for c in ["Definitive", "Strong", "Moderate"]),
    },
}
(OUT / "clingen_manifest.json").write_text(json.dumps(manifest, indent=2))
print(json.dumps(manifest, indent=2))
