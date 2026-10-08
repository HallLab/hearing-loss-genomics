#!/usr/bin/env python3
"""
Phase 5, step 1 -- give every tested gene a chromosome and a position.

A Manhattan plot needs a coordinate per gene; the SAIGE output has only the
symbol. Two sources, in order:

  1. The release's own Ensembl metadata,
     Exome/metadata/Homo_sapiens.GRCh38.113.ENSG_locations_symbols.tsv
     -- gene_id, chromosome, seq_region_start/end, gene_symbol. Authoritative,
     and an institutional release file, which keeps Phase 5 inside this
     folder's isolation rule.

  2. Our own masks, for any symbol the metadata does not carry. Each mask line
     holds the variant IDs actually tested, so the midpoint of a gene's
     variants is a build-matched coordinate that cannot disagree with what was
     tested. This is the approach Elena's step9_2_genemap.R takes, applied to
     our masks instead of her group files.

Source 2 is not only a fallback: it is also a CHECK on source 1. Where both
have a gene, the metadata interval should contain the variants we tested. A
gene where it does not is either a symbol collision or a build mismatch, and
the report counts them.

Input : phase_4/results/burden_all_cohorts.tsv
        phase_2/results/masks_v2/*.txt
Output: phase_5/results/gene_positions.tsv
        phase_5/results/01_gene_positions.json
"""
import json
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parents[1]
REPL = HERE.parent
META = Path("/static/PMBB/PMBB-Release-2026-4.0/Exome/metadata/"
            "Homo_sapiens.GRCh38.113.ENSG_locations_symbols.tsv")
OUT = HERE / "results"
OUT.mkdir(parents=True, exist_ok=True)

AUTOSOMES = {str(c) for c in range(1, 23)}

# ---- genes we actually tested ----
burden = pd.read_csv(HERE / "results/burden_all_cohorts.tsv", sep="\t",
                     usecols=["Region", "CHR"])
tested = burden.drop_duplicates("Region").set_index("Region")["CHR"].astype(int)
print(f"genes tested: {len(tested):,}")

# ---- source 2 first, because it is also the check ----
from_masks = {}
for mask in ("pLOF", "pDM", "pLOF_pDM"):
    with (REPL / f"phase_2/results/masks_v2/{mask}.txt").open() as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 3 or f[1] != "var":
                continue
            gene = f[0]
            pos = [int(v.split(":")[1]) for v in f[2:] if v]
            if not pos:
                continue
            lo, hi = min(pos), max(pos)
            if gene in from_masks:
                plo, phi = from_masks[gene][1:]
                from_masks[gene] = (f[2].split(":")[0], min(lo, plo), max(hi, phi))
            else:
                from_masks[gene] = (f[2].split(":")[0], lo, hi)
print(f"genes with coordinates from our masks: {len(from_masks):,}")

# ---- source 1 ----
meta = pd.read_csv(META, sep="\t", dtype=str)
meta = meta[meta.chromosome.isin(AUTOSOMES)]
meta["start"] = pd.to_numeric(meta.seq_region_start)
meta["end"] = pd.to_numeric(meta.seq_region_end)
# a symbol can appear more than once; keep the longest span, which is the
# gene rather than a fragment of it
meta["span"] = meta.end - meta.start
meta = (meta.sort_values("span", ascending=False)
            .drop_duplicates("gene_symbol")
            .set_index("gene_symbol"))
print(f"autosomal symbols in release metadata: {len(meta):,}")

rows, report = [], {"source": {"metadata": 0, "masks": 0, "neither": 0},
                    "checks": {"chrom_disagreement": 0,
                               "variants_outside_metadata_interval": 0,
                               "examples_chrom": [], "examples_outside": []}}

for gene, chr_tested in tested.items():
    m = meta.loc[gene] if gene in meta.index else None
    k = from_masks.get(gene)

    if m is not None and k is not None:
        if m.chromosome != k[0]:
            report["checks"]["chrom_disagreement"] += 1
            if len(report["checks"]["examples_chrom"]) < 8:
                report["checks"]["examples_chrom"].append(
                    {"gene": gene, "metadata": m.chromosome, "our_variants": k[0]})
        elif not (m.start <= k[1] and k[2] <= m.end):
            report["checks"]["variants_outside_metadata_interval"] += 1
            if len(report["checks"]["examples_outside"]) < 8:
                report["checks"]["examples_outside"].append(
                    {"gene": gene, "metadata": f"{m.start}-{m.end}",
                     "our_variants": f"{k[1]}-{k[2]}"})

    # the coordinate we use: metadata when it agrees on the chromosome we
    # tested, our masks otherwise
    if m is not None and m.chromosome == str(chr_tested):
        rows.append((gene, int(chr_tested), int((m.start + m.end) // 2), "metadata"))
        report["source"]["metadata"] += 1
    elif k is not None and k[0] == str(chr_tested):
        rows.append((gene, int(chr_tested), int((k[1] + k[2]) // 2), "masks"))
        report["source"]["masks"] += 1
    else:
        report["source"]["neither"] += 1

pos = pd.DataFrame(rows, columns=["Region", "CHR", "POS", "coord_source"])
pos.to_csv(OUT / "gene_positions.tsv", sep="\t", index=False)
report["genes_tested"] = int(len(tested))
report["genes_placed"] = int(len(pos))
(OUT / "01_gene_positions.json").write_text(json.dumps(report, indent=2))

print(f"\nplaced {len(pos):,} of {len(tested):,} genes")
for k, v in report["source"].items():
    print(f"  from {k:9s}: {v:,}")
c = report["checks"]
print(f"\nchecks against our own masks:")
print(f"  chromosome disagreements            : {c['chrom_disagreement']:,}")
print(f"  variants outside metadata interval  : {c['variants_outside_metadata_interval']:,}")
