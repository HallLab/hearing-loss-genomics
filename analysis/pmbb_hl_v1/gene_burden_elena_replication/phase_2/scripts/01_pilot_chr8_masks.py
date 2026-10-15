#!/usr/bin/env python3
"""
Phase 2 / check 01 — pilot: compare the pipeline's masks against the release's
group files, on chr8 only.

A pilot, deliberately. The two sides disagree on both identifier conventions:

    gene     pipeline = HGNC symbol      release = Ensembl gene id
    variant  pipeline = 11:2147642:C:G   release = chr11_2147642_C_G

Any comparison therefore runs through a translation layer, which is where a
mis-specified diff would silently produce a wrong answer. Validating that layer on
one chromosome costs minutes; validating it after a genome-wide run costs hours.

Which masks the run actually used was established by resolving the Nextflow staged
symlinks, not by reading a directory name:
    rarevariant_geneburden/masks_deduplicated/{ALL,pDM,pLOF,pLOF_pDM}.txt

Writes only under phase_2/results/.
"""
import json, re
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent
REL  = Path("/static/PMBB/PMBB-Release-2026-4.0/Exome")
MASK = Path("/project/hall/analysis/hearing-loss-genomics/analysis/elena/"
            "rarevariant_geneburden/masks_deduplicated")
OUT  = HERE / "results"; OUT.mkdir(parents=True, exist_ok=True)

CHROM = "8"

# ---------------------------------------------------------------- release side
grp = REL / f"group_file_annotations/PMBB-Release-2026-4.0_genetic_exome.groupfiles.chr{CHROM}.txt"
rel_gene = defaultdict(dict)                      # ensg -> {variant -> anno}
with grp.open() as fh:
    pending = {}
    for line in fh:
        f = line.rstrip("\n").split("\t")
        gene, kind, rest = f[0], f[1], f[2:]
        if kind == "var":
            pending[gene] = rest
        elif kind == "anno":
            for v, a in zip(pending.pop(gene, []), rest):
                rel_gene[gene][v] = a

# ensg -> symbol, from the release's own VEP output
vep = REL / f"vep_annotations/PMBB-Release-2026-4.0_genetic_exome.vep_annotations.chr{CHROM}.tsv"
ensg2sym, sym2ensg = {}, defaultdict(set)
with vep.open() as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    gi, si = hdr.index("Gene"), hdr.index("SYMBOL")
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) > si and f[gi] and f[si] and f[si] != "-":
            ensg2sym[f[gi]] = f[si]
            sym2ensg[f[si]].add(f[gi])

# ---------------------------------------------------------------- pipeline side
def norm(v):
    """11:2147642:C:G -> chr11_2147642_C_G"""
    p = v.split(":")
    return f"chr{p[0]}_{p[1]}_{p[2]}_{p[3]}" if len(p) == 4 else None

pipe = {}
for mask in ["pLOF", "pDM", "pLOF_pDM", "ALL"]:
    per_gene = defaultdict(set)
    with (MASK / f"{mask}.txt").open() as fh:
        for line in fh:
            f = line.rstrip("\n").split()
            if len(f) < 3 or f[1] != "var":
                continue
            for v in f[2:]:
                if v.startswith(f"{CHROM}:"):
                    n = norm(v)
                    if n:
                        per_gene[f[0]].add(n)
    pipe[mask] = {g: s for g, s in per_gene.items() if s}

# ---------------------------------------------------------------- compare
rel_by_anno = defaultdict(lambda: defaultdict(set))   # anno -> symbol -> variants
unmapped_ensg = set()
for ensg, vmap in rel_gene.items():
    sym = ensg2sym.get(ensg)
    if sym is None:
        unmapped_ensg.add(ensg); continue
    for v, a in vmap.items():
        rel_by_anno[a][sym].add(v)

PAIRS = [("pLOF", "pLoF"), ("pDM", "damaging_missense")]
report = {
    "chromosome": CHROM,
    "release": {
        "genes": len(rel_gene),
        "variants": sum(len(v) for v in rel_gene.values()),
        "ensg_unmapped_to_symbol": len(unmapped_ensg),
        "by_annotation": {a: sum(len(s) for s in g.values()) for a, g in rel_by_anno.items()},
    },
    "pipeline": {m: {"genes": len(g), "variants": sum(len(s) for s in g.values())}
                 for m, g in pipe.items()},
    "comparison": {},
}

for pmask, ranno in PAIRS:
    P, R = pipe.get(pmask, {}), rel_by_anno.get(ranno, {})
    shared = set(P) & set(R)
    agree = sum(1 for g in shared if P[g] == R[g])
    pv = set().union(*P.values()) if P else set()
    rv = set().union(*R.values()) if R else set()
    report["comparison"][f"{pmask}_vs_{ranno}"] = {
        "genes_pipeline_only": len(set(P) - set(R)),
        "genes_release_only":  len(set(R) - set(P)),
        "genes_shared":        len(shared),
        "genes_identical_variant_set": agree,
        "variants_pipeline_only": len(pv - rv),
        "variants_release_only":  len(rv - pv),
        "variants_shared":        len(pv & rv),
    }

(OUT / "01_pilot_chr8.json").write_text(json.dumps(report, indent=2, default=str))
print(json.dumps(report, indent=2, default=str))
