#!/usr/bin/env python3
"""
Phase 2 / check 02 — is the pLOF divergence a universe artifact, a transcript
policy, or a mis-assignment?

Check 01 found the pipeline's pLOF mask carries ~3.7x more chr8 variants than the
release's pLoF annotation. Three explanations had to be separated before the gap
could be called anything:

  1. different variant universes   -> the release group files cover fewer variants
  2. transcript policy             -> pLOF on a non-canonical transcript is
                                      legitimately pLOF under an any-transcript rule
  3. mis-assignment                -> variants that are not loss-of-function at all

The pipeline's own documented definition, from the analysis plan in
elena_publishes/ ("Hall Lab_ Analysis Plan Draft - Rare Variant ExWAS.docx"):

    pLOF (predicted loss-of-function)
    Frameshift, stop-gained, start-lost, stop-lost,
    or splice-site variants (SpliceAI >= 0.2)
    -- "VEP Consequence + SpliceAI"

That definition is broader than a consequence-only rule, because SpliceAI scores every
variant including synonymous ones. It is therefore tested here explicitly rather than
assumed to be the explanation.

This check quantifies each. Release VEP is used as the reference for consequences
across every transcript, not just the canonical one.

Scope: chr8 pilot. Genome-wide is check 03.
Writes only under phase_2/results/.
"""
import json
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent
REL  = Path("/static/PMBB/PMBB-Release-2026-4.0/Exome")
MASK = Path("/project/hall/analysis/hearing-loss-genomics/analysis/elena/"
            "rarevariant_geneburden/masks_deduplicated")
OUT  = HERE / "results"; OUT.mkdir(parents=True, exist_ok=True)
CHROM = "8"

LOF_CSQ = {"stop_gained", "frameshift_variant", "splice_acceptor_variant",
           "splice_donor_variant", "start_lost", "stop_lost", "transcript_ablation"}

# ---- release: canonical annotation per variant, from the group file
rel, pend = {}, {}
for line in (REL / f"group_file_annotations/PMBB-Release-2026-4.0_genetic_exome.groupfiles.chr{CHROM}.txt").open():
    f = line.rstrip("\n").split("\t"); g, k, rest = f[0], f[1], f[2:]
    if k == "var":
        pend[g] = rest
    else:
        for v, a in zip(pend.pop(g, []), rest):
            rel[v] = a

# ---- pipeline mask, normalised to release variant ids
def norm(v):
    p = v.split(":")
    return f"chr{p[0]}_{p[1]}_{p[2]}_{p[3]}" if len(p) == 4 else None

def load(mask):
    s = set()
    for line in (MASK / f"{mask}.txt").open():
        f = line.rstrip("\n").split()
        if len(f) < 3 or f[1] != "var":
            continue
        for v in f[2:]:
            if v.startswith(CHROM + ":"):
                n = norm(v)
                if n:
                    s.add(n)
    return s

# ---- structural check: do var and anno rows pair up?
misaligned = {}
for mask in ["pLOF", "pDM", "pLOF_pDM", "ALL"]:
    nv, na = {}, {}
    for line in (MASK / f"{mask}.txt").open():
        f = line.rstrip("\n").split()
        if len(f) < 2: continue
        (nv if f[1] == "var" else na)[f[0]] = len(f) - 2
    misaligned[mask] = sum(1 for g in nv if nv.get(g) != na.get(g))

# ---- consequences across ALL transcripts, from release VEP
def csq_any(variants):
    """Consequences on every transcript, plus max SpliceAI delta score."""
    seen, spl = defaultdict(set), defaultdict(float)
    vep = REL / f"vep_annotations/PMBB-Release-2026-4.0_genetic_exome.vep_annotations.chr{CHROM}.tsv"
    with vep.open() as fh:
        h = fh.readline().rstrip("\n").split("\t")
        iu, ic = h.index("#Uploaded_variation"), h.index("Consequence")
        isp = h.index("SpliceAI_pred")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            v = f[iu]
            if v in variants:
                seen[v].update(f[ic].split("&"))
                if len(f) > isp and f[isp] not in ("", "-"):
                    for pred in f[isp].split(","):
                        q = pred.split("|")
                        if len(q) >= 5:
                            try:
                                spl[v] = max(spl[v], max(float(x) for x in q[1:5]))
                            except ValueError:
                                pass
    return seen, spl

report = {"chromosome": CHROM,
          "release_universe": len(rel),
          "var_anno_misaligned_genes": misaligned}

for pmask, ranno in [("pLOF", "pLoF"), ("pDM", "damaging_missense")]:
    P = load(pmask)
    R = {v for v, a in rel.items() if a == ranno}
    only = P - R
    in_universe = {v for v in only if v in rel}
    by_anno = defaultdict(int)
    for v in in_universe:
        by_anno[rel[v]] += 1

    seen, spl = csq_any(in_universe)
    lof_somewhere = {v for v in in_universe if seen.get(v, set()) & LOF_CSQ}
    no_lof = in_universe - lof_somewhere
    residual = {v for v in no_lof if spl.get(v, 0.0) < 0.2}
    resid_by_anno = defaultdict(int)
    for v in residual:
        resid_by_anno[rel[v]] += 1

    report[pmask] = {
        "pipeline": len(P), "release": len(R), "shared": len(P & R),
        "pipeline_only": len(only),
        "  absent_from_release_universe": len(only - in_universe),
        "  present_with_other_annotation": len(in_universe),
        "    by_release_annotation": dict(by_anno),
        "    LoF_on_some_transcript": len(lof_somewhere),
        "    LoF_on_NO_transcript": len(no_lof),
        "      of those, SpliceAI >= 0.2": len(no_lof) - len(residual),
        "      SpliceAI at 0.5": sum(1 for v in no_lof if spl.get(v, 0.0) >= 0.5),
        "      SpliceAI at 0.8": sum(1 for v in no_lof if spl.get(v, 0.0) >= 0.8),
        "    RESIDUAL_meets_neither_criterion": len(residual),
        "      pct_of_pipeline_mask": round(100 * len(residual) / max(len(P), 1), 1),
        "      by_release_annotation": dict(resid_by_anno),
        "release_only": len(R - P),
    }

(OUT / "02_plof_diagnosis.json").write_text(json.dumps(report, indent=2))
print(json.dumps(report, indent=2))
