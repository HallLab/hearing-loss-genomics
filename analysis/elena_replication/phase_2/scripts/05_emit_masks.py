#!/usr/bin/env python3
"""
Phase 2 / check 05 — emit the two mask sets, and close the phase.

Checks 01-04 established what is wrong with the pLOF mask. This writes down the
deliverable, mirroring Phase 1's two arms:

  reproduction  the pipeline's masks exactly as they ran
  corrected     pLOF and pLOF_pDM rebuilt with the documented rule applied

The corrected rule, from the analysis plan in elena_publishes/:

    pLOF = frameshift, stop_gained, start_lost, stop_lost,
           splice_acceptor_variant, splice_donor_variant
         OR a splice-site annotation WITH SpliceAI >= 0.2

Two changes from what was implemented. Terms are matched exactly against the
comma-separated Consequence list rather than as substrings, and the three VEP
IMPACT=LOW splice terms are admitted only through the SpliceAI gate the plan
specifies -- they are no longer unconditional.

ALL and pDM are untouched: pDM is built from is_AlphaMissense_DM / is_REVEL_DM,
which the defect does not involve, and ALL is every variant. They are symlinked
rather than copied, because ALL.txt is 425 MB and duplicating it buys nothing.

pLOF_pDM needs care. The builder gives the pLOF annotation precedence, so a variant
labelled pLOF there may also be pDM. Dropping it would silently lose a legitimate
pDM variant, so a variant that fails corrected-pLOF but is pDM is RE-ANNOTATED
rather than removed.

Writes only under phase_2/results/.
"""
import json, os
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent
CLS  = Path("/project/hall/analysis/hearing-loss-genomics/analysis/elena/"
            "rarevariantExWAS/variant_categories")
SRC  = Path("/project/hall/analysis/hearing-loss-genomics/analysis/elena/"
            "rarevariant_geneburden/masks_deduplicated")
OUT  = HERE / "results" / "masks"; OUT.mkdir(parents=True, exist_ok=True)

REAL_LOF = {"stop_gained", "frameshift_variant", "splice_acceptor_variant",
            "splice_donor_variant", "start_lost", "stop_lost", "transcript_ablation"}

def to_mask_id(v):
    """chr8_100002105_G_T -> 8:100002105:G:T"""
    p = v.split("_")
    return f"{p[0][3:]}:{p[1]}:{p[2]}:{p[3]}" if len(p) == 4 and p[0].startswith("chr") else None

# ---------------------------------------------- per-variant corrected decisions
corrected_plof, is_pdm = set(), set()
stats = {"rows": 0, "variants": 0}
for c in range(1, 23):
    with (CLS / f"chr{c}.classified.tsv").open() as fh:
        h = fh.readline().rstrip("\n").split("\t")
        iv, ic = h.index("#Uploaded_variation"), h.index("Consequence")
        isp, iam, ire = h.index("SpliceAI_max"), h.index("is_AlphaMissense_DM"), h.index("is_REVEL_DM")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            stats["rows"] += 1
            vid = to_mask_id(f[iv])
            if vid is None:
                continue
            terms = set(f[ic].split(","))
            if terms & REAL_LOF:
                corrected_plof.add(vid)
            else:
                try:
                    sai = float(f[isp])
                except (ValueError, IndexError):
                    sai = 0.0
                if sai >= 0.2 and any("splice" in t for t in terms):
                    corrected_plof.add(vid)
            if f[iam].lower() == "true" or f[ire].lower() == "true":
                is_pdm.add(vid)

# ---------------------------------------------- rewrite the two affected masks
def rewrite(name, keep_plof, allow_pdm):
    """Filter a mask file, preserving gene order and formatting."""
    src, dst = SRC / f"{name}.txt", OUT / f"corrected_{name}.txt"
    kept = dropped = reannotated = 0
    with src.open() as fh, dst.open("w") as out:
        pending = {}
        for line in fh:
            f = line.rstrip("\n").split()
            gene, kind, items = f[0], f[1], f[2:]
            if kind == "var":
                pending[gene] = items
                continue
            varlist, annos, keep_v, keep_a = pending.pop(gene, []), items, [], []
            for v, a in zip(varlist, annos):
                if a == "pLOF":
                    if v in keep_plof:
                        keep_v.append(v); keep_a.append("pLOF"); kept += 1
                    elif allow_pdm and v in is_pdm:
                        keep_v.append(v); keep_a.append("pDM"); reannotated += 1
                    else:
                        dropped += 1
                else:
                    keep_v.append(v); keep_a.append(a); kept += 1
            if keep_v:
                out.write("\t".join([gene, "var"] + keep_v) + "\n")
                out.write("\t".join([gene, "anno"] + keep_a) + "\n")
    return {"kept": kept, "dropped": dropped, "reannotated_as_pDM": reannotated,
            "file": dst.name}

manifest = {
    "corrected_rule": "exact term match on real LoF, OR any splice annotation with SpliceAI >= 0.2",
    "annotation_rows_read": stats["rows"],
    "variants_passing_corrected_pLOF": len(corrected_plof),
    "variants_that_are_pDM": len(is_pdm),
    "rewritten": {}, "unchanged": {},
}
manifest["rewritten"]["pLOF"] = rewrite("pLOF", corrected_plof, allow_pdm=False)
manifest["rewritten"]["pLOF_pDM"] = rewrite("pLOF_pDM", corrected_plof, allow_pdm=True)

for name in ["ALL", "pDM"]:
    link = OUT / f"corrected_{name}.txt"
    if link.is_symlink() or link.exists():
        link.unlink()
    link.symlink_to(SRC / f"{name}.txt")
    manifest["unchanged"][name] = {"symlink_to": str(SRC / f'{name}.txt'),
                                   "reason": "the defect does not involve this mask"}

manifest["reproduction_arm"] = {
    "location": str(SRC),
    "note": "the pipeline's masks as they ran; not copied, they are a stable input",
}
(HERE / "results" / "05_masks_manifest.json").write_text(json.dumps(manifest, indent=2))
print(json.dumps(manifest, indent=2))
