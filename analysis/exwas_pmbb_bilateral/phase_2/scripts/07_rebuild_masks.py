#!/usr/bin/env python3
"""
Phase 2 / check 07 — rebuild all four masks from the release annotation.

Supersedes check 05, which filtered the pipeline's masks down. That worked while the
only defect was over-inclusion in pLOF. It cannot work now: the REVEL parsing bug
means pDM is missing roughly three quarters of its variants, and nothing can be
*added* by filtering. So the masks are rebuilt from source.

SOURCE. The release's own VEP (Exome/vep_annotations/), not the pipeline's
variant_categories/, because the latter already carries the damaged REVEL -- its
REVEL_score column holds no multi-valued entries at all, the lists having been
coerced to NaN before it was written.

THE FOUR RULES

  pLOF      exact term match on {frameshift, stop_gained, start_lost, stop_lost,
            splice_acceptor, splice_donor, transcript_ablation}
            OR any other splice annotation WITH SpliceAI >= 0.2
            AND BIOTYPE == protein_coding

  pDM       am_class in {pathogenic, likely_pathogenic}
            OR max REVEL across transcripts >= 0.5
            AND BIOTYPE == protein_coding

  pLOF_pDM  the union, pLOF annotation taking precedence (as the pipeline does)
  ALL       every variant present in the cohort PLINK set
            AND BIOTYPE == protein_coding   <- a deliberate choice, see below

Three fixes against what ran, one per defect:
  - terms matched exactly rather than as substrings, and the three IMPACT=LOW splice
    terms admitted only through the SpliceAI gate the analysis plan specifies
  - REVEL read as the per-transcript list VEP writes, maximum taken
  - non-coding genes excluded

Terms are matched against the consequence list, NOT against VEP's IMPACT field, which
is per row: a row reading `stop_gained,splice_polypyrimidine_tract_variant` is HIGH as
a whole, so filtering on IMPACT would readmit exactly what is being removed.

THE BIOTYPE RESTRICTION APPLIES TO `ALL` TOO, decided by Andre on 2026-10-05.

`ALL` had no defect: it is the unrestricted baseline, every rare variant with no damage
filter. Applying the biotype restriction to it removes 6,532 genes (24,570 -> 18,038)
and is therefore a second difference between the arms, beyond the three fixes.

Recorded because it was argued both ways and the quieter option was not taken.

  for  -- a burden test on a lncRNA is meaningless in any mask, not only in pLOF. A
         gene with no protein has no function to lose whichever variant set is used,
         and 6,532 fewer genes is 6,532 fewer tests carrying multiple-testing
         correction.
  against -- a replication should change only what is defective, and the non-coding
         problem was raised about the pLOF masks. A second difference makes any later
         divergence harder to attribute.

The decision is for. Later phases comparing the arms must therefore hold TWO
differences constant in `ALL`: the three damage fixes, and this gene-set restriction.

SCOPE. Variants are kept only if present in the cohort's PLINK set, as the pipeline
does -- the mask is cohort-specific by construction.

Writes only under phase_2/results/.
"""
import json, sys
from collections import defaultdict
from pathlib import Path

HERE  = Path(__file__).resolve().parent.parent
REL   = Path("/static/PMBB/PMBB-Release-2026-4.0/Exome/vep_annotations")
PLINK = Path("/project/hall/analysis/hearing-loss-genomics/analysis/elena/"
             "rarevariant_geneburden/plink_deduplicated")
OUT   = HERE / "results" / "masks_v2"; OUT.mkdir(parents=True, exist_ok=True)

REAL_LOF = {"stop_gained", "frameshift_variant", "splice_acceptor_variant",
            "splice_donor_variant", "start_lost", "stop_lost", "transcript_ablation"}
AM_DM    = {"pathogenic", "likely_pathogenic"}

def max_revel(s):
    best = None
    for x in s.split(","):
        if x in (".", "", "-"):
            continue
        try:
            v = float(x)
        except ValueError:
            continue
        best = v if best is None else max(best, v)
    return best

def max_spliceai(s):
    best = 0.0
    for pred in s.split(","):
        q = pred.split("|")
        if len(q) >= 5:
            for x in q[1:5]:
                try:
                    best = max(best, float(x))
                except ValueError:
                    pass
    return best

stats = defaultdict(int)
masks = {m: defaultdict(set) for m in ["pLOF", "pDM", "pLOF_pDM", "ALL"]}
annot = {}

for c in range(1, 23):
    # variant IDs present in the cohort, and the id form the mask uses
    bim_ids = {}
    with (PLINK / f"chr{c}_deduplicated.bim").open() as fh:
        for line in fh:
            f = line.split()
            p = f[1].split(":")
            if len(p) == 4:
                bim_ids[f"chr{p[0]}_{p[1]}_{p[2]}_{p[3]}"] = f[1]

    per_var = defaultdict(lambda: {"plof": False, "pdm": False, "sym": None})
    with (REL / f"PMBB-Release-2026-4.0_genetic_exome.vep_annotations.chr{c}.tsv").open() as fh:
        h = fh.readline().rstrip("\n").split("\t")
        iu, ic, iam = h.index("#Uploaded_variation"), h.index("Consequence"), h.index("am_class")
        irv, isp = h.index("REVEL_score"), h.index("SpliceAI_pred")
        isy, ibt = h.index("SYMBOL"), h.index("BIOTYPE")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            vid = f[iu]
            if vid not in bim_ids:
                continue
            stats["rows_in_cohort"] += 1
            if f[ibt] != "protein_coding":          # the biotype fix
                continue
            sym = f[isy]
            if not sym or sym == "-":
                continue
            r = per_var[(vid, sym)]
            r["sym"] = sym
            terms = set(f[ic].split(","))
            if terms & REAL_LOF:
                r["plof"] = True
            elif any("splice" in t for t in terms) and max_spliceai(f[isp]) >= 0.2:
                r["plof"] = True
            rev = max_revel(f[irv])                  # the REVEL fix
            if f[iam] in AM_DM or (rev is not None and rev >= 0.5):
                r["pdm"] = True

    for (vid, sym), r in per_var.items():
        mid = bim_ids[vid]
        masks["ALL"][sym].add((mid, "ALL"))
        if r["plof"]:
            masks["pLOF"][sym].add((mid, "pLOF"))
        if r["pdm"]:
            masks["pDM"][sym].add((mid, "pDM"))
        if r["plof"] or r["pdm"]:
            masks["pLOF_pDM"][sym].add((mid, "pLOF" if r["plof"] else "pDM"))
    print(f"  chr{c} done", file=sys.stderr, flush=True)

manifest = {"source": "release VEP, rebuilt from scratch", "masks": {}}
for m, genes in masks.items():
    path = OUT / f"{m}.txt"
    n = 0
    with path.open("w") as out:
        for g in sorted(genes):
            items = sorted(genes[g])
            out.write("\t".join([g, "var"] + [v for v, _ in items]) + "\n")
            out.write("\t".join([g, "anno"] + [a for _, a in items]) + "\n")
            n += len(items)
    manifest["masks"][m] = {"genes": len(genes), "entries": n, "file": path.name}

(HERE / "results" / "07_masks_v2_manifest.json").write_text(json.dumps(manifest, indent=2))
print(json.dumps(manifest, indent=2))
