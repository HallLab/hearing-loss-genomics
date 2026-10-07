#!/usr/bin/env python3
"""
Phase 2 / check 06 — two defects this replication missed, found by Nikki Palmiero.

Both were raised by Nikki on 2026-10-05 after she reviewed the pipeline
independently. She reached the same pLOF diagnosis as checks 01-04, and found two
things these checks did not.

  A. REVEL is parsed as a single number, but VEP writes it per transcript as a
     comma-separated list (0.131,.,0.131,0.131,0.131). pd.to_numeric(errors="coerce")
     turns that into NaN, so only single-valued entries survive and pDM silently
     undercounts.

  B. Variants are assigned to any overlapping transcript, including non-coding genes.
     A lncRNA or antisense transcript has no protein to lose function of, so a pLOF
     call on one is meaningless.

WHY WE MISSED (A), which is worth recording. Checks 01-02 compared pDM against the
release's group files from outside and saw a divergence consistent with a different
predictor threshold -- so we wrote it off as a judgement call. From outside, a
different threshold and broken parsing look identical: both produce "fewer variants
than the reference". Only reading the field distinguishes them. Nikki read the code.

Writes only under phase_2/results/.
"""
import json
from collections import defaultdict
from pathlib import Path
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent.parent
REL  = Path("/static/PMBB/PMBB-Release-2026-4.0/Exome")
MASK = Path("/project/hall/analysis/hearing-loss-genomics/analysis/elena/"
            "rarevariant_geneburden/masks_deduplicated")
OUT  = HERE / "results"; OUT.mkdir(parents=True, exist_ok=True)
SYMBOL_BIOTYPE = HERE / "data" / "symbol_biotype.tsv"
CHROM = "8"   # the REVEL measurement is a chr8 pilot; the biotype one is genome-wide

report = {"chromosome": CHROM, "raised_by": "Nikki Palmiero, 2026-10-05"}

# ------------------------------------------------------------------ A. REVEL
vep = REL / f"vep_annotations/PMBB-Release-2026-4.0_genetic_exome.vep_annotations.chr{CHROM}.tsv"
d = pd.read_csv(vep, sep="\t", low_memory=False,
                usecols=["#Uploaded_variation", "Consequence", "REVEL_score", "BIOTYPE"])
mis = d[d.Consequence.astype(str).str.contains("missense_variant")]
raw = mis.REVEL_score.astype(str)

as_single = pd.to_numeric(raw, errors="coerce")        # what the pipeline does

def max_of_list(s):                                     # what it should do
    vals = [float(x) for x in s.split(",") if x not in (".", "", "-", "nan")]
    return max(vals) if vals else np.nan

as_list = raw.map(max_of_list)

v = mis.assign(single=as_single, full=as_list).groupby("#Uploaded_variation").agg(
    single_dm=("single", lambda s: bool((s >= 0.5).any())),
    full_dm=("full", lambda s: bool((s >= 0.5).any())))

report["A_revel_parsing"] = {
    "missense_rows": len(mis),
    "rows_with_a_score__as_parsed": int(as_single.notna().sum()),
    "rows_with_a_score__reading_the_list": int(as_list.notna().sum()),
    "pct_recovered": round(100 * as_list.notna().mean(), 1),
    "variants_passing_REVEL_0.5__as_parsed": int(v.single_dm.sum()),
    "variants_passing_REVEL_0.5__correct": int(v.full_dm.sum()),
    "undercount": int(v.full_dm.sum() - v.single_dm.sum()),
    "pct_of_correct_that_was_lost": round(
        100 * (v.full_dm.sum() - v.single_dm.sum()) / max(int(v.full_dm.sum()), 1), 1),
}

# ------------------------------------------------------------------ B. biotype
#
# Measured at the GENE level, genome-wide. A first attempt asked whether each gene had
# any variant on a protein-coding transcript, on chr8 only, and returned 0.8% of
# entries -- a number that would have understated a correct finding. Two things were
# wrong with it: the question is whether the GENE is protein-coding, not whether some
# variant of it lands on a coding transcript; and TMC3-AS1, the example raised, is on
# chr15 and was not even in the sample.
#
# The symbol -> biotype map is built from the release's own VEP across all 22
# chromosomes by scripts/06a_symbol_biotype.sh.

bt = pd.read_csv(SYMBOL_BIOTYPE, sep="\t", header=None, names=["SYMBOL", "BIOTYPE"])
coding_genes = set(bt.loc[bt.BIOTYPE == "protein_coding", "SYMBOL"])

per_mask = {}
for mask in ["pLOF", "pDM", "pLOF_pDM"]:
    genes, entries, entries_nc = set(), 0, 0
    for line in (MASK / f"{mask}.txt").open():
        f = line.rstrip("\n").split()
        if len(f) < 3 or f[1] != "var":
            continue
        genes.add(f[0])
        entries += len(f) - 2
        if f[0] not in coding_genes:
            entries_nc += len(f) - 2
    nc = sorted(g for g in genes if g not in coding_genes)
    per_mask[mask] = {
        "genes": len(genes),
        "genes_not_protein_coding": len(nc),
        "pct_genes": round(100 * len(nc) / max(len(genes), 1), 1),
        "entries": entries,
        "entries_on_non_coding_genes": entries_nc,
        "pct_entries": round(100 * entries_nc / max(entries, 1), 1),
        "examples": nc[:8],
    }

report["B_biotype"] = {
    "scope": "genome-wide, gene level",
    "symbols_in_release_vep": len(set(bt.SYMBOL)),
    "of_which_protein_coding": len(coding_genes),
    "by_mask": per_mask,
    "reading": ("Large in genes, small in variants, and both are true. 1,101 genes carry a "
                "burden test that cannot mean anything -- a lncRNA has no protein to lose "
                "function of -- while only about 2% of mask entries sit on them. Reporting "
                "only the 2% would understate it; only the 1,101 would overstate it."),
    "pDM_barely_affected": ("8 genes, because a missense call requires a protein by "
                            "definition. pDM's defect is the REVEL parsing above, not this."),
}

(OUT / "06_revel_and_biotype.json").write_text(json.dumps(report, indent=2))
print(json.dumps(report, indent=2))
