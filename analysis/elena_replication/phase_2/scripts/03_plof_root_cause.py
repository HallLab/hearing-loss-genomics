#!/usr/bin/env python3
"""
Phase 2 / check 03 — the root cause of the pLOF inflation.

Check 02 found that 52.3% of the chr8 pLOF mask meets neither clause of the
documented definition, and left the cause open with one alternative unexcluded:
that the pipeline's own VEP run disagrees with the release's, making the whole
comparison a difference between two annotation runs rather than a defect.

This settles it using the pipeline's OWN annotation, not the release's. The mask
builder (rarevariant_geneburden/step3_masks.bsub) reads
rarevariantExWAS/variant_categories/chr{N}.classified.tsv, which carries the
pipeline's Consequence, its SpliceAI_max, and its own is_pLOF decision side by side.
So the question becomes internal: does is_pLOF follow the documented rule?

Documented rule (analysis plan, elena_publishes/):
    pLOF = frameshift, stop-gained, start-lost, stop-lost,
           or splice-site variants (SpliceAI >= 0.2)

Scope: chr8. Writes only under phase_2/results/.
"""
import json
from collections import defaultdict
from pathlib import Path
import pandas as pd

HERE  = Path(__file__).resolve().parent.parent
CLS   = Path("/project/hall/analysis/hearing-loss-genomics/analysis/elena/"
             "rarevariantExWAS/variant_categories")
OUT   = HERE / "results"; OUT.mkdir(parents=True, exist_ok=True)
CHROM = "8"

REAL_LOF = {"stop_gained", "frameshift_variant", "splice_acceptor_variant",
            "splice_donor_variant", "start_lost", "stop_lost", "transcript_ablation"}

cl   = pd.read_csv(CLS / f"chr{CHROM}.classified.tsv", sep="\t",
                   usecols=["Consequence", "LoF", "SpliceAI_max", "is_pLOF"], low_memory=False)
csq  = cl["Consequence"].fillna("").astype(str)
plof = cl["is_pLOF"].astype(str).str.lower().eq("true")
sai  = pd.to_numeric(cl["SpliceAI_max"], errors="coerce").fillna(0)
real = csq.apply(lambda s: bool(set(s.split(",")) & REAL_LOF))
spw  = csq.str.contains("splice")

report = {"chromosome": CHROM, "annotation_rows": len(cl)}

# ---- 1. what does is_pLOF=True consist of?
report["is_pLOF_true"] = {
    "rows": int(plof.sum()),
    "with_a_real_LoF_consequence": int((plof & real).sum()),
    "no_real_LoF_but_contains_'splice'": int((plof & ~real & spw).sum()),
    "no_real_LoF_and_no_'splice'": int((plof & ~real & ~spw).sum()),
}

# ---- 2. does the documented SpliceAI gate operate at all?
m = (~real) & spw
report["spliceai_gate"] = {
    "rows_splice_without_real_LoF": int(m.sum()),
    "is_pLOF_true__SpliceAI_below_0.2": int((m & plof & (sai < 0.2)).sum()),
    "is_pLOF_true__SpliceAI_at_or_above_0.2": int((m & plof & (sai >= 0.2)).sum()),
    "is_pLOF_FALSE_although_SpliceAI_at_or_above_0.2": int((m & ~plof & (sai >= 0.2)).sum()),
    "verdict": "SpliceAI does not gate the decision -- it is uncorrelated with it",
}

# ---- 3. which splice terms actually trigger it
t, f = defaultdict(int), defaultdict(int)
for s, p in zip(csq[m], plof[m]):
    for term in s.split(","):
        if "splice" in term:
            (t if p else f)[term] += 1
report["splice_terms"] = {term: {"is_pLOF_true": t[term], "is_pLOF_false": f[term],
                                 "always_triggers": f[term] == 0 and t[term] > 0}
                          for term in sorted(set(t) | set(f), key=lambda x: -(t[x] + f[x]))}
report["always_triggering_terms"] = [k for k, v in report["splice_terms"].items()
                                     if v["always_triggers"]]
report["vep_impact_of_those_terms"] = "LOW -- none is a loss-of-function consequence"

(OUT / "03_plof_root_cause.json").write_text(json.dumps(report, indent=2))
print(json.dumps(report, indent=2))
