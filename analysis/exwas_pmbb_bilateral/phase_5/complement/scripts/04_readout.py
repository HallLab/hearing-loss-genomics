#!/usr/bin/env python3
"""
Complement run, the readout -- where does the ClinGen enrichment live?

Sets the complement beside everything already measured, so the comparison is
made against a full set of reference points rather than against the one that
suits a story.

Decision rule, written before the result exists:

  complement near or above the power replicates (median 2 of the top 50)
      The signal lives in the cases the restriction DISCARDS. Since those are
      unilateral, conductive, mixed and -- mostly -- unspecified hearing loss,
      the likely cause is ICD coding imprecision rather than biology: people
      with genuine bilateral sensorineural loss coded as something vaguer. The
      restriction then costs real cases, and the phenotype-by-code approach has
      a measured ceiling.

  complement near zero, like the restricted arm
      The signal lives in the cases the restriction KEEPS, and the restricted
      arm's null is lost power. The phenotype decision comes out clean, and
      what is left is a sample-size problem rather than a definition problem.

  both near zero and the broad arm's enrichment unexplained
      Would mean the enrichment needs the full 6,752 and neither half carries
      it, which is itself an answer about power rather than definition.

The complement has 3,588 cases against 3,164 for every comparison point --
13% more, which biases toward finding enrichment. A null here is therefore
stronger than a size-matched null would have been.

Output: phase_5/complement/results/04_readout.json
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hypergeom

HERE = Path(__file__).resolve().parents[1]
EX = HERE.parent.parent
PC = EX / "phase_5/power_control"
BROAD = Path("/project/hall/analysis/hearing-loss-genomics/analysis/"
             "elena_replication/phase_4/results")

cg = pd.read_csv(EX / "phase_5/data/clingen_hl_genes.tsv", sep="\t")
GENES = set(cg[cg.best_classification.isin(["Definitive", "Strong"])].gene)


def readout(d, label, n_cases):
    d = d.dropna(subset=["omnibus_p"]).sort_values("omnibus_p").reset_index(drop=True)
    d["rank"] = np.arange(1, len(d) + 1)
    inset = d.Region.isin(GENES)
    n_in = int(inset.sum())
    out = {"label": label, "cases": n_cases, "clingen_tested": n_in,
           "median_rank_clingen": int(d["rank"][inset].median()),
           "genes_tested": int(len(d)), "top_n": {}}
    for N in (50, 100, 250):
        hit = int((d["rank"][inset] <= N).sum())
        out["top_n"][N] = {"observed": hit,
                           "expected": round(N * n_in / len(d), 2),
                           "p": float(hypergeom.sf(hit - 1, len(d), n_in, N))}
    out["top50_genes"] = sorted(d.Region[inset & (d["rank"] <= 50)].tolist())
    return out


def omnibus_from(dirpath, stem):
    frames = []
    for chrom in range(1, 23):
        f = dirpath / f"{stem}_chr{chrom}.txt"
        if not f.is_file() or f.stat().st_size == 0:
            raise SystemExit(f"missing or empty: {f}")
        frames.append(pd.read_csv(f, sep="\t", dtype={"Region": str, "Group": str}))
    a = pd.concat(frames, ignore_index=True)
    o = a[a.Group == "Cauchy"][["Region", "Pvalue"]].rename(
        columns={"Pvalue": "omnibus_p"})
    o["omnibus_p"] = pd.to_numeric(o.omnibus_p)
    if o.Region.duplicated().any():
        raise SystemExit(f"{stem}: more than one Cauchy row for some gene")
    return o


res = {"rows": []}
res["rows"].append(readout(
    pd.read_csv(BROAD / "omnibus_combined.tsv", sep="\t")[["Region", "omnibus_p"]],
    "broad, all cases", 6752))
for r in range(1, 6):
    res["rows"].append(readout(
        omnibus_from(PC / f"results/step2/rep{r}", f"rep{r}"),
        f"broad subsampled, draw {r}", 3164))
res["rows"].append(readout(
    pd.read_csv(EX / "phase_5/results/omnibus_combined.tsv", sep="\t")[
        ["Region", "omnibus_p"]],
    "restricted (bilateral SN)", 3164))
res["rows"].append(readout(
    omnibus_from(HERE / "results/step2", "complement"),
    "COMPLEMENT (what the restriction drops)", 3588))

(HERE / "results" / "04_readout.json").write_text(json.dumps(res, indent=2))

print(f"{'':42s} {'cases':>6s} {'top50':>6s} {'exp':>5s} {'p':>9s}  genes")
for v in res["rows"]:
    t = v["top_n"][50]
    print(f"{v['label']:42s} {v['cases']:>6,} {t['observed']:>6} {t['expected']:>5} "
          f"{t['p']:>9.2g}  {', '.join(v['top50_genes'])}")
