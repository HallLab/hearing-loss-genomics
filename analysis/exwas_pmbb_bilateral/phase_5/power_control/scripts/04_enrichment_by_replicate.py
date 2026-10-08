#!/usr/bin/env python3
"""
Power control, the readout -- does the ClinGen enrichment survive at matched N?

Runs the same top-N test as phase_5/scripts/06_clingen_enrichment.py on each of
the five size-matched replicates, and sets them beside the two arms it is
adjudicating between.

HOW TO READ IT, decided before seeing the result so that it cannot be decided
after. The broad arm at full size has 5 ClinGen Definitive/Strong genes in its
top 50 against 0.26 expected; the restricted arm has 0.

  replicates cluster near 5  ->  the restriction removed SIGNAL. The phenotype
                                 decision costs real biology.
  replicates cluster near 0  ->  the restriction removed only CASES. The null is
                                 power, the phenotype decision is exonerated,
                                 and the broad arm's enrichment is a
                                 sample-size effect rather than a phenotype one.
  replicates in between      ->  both, in proportion. Report the median and do
                                 not round it to whichever story is tidier.

Output: phase_5/power_control/results/04_enrichment.json
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hypergeom

HERE = Path(__file__).resolve().parents[1]
EX = HERE.parent.parent
BROAD = Path("/project/hall/analysis/hearing-loss-genomics/analysis/"
             "elena_replication/phase_4/results")

cg = pd.read_csv(EX / "phase_5/data/clingen_hl_genes.tsv", sep="\t")
GENES = set(cg[cg.best_classification.isin(["Definitive", "Strong"])].gene)


def readout(df, label):
    d = df.dropna(subset=["omnibus_p"]).sort_values("omnibus_p").reset_index(drop=True)
    d["rank"] = np.arange(1, len(d) + 1)
    inset = d.Region.isin(GENES)
    n_in = int(inset.sum())
    out = {"label": label, "genes_tested": int(len(d)), "clingen_tested": n_in,
           "median_rank_clingen": int(d["rank"][inset].median()), "top_n": {}}
    for N in (50, 100, 250):
        hit = int((d["rank"][inset] <= N).sum())
        exp = N * n_in / len(d)
        out["top_n"][N] = {"observed": hit, "expected": round(exp, 2),
                           "p": float(hypergeom.sf(hit - 1, len(d), n_in, N))}
    out["top50_genes"] = sorted(d.Region[inset & (d["rank"] <= 50)].tolist())
    return out


def merge_replicate(r):
    frames = []
    for chrom in range(1, 23):
        f = HERE / f"results/step2/rep{r}/rep{r}_chr{chrom}.txt"
        if not f.is_file() or f.stat().st_size == 0:
            raise SystemExit(f"missing or empty: {f}")
        frames.append(pd.read_csv(f, sep="\t", dtype={"Region": str, "Group": str}))
    a = pd.concat(frames, ignore_index=True)
    o = a[a.Group == "Cauchy"][["Region", "Pvalue"]].rename(
        columns={"Pvalue": "omnibus_p"})
    o["omnibus_p"] = pd.to_numeric(o.omnibus_p)
    if o.Region.duplicated().any():
        raise SystemExit(f"rep{r}: more than one Cauchy row for some gene")
    return o


res = {"gene_set": "ClinGen HL GCEP, Definitive+Strong",
       "n_gene_set": len(GENES), "arms": {}, "replicates": {}}

res["arms"]["broad_full_6752_cases"] = readout(
    pd.read_csv(BROAD / "omnibus_combined.tsv", sep="\t")[["Region", "omnibus_p"]],
    "broad phenotype, 6,752 cases")
res["arms"]["restricted_3164_cases"] = readout(
    pd.read_csv(EX / "phase_5/results/omnibus_combined.tsv", sep="\t")[
        ["Region", "omnibus_p"]],
    "bilateral sensorineural, 3,164 cases")

counts = []
for r in range(1, 6):
    res["replicates"][f"rep{r}"] = readout(
        merge_replicate(r), f"broad phenotype subsampled to 3,164 cases, draw {r}")
    counts.append(res["replicates"][f"rep{r}"]["top_n"][50]["observed"])

res["verdict_inputs"] = {
    "broad_full_top50": res["arms"]["broad_full_6752_cases"]["top_n"][50]["observed"],
    "restricted_top50": res["arms"]["restricted_3164_cases"]["top_n"][50]["observed"],
    "replicate_top50": counts,
    "replicate_median": float(np.median(counts)),
}

(HERE / "results" / "04_enrichment.json").write_text(json.dumps(res, indent=2))

print(f"ClinGen Definitive+Strong in the top 50, combined cohort\n")
print(f"{'':44s} {'top50':>6s} {'exp':>5s} {'p':>9s}   genes")
for k in ["broad_full_6752_cases", "restricted_3164_cases"]:
    v = res["arms"][k]; t = v["top_n"][50]
    print(f"{v['label']:44s} {t['observed']:>6} {t['expected']:>5} {t['p']:>9.2g}   "
          f"{', '.join(v['top50_genes'])}")
print()
for r in range(1, 6):
    v = res["replicates"][f"rep{r}"]; t = v["top_n"][50]
    print(f"{v['label']:44s} {t['observed']:>6} {t['expected']:>5} {t['p']:>9.2g}   "
          f"{', '.join(v['top50_genes'])}")
print(f"\nreplicate median: {res['verdict_inputs']['replicate_median']}  "
      f"(broad {res['verdict_inputs']['broad_full_top50']}, "
      f"restricted {res['verdict_inputs']['restricted_top50']})")
