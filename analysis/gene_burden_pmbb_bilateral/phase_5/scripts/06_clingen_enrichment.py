#!/usr/bin/env python3
"""
Phase 5, step 5.5b -- are ClinGen hearing-loss genes nearer the top than chance?

THE QUESTION THIS SETTLES. Four known deafness genes moved the wrong way when
the phenotype was restricted -- SIX1 from rank 4 to 6,703. Two readings, and
picking between them by eye is not possible:

  (a) the restriction removed the cases that carried real signal
  (b) neither arm has signal, so ranks wander and four genes picked by hand
      moving together means nothing

Four genes cannot distinguish these. A hundred can. The test runs on BOTH arms,
because the comparison is the point: if the broad phenotype shows enrichment and
the restricted one does not, (a) gains support. If neither does, (b) does.

TWO TESTS, because they fail differently.

  Mann-Whitney U on -log10(p), one-sided. Uses every gene, so a broad weak shift
  shows up. Blind to a handful of genes at the very top, which is what a
  rare-variant study would actually produce.

  Top-N counts against the hypergeometric. Sees exactly that -- genes at the
  top -- and is blind to a broad shift. N is reported at several values rather
  than one, because picking the N that looks best after seeing the data is how
  enrichment gets manufactured.

Gene set: ClinGen Hearing Loss GCEP, each gene at its best classification.
Primary Definitive+Strong, sensitivity +Moderate -- the tiers Andre set
provisionally on 2026-08-26.

Input : phase_5/data/clingen_hl_genes.tsv, both arms' omnibus tables
Output: phase_5/results/06_clingen_enrichment.json
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, hypergeom

HERE = Path(__file__).resolve().parents[1]
BROAD = Path("/project/hall/analysis/hearing-loss-genomics/analysis/"
             "elena_replication/phase_4/results")
OUT = HERE / "results"

cg = pd.read_csv(HERE / "data/clingen_hl_genes.tsv", sep="\t")
TIERS = {
    "definitive_strong": set(cg[cg.best_classification.isin(["Definitive", "Strong"])].gene),
    "plus_moderate": set(cg[cg.best_classification.isin(
        ["Definitive", "Strong", "Moderate"])].gene),
}

ARMS = {
    "restricted_bilateral_sn": lambda c: pd.read_csv(
        HERE / f"results/omnibus_{c}.tsv", sep="\t")[["Region", "omnibus_p"]],
    "broad_SO_396": lambda c: pd.read_csv(
        BROAD / f"omnibus_{c}.tsv", sep="\t")[["Region", "omnibus_p"]],
}

report = {"gene_set_sizes": {k: len(v) for k, v in TIERS.items()}, "arms": {}}

for arm, load in ARMS.items():
    report["arms"][arm] = {}
    for cohort in ["combined", "EUR", "AFR"]:
        d = load(cohort).dropna()
        d["score"] = -np.log10(pd.to_numeric(d.omnibus_p).clip(lower=1e-300))
        d = d.sort_values("omnibus_p").reset_index(drop=True)
        d["rank"] = np.arange(1, len(d) + 1)
        res = {}
        for tier, genes in TIERS.items():
            inset = d.Region.isin(genes)
            n_in = int(inset.sum())
            if n_in < 10:
                res[tier] = {"tested": n_in, "note": "too few tested to test"}
                continue
            u, p = mannwhitneyu(d.score[inset], d.score[~inset], alternative="greater")
            topn = {}
            for N in (50, 100, 250, 500, 1000):
                hit = int((d["rank"][inset] <= N).sum())
                exp = N * n_in / len(d)
                topn[N] = {
                    "observed": hit, "expected": round(exp, 2),
                    "p_hypergeom": float(hypergeom.sf(hit - 1, len(d), n_in, N)),
                }
            res[tier] = {
                "tested": n_in, "of_set": len(genes),
                "median_rank_in_set": int(d["rank"][inset].median()),
                "median_rank_overall": int(d["rank"].median()),
                "mannwhitney_p": float(p),
                "top_n": topn,
            }
        report["arms"][arm][cohort] = res

OUT.mkdir(parents=True, exist_ok=True)
(OUT / "06_clingen_enrichment.json").write_text(json.dumps(report, indent=2))

print(f"ClinGen HL: {len(TIERS['definitive_strong'])} Definitive+Strong, "
      f"{len(TIERS['plus_moderate'])} +Moderate\n")
for tier in ["definitive_strong", "plus_moderate"]:
    print(f"=== {tier} ===")
    print(f"{'arm':26s} {'cohort':9s} {'tested':>6s} {'med rank':>9s} "
          f"{'of':>6s} {'MW p':>8s} {'top100 obs/exp':>15s}")
    for arm in ARMS:
        for cohort in ["combined", "EUR", "AFR"]:
            r = report["arms"][arm][cohort][tier]
            if "mannwhitney_p" not in r:
                continue
            t = r["top_n"][100]
            print(f"{arm:26s} {cohort:9s} {r['tested']:>6} "
                  f"{r['median_rank_in_set']:>9,} {r['median_rank_overall']:>6,} "
                  f"{r['mannwhitney_p']:>8.3f} {t['observed']:>7}/{t['expected']:<7}")
    print()
