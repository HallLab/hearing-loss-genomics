#!/usr/bin/env python3
"""
Add Benjamini-Hochberg FDR to the gene burden arms, as asked in the 2026-10-14 meeting.

WHAT WAS ASKED, and the scope, which is deliberate:

  "let's keep Bonferroni, but also include FDR [...] add FDR P of less than 0.1"
  "you don't have to do FDR and ExWAS. ExWAS just stays as that p-value, but FDR
   you can add to gene burden. You do not do FDR for ExWAS."
                                                        -- Nikki, 2026-10-14

So: both bars on the burden arms, nothing on the ExWAS arms. That split is principled
rather than arbitrary. FDR is a threshold choice, and the two analyses genuinely use
different thresholds -- the burden corrects across genes, while the ExWAS is judged
against the 5e-8 convention, which is not a count of anything we did.

WHAT THE TWO BARS CONTROL, since they answer different questions:

  Bonferroni   the chance of even ONE false positive across the whole analysis
  FDR          the expected PROPORTION of false positives among whatever is called

FDR is the looser and the more useful of the two when the point is to see whether
anything is worth following up. Hence q < 0.1 rather than 0.05.

Written into each arm's own results folder so the arms stay self-contained.
"""
import json
from pathlib import Path

import pandas as pd
from statsmodels.stats.multitest import multipletests

ROOT = Path("/project/hall/analysis/hearing-loss-genomics/analysis")
Q_BAR = 0.10
ARMS = {
    "Replication": ROOT / "gene_burden_elena_replication/phase_4/results",
    "Bilateral": ROOT / "gene_burden_pmbb_bilateral/phase_5/results",
}
COHORTS = ["combined", "EUR", "AFR"]


def main():
    summary = {"q_bar": Q_BAR, "arms": {}}
    print(f"Benjamini-Hochberg on the gene omnibus, q < {Q_BAR}\n")

    for arm, base in ARMS.items():
        summary["arms"][arm] = {}
        for c in COHORTS:
            src = base / f"omnibus_{c}.tsv"
            d = pd.read_csv(src, sep="\t")
            p = pd.to_numeric(d["omnibus_p"], errors="coerce")
            ok = p.notna() & (p > 0) & (p <= 1)
            if not ok.all():
                print(f"  NOTE {arm}/{c}: {int((~ok).sum())} rows without a usable "
                      f"p-value, excluded from the FDR")

            rej, q, _, _ = multipletests(p[ok].values, alpha=Q_BAR, method="fdr_bh")
            d.loc[ok, "fdr_q"] = q
            d.loc[ok, "fdr_significant"] = rej
            bar = float(d["bonferroni_bar"].iloc[0])
            d["bonferroni_significant"] = p < bar

            dest = base / f"omnibus_{c}_with_fdr.tsv"
            d.sort_values("omnibus_p").to_csv(dest, sep="\t", index=False)

            top = d.loc[ok].nsmallest(5, "omnibus_p")[
                ["Region", "omnibus_p", "fdr_q"]]
            summary["arms"][arm][c] = {
                "genes": int(ok.sum()),
                "min_p": float(p[ok].min()),
                "min_q": float(q.min()),
                "bonferroni_bar": bar,
                "bonferroni_hits": int((p[ok] < bar).sum()),
                "fdr_hits": int(rej.sum()),
                "top5": [{"gene": r.Region, "p": float(r.omnibus_p),
                          "q": float(r.fdr_q)} for r in top.itertuples()],
            }
            s = summary["arms"][arm][c]
            print(f"  {arm:<12}{c:<10} genes {s['genes']:>7,}   "
                  f"min p {s['min_p']:.2e}   min q {s['min_q']:.4f}   "
                  f"Bonf {s['bonferroni_hits']}   FDR {s['fdr_hits']}")
        print()

    print("Top 5 genes per cohort, with both measures\n")
    for arm in ARMS:
        for c in COHORTS:
            s = summary["arms"][arm][c]
            print(f"  {arm} / {c}")
            for t in s["top5"]:
                print(f"      {t['gene']:<14} p = {t['p']:.2e}   q = {t['q']:.4f}")
        print()

    dest = ROOT / "synthesis/results/fdr.json"
    dest.write_text(json.dumps(summary, indent=2) + "\n")
    print(f"written: {dest.relative_to(ROOT)}  "
          f"and omnibus_<cohort>_with_fdr.tsv in each arm")


if __name__ == "__main__":
    main()
