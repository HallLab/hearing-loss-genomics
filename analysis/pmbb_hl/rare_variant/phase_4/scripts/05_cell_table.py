#!/usr/bin/env python3
"""
Phase 4 -- one table with every cell of every arm, and the omnibus beside it.

WHY THIS EXISTS. The per-cell results are spread across nine `cells_*.tsv` and the
gene-level verdict across nine `omnibus_*.tsv`. Answering "where did this gene's signal
come from?" meant joining them by hand, which is how the FBXO25 question got asked. One
long table answers it by reading one file.

ROW = one (phenotype, cohort, gene, mask, max_MAF) cell. 1.4M rows, every p-value the
burden arm produced, nothing filtered.

BOTH OMNIBUS VARIANTS ARE CARRIED, deliberately:

  p_omnibus            what phase 4 ships and what every FINDINGS quotes: ACAT over
                       the cells left after the EXISTING deduplication, which is
                       flagged here as `in_omnibus`
  p_omnibus_distinct   ACAT after also collapsing cells that are the same test at a
                       looser frequency cutoff (`is_representative`)

EVERY ROW IS KEPT, including the cells the omnibus excluded, because the point of the
table is to show what was scored -- `in_omnibus` says what was then combined. The first
version of this script recomputed ACAT over all rows and did not reproduce the shipped
p-value; the check below now refuses rather than shipping a table that disagrees with
the files it summarises.

They differ because a MAF cutoff above every variant in a gene produces an identical
result at every looser cutoff. FBXO25 has 16 pLOF variants, all with MAC <= 8, so its
three pLOF cells are one test counted three times; its best p enters the combination
with weight 3/9 instead of 1/5.

Neither is invalid -- ACAT admits arbitrary non-negative weights -- but they are not the
same quantity, and the implicit weighting in the first is "how many cutoffs did this
test survive", which has no stated justification and systematically favours genes whose
signal is ultra-rare. Carrying both lets the choice be made on evidence rather than
inherited. Measured effect: no arm crosses any bar under either, in any cohort; the
largest single change is FBXO25 in bilateral EUR, q 0.113 -> 0.188.

The existing phase 4 deduplication is a DIFFERENT axis and is already applied upstream:
it removes a mask whose cell is identical to another mask's AT THE SAME max_MAF, which
is the case Nikki raised (pLOF == pLOF;pDM when a gene has no pDM variants). This one is
about the same mask across frequency cutoffs. Both can apply to one gene.

SORTED BY OMNIBUS RANK within each arm, not alphabetically, so the genes worth looking
at are at the top of the file.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[1]
RESULTS = HERE / "results"
OUT = RESULTS / "step2"
PHENOTYPES = ["broad", "bilateral", "tinnitus"]
COHORTS = ["combined", "EUR", "AFR"]

# The columns that identify one TEST. Two cells agreeing on all three are the same
# variants scored the same way -- not two pieces of evidence.
TEST_KEY = ["Group", "MAC", "Pvalue"]
# Phase 4's existing dedup key, replicated EXACTLY from 04_merge_results.py. It is a
# different axis: same max_MAF, different mask, identical result -- the pLOF ==
# pLOF;pDM collision. Rounding and dtype have to match or the flag drifts from the
# omnibus it is meant to explain.
MERGE_KEY = ["Region", "_k"]

RENAME = {
    "Region": "gene", "CHR": "chr", "Group": "mask", "max_MAF": "max_maf",
    "Pvalue": "p_cell", "Pvalue_Burden": "p_burden", "Pvalue_SKAT": "p_skat",
    "BETA_Burden": "beta_burden", "SE_Burden": "se_burden",
    "MAC": "mac", "MAC_case": "mac_case", "MAC_control": "mac_control",
    "Number_rare": "n_rare", "Number_ultra_rare": "n_ultra_rare",
}
ORDER = ["phenotype", "cohort", "gene", "chr", "mask", "max_maf",
         "p_cell", "p_burden", "p_skat", "beta_burden", "se_burden",
         "mac", "mac_case", "mac_control", "n_rare", "n_ultra_rare",
         "in_omnibus", "n_cells_identical", "is_representative",
         "p_omnibus", "p_omnibus_distinct", "omnibus_rank",
         "bonferroni_bar", "fdr_q", "fdr_q_distinct"]


def acat(p):
    """Cauchy combination. A MEAN of tan-transformed p, so weights are implicit in
    how many times a value appears."""
    p = np.clip(np.asarray(p, float), 1e-300, 1 - 1e-16)
    return float(0.5 - np.arctan(np.tan((0.5 - p) * np.pi).mean()) / np.pi)


def acat_by_gene(df):
    p = np.clip(df.Pvalue.to_numpy(float), 1e-300, 1 - 1e-16)
    m = df.assign(_t=np.tan((0.5 - p) * np.pi)).groupby("Region")._t.mean()
    return 0.5 - np.arctan(m) / np.pi


def bh(p):
    """Benjamini-Hochberg step-up. The running minimum propagates from the largest
    rank downward, so a gene's q can be set by a worse-ranked gene."""
    o = np.argsort(p)
    ps = np.asarray(p, float)[o]
    n = len(ps)
    q = np.minimum.accumulate((ps * n / np.arange(1, n + 1))[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.clip(q, 0, 1)
    return out


def one_arm(pheno, cohort):
    tag = f"{pheno}_{cohort}"
    cells = pd.read_csv(RESULTS / f"cells_{tag}.tsv", sep="\t")
    om = pd.read_csv(RESULTS / f"omnibus_{tag}.tsv", sep="\t")

    # The join must be total in both directions. A gene in one file and not the other
    # means a cell was dropped or an omnibus was invented, and either is a defect that
    # a merged table would hide rather than show.
    cg, og = set(cells.Region), set(om.Region)
    if cg != og:
        sys.exit(f"REFUSING: {tag} has {len(cg - og):,} genes with cells but no "
                 f"omnibus and {len(og - cg):,} with an omnibus but no cells. "
                 f"The two files disagree about which genes were tested.")

    # Which cells the shipped omnibus combined.
    cells["_k"] = list(zip(cells.max_MAF.astype(str),
                           cells.MAC.round(10), cells.Pvalue.round(12)))
    cells["in_omnibus"] = ~cells.duplicated(subset=MERGE_KEY, keep="first")

    # Identical tests across frequency cutoffs, among the cells that were combined.
    # The representative is the TIGHTEST cutoff producing the result -- the looser
    # ones add no variants.
    kept = cells[cells.in_omnibus]
    grp = kept.groupby(["Region"] + TEST_KEY, sort=False)
    cells["n_cells_identical"] = grp.max_MAF.transform("size").reindex(cells.index)
    cells["is_representative"] = (
        cells.in_omnibus & (cells.max_MAF == grp.max_MAF.transform("min")
                            .reindex(cells.index)))
    cells["n_cells_identical"] = cells.n_cells_identical.fillna(0).astype(int)

    p_dist = acat_by_gene(cells[cells.is_representative]).rename("p_omnibus_distinct")
    if len(p_dist) != len(om):
        sys.exit(f"REFUSING: {tag} collapsed to {len(p_dist):,} genes from "
                 f"{len(om):,}. Collapsing identical cells must not remove a gene.")

    # The flag has to reproduce the shipped omnibus exactly, or it is describing
    # some other combination than the one the project reports.
    rec = acat_by_gene(cells[cells.in_omnibus]).reindex(om.Region.values)
    rel = (np.abs(rec.to_numpy() - om.omnibus_p.to_numpy())
           / om.omnibus_p.to_numpy()).max()
    if not rel < 1e-12:
        sys.exit(f"REFUSING: {tag} -- ACAT over the cells flagged `in_omnibus` does "
                 f"not reproduce omnibus_{tag}.tsv (max relative error {rel:.2e}). "
                 f"The dedup key here has drifted from 04_merge_results.py.")
    used = cells.groupby("Region").in_omnibus.sum().reindex(om.Region.values)
    if not (used.to_numpy() == om.cells_used.to_numpy()).all():
        sys.exit(f"REFUSING: {tag} -- `in_omnibus` counts disagree with the "
                 f"cells_used column of omnibus_{tag}.tsv.")

    om = om.merge(p_dist, left_on="Region", right_index=True, how="left")
    om["fdr_q_distinct"] = bh(om.p_omnibus_distinct.to_numpy())
    om["omnibus_rank"] = om.omnibus_p.rank(method="min").astype(int)

    d = cells.merge(
        om[["Region", "omnibus_p", "p_omnibus_distinct", "omnibus_rank",
            "bonferroni_bar", "fdr_q", "fdr_q_distinct"]],
        on="Region", how="left", validate="many_to_one")
    d = d.rename(columns={**RENAME, "omnibus_p": "p_omnibus"})
    d["phenotype"], d["cohort"] = pheno, cohort
    # SAIGE writes allele and variant counts as floats, so they print as "18.0".
    # They are counts; printing them as counts costs nothing and makes the file
    # both smaller and readable.
    for c in ("mac", "mac_case", "mac_control", "n_rare", "n_ultra_rare", "chr"):
        d[c] = d[c].astype("Int64")
    return d[ORDER].sort_values(["omnibus_rank", "mask", "max_maf"])


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    frames = []
    print(f"  {'arm':<22}{'rows':>12}{'genes':>9}{'cells w/ a twin':>18}"
          f"{'min q':>9}{'min q dist':>12}")
    for pheno in PHENOTYPES:
        for cohort in COHORTS:
            d = one_arm(pheno, cohort)
            frames.append(d)
            twins = int((d.n_cells_identical > 1).sum())
            print(f"  {pheno + '_' + cohort:<22}{len(d):>12,}{d.gene.nunique():>9,}"
                  f"{twins:>13,} ({twins / len(d):.0%}){d.fdr_q.min():>9.3f}"
                  f"{d.fdr_q_distinct.min():>12.3f}")

    all_cells = pd.concat(frames, ignore_index=True)
    out = OUT / "all_cells.tsv"
    # NO float_format. At 6 significant digits the p-values no longer round-trip and
    # the table disagreed with omnibus_*.tsv at 5e-06 relative -- small, but this file
    # exists to be the record.
    all_cells.to_csv(out, sep="\t", index=False)
    print(f"\n  {out.relative_to(HERE)}  "
          f"{len(all_cells):,} rows, {out.stat().st_size / 1e6:.0f} MB")

    # A gene-level companion, one row per (arm, gene), for anyone who wants the
    # verdict without the 1.4M cells.
    genes = (all_cells.groupby(["phenotype", "cohort", "gene"], as_index=False)
             .agg(chr=("chr", "first"), omnibus_rank=("omnibus_rank", "first"),
                  p_omnibus=("p_omnibus", "first"),
                  p_omnibus_distinct=("p_omnibus_distinct", "first"),
                  fdr_q=("fdr_q", "first"), fdr_q_distinct=("fdr_q_distinct", "first"),
                  bonferroni_bar=("bonferroni_bar", "first"),
                  best_cell_p=("p_cell", "min"), cells=("p_cell", "size"),
                  cells_used=("in_omnibus", "sum"),
                  distinct_tests=("is_representative", "sum"),
                  mac_widest=("mac", "max")))
    # Ties are broken toward the TIGHTEST cutoff. A gene whose variants are all
    # ultra-rare has the same p at every cutoff, and reporting the loosest of them as
    # "where the signal is" points at a frequency band that contributed nothing.
    best = (all_cells.sort_values(["p_cell", "max_maf"])
            .groupby(["phenotype", "cohort", "gene"])
            .agg(best_mask=("mask", "first"), best_max_maf=("max_maf", "first")))
    genes = genes.merge(best, on=["phenotype", "cohort", "gene"], how="left")
    genes = genes.sort_values(["phenotype", "cohort", "omnibus_rank"])
    gout = OUT / "genes.tsv"
    genes.to_csv(gout, sep="\t", index=False)
    print(f"  {gout.relative_to(HERE)}  {len(genes):,} rows, "
          f"{gout.stat().st_size / 1e6:.0f} MB")


if __name__ == "__main__":
    main()
