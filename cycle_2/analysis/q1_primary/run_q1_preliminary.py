#!/usr/bin/env python3
"""
Cycle 2 / Q1 preliminary — ClinGen primary gene set applied to Elena's existing PMBB v4
SAIGE-GENE+ results. No new compute: the genes were already tested exome-wide.

Pre-declared before looking at any result (see cycle_2/README.md §3):
  - primary set   : ClinGen HL GCEP Definitive+Strong mapped to phecode 389
  - primary p     : SAIGE-GENE+ Cauchy omnibus (falls back to the single group when absent)
  - FDR           : Benjamini-Hochberg WITHIN each (ancestry, mask, MAF) stratum
  - MOI split     : AD vs AR declared in advance; XL reported separately (chr1-22 pipeline)
  - calibration   : lambda_GC computed exome-wide, per stratum
"""
import json, pathlib, numpy as np, pandas as pd
from scipy.stats import chi2

ROOT = pathlib.Path("/project/hall/analysis/hearing-loss-genomics")
CLIN = ROOT / "cycle_2/data/clingen"
RES  = ROOT / "analysis/elena/HL_only_rarevariant/newmasks_combined_gene_results.csv"
OUT  = ROOT / "cycle_2/results/q1_preliminary"
OUT.mkdir(parents=True, exist_ok=True)

# ---------- 1. primary gene set -------------------------------------------------
summary = pd.read_csv(CLIN / "hl_phecode_summary.tsv", sep="\t")
row389 = summary.loc[summary.PHECODE.astype(float) == 389.0].iloc[0]
primary = sorted(set(row389.GENES.split(";")))

gcep = pd.read_csv(CLIN / "clingen_hl_gcep.tsv", sep="\t")
gcep.columns = [c.strip() for c in gcep.columns]
ds = gcep[gcep["CLASSIFICATION"].isin(["Definitive", "Strong"])]
moi = (ds.groupby("GENE SYMBOL")["MOI"]
         .agg(lambda s: "/".join(sorted(set(s)))).to_dict())

gene_meta = pd.DataFrame({"Region": primary})
gene_meta["MOI"] = gene_meta.Region.map(moi).fillna("?")
gene_meta["MOI_class"] = np.where(gene_meta.MOI.str.contains("XL"), "XL",
                          np.where(gene_meta.MOI.str.contains("AD"), "AD",
                          np.where(gene_meta.MOI.str.contains("AR"), "AR", "?")))

# ---------- 2. results, with the malformed-header fix ---------------------------
hdr = open(RES).readline().rstrip("\n").split(",")
COLS = hdr[:16]                      # fields 17-41 are junk from a bad concatenation
lost = dict(zip(COLS, hdr[16:32]))   # one real data row got glued into the header
df = pd.read_csv(RES, skiprows=1, header=None, usecols=range(16),
                 names=COLS, low_memory=False)
n_raw = len(df)
df = df[df.Region != "-"]
df = df[df.Region != "Region"]
df["Pvalue"] = pd.to_numeric(df.Pvalue, errors="coerce")
df = df.dropna(subset=["Pvalue"])
df = df[(df.Pvalue > 0) & (df.Pvalue <= 1)]
df = df.drop_duplicates(subset=["Region", "Group", "max_MAF", "Ancestry", "Mask"])

STRATA = ["Ancestry", "Mask", "max_MAF"]

# one row per (stratum, gene): prefer the SAIGE-GENE+ Cauchy omnibus when present
df["_rank"] = np.where(df.Group == "Cauchy", 0, 1)
omni = (df.sort_values(STRATA + ["Region", "_rank"])
          .drop_duplicates(subset=STRATA + ["Region"], keep="first")
          .drop(columns="_rank"))

# ---------- 3. calibration, exome-wide -----------------------------------------
def lam(p):
    p = np.asarray(p, float)
    p = p[(p > 0) & (p <= 1)]
    return np.median(chi2.isf(p, 1)) / chi2.ppf(0.5, 1) if len(p) else np.nan

cal = (omni.groupby(STRATA)
           .agg(n_genes_exomewide=("Region", "nunique"),
                lambda_GC=("Pvalue", lam))
           .reset_index().sort_values(STRATA))

# ---------- 4. primary set + BH within stratum ---------------------------------
prim = omni[omni.Region.isin(primary)].merge(gene_meta, on="Region", how="left")

def bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p); n = len(p)
    q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.minimum(q, 1.0), n

parts = []
for _, g in prim.groupby(STRATA, sort=False):
    g = g.copy()
    g["FDR"], g["n_in_stratum"] = bh(g.Pvalue.values)[0], len(g)
    parts.append(g)
prim = pd.concat(parts).sort_values("Pvalue")

# ---------- 5. write ------------------------------------------------------------
prim.to_csv(OUT / "q1_primary_full.csv", index=False)
cal.to_csv(OUT / "calibration_lambda_gc.csv", index=False)
gene_meta.to_csv(OUT / "primary_gene_set_61.csv", index=False)

hits = prim[prim.FDR < 0.05]
hits.to_csv(OUT / "q1_primary_fdr05_hits.csv", index=False)

manifest = {
    "results_file": str(RES),
    "rows_raw": int(n_raw),
    "rows_used_after_qc": int(len(df)),
    "malformed_header_note": "file has 41 fields; only first 16 are real. One data row was "
                             "glued into the header and is therefore missing from the table.",
    "row_lost_in_header": lost,
    "primary_set": {"source": "ClinGen HL GCEP Definitive+Strong -> phecode 389",
                    "n_genes": len(primary)},
    "moi_breakdown": gene_meta.MOI_class.value_counts().to_dict(),
    "genes_absent_from_results": sorted(set(primary) - set(omni.Region.unique())),
    "n_strata": int(len(cal)),
    "n_hits_fdr05": int(len(hits)),
}
(OUT / "manifest.json").write_text(json.dumps(manifest, indent=2))

pd.set_option("display.width", 200)
print("=== calibration (exome-wide lambda_GC) ===")
print(cal.to_string(index=False))
print("\n=== primary set MOI ===");  print(gene_meta.MOI_class.value_counts().to_string())
print("\ngenes in primary set absent from results:", manifest["genes_absent_from_results"])
print("\n=== top 15 primary-set results (by p) ===")
print(prim.head(15)[["Region","MOI_class","Ancestry","Mask","max_MAF","Pvalue",
                     "FDR","BETA_Burden","MAC_case","MAC_control"]].to_string(index=False))
print(f"\nhits at FDR<0.05 within stratum: {len(hits)}")
if len(hits):
    print(hits[["Region","MOI_class","Ancestry","Mask","max_MAF","Pvalue","FDR","BETA_Burden"]].to_string(index=False))
