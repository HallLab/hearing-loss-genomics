#!/usr/bin/env python3
"""
Cycle 2 · step 01 — hearing-loss / tinnitus / combined phenotypes from PMBB v4.

Built independently inside cycle_2, reading only the institutional release at
/static/PMBB/PMBB-Release-2026-4.0/. No dependency on any other analyst's outputs.

CREDIT: the fact that PMBB v4 relocated standard tinnitus ICD codes into the OMOP
`observation` table was established by Nikki Palmiero and Elena's notebook
`PMBB_4_PhecodeX_Hearing_Tinnitus.ipynb`. That discovery is load-bearing here and is
documented nowhere else. Verified independently before adopting (2026-08-26):

    tinnitus events in conditions_phecode_12.txt (389.4) :  3,016
    tinnitus events in conditions_phecode_x.txt (SO_397.1):  3,016   <- pulsatile only
    tinnitus events in observation.txt (388.3x + H93.1x)  : 25,094   <- the real bulk

Neither phecode file alone is usable for tinnitus in v4. Hearing impairment is fine in
either (SO_396 = 134,917 events).

DIFFERENCES from that notebook, deliberate:
  1. adds the combined HL-and/or-tinnitus phenotype (2026-07-01 decision, never executed);
  2. intersects with the EXOME sample list, not the imputed LD-pruned .fam. The imputed
     list drops 517 exome-sequenced participants and admits 85 without exome data;
  3. versioned script with a manifest instead of a notebook.

Case/control rules (unchanged, they are sound):
  case                            : qualifying evidence on >= 2 distinct dates
  excluded_one_date               : evidence on exactly 1 date
  excluded_related_ear_phenotype  : no target evidence but some other ear-family evidence
  control                         : neither
"""
import json, subprocess, sys
from pathlib import Path
import numpy as np, pandas as pd

REL   = Path("/static/PMBB/PMBB-Release-2026-4.0")
PHENO = REL / "Phenotype/4.0"
PFX   = "PMBB-Release-2026-4.0_phenotype_"
OUT   = Path("/project/hall/analysis/hearing-loss-genomics/cycle_2/data/phenotype")
OUT.mkdir(parents=True, exist_ok=True)

PHECODEX = PHENO / f"{PFX}conditions_phecode_x.txt"
OBS      = PHENO / f"{PFX}observation.txt"
COV      = PHENO / f"{PFX}covariates.txt"
EXOME_ANC = REL / "Exome/PCA/combined/PMBB-Release-2026-4.0_genetic_exome.commonsnps.samples_ancestries.tsv"

TINNITUS_OBS_CODES = {"388.3","388.30","388.31","388.32",
                      "H93.1","H93.11","H93.12","H93.13","H93.19"}

# ---------------------------------------------------------------- extraction (cached)
AWK_EAR = r"""
BEGIN{FS=OFS="\t"}
NR==1{for(i=1;i<=NF;i++) h[$i]=i;
      print "condition_occurrence_id","person_id","condition_start_date","condition_source_value"; next}
$h["condition_source_value"] ~ /^SO_39[0-9](\.|$)/ {
      print $h["condition_occurrence_id"],$h["person_id"],$h["condition_start_date"],$h["condition_source_value"]}
"""

AWK_TIN = r"""
BEGIN{FS=OFS="\t"}
NR==1{for(i=1;i<=NF;i++) h[$i]=i;
      print "observation_id","person_id","observation_date","observation_source_value"; next}
{v=$h["observation_source_value"]}
v=="388.3"||v=="388.30"||v=="388.31"||v=="388.32"||v=="H93.1"||v=="H93.11"||v=="H93.12"||v=="H93.13"||v=="H93.19"{
      print $h["observation_id"],$h["person_id"],$h["observation_date"],v}
"""

def awk_extract(src, dst, prog):
    """Pre-filter a very large TSV with awk; pandas on 12 GB is not worth it."""
    if dst.exists() and dst.stat().st_size > 0:
        print(f"  [cache] {dst.name}"); return
    print(f"  [awk ] {src.name} -> {dst.name}", flush=True)
    with open(dst, "w") as fh:
        subprocess.run(["awk", prog, str(src)], stdout=fh, check=True)

print("== extracting ==", flush=True)
ear_raw = OUT / "_ear_family_phecodex.tsv"
awk_extract(PHECODEX, ear_raw, AWK_EAR)
tin_raw = OUT / "_tinnitus_observation.tsv"
awk_extract(OBS, tin_raw, AWK_TIN)

# ---------------------------------------------------------------- load
print("== building ==")
ear = pd.read_csv(ear_raw, sep="\t", dtype="string")
ear["date"] = pd.to_datetime(ear.condition_start_date, errors="coerce").dt.date
ear["code"] = ear.condition_source_value.str.strip().str.upper()
ear = ear.dropna(subset=["person_id","date"])

tin = pd.read_csv(tin_raw, sep="\t", dtype="string")
tin["date"] = pd.to_datetime(tin.observation_date, errors="coerce").dt.date
tin["code"] = tin.observation_source_value.str.strip().str.upper()
tin = tin.dropna(subset=["person_id","date"])

def dates(df):
    return df[["person_id","date"]].drop_duplicates()

hl_d  = dates(ear[ear.code.str.match(r"^SO_396(\.|$)")])
tin_d = pd.concat([dates(ear[ear.code.str.match(r"^SO_397(\.|$)")]), dates(tin)]).drop_duplicates()
both_d = pd.concat([hl_d, tin_d]).drop_duplicates()
related_d = pd.concat([dates(ear), dates(tin)]).drop_duplicates()   # any ear-family evidence

TARGETS = {"hearing_impairment": hl_d, "tinnitus": tin_d, "hl_or_tinnitus": both_d}

# ---------------------------------------------------------------- sample frame
anc = pd.read_csv(EXOME_ANC, sep="\t", dtype="string")
id_col  = next(c for c in anc.columns if c.upper() in {"IID","PERSON_ID","SAMPLE_ID","#IID"})
anc_col = next(c for c in anc.columns if "ancest" in c.lower() or c.lower() in {"pop","class"})
anc = anc.rename(columns={id_col:"person_id", anc_col:"ancestry"})[["person_id","ancestry"]]
anc["person_id"] = anc.person_id.str.strip()

cov = pd.read_csv(COV, sep="\t", low_memory=False)
cov.columns = cov.columns.str.strip()
cov["person_id"] = cov.person_id.astype(str).str.strip()
cov = cov.drop_duplicates("person_id")

base = anc.merge(cov, on="person_id", how="inner")
rel_counts = related_d.groupby("person_id").date.nunique().rename("n_related_ear_dates")

rows = []
for name, tdates in TARGETS.items():
    tc = tdates.groupby("person_id").date.nunique().rename("n_target_dates")
    r = base[["person_id","ancestry"]].copy()
    r["phenotype"] = name
    r = r.join(tc, on="person_id").join(rel_counts, on="person_id")
    r[["n_target_dates","n_related_ear_dates"]] = r[["n_target_dates","n_related_ear_dates"]].fillna(0).astype(int)
    r["status"] = np.select(
        [r.n_target_dates >= 2, r.n_target_dates == 1, r.n_related_ear_dates >= 1],
        ["case", "excluded_one_date", "excluded_related_ear_phenotype"], default="control")
    rows.append(r)
status = pd.concat(rows, ignore_index=True)

# ---------------------------------------------------------------- QC (hard assertions)
assert (status.loc[status.status=="case","n_target_dates"] >= 2).all()
assert (status.loc[status.status=="excluded_one_date","n_target_dates"] == 1).all()
assert (status.loc[status.status=="control","n_target_dates"] == 0).all()
assert (status.loc[status.status=="control","n_related_ear_dates"] == 0).all()
for p in TARGETS:
    s = status[status.phenotype==p]
    assert len(s) == s.person_id.nunique(), f"duplicate person in {p}"
# the combined phenotype must contain every HL case and every tinnitus case
c = {p: set(status[(status.phenotype==p)&(status.status=="case")].person_id) for p in TARGETS}
assert c["hearing_impairment"] <= c["hl_or_tinnitus"] and c["tinnitus"] <= c["hl_or_tinnitus"]

# ---------------------------------------------------------------- write
status.to_csv(OUT / "phenotype_status.csv.gz", index=False)
summary = (status.groupby(["phenotype","status"]).size().unstack(fill_value=0)
                 .reindex(columns=["case","control","excluded_one_date",
                                   "excluded_related_ear_phenotype"], fill_value=0))
summary["analysis_N"] = summary.case + summary.control
summary["case_rate_%"] = (100*summary.case/summary.analysis_N).round(2)
summary.to_csv(OUT / "phenotype_summary.csv")

by_anc = (status[status.status.isin(["case","control"])]
          .groupby(["phenotype","ancestry"])
          .agg(N=("person_id","size"), cases=("status", lambda s:(s=="case").sum()))
          .reset_index())
by_anc["case_rate_%"] = (100*by_anc.cases/by_anc.N).round(2)
by_anc.to_csv(OUT / "phenotype_by_ancestry.csv", index=False)

(OUT / "manifest.json").write_text(json.dumps({
    "release": str(REL),
    "sources": {"phecodex": str(PHECODEX), "observation": str(OBS),
                "covariates": str(COV), "exome_samples": str(EXOME_ANC)},
    "tinnitus_observation_codes": sorted(TINNITUS_OBS_CODES),
    "rules": "case >=2 distinct dates; 1 date excluded; related ear-family evidence excluded from controls",
    "sample_frame": "exome sample list (not the imputed LD-pruned .fam)",
    "n_exome_samples": int(len(base)),
    "summary": json.loads(summary.to_json(orient="index")),
}, indent=2))

pd.set_option("display.width", 200)
print("\n=== phenotype summary ===");  print(summary.to_string())
print("\n=== by ancestry (case/control only) ==="); print(by_anc.to_string(index=False))
print(f"\nwritten to {OUT}")
