#!/usr/bin/env python3
"""
Phase 3 — principal components: how many, and why.

Shared code for the notebook. Kept as a module so the notebook stays readable and
the numbers in the write-ups come from one place.

Three eigenvalue sources, because the pipeline used three different PCAs:
  combined   the release's own exome PCA   (Exome/PCA/combined/...eigenvalues.tsv)
  EUR        a within-ancestry PCA Elena ran  (EUR_variance_explained.tsv)
  AFR        a within-ancestry PCA Elena ran  (AFR_variance_explained.tsv)

Palette: the dataviz reference instance. Scatter plots are an all-pairs form, where
only the first three categorical slots validate, so grouped scatters cap at three
series and everything else is drawn as small multiples with one hue per panel.
"""
from pathlib import Path
import pandas as pd

REL  = Path("/static/PMBB/PMBB-Release-2026-4.0/Exome")
PIPE = Path("/project/hall/analysis/hearing-loss-genomics/analysis/elena/rarevariantExWAS")

# --- dataviz reference palette, light mode
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"]
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#8a8984", "#e4e3df"
CONTEXT = "#d8d7d2"          # background points in small multiples
WHAT_RAN = {"combined": 5, "EUR": 9, "AFR": 10}


SOURCES = {
    "combined": REL / "PCA/combined/PMBB-Release-2026-4.0_genetic_exome.no1KG.eigenvalues.tsv",
    "EUR":      PIPE / "EUR_variance_explained.tsv",
    "AFR":      PIPE / "AFR_variance_explained.tsv",
    "pcs_per_person": REL / "PCA/combined/PMBB-Release-2026-4.0_genetic_exome.commonsnps.samples_ancestries.tsv",
}


def sources():
    """Where every number in this notebook comes from, with owner and date.

    Printed in the notebook rather than left inside this module: the notebook's job
    is to let someone verify the claim, and a reader who cannot see which files were
    read cannot do that.

    Note that no PCA is computed anywhere here. These are eigenvalues that already
    existed -- one table published by PMBB, two produced by the pipeline itself. A
    scree plot is a view of a table, not a calculation.
    """
    import pwd, datetime
    rows = []
    for name, path in SOURCES.items():
        st = path.stat()
        rows.append({
            "what": name,
            "owner": pwd.getpwuid(st.st_uid).pw_name,
            "modified": datetime.datetime.fromtimestamp(st.st_mtime).strftime("%Y-%m-%d"),
            "path": str(path),
        })
    return pd.DataFrame(rows)


def eigenvalues():
    """PC -> eigenvalue, for each of the three PCAs. Returns tidy long form."""
    rows = []
    rel = pd.read_csv(SOURCES["combined"], sep="\t").iloc[0].astype(float)
    for i, v in enumerate(rel.values, 1):
        rows.append({"cohort": "combined", "PC": i, "eigenvalue": v})
    for grp in ["EUR", "AFR"]:
        d = pd.read_csv(SOURCES[grp], sep="\t")
        for _, r in d.iterrows():
            rows.append({"cohort": grp, "PC": int(r.PC), "eigenvalue": float(r.Eigenvalue)})
    df = pd.DataFrame(rows)
    df["variance_pct"] = df.groupby("cohort").eigenvalue.transform(lambda s: 100 * s / s.sum())
    df["drop_pct"] = df.groupby("cohort").eigenvalue.transform(lambda s: -100 * s.pct_change())
    return df


def elbow(df, cohort, flat_threshold=10.0):
    """First PC whose drop from the previous one falls below the threshold.

    The elbow is where the curve stops falling. Everything from that PC onward
    explains about as much as its neighbour, which is the signature of noise rather
    than population structure. Reported as 'keep up to the PC before it'.
    """
    d = df[df.cohort == cohort].sort_values("PC")
    for _, r in d.iterrows():
        if pd.notna(r.drop_pct) and r.drop_pct < flat_threshold:
            return int(r.PC) - 1
    return int(d.PC.max())


def pcs_with_ancestry():
    """Release exome PCs 1-4 with the ancestry label, for the scatter plots."""
    d = pd.read_csv(SOURCES["pcs_per_person"], sep="\t",
                    usecols=["IID", "Class", "PC1", "PC2", "PC3", "PC4"])
    d["Class"] = d.Class.fillna("UNKNOWN").str.replace(r"UNKNOWN\d", "UNKNOWN", regex=True)
    return d


def style(ax):
    """Recessive axes and grid; the data carries the ink."""
    ax.set_facecolor("#fcfcfb")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=9, length=3, color=GRID)
    ax.grid(True, color=GRID, linewidth=0.6, alpha=0.9)
    ax.set_axisbelow(True)
    return ax
