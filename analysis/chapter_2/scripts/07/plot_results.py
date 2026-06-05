#!/usr/bin/env python3
# Figure for the ZNF175 v1-vs-v2 burden result (NB 07). Saves PNG to results/07/.
from pathlib import Path
import numpy as np, pandas as pd, statsmodels.api as sm
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

R07 = Path("/project/hall/analysis/hearing-loss-genomics/analysis/chapter_2/results/07")
PCS = [f"PC{i}" for i in range(1, 11)]
d1 = pd.read_csv(R07 / "v1_analysis.csv")
d2 = pd.read_csv(R07 / "v2_analysis.csv")

def fit_or(d):
    dd = d.dropna(subset=["carrier", "age", "sex_male", "tinnitus"] + PCS).copy()
    dd["age2"] = dd["age"] ** 2
    X = sm.add_constant(dd[["carrier", "age", "age2", "sex_male"] + PCS])
    r = sm.Logit(dd["tinnitus"], X).fit(disp=0, maxiter=200)
    b, se = r.params["carrier"], r.bse["carrier"]
    return np.exp(b), np.exp(b - 1.96 * se), np.exp(b + 1.96 * se), r.pvalues["carrier"]

def rates(d):
    c = d[d.carrier == 1]; n = d[d.carrier == 0]
    return dict(cr=c.tinnitus.mean()*100, cc=int(c.tinnitus.sum()), cn=len(c),
                nr=n.tinnitus.mean()*100, nc=int(n.tinnitus.sum()), nn=len(n))

r1, r2 = rates(d1), rates(d2)
or1 = fit_or(d1); or2 = fit_or(d2)

fig, ax = plt.subplots(1, 3, figsize=(15.5, 4.8))
C_CAR, C_NON = "#c0392b", "#7f8c8d"

# --- Panel A: tinnitus rate among carriers vs non-carriers ---
a = ax[0]
x = np.arange(2); w = 0.38
a.bar(x - w/2, [r1["cr"], r2["cr"]], w, color=C_CAR, label="ZNF175 carrier")
a.bar(x + w/2, [r1["nr"], r2["nr"]], w, color=C_NON, label="non-carrier")
for i, r in enumerate([r1, r2]):
    a.text(i - w/2, r["cr"], f'{r["cr"]:.1f}%\n({r["cc"]}/{r["cn"]})', ha="center", va="bottom", fontsize=9, color=C_CAR)
    a.text(i + w/2, r["nr"], f'{r["nr"]:.1f}%\n({r["nc"]}/{r["nn"]})', ha="center", va="bottom", fontsize=8, color=C_NON)
a.set_xticks(x); a.set_xticklabels(["v1 (~11K)", "v2 (~44K)"])
a.set_ylabel("tinnitus rate (%)"); a.set_ylim(0, 19)
a.set_title("A. Tinnitus enrichment among carriers"); a.legend(frameon=False, fontsize=9)

# --- Panel B: forest plot of carrier OR (adjusted logistic, 95% CI) ---
b = ax[1]
labels = [f"v2 (~44K)\nOR={or2[0]:.1f}, p={or2[3]:.1e}", f"v1 (~11K)\nOR={or1[0]:.1f}, p={or1[3]:.1e}"]
ors = [or2, or1]; ys = [0, 1]
for y, o in zip(ys, ors):
    b.plot([o[1], o[2]], [y, y], "-", color="#2c3e50", lw=2)
    b.plot(o[0], y, "o", color="#c0392b", ms=9)
b.axvline(1, color="gray", ls="--", lw=1)
b.set_yticks(ys); b.set_yticklabels(labels, fontsize=9)
b.set_xscale("log"); b.set_xlabel("carrier odds ratio (95% CI, log scale)")
b.set_ylim(-0.6, 1.6); b.set_title("B. Adjusted OR decays 11K → 44K")
b.text(1.02, -0.45, "OR=1\n(no effect)", color="gray", fontsize=8, va="center")

# --- Panel C: dilution — carrier pool grows, carrier-cases stay flat ---
c = ax[2]
x = np.arange(2)
c.bar(x, [r1["cn"], r2["cn"]], 0.5, color="#bdc3c7", label="all carriers")
c.bar(x, [r1["cc"], r2["cc"]], 0.5, color="#c0392b", label="carriers WITH tinnitus")
for i, r in enumerate([r1, r2]):
    c.text(i, r["cn"], f'{r["cn"]}', ha="center", va="bottom", fontsize=10)
    c.text(i, r["cc"], f'{r["cc"]}', ha="center", va="bottom", fontsize=10, color="white", fontweight="bold")
c.set_xticks(x); c.set_xticklabels(["v1 (~11K)", "v2 (~44K)"])
c.set_ylabel("number of carriers"); c.set_title("C. Carrier pool grows, cases stay ~4")
c.legend(frameon=False, fontsize=9, loc="upper left")

fig.suptitle("ZNF175 pLOF burden × tinnitus — same (Park) pipeline, v1 vs v2", fontsize=13, y=1.02)
fig.tight_layout()
fig.savefig(R07 / "fig_znf175_v1_v2.png", dpi=130, bbox_inches="tight")
print("saved -> results/07/fig_znf175_v1_v2.png")
print(f"v1 OR={or1[0]:.1f} CI[{or1[1]:.1f},{or1[2]:.1f}] p={or1[3]:.1e}")
print(f"v2 OR={or2[0]:.1f} CI[{or2[1]:.1f},{or2[2]:.1f}] p={or2[3]:.1e}")
