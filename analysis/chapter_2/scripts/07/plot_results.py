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

def _fit(dd):
    dd = dd.copy(); dd["age2"] = dd["age"] ** 2
    X = sm.add_constant(dd[["carrier", "age", "age2", "sex_male"] + PCS])
    r = sm.Logit(dd["tinnitus"], X).fit(disp=0, maxiter=200)
    b, se = r.params["carrier"], r.bse["carrier"]
    return np.exp(b), np.exp(b - 1.96 * se), np.exp(b + 1.96 * se), r.pvalues["carrier"]

def fit_or(d):  # pooled adjusted logistic
    return _fit(d.dropna(subset=["carrier", "age", "sex_male", "tinnitus"] + PCS))

def fit_meta(d):  # EUR/AFR-stratified + IVW meta (AFR has 0 carrier-cases -> meta = EUR stratum)
    rows = []
    for anc in ["EUR", "AFR"]:
        sub = d[(d["anc"] == anc)].dropna(subset=["carrier", "age", "sex_male", "tinnitus"] + PCS)
        if sub["tinnitus"].sum() == 0 or sub["carrier"].sum() == 0:
            continue
        try:
            dd = sub.copy(); dd["age2"] = dd["age"] ** 2
            X = sm.add_constant(dd[["carrier", "age", "age2", "sex_male"] + PCS])
            r = sm.Logit(dd["tinnitus"], X).fit(disp=0, maxiter=200)
            rows.append((r.params["carrier"], r.bse["carrier"]))
        except Exception:
            pass
    w = [1/se**2 for _, se in rows if np.isfinite(se) and se > 0]
    bs = [b for (b, se) in rows if np.isfinite(se) and se > 0]
    bm = sum(b*wi for b, wi in zip(bs, w)) / sum(w)
    sem = np.sqrt(1/sum(w))
    from scipy.stats import norm
    p = 2*(1-norm.cdf(abs(bm/sem)))
    return np.exp(bm), np.exp(bm-1.96*sem), np.exp(bm+1.96*sem), p

def rates(d):
    c = d[d.carrier == 1]; n = d[d.carrier == 0]
    return dict(cr=c.tinnitus.mean()*100, cc=int(c.tinnitus.sum()), cn=len(c),
                nr=n.tinnitus.mean()*100, nc=int(n.tinnitus.sum()), nn=len(n))

r1, r2 = rates(d1), rates(d2)
or1 = fit_or(d1); or2 = fit_or(d2)
m1 = fit_meta(d1); m2 = fit_meta(d2)

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

# --- Panel B: forest plot — adjusted logistic (robust) + EUR/AFR meta (Park design), per cohort ---
b = ax[1]
# rows top->bottom: v1 adj, v1 meta, v2 adj, v2 meta
entries = [  # (y, OR-tuple, method, label)
    (3, or1, "adj",  f"v1 adjusted   OR={or1[0]:.1f}, p={or1[3]:.0e}"),
    (2, m1,  "meta", f"v1 EUR/AFR meta   OR={m1[0]:.1f}, p={m1[3]:.0e}"),
    (1, or2, "adj",  f"v2 adjusted   OR={or2[0]:.1f}, p={or2[3]:.0e}"),
    (0, m2,  "meta", f"v2 EUR/AFR meta   OR={m2[0]:.1f}, p={m2[3]:.0e}"),
]
sty = {"adj": dict(color="#c0392b", marker="o"), "meta": dict(color="#2c3e50", marker="D")}
for y, o, meth, _ in entries:
    b.plot([o[1], o[2]], [y, y], "-", color=sty[meth]["color"], lw=2, alpha=0.8)
    b.plot(o[0], y, sty[meth]["marker"], color=sty[meth]["color"], ms=9)
b.axvline(1, color="gray", ls="--", lw=1)
b.set_yticks([e[0] for e in entries]); b.set_yticklabels([e[3] for e in entries], fontsize=8)
b.set_xscale("log"); b.set_xlabel("carrier odds ratio (95% CI, log scale)")
b.set_ylim(-0.6, 3.6); b.set_title("B. Carrier OR decays 11K → 44K")
from matplotlib.lines import Line2D
b.legend(handles=[Line2D([0],[0],color="#c0392b",marker="o",ls="-",label="adjusted logistic (robust)"),
                  Line2D([0],[0],color="#2c3e50",marker="D",ls="-",label="EUR/AFR meta (Park; wide CI)")],
         fontsize=7.5, loc="lower right", frameon=False)

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
