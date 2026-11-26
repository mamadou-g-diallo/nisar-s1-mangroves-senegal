"""Replication: AUC of mangroves against dryland forest, Saloum and Casamance (labels in French).

Inputs: WORK/saloum_full/separability_t95.csv, WORK/casamance/separability.csv.
Outputs: results/f6_replication.csv, figures/f6_replication.png. Usage: python scripts/11_replication_casamance.py
"""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import CASAMANCE, FIGURES, RESULTS, SALOUM

L, L2, C, C2 = "#1f5fae", "#6fa3e0", "#d4541f", "#f0a074"
plt.rcParams.update({"font.size": 15, "font.family": "DejaVu Sans"})
sites = [("Saloum", SALOUM / "separability_t95.csv"), ("Casamance", CASAMANCE / "separability.csv")]
vars_ = [("L HV", L), ("L HH", L2), ("C VH", C), ("C VV", C2)]
rows = []
for site, f in sites:
    s = pd.read_csv(f); s["var"] = s.variable.str.replace(r" \(dB\) \[n=\d+\]", "", regex=True)
    s = s[s.class_b_name == "Tree cover"]
    for v, _ in vars_:
        rows.append(dict(site=site, var=v, auc=float(s[s["var"] == v].auc_sep.iloc[0])))
d = pd.DataFrame(rows); d.to_csv(RESULTS / "f6_replication.csv", index=False)
fig, ax = plt.subplots(figsize=(8, 3.9), constrained_layout=True)
w = 0.19
for k, (v, col) in enumerate(vars_):
    vals = [d[(d.site == s) & (d["var"] == v)].auc.iloc[0] for s, _ in sites]
    x = np.arange(len(sites)) + (k - 1.5) * (w + 0.02)
    ax.bar(x, np.array(vals) - 0.5, w, bottom=0.5, color=col, label=v)
    for xi, vi in zip(x, vals, strict=True):
        ax.text(xi, vi + 0.008, f"{vi:.2f}".replace(".", ","), ha="center", fontsize=13)
ax.set_xticks(range(len(sites)), [s for s, _ in sites]); ax.set_ylim(0.5, 1.0)
ax.set_ylabel("AUC mangrove / forêt continentale")
ax.legend(frameon=False, ncol=4, loc="upper center", bbox_to_anchor=(0.5, 1.13), fontsize=13)
ax.grid(axis="y", color="#e6e6e6"); ax.set_axisbelow(True)
for sp in ("top", "right"): ax.spines[sp].set_visible(False)
fig.savefig(FIGURES / "f6_replication.png", dpi=200)
print(d.pivot(index="site", columns="var", values="auc").round(2))
