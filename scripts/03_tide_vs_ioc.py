"""Predicted tide checked against the IOC Dakar observations (raw, sensor ra2), June-November 2026.

Inputs: data/ioc/dk_*.json (00_fetch_inputs.py), results/tide_prediction_2026.csv, results/passes_tide.csv.
Outputs: results/passes_tide_vs_ioc.csv, figures/f3_maree.png. Usage: python scripts/03_tide_vs_ioc.py
"""
import glob
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from common import DATA, FIGURES, RESULTS

L, C, INK = "#1f5fae", "#d4541f", "#1d1d1d"
plt.rcParams.update({"font.size": 23, "xtick.labelsize": 21, "ytick.labelsize": 21, "font.family": "DejaVu Sans"})
rows = []
for f in sorted(glob.glob(str(DATA / "ioc" / "dk_*.json"))):
    rows += json.load(open(f))
d = pd.DataFrame(rows); d = d[d.sensor == "ra2"]; d["t"] = pd.to_datetime(d.stime, utc=True)
g = d.drop_duplicates("t").set_index("t").slevel.sort_index()
g = g[(g - g.rolling("2h", center=True).median()).abs() < 0.25]          # spikes out
pred = pd.read_csv(RESULTS / "tide_prediction_2026.csv", parse_dates=["time"]).set_index("time").tide_m
g15 = g.resample("15min").median().dropna(); m = pred.reindex(g15.index).dropna(); o = g15.reindex(m.index)
off = (o - m).mean(); obs = g - off
r = np.corrcoef(o, m)[0, 1]; gain = np.polyfit(m, o, 1)[0]
ps = pd.read_csv(RESULTS / "passes_tide.csv", parse_dates=["time"])
q = []
for _, p in ps.iterrows():
    w = obs[(obs.index >= p.time - pd.Timedelta("8min")) & (obs.index <= p.time + pd.Timedelta("8min"))]
    q.append(dict(sensor=p.sensor, time=p.time, model_m=p.tide_m, obs_m=w.median() if len(w) >= 3 else np.nan))
q = pd.DataFrame(q); q["diff_cm"] = (q.obs_m - q.model_m) * 100; q.to_csv(RESULTS / "passes_tide_vs_ioc.csv", index=False)
summ = {}
for key, lab in (("NISAR", "NISAR"), ("R133A", "S1")):
    h = q[q.sensor.str.contains(key)].dropna()
    summ[lab] = dict(n=len(h), rmse=np.sqrt(np.mean(h.diff_cm ** 2)), rho=spearmanr(h.model_m, h.obs_m)[0],
                     rng_obs=np.ptp(h.obs_m), rng_mod=np.ptp(h.model_m))
print(f"15-min series: r={r:.3f} gain={gain:.3f}"); print(pd.DataFrame(summ).round(3))

# Figure (labels in French, as on the poster)
pr = pred[(pred.index >= "2026-06-01") & (pred.index <= "2026-11-20")]
fig, ax = plt.subplots(figsize=(15.3, 3.6), constrained_layout=True)
ax.plot(pr.index, pr.values, color="#c4c4c4", lw=0.5)
qq = q[q.time >= pd.Timestamp("2026-06-01", tz="UTC")]
n = qq[qq.sensor.str.startswith("NISAR")]; s1 = qq[qq.sensor.str.contains("R133A")]
ax.scatter(s1.time, s1.model_m, s=150, color=C, edgecolor="white", lw=1.5, zorder=3, label="Sentinel-1 (C), 19 h 17 UTC")
ax.scatter(n.time, n.model_m, s=190, marker="D", color=L, edgecolor="white", lw=1.5, zorder=4, label="NISAR (L), 19 h 34 UTC")
nn = n.dropna(subset=["obs_m"])
ax.scatter(nn.time, nn.obs_m, s=330, marker="D", facecolor="none", edgecolor=INK, lw=1.6, zorder=5,
           label="NISAR, hauteur mesurée (IOC Dakar)")
ax.set_ylabel("Marée (m)")
MOIS = ["janv.", "févr.", "mars", "avr.", "mai", "juin", "juil.", "août", "sept.", "oct.", "nov.", "déc."]
ax.xaxis.set_major_locator(mdates.MonthLocator())
ax.xaxis.set_major_formatter(mticker.FuncFormatter(
    lambda v, _: MOIS[mdates.num2date(v).month - 1] + (" 2026" if mdates.num2date(v).month == 6 else "")))
ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:g}".replace(".", ",").replace("-", "−")))
ax.set_ylim(-1.0, 1.55)
ax.legend(frameon=False, loc="upper left", ncol=2, fontsize=19, columnspacing=1.0, handletextpad=0.2)
ax.grid(color="#e6e6e6"); ax.set_axisbelow(True)
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)
fig.savefig(FIGURES / "f3_maree.png", dpi=200)
