"""Median of every layer (NISAR and Sentinel-1 gamma0, NISAR 12-day coherence) per WorldCover class and date,
with the tide at the pass (gamma0) or the tide difference |dh| between the two dates (coherence).

Inputs: the nisar-s1 run output of the Saloum (WORK/saloum_full), results/tide_prediction_2026.csv.
Outputs: results/class_medians_per_date.csv, WORK/classes.npy. Usage: python scripts/04_class_medians.py
"""
import os

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from nisar_s1.landcover import erode_classes, worldcover_like
from nisar_s1.raster import load_raster

from common import CLASSES, PASS_TIME, RESULTS, SALOUM

idx = pd.read_csv(SALOUM / "files.csv")
idx["path"] = [p if os.path.isabs(p) else str(SALOUM / os.path.basename(p)) for p in idx.path]
ref = load_raster(idx.path.iloc[0])
print("WorldCover on the grid...", flush=True)
cls = erode_classes(worldcover_like(ref), 1)
np.save(CLASSES, cls)
tide = pd.read_csv(RESULTS / "tide_prediction_2026.csv", parse_dates=["time"]).set_index("time").tide_m


def h(date, hhmm):
    t = pd.Timestamp(f"{date} {hhmm}", tz="UTC")
    return float(tide.iloc[tide.index.get_indexer([t], method="nearest")[0]])


names = {10: "Tree cover", 95: "Mangroves", 60: "Bare/tannes", 90: "Herb. wetland", 40: "Cropland", 80: "Water"}
rows = []
for i, r in enumerate(idx.itertuples(), 1):
    print(f"  [{i}/{len(idx)}] {r.variable} {r.date}", flush=True)
    a = load_raster(r.path).data
    if r.product == "GUNW":
        x = abs(h(r.date2, PASS_TIME["NISAR"]) - h(r.date, PASS_TIME["NISAR"]))
    else:
        x = h(r.date, PASS_TIME[r.sensor])
    for k, nm in names.items():
        v = a[(cls == k) & np.isfinite(a)]
        if len(v) > 200:
            rows.append(dict(variable=r.variable, date=r.date, cls=nm, x=x, med=np.median(v), n=len(v)))
o = pd.DataFrame(rows); o.to_csv(RESULTS / "class_medians_per_date.csv", index=False)
out = []
for (var, c), g in o.groupby(["variable", "cls"]):
    if len(g) >= 5:
        rho, p = spearmanr(g.x, g.med); sl = np.polyfit(g.x, g.med, 1)[0]
        out.append(dict(variable=var, cls=c, n=len(g), rho=round(rho, 2), p=round(p, 3), slope_per_m=round(sl, 2),
                        range=round(np.ptp(g.med), 2)))
pd.set_option("display.width", 160)
print(pd.DataFrame(out).sort_values(["variable", "cls"]).to_string(index=False))
