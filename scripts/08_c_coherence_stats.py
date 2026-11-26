"""C-band coherence statistics on the common 20 m grid, compared with L band.

Inputs: WORK/s1_slc/processed/t133/coh_*.dim (07_s1_coherence.py), WORK/saloum_full, WORK/classes.npy,
results/tide_prediction_2026.csv, results/class_medians_per_date.csv.
Outputs: results/coh_vs_tide_LC.csv, results/c_coh_separability.csv, results/c_coh_ci.csv,
WORK/c_coh_median_VV.tif (median of the 12-day pairs). Usage: python scripts/08_c_coherence_stats.py
"""
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from nisar_s1 import workflow
from nisar_s1.grid import read_geotiff, to_grid
from nisar_s1.landcover import WORLDCOVER_CLASSES, class_separability
from nisar_s1.raster import Raster, load_raster, save_raster
from nisar_s1.snap import dimap_band
from nisar_s1.stats import auc, block_bootstrap_auc

from common import C_COH_MEDIAN, CLASSES, PASS_TIME, RESULTS, S1_SLC, SALOUM

COH = S1_SLC / "processed" / "t133"
ref = load_raster(SALOUM / "nisar_GCOV_T124D_20260902_HV_gamma0_dB.tif")
cls = np.load(CLASSES); assert cls.shape == ref.data.shape
pred = pd.read_csv(RESULTS / "tide_prediction_2026.csv", parse_dates=["time"]).set_index("time").tide_m


def tide(day, hhmm):
    t = pd.Timestamp(f"{day} {hhmm}", tz="UTC")
    return float(pred.iloc[pred.index.get_indexer([t], method="nearest")[0]])


dims = sorted(p for p in COH.glob("coh_*.dim") if "_part" not in p.name and "_radar" not in p.name)
print(f"{len(dims)} C-band pairs", flush=True)
stack = []; rows = []
for k, d in enumerate(dims, 1):
    a, b = d.stem.split("_")[1:3]; dt = (pd.Timestamp(b) - pd.Timestamp(a)).days
    print(f"  [{k}/{len(dims)}] {a} -> {b} ({dt} d)", flush=True)
    dh = abs(tide(b, PASS_TIME["Sentinel-1"]) - tide(a, PASS_TIME["Sentinel-1"]))
    r = to_grid(read_geotiff(dimap_band(d, "coh*VV")), ref)
    r.data[r.data <= 0] = np.nan   # 0 = no data (sea) in the SNAP output
    if dt == 12:
        stack.append(r.data)
    for c, nm in ((95, "Mangroves"), (60, "Bare/tannes"), (10, "Tree cover")):
        v = r.data[(cls == c) & np.isfinite(r.data)]
        rows.append(dict(band="C", pol="VV", ref=a, sec=b, dt=dt, dh=dh, cls=nm, med=np.median(v), n=v.size))
# L band, same table
fl = pd.read_csv(RESULTS / "class_medians_per_date.csv").drop_duplicates(["variable", "date", "cls"])
for r in fl[fl.variable == "L coh20"].itertuples():
    rows.append(dict(band="L", pol="HH", ref=r.date, sec="", dt=12, dh=r.x, cls=r.cls, med=r.med, n=r.n))
t = pd.DataFrame(rows); t.to_csv(RESULTS / "coh_vs_tide_LC.csv", index=False)
out = []
for (band, pol, dt, c), g in t.groupby(["band", "pol", "dt", "cls"]):
    if len(g) >= 5:
        rho, p = spearmanr(g.dh, g.med)
        out.append(dict(band=band, pol=pol, dt=dt, cls=c, n=len(g), median=g.med.median(), rho=rho, p=p))
print(pd.DataFrame(out).round(3).to_string(index=False))
m = np.nanmedian(np.stack(stack), axis=0).astype("float32")
layers = {f"C coh VV [n={len(stack)}]": Raster(m, ref.transform, ref.crs)}
save_raster(Raster(m, ref.transform, ref.crs), C_COH_MEDIAN)
med = workflow.temporal_medians(SALOUM)
layers.update({k: v for k, v in med.items() if k.startswith("L coh")})
seps = [class_separability(layers, cls, 95, int(c)) for c in np.unique(cls) if c and c != 95 and (cls == c).sum() >= 50]
sep = pd.concat(seps, ignore_index=True); sep["class_b_name"] = sep.class_b.map(WORLDCOVER_CLASSES)
sep.to_csv(RESULTS / "c_coh_separability.csv", index=False)
print(sep.pivot_table(index="class_b_name", columns="variable", values="auc_sep").round(2).to_string())
ci = []
for name, r in layers.items():
    for c in (10, 60):
        u = auc(r.data[(cls == 95) & np.isfinite(r.data)], r.data[(cls == c) & np.isfinite(r.data)])
        lo, hi = block_bootstrap_auc(r.data, cls == 95, cls == c, block=25, n=300)
        ci.append(dict(var=name, cls=c, auc=round(max(u, 1 - u), 3), lo=lo, hi=hi))
print(pd.DataFrame(ci).to_string(index=False)); pd.DataFrame(ci).to_csv(RESULTS / "c_coh_ci.csv", index=False)
