"""95 % block-bootstrap intervals of the mangrove separation AUC (Saloum and Casamance).

AUC is reported as max(AUC, 1 - AUC); blocks of 200 m and 500 m (10 and 25 pixels of 20 m), 300 draws.
Inputs: WORK/saloum_full, WORK/casamance. Output: results/auc_ci.csv. Usage: python scripts/06_auc_ci.py
"""
import numpy as np
import pandas as pd

from nisar_s1 import workflow
from nisar_s1.landcover import erode_classes, worldcover_like
from nisar_s1.stats import auc, block_bootstrap_auc

from common import CASAMANCE, RESULTS, SALOUM

SITES = {"saloum_full": SALOUM, "casamance": CASAMANCE}
VARS = ("L HV", "L HH", "C VH", "C VV")
n_jobs = len(SITES) * len(VARS) * 2 * 2
rows, k = [], 0
for site, folder in SITES.items():
    med = workflow.temporal_medians(folder)
    cls = erode_classes(worldcover_like(next(iter(med.values()))), 1)
    for v in VARS:
        img = next(r for name, r in med.items() if name.startswith(v)).data
        for c in (10, 20):            # tree cover, shrubland
            for b in (10, 25):
                k += 1
                print(f"[{k}/{n_jobs}] {site} {v} mangrove/{c} block {b * 20} m", flush=True)
                u = auc(img[(cls == 95) & np.isfinite(img)], img[(cls == c) & np.isfinite(img)])
                lo, hi = block_bootstrap_auc(img, cls == 95, cls == c, block=b, n=300)
                rows.append(dict(site=site, var=v, cls=c, block_m=b * 20, auc=round(max(u, 1 - u), 3), lo=lo, hi=hi))
d = pd.DataFrame(rows); d.to_csv(RESULTS / "auc_ci.csv", index=False)
print(d.to_string(index=False))
