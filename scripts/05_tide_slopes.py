"""Slope of the mangrove gamma0 against the tide at the pass, with its 95 % confidence interval.

Ordinary least squares on the raw per-date medians (no seasonal correction); the interval is
slope +/- t(0.975, n-2) x standard error, and p tests a zero slope.

Input: results/class_medians_per_date.csv. Output: results/tide_slopes_ci.csv.
Usage: python scripts/05_tide_slopes.py
"""
import pandas as pd
from scipy.stats import linregress, t

from common import RESULTS

o = pd.read_csv(RESULTS / "class_medians_per_date.csv").drop_duplicates(["variable", "date", "cls"])
rows = []
for v in ("L HV (dB)", "L HH (dB)", "C VH (dB)", "C VV (dB)"):
    g = o[(o.variable == v) & (o.cls == "Mangroves")]
    r = linregress(g.x, g.med)
    half = t.ppf(0.975, len(g) - 2) * r.stderr
    rows.append(dict(var=v, n=len(g), slope=round(r.slope, 2), lo=round(r.slope - half, 2),
                     hi=round(r.slope + half, 2), p=round(r.pvalue, 2)))
d = pd.DataFrame(rows); d.to_csv(RESULTS / "tide_slopes_ci.csv", index=False)
print(d.to_string(index=False))
