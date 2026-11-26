"""Harmonic tide model fitted on the Dakar tide gauge (UHSLC 223, hourly), predicted at the satellite passes.

h(t) = Z0 + sum_k [a_k cos(w_k t) + b_k sin(w_k t)], 22 constituents, least squares. Fitted on 2015-2020 and
checked on 2021-2022 (not used in the fit), then refitted on 2015-2022.

Inputs: data/tide/dakar_global_hourly_rqds.csv (00_fetch_inputs.py), data/nisar_times.csv (01_satellite_passes.py).
Outputs: results/passes_tide.csv, results/tide_prediction_2026.csv. Usage: python scripts/02_tide_model.py
"""
import numpy as np
import pandas as pd

from nisar_s1.cband import search_s1_rtc

from common import BBOX, DATA, RESULTS

d = pd.read_csv(DATA / "tide" / "dakar_global_hourly_rqds.csv", skiprows=[1])   # row 2 holds the units
d["time"] = pd.to_datetime(d["time"]); d = d.dropna()
d = d[(d.sea_level > -5000)]
med = d.sea_level.median(); mad = (d.sea_level - med).abs().median()
d = d[(d.sea_level - med).abs() < 10 * mad]
print(f"Dakar: {len(d)} valid hourly values, {d.time.min().date()} -> {d.time.max().date()}")
# constituent periods (hours)
P = {"M2": 12.4206012, "S2": 12.0, "N2": 12.65834751, "K2": 11.96723606, "2N2": 12.90537297, "MU2": 12.8717576,
     "NU2": 12.62600509, "L2": 12.19162085, "T2": 12.01644934, "K1": 23.93447213, "O1": 25.81933871, "P1": 24.06588766,
     "Q1": 26.868350, "J1": 23.09848146, "M4": 6.210300601, "MS4": 6.103339275, "MN4": 6.269173724, "M6": 4.140200401,
     "MF": 327.8599387, "MM": 661.3111655, "SSA": 4383.076325, "SA": 8766.15265}
T0 = pd.Timestamp("2020-01-01", tz="UTC")


def design(t):
    h = (pd.DatetimeIndex(t) - T0).total_seconds().values / 3600
    cols = [np.ones_like(h)]
    for p in P.values():
        w = 2 * np.pi / p; cols += [np.cos(w * h), np.sin(w * h)]
    return np.column_stack(cols)


def fit(df):
    c, *_ = np.linalg.lstsq(design(df.time), df.sea_level.values / 1000, rcond=None)
    return c


# validation: fit 2015-2020, predict 2021-2022
tr, te = d[d.time < "2021-01-01"], d[d.time >= "2021-01-01"]
c = fit(tr); res = te.sea_level.values / 1000 - design(te.time) @ c
print(f"Validation 2021-2022 (not used in the fit): RMSE {np.sqrt(np.mean(res**2)) * 100:.1f} cm, "
      f"tidal range std {te.sea_level.std() / 10:.1f} cm, explained {1 - np.var(res) / np.var(te.sea_level / 1000):.3f}")
c = fit(d)
amp = {k: np.hypot(c[1 + 2 * i], c[2 + 2 * i]) for i, k in enumerate(P)}
print("Main amplitudes (cm):", {k: round(v * 100, 1) for k, v in sorted(amp.items(), key=lambda x: -x[1])[:6]})


def pred(t):
    return design(t) @ c - c[0]


# satellite passes
n = pd.read_csv(DATA / "nisar_times.csv", parse_dates=["startTime"])
n["sensor"] = "NISAR L (T124D)"; n = n.rename(columns={"startTime": "time"})[["sensor", "time"]]
s = search_s1_rtc(BBOX["saloum"], "2026-05-01", "2026-11-22")
s["sensor"] = "S1 C (R" + s.relative_orbit.astype(str) + s.orbit_state.str[0].str.upper() + ")"
s = s.rename(columns={"datetime": "time"})[["sensor", "time"]]
a = pd.concat([n, s]); a["time"] = pd.to_datetime(a.time, utc=True)
a["tide_m"] = pred(a.time)
a["dhdt_m_per_h"] = (pred(a.time + pd.Timedelta("10min")) - pred(a.time - pd.Timedelta("10min"))) * 3
a["range_day_m"] = [np.ptp(pred(pd.date_range(t - pd.Timedelta("12h"), t + pd.Timedelta("12h"), freq="10min")))
                    for t in a.time]
a = a.sort_values(["sensor", "time"])
a.to_csv(RESULTS / "passes_tide.csv", index=False)
pd.set_option("display.width", 160)
print(a.assign(time=a.time.dt.strftime("%Y-%m-%d %H:%M")).round(2).to_string(index=False))
for k, g in a.groupby("sensor"):
    print(f"{k}: n={len(g)}, tide min {g.tide_m.min():.2f} max {g.tide_m.max():.2f} m, span {np.ptp(g.tide_m):.2f} m, "
          f"std {g.tide_m.std():.2f} m; typical daily range {g.range_day_m.median():.2f} m")
t = pd.date_range("2026-05-01", "2026-11-20", freq="15min", tz="UTC")
pd.DataFrame({"time": t, "tide_m": pred(t)}).to_csv(RESULTS / "tide_prediction_2026.csv", index=False)
