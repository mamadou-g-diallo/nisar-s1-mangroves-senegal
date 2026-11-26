"""Acquisition times of NISAR GCOV and Sentinel-1 RTC over the Saloum window.

Writes data/nisar_times.csv and data/s1_times.csv. Usage: python scripts/01_satellite_passes.py
"""
import pandas as pd

from nisar_s1.cband import search_s1_rtc
from nisar_s1.workflow import nisar_catalog

from common import BBOX, DATA

bb = BBOX["saloum"]
df = nisar_catalog(bb, "2025-10-01", "2026-11-22", products=("GCOV",))
print(df[["track", "direction", "bw_code", "pol_code", "startTime"]].to_string())
df[["track", "direction", "bw_code", "pol_code", "startTime", "sceneName"]].to_csv(DATA / "nisar_times.csv", index=False)

s = search_s1_rtc(bb, "2026-01-01", "2026-11-22")
s = s.rename(columns={"relative_orbit": "orbit", "orbit_state": "state", "datetime": "time"})
s = s[["orbit", "state", "time"]].sort_values("time")
print(s.groupby(["orbit", "state"]).size())
s.to_csv(DATA / "s1_times.csv", index=False)
