"""Download the small public inputs of the study into data/ (about 25 MB, not versioned).

- Dakar tide gauge, hourly sea level, UHSLC station 223: research-quality (2015-2022, used for the
  harmonic fit) and fast-delivery series (Caldwell et al., UHSLC v5.5).
- Dakar tide gauge, IOC Sea Level Station Monitoring Facility, raw 1-min data, June-November 2026
  (independent check of the predicted tide at the satellite passes).
- Natural Earth 1:50m countries (map of the study area).

Usage: python scripts/00_fetch_inputs.py
"""
import time
import urllib.request

import pandas as pd

from common import DATA

UA = {"User-Agent": "nisar-s1-saloum-mangroves/1.0"}


def get(url, out, label):
    out.parent.mkdir(parents=True, exist_ok=True)
    if out.exists() and out.stat().st_size > 0:
        print(f"{label}: already present", flush=True)
        return
    t0 = time.time()
    with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=300) as r:
        out.write_bytes(r.read())
    print(f"{label}: {out.stat().st_size / 1e6:.1f} MB in {time.time() - t0:.0f} s", flush=True)


uhslc = "https://uhslc.soest.hawaii.edu/erddap/tabledap/{ds}.csv?time,sea_level&uhslc_id=223&time>=2015-01-01"
jobs = [(uhslc.format(ds=ds), DATA / "tide" / f"dakar_{ds}.csv", f"UHSLC {ds}")
        for ds in ("global_hourly_rqds", "global_hourly_fast")]
starts = pd.date_range("2026-06-01", "2026-11-16", freq="SMS")   # 1st and 15th of each month
for a in starts:
    b = a + (pd.offsets.SemiMonthBegin())
    jobs.append(("https://www.ioc-sealevelmonitoring.org/service.php?query=data&code=dakar"
                 f"&timestart={a:%Y-%m-%d}&timestop={b:%Y-%m-%d}&format=json",
                 DATA / "ioc" / f"dk_{a:%m_%d}.json", f"IOC Dakar {a:%Y-%m-%d}"))
jobs.append(("https://naciscdn.org/naturalearth/50m/cultural/ne_50m_admin_0_countries.zip",
             DATA / "naturalearth" / "ne_50m_admin_0_countries.zip", "Natural Earth countries"))

for i, (url, out, label) in enumerate(jobs, 1):
    get(url, out, f"[{i}/{len(jobs)}] {label}")
