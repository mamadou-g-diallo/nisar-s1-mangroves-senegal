"""Sentinel-1 C-band coherence over the Saloum: orbit 133 A, sub-swath IW2, VV, with ESA SNAP.

  python scripts/07_s1_coherence.py download [i/n]   # reduced SAFE (IW2 only) of the 16 SLC of data/s1_slc_t133.csv
  python scripts/07_s1_coherence.py process          # coherence: 13 pairs at 12 days + 2 pairs at 6 days
  python scripts/07_s1_coherence.py watch            # process each pair as soon as its two SAFE are present

Needs Earthdata credentials (~/.netrc) for the download and SNAP's gpt (set GPT, default ~/esa-snap/bin/gpt).
Outputs: WORK/s1_slc/sentinel1/t133/*.SAFE ; WORK/s1_slc/processed/t133/coh_<ref>_<sec>.dim (EPSG:32628, 20 m).
"""
import os
import sys
import time
from pathlib import Path

import geopandas as gpd
import pandas as pd
from shapely.geometry import box

from nisar_s1.download import earthdata_session, extract_safe_subset
from nisar_s1.snap import s1_coherence_two_step as s1_coherence

from common import BBOX, DATA, S1_SLC

SEL = DATA / "s1_slc_t133.csv"
CRS, SWATH = "EPSG:32628", "iw2"
GPT = os.environ.get("GPT", str(Path.home() / "esa-snap/bin/gpt"))
RAW, OUT = S1_SLC / "sentinel1" / "t133", S1_SLC / "processed" / "t133"


def dates():
    s = pd.read_csv(SEL)
    s["day"] = pd.to_datetime(s.startTime).dt.strftime("%Y-%m-%d")
    return s.sort_values("day").reset_index(drop=True)


def pairs(s):
    d = list(s.day)
    out = []
    for a in d:
        for b in d:
            dt = (pd.Timestamp(b) - pd.Timestamp(a)).days
            if dt in (12, 6):
                out.append((a, b, dt))
    return out


def download(shard="0/1"):
    i, n = map(int, shard.split("/"))
    s = dates().iloc[i::n]; sess = earthdata_session(); t0 = time.time(); RAW.mkdir(parents=True, exist_ok=True)
    for k, r in enumerate(s.itertuples(), 1):
        tag = f"[stream {i}: {k}/{len(s)}] {r.day}"
        if (RAW / f"{r.sceneName}.SAFE").exists():
            print(f"{tag} already present", flush=True); continue
        try:
            extract_safe_subset(r.url, RAW, SWATH, session=sess, label=tag)
        except Exception as e:
            print(f"{tag} FAILED {type(e).__name__}: {e}", flush=True)
        print(f"{tag} {(time.time() - t0) / 60:.1f} min elapsed", flush=True)


def process():
    s = dates(); safe = {r.day: RAW / f"{r.sceneName}.SAFE" for r in s.itertuples()}
    win = gpd.GeoSeries([box(*BBOX["saloum"])], crs=4326)
    x0, y0, x1, y1 = win.to_crs(CRS).total_bounds
    region = gpd.GeoSeries([box(x0 - 1000, y0 - 1000, x1 + 1000, y1 + 1000)], crs=CRS).to_crs(4326).iloc[0]
    OUT.mkdir(parents=True, exist_ok=True); t0 = time.time()
    jobs = sorted(pairs(s), key=lambda j: j[2] != 12)   # 12-day pairs first, 6-day pairs last
    for i, (a, b, dt) in enumerate(jobs, 1):
        tag = f"[coherence {i}/{len(jobs)}] {a} -> {b} ({dt} d)"
        if not (safe[a].exists() and safe[b].exists()):
            print(f"{tag} SAFE missing, skipped", flush=True); continue
        try:
            s1_coherence(safe[a], safe[b], OUT / f"coh_{a}_{b}.dim", SWATH, win.iloc[0].wkt, region.wkt, GPT,
                         pixel_m=20.0, crs=CRS, label=tag, threads=int(os.environ.get("SNAP_THREADS", 1)), pols="VV")
        except Exception as e:
            print(f"{tag} FAILED {e}", flush=True)
        print(f"{tag} {(time.time() - t0) / 60:.1f} min elapsed", flush=True)


def watch():
    """Process each pair as soon as its two SAFE are present; stop when all pairs exist."""
    jobs = pairs(dates())
    while True:
        process()
        done = sum((OUT / f"coh_{a}_{b}.dim").exists() for a, b, _ in jobs)
        print(f"[watch] {done}/{len(jobs)} pairs done", flush=True)
        if done == len(jobs):
            break
        time.sleep(300)


if __name__ == "__main__":
    {"download": download, "process": process, "watch": watch}[sys.argv[1]](*sys.argv[2:])
