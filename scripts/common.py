"""Paths and constants shared by the scripts of the Saloum study.

Small inputs and all tables and figures live in the repository. Heavy products (the GeoTIFFs written by
``nisar-s1 run``, the Sentinel-1 SLC and the SNAP coherence products) live in a work folder outside it,
set with the ``SALOUM_WORK`` environment variable (default: ``./work``, ignored by Git).
"""
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"          # small inputs: committed (footprints, SLC list) or fetched by 00_fetch_inputs.py
RESULTS = ROOT / "results"    # tables written by the scripts
FIGURES = ROOT / "figures"    # figures of the poster
WORK = Path(os.environ.get("SALOUM_WORK", ROOT / "work")).expanduser()

SALOUM = WORK / "saloum_full"   # nisar-s1 run output, Saloum (see README, step 1)
CASAMANCE = WORK / "casamance"  # nisar-s1 run output, Casamance (replication site)
S1_SLC = WORK / "s1_slc"        # Sentinel-1 SLC (reduced SAFE) and SNAP coherence products
CLASSES = WORK / "classes.npy"  # WorldCover classes on the Saloum grid, eroded by one pixel
C_COH_MEDIAN = WORK / "c_coh_median_VV.tif"   # median C-band 12-day coherence on the same grid

# Optional: GADM 4.1 Senegal boundaries (gadm41_SEN_0.shp, gadm41_SEN_1.shp), not redistributable.
# Without it, the study-area map uses the Natural Earth country outline and no regions.
GADM = Path(os.environ["SALOUM_GADM"]).expanduser() if os.environ.get("SALOUM_GADM") else None

BBOX = {"saloum": [-16.75, 13.75, -16.45, 14.05], "casamance": [-16.40, 12.50, -16.20, 12.70]}
PASS_TIME = {"NISAR": "19:34", "Sentinel-1": "19:17"}   # UTC, NISAR track 124 D and Sentinel-1 orbit 133 A

for d in (RESULTS, FIGURES):
    d.mkdir(exist_ok=True)
