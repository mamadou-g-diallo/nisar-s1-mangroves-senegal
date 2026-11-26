"""Study-area figure (labels in French): Senegal, radar footprints, colour composites of the Saloum window.

Inputs: data/naturalearth (00_fetch_inputs.py), data/footprints, WORK/saloum_full; optional GADM boundaries
(SALOUM_GADM). Output: figures/f0_localisation.png. Usage: python scripts/10_study_area.py
"""
import json
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.patches import Rectangle, Patch
from matplotlib.lines import Line2D
import numpy as np
import geopandas as gpd
from shapely.geometry import shape
from nisar_s1 import workflow

from common import DATA as D, FIGURES as F, GADM, SALOUM
L, C, SITE, INK = "#1f5fae", "#d4541f", "#7b2cbf", "#1d1d1d"
plt.rcParams.update({"font.size": 14, "font.family": "DejaVu Sans"})
world = gpd.read_file(f"zip://{D}/naturalearth/ne_50m_admin_0_countries.zip")
if GADM:   # finer coastline and regions
    sen = gpd.read_file(GADM / "gadm41_SEN_0.shp"); reg = gpd.read_file(GADM / "gadm41_SEN_1.shp")
else:
    sen = world[world.ADMIN == "Senegal"]; reg = None
nis = json.load(open(D / "footprints/nisar_footprints.json")); s1 = json.load(open(D / "footprints/s1_footprints.json"))
sites = {"saloum": ([-16.75, 13.75, -16.45, 14.05], "Delta du Saloum")}

fig = plt.figure(figsize=(16.5, 8.8), constrained_layout=True)
gs = fig.add_gridspec(2, 2, width_ratios=[2.55, 1])
ax = fig.add_subplot(gs[:, 0])
sen.plot(ax=ax, color="#f6ecd2", edgecolor="#3a3a3a", lw=1.6, zorder=2)
if reg is not None:
    reg.boundary.plot(ax=ax, color="#b9a77a", lw=0.7, zorder=2)
for k, ls in (("saloum", "-"),):
    g = shape(nis[k]["geom"]); x, y = g.exterior.xy; ax.plot(x, y, color=L, lw=2.2, ls=ls)
    g = shape(s1[k]["geom"]); x, y = g.exterior.xy; ax.plot(x, y, color=C, lw=1.6, ls=":")
for k, (b, lab) in sites.items():
    ax.add_patch(Rectangle((b[0], b[1]), b[2] - b[0], b[3] - b[1], fc=SITE, ec=SITE, lw=1.5, alpha=0.85, zorder=5))
    ax.annotate(lab, (b[2], (b[1] + b[3]) / 2), xytext=(b[2] + 0.25, b[1] - 0.05),
                fontsize=15, fontweight="bold", color=SITE, arrowprops=dict(arrowstyle="-", color=SITE), zorder=6,
                path_effects=[pe.withStroke(linewidth=4, foreground="white")])
for name, lon, lat, dx, dy in (("Dakar", -17.45, 14.69, 0.12, 0.12), ("Kaolack", -16.07, 14.15, 0.1, 0.12),
                              ("Saint-Louis", -16.49, 16.03, 0.1, 0.1), ("Tambacounda", -13.67, 13.77, 0.1, 0.1),
                              ("Ziguinchor", -16.27, 12.56, 0.1, 0.1)):
    ax.plot(lon, lat, "o", color=INK, ms=5, zorder=6); ax.text(lon + dx, lat + dy, name, fontsize=13, color=INK, zorder=6)
ax.plot(-17.42, 14.676, marker="^", color="#2b8a3e", ms=13, mec="white", zorder=7)
ax.text(-14.9, 15.0, "SÉNÉGAL", fontsize=24, color="#7a6a3a", fontweight="bold", alpha=0.8)
ax.text(-15.6, 13.42, "Gambie", fontsize=11.5, color="#777777", style="italic")
ax.set_xlim(-18.2, -11.1); ax.set_ylim(11.45, 16.85); ax.set_aspect("equal")
ax.set_xlabel("Longitude (°)"); ax.set_ylabel("Latitude (°)")
ax.grid(color="#dddddd", lw=0.6)
ax.legend(handles=[Line2D([], [], color=L, lw=2.2, label="NISAR GCOV, trace 124 D"),
                   Line2D([], [], color=C, lw=1.6, ls=":", label="Sentinel-1, orbite 133 A"),
                   Patch(color=SITE, label="Fenêtre d'étude (33 × 33 km)"),
                   Line2D([], [], marker="^", color="#2b8a3e", ls="", ms=11, label="Marégraphe de Dakar")],
          loc="lower right", fontsize=12.5, framealpha=0.95)
ax.set_title("a)", loc="left", fontsize=17, fontweight="bold")
# inset Africa
ins = ax.inset_axes([0.795, 0.715, 0.195, 0.275])
world[world.CONTINENT == "Africa"].plot(ax=ins, color="#dddddd", edgecolor="white", lw=0.3)
sen.plot(ax=ins, color="#c0002a")
ins.set_xlim(-19, 52); ins.set_ylim(-36, 38); ins.set_xticks([]); ins.set_yticks([]); ins.set_xlabel(""); ins.set_ylabel("")
med = workflow.temporal_medians(SALOUM)
def get(v):
    return next(r for n, r in med.items() if n.startswith(v)).data
def stretch(a, land):
    lo, hi = np.nanpercentile(a[land], [2, 98]); return np.clip((a - lo) / (hi - lo), 0, 1)
for pos, (a1, a2), lab in ((gs[0, 1], ("L HH", "L HV"), "b)"),
                           (gs[1, 1], ("C VV", "C VH"), "c)")):
    a = fig.add_subplot(pos)
    x, y = get(a1), get(a2)
    land = np.isfinite(x) & np.isfinite(y) & (y > np.nanpercentile(y, 50))   # stretch on land, not water
    r, g = stretch(x, land), stretch(y, land)
    rgb = np.dstack([r, g, stretch(x - y, land) * r])   # ratio weighted by HH: water stays dark
    rgb[~np.isfinite(x) | ~np.isfinite(y)] = 0
    a.imshow(rgb)
    ny, nx = x.shape
    a.plot([nx * 0.05, nx * 0.05 + 5000 / 20], [ny * 0.93] * 2, color="white", lw=4)
    a.text(nx * 0.05, ny * 0.89, "5 km", color="white", fontsize=12)
    a.set_xticks([]); a.set_yticks([]); a.set_title(lab, loc="left", fontsize=17, fontweight="bold")
    for sp in a.spines.values(): sp.set_color(SITE); sp.set_linewidth(2.5)
fig.savefig(F / "f0_localisation.png", dpi=200)
print("ok")
