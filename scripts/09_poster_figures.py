"""Figures 1, 2, 4 and 5 of the Saloum poster, and their tables (labels in French, as on the poster).

Inputs: WORK/saloum_full, WORK/classes.npy, WORK/c_coh_median_VV.tif and the tables of results/.
Usage: python scripts/09_poster_figures.py
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from nisar_s1 import workflow

from common import C_COH_MEDIAN, CLASSES, FIGURES as F, RESULTS, SALOUM as SAL

L, L2, C, C2 = "#1f5fae", "#6fa3e0", "#d4541f", "#f0a074"
INK, MUTED, GRID = "#1d1d1d", "#6b6b6b", "#e6e6e6"
# Figures are drawn at their printed size (column width 15.3 in), so 22 pt here is about 22 pt on the A0 sheet
W = 15.3
plt.rcParams.update({"font.size": 23, "axes.titlesize": 23, "axes.labelsize": 23, "xtick.labelsize": 21, "ytick.labelsize": 21, "legend.fontsize": 20, "font.family": "DejaVu Sans"})


MOIS = ["janv.", "févr.", "mars", "avr.", "mai", "juin", "juil.", "août", "sept.", "oct.", "nov.", "déc."]
import matplotlib.ticker as mticker
import matplotlib.dates as mdates
VIRG = mticker.FuncFormatter(lambda v, _: f"{v:g}".replace(".", ",").replace("-", "−"))
MOISF = mticker.FuncFormatter(lambda v, _: MOIS[mdates.num2date(v).month - 1])


def tidy(ax):
    ax.grid(color=GRID, lw=0.9); ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_color("#aaaaaa")
    ax.tick_params(colors="#444444")


# 1. Maps: temporal medians of L HV and C VH, mangrove outline ---------------------------------
print("[1/4] maps", flush=True)
med = workflow.temporal_medians(SAL)
hv = next(v for k, v in med.items() if k.startswith("L HV")); vh = next(v for k, v in med.items() if k.startswith("C VH"))
cls = np.load(CLASSES)
fig, axs = plt.subplots(1, 2, figsize=(14, 6.6), constrained_layout=True)
for ax, r, title in ((axs[0], hv, "NISAR bande L, HV (médiane de 11 dates)"),
                     (axs[1], vh, "Sentinel-1 bande C, VH (médiane de 16 dates)")):
    a = r.data
    lo, hi = np.nanpercentile(a, [2, 98])
    im = ax.imshow(a, cmap="gray", vmin=lo, vmax=hi)
    ax.contour(cls == 95, levels=[0.5], colors=["#2fbf71"], linewidths=0.7)
    ax.set_title(title, loc="left", color=INK); ax.set_xticks([]); ax.set_yticks([])
    cb = fig.colorbar(im, ax=ax, shrink=0.8, pad=0.01); cb.set_label("γ⁰ (dB)")
ny, nx = hv.data.shape
axs[0].plot([nx * 0.06, nx * 0.06 + 5000 / 20], [ny * 0.95] * 2, color="white", lw=4)
axs[0].text(nx * 0.06, ny * 0.93, "5 km", color="white", fontsize=14)
fig.savefig(F / "f1_cartes.png", dpi=200); plt.close(fig)

# 2. Separability of mangroves vs each class -----------------------------------------------------
print("[2/4] separability", flush=True)
names = {"Tree cover": "Forêt continentale", "Shrubland": "Arbustes", "Built-up": "Bâti",
         "Cropland": "Cultures", "Grassland": "Savane herbeuse", "Herbaceous wetland": "Prairie humide",
         "Bare / sparse vegetation": "Tannes, sol nu", "Permanent water bodies": "Eau"}
s = pd.read_csv(SAL / "separability_t95.csv")
cc = RESULTS / "c_coh_separability.csv"
if cc.exists():   # C-band coherence (Sentinel-1 SLC, VV, 12-day pairs), same classes and grid
    c2 = pd.read_csv(cc); c2 = c2[c2.variable.str.startswith("C coh")]
    s = pd.concat([s, c2], ignore_index=True)
s["var"] = s["variable"].str.replace(r" \[n=\d+\]", "", regex=True)
piv = s.pivot_table(index="class_b_name", columns="var", values="auc_sep").reindex(list(names))
vars_ = [("L HV (dB)", L, "o", "L HV"), ("L HH (dB)", L2, "s", "L HH"), ("L coh20", "#7a4fb3", "^", "L cohérence 12 j"),
         ("C VH (dB)", C, "o", "C VH"), ("C VV (dB)", C2, "s", "C VV"),
         ("C coh VV", "#8c2d04", "v", "C cohérence 12 j")]
fig, ax = plt.subplots(figsize=(W, 4.3), constrained_layout=True)
y = np.arange(len(names))[::-1]
for k, (v, col, m, lab) in enumerate(vars_):
    ax.scatter(piv[v], y + (k - 2.5) * 0.13, color=col, marker=m, s=260, edgecolor="white", lw=1.4, label=lab, zorder=3)
ax.set_yticks(y, [names[n] for n in piv.index]); ax.set_xlim(0.5, 1.01); ax.set_ylim(y[-1] - 0.6, y[0] + 0.7)
ax.axvline(0.8, color="#999999", ls="--", lw=1.2); ax.text(0.805, y[0] + 0.35, "AUC 0,8", color=MUTED, fontsize=20)
ax.xaxis.set_major_formatter(VIRG)   # axis described in the poster caption
ax.legend(frameon=False, loc="center left", bbox_to_anchor=(1.0, 0.5), handletextpad=0.3); tidy(ax)
fig.savefig(F / "f2_separabilite.png", dpi=200); plt.close(fig)

# 4. Mangrove stability through time and tide -----------------------------------------------------
print("[3/4] stability", flush=True)
o = pd.read_csv(RESULTS / "class_medians_per_date.csv").drop_duplicates(["variable", "date", "cls"])
fig, axs = plt.subplots(1, 2, figsize=(W, 3.6), constrained_layout=True)
rows = []
for v, col, m, lab in vars_[:2] + vars_[3:5]:
    g = o[(o.variable == v) & (o.cls == "Mangroves")].sort_values("date")
    an = g.med - g.med.mean()
    axs[0].plot(pd.to_datetime(g.date), an, color=col, marker=m, ms=12, lw=2.5, label=lab)
    k_ = v.replace(" (dB)", "")
    ci = pd.read_csv(RESULTS / "tide_slopes_ci.csv").set_index("var").loc[v]
    f = lambda z: f"{z:+.2f}".replace(".", ",").replace("-", "−")
    axs[1].scatter(g.x, an, color=col, marker=m, s=170, edgecolor="white", lw=1.2, zorder=3, label=lab)
    rho, p = spearmanr(g.x, g.med); sl = np.polyfit(g.x, g.med, 1)[0]
    rows.append(dict(var=lab, n=len(g), range_db=np.ptp(g.med), slope=sl, rho=rho, p=p))
axs[0].set_ylabel("Écart (dB)"); axs[0].set_title("a)", loc="left", fontweight="bold")
axs[1].set_xlabel("Marée au passage (m)"); axs[1].set_title("b)", loc="left", fontweight="bold")
axs[0].legend(frameon=False, ncol=4, loc="lower center", handletextpad=0.2, columnspacing=0.8, borderaxespad=0.1)
for a in axs:
    tidy(a); a.set_ylim(-1.3, 0.75)
axs[0].xaxis.set_major_locator(mdates.MonthLocator()); axs[0].xaxis.set_major_formatter(MOISF)
for a_ in axs: a_.yaxis.set_major_formatter(VIRG)
axs[1].xaxis.set_major_formatter(VIRG)
fig.savefig(F / "f4_stabilite.png", dpi=200); plt.close(fig)
pd.DataFrame(rows).round(3).to_csv(RESULTS / "f4_stabilite.csv", index=False)

# 5. Coherence vs tide difference, L and C ------------------------------------------------------
print("[4/4] coherence", flush=True)
lc = pd.read_csv(RESULTS / "coh_vs_tide_LC.csv")
lab = {"Mangroves": ("Mangrove", "#2fbf71", "o"), "Bare/tannes": ("Tannes, vasières", "#a07a4a", "s"),
       "Tree cover": ("Forêt continentale", "#555555", "^")}
# a, b: maps of the temporal median coherence (same grid); c, d: coherence vs tide difference
import rasterio
fig = plt.figure(figsize=(W, 5.1), layout="constrained")
sf_map, sf_sc = fig.subfigures(1, 2, width_ratios=[1, 1.08], wspace=0.02)
lco = next(v for k, v in med.items() if k.startswith("L coh")).data
cco = rasterio.open(C_COH_MEDIAN).read(1).astype(float); cco[cco <= 0] = np.nan
water = cls == 80
cm = plt.get_cmap("viridis").copy(); cm.set_bad("#d9d9d9")
maxs = sf_map.subplots(1, 2)
for a, img, t in ((maxs[0], lco, "a)"), (maxs[1], cco, "b)")):
    im = a.imshow(np.where(water, np.nan, img), cmap=cm, vmin=0, vmax=1, interpolation="nearest")
    a.set_xticks([]); a.set_yticks([]); a.set_title(t, loc="left", fontweight="bold")
    ny, nx = img.shape; a.plot([nx * 0.06, nx * 0.06 + 250], [ny * 0.94] * 2, color="black", lw=3)
    a.text(nx * 0.06, ny * 0.89, "5 km", color="black", fontsize=17)
cb = sf_map.colorbar(im, ax=maxs, orientation="horizontal", shrink=0.85, aspect=25, pad=0.02)
cb.set_label("Cohérence 12 j (médiane)"); cb.ax.xaxis.set_major_formatter(VIRG)
axs = sf_sc.subplots(1, 2, sharey=True)
out = []
for ax, (band, pol, title) in zip(axs, (("L", "HH", "c)"), ("C", "VV", "d)")), strict=True):
    for k, (nm, col, m) in lab.items():
        g = lc[(lc.band == band) & (lc.pol == pol) & (lc.dt == 12) & (lc.cls == k)]
        rho, p = spearmanr(g.dh, g.med)
        ax.scatter(g.dh, g.med, color=col, marker=m, s=150, edgecolor="white", lw=1.0, label=nm)
        out.append(dict(band=band, cls=k, n=len(g), median=g.med.median(), rho=rho, p=p))
    ax.set_title(title, loc="left", fontweight="bold"); ax.set_xlabel("|Δh| (m)")
    ax.set_ylim(0, 0.92); tidy(ax); ax.xaxis.set_major_formatter(VIRG); ax.yaxis.set_major_formatter(VIRG)
axs[0].set_ylabel("Cohérence 12 j (médiane)")
sf_sc.legend(*axs[0].get_legend_handles_labels(), loc="outside lower center", ncol=3, frameon=False,
             handletextpad=0.1, columnspacing=0.6, fontsize=18)
fig.savefig(F / "f5_coherence.png", dpi=200); plt.close(fig)
pd.DataFrame(out).round(3).to_csv(RESULTS / "f5_coherence.csv", index=False)
print(pd.DataFrame(rows).round(2).to_string(index=False)); print(pd.DataFrame(out).round(2).to_string(index=False))
