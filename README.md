# Saloum delta mangroves seen by NISAR — L band vs C band, and the tide at the satellite pass

Independent study comparing **NISAR (L band, HH/HV)** and **Sentinel-1 (C band, VV/VH)** over the mangroves of the
Sine-Saloum delta (Senegal), June–November 2026, with the tide height at every satellite pass.

The processing code is the [`nisar_s1`](https://github.com/mamadou-g-diallo/nisar-s1) package (pinned release
v0.5.0). This repository holds the study itself: configuration, analysis scripts, result tables, figures and the
poster. A sister study applies the same package to forests in France:
[nisar-s1-forest-france](https://github.com/mamadou-g-diallo/nisar-s1-forest-france).

**Poster (French):** [poster/poster_saloum_fr.pdf](poster/poster_saloum_fr.pdf)

## Questions

1. Which band and polarisation best separate mangroves from the other land-cover classes, in particular from
   dryland forest, the class most often confused with them?
2. Does the tide at the time of the pass change the mangrove signal (backscatter, 12-day coherence)?

## Main results

Study window 33 × 33 km, 20 m common grid, ESA WorldCover 2021 classes eroded by one pixel. AUC is reported as
max(AUC, 1 − AUC) (0.5 = no separation, 1 = perfect); 95 % intervals by spatial block bootstrap (500 m blocks).

| Mangroves vs dryland forest | AUC [95 % CI] |
|---|---|
| NISAR L HV, γ⁰ | **0.87** [0.86–0.89] |
| NISAR L HH, γ⁰ | 0.78 [0.76–0.79] |
| Sentinel-1 C VV, γ⁰ | 0.62 [0.59–0.64] |
| Sentinel-1 C VH, γ⁰ | 0.54 [0.51–0.57] |
| NISAR L 12-day coherence (HH) | 0.77 [0.75–0.79] |
| Sentinel-1 C 12-day coherence (VV) | 0.63 [0.62–0.65] |

- **L band separates mangroves from dryland forest; C band hardly does.** Mangroves are darker in L HV
  (−14.2 dB against −11.6 dB). Against water, built-up areas and bare flats (*tannes*), both bands work (AUC ≥ 0.9);
  against shrubland, no single variable suffices.
- **Coherence:** median 12-day coherence over mangroves is 0.72 in L band and 0.27 in C band: the canopy
  decorrelates C band, while L band stays coherent.
- **Tide:** the 11 NISAR passes span 1.00 m of predicted tide (0.90 m measured). No significant effect of the tide
  on mangrove backscatter in either band (all 95 % CIs of the slopes contain 0), nor on coherence against the tide
  difference |Δh| between the two dates (all p ≥ 0.07).
- **Tide model check:** the predicted tide matches the independent IOC Dakar measurements at the passes within
  4.5 cm RMS (NISAR, ρ = 0.89) and 5.0 cm (Sentinel-1).
- **Replication:** the L-band advantage against dryland forest holds in the lower Casamance (L HV 0.94, C VH 0.78).

These are results of one dry-to-wet season, one site and land-cover classes from WorldCover, without field
validation. The tide is the Dakar gauge, about 120 km from the delta, not the tide inside the delta.

## Data

| Data | Use | Source |
|---|---|---|
| NISAR L2 GCOV (γ⁰ HH, HV) and GUNW (12-day coherence), track 124 D, 11 dates | L band | NASA/ISRO, ASF DAAC; products generated with ISCE3 |
| Sentinel-1 RTC (γ⁰ VV, VH), orbit 133 A, 16 dates | C band backscatter | ESA Copernicus, via Microsoft Planetary Computer |
| Sentinel-1 IW SLC, orbit 133 A, sub-swath IW2 | C-band coherence (SNAP), 13 pairs at 12 days | ESA Copernicus, ASF |
| ESA WorldCover 10 m 2021 v200 | Land-cover classes | Zanaga et al. (2022), CC BY 4.0 |
| Dakar tide gauge, hourly, 2015–2022 | Harmonic tide model | UHSLC station 223 |
| Dakar tide gauge, raw, June–November 2026 | Independent check of the model | IOC Sea Level Station Monitoring Facility |

Sentinel-1 IW acquires VV+VH over land, so C-band coherence is computed in VV, the co-polarised counterpart of the
NISAR GUNW HH coherence.

## Method

1. **Extraction** with `nisar-s1 run`: NISAR GCOV/GUNW read remotely by window and Sentinel-1 RTC, all resampled on a
   common UTM 20 m grid; per-class statistics with mangroves (WorldCover 95) as the target class.
2. **Tide**: harmonic model with 22 constituents,
   h(t) = Z₀ + Σₖ [aₖ cos(ωₖt) + bₖ sin(ωₖt)], fitted by least squares on 2015–2020 and checked on 2021–2022
   (6.7 cm RMS, 97 % of variance explained), then refitted on 2015–2022 and predicted at every pass.
3. **C-band coherence** of Sentinel-1 SLC pairs with ESA SNAP, in two `gpt` runs (coherence in radar geometry,
   then terrain correction) to fit in 8 GB of memory.
4. **Statistics**: separation AUC with block-bootstrap intervals; per-date class medians against the tide at the
   pass (OLS slope with 95 % CI); coherence against |Δh| (Spearman ρ).

## Reproduce

```bash
pip install -r requirements.txt          # nisar_s1 v0.5.0 and plotting libraries
export SALOUM_WORK=/path/to/work         # heavy products, outside the repository (default ./work)
python scripts/00_fetch_inputs.py        # tide gauges and Natural Earth into data/ (about 25 MB)
bash scripts/run_extraction.sh           # step 1, needs NASA Earthdata credentials in ~/.netrc
python scripts/01_satellite_passes.py    # acquisition times
python scripts/02_tide_model.py          # tide model, tide at every pass
python scripts/03_tide_vs_ioc.py         # check against IOC measurements, figure 3
python scripts/04_class_medians.py       # per-date class medians, WorldCover classes on the grid
python scripts/05_tide_slopes.py         # slopes against the tide with 95 % CI
python scripts/06_auc_ci.py              # AUC intervals (block bootstrap)
python scripts/07_s1_coherence.py download && python scripts/07_s1_coherence.py process   # SNAP, long
python scripts/08_c_coherence_stats.py   # C-band coherence statistics
python scripts/09_poster_figures.py      # figures 2, 4, 5
python scripts/10_study_area.py          # figure 1 (optional GADM boundaries: export SALOUM_GADM=...)
python scripts/11_replication_casamance.py
cd poster && lualatex poster_saloum_fr.tex
```

Every script prints its progress and writes its tables to `results/` and its figures to `figures/`. The tables
committed in `results/` are those behind the poster, so steps 9–11 can be rerun from them once the work folder
exists. Figure labels are in French, as on the poster.

## Layout

```
config/sites_senegal.yaml   study windows (Saloum, Casamance) for nisar-s1
data/                       small inputs: footprints, Sentinel-1 SLC list, acquisition times (+ fetched inputs)
scripts/                    00–11 analysis steps, common.py (paths), run_extraction.sh
results/                    tables (CSV)
figures/                    poster figures (PNG)
poster/                     poster source (LaTeX, tikzposter), PDF and UCAD logo
```

## References

- Lucas, R. M. et al. (2007). The potential of L-band SAR for quantifying mangrove characteristics and change:
  case studies from the tropics. *Aquatic Conservation* 17(3), 245–264.
- Bunting, P. et al. (2022). Global Mangrove Extent Change 1996–2020: Global Mangrove Watch Version 3.0.
  *Remote Sensing* 14(15), 3657.
- Simard, M. et al. (2019). Mangrove canopy height globally related to precipitation, temperature and cyclone
  frequency. *Nature Geoscience* 12, 40–45.
- Lombard, F. & Andrieu, J. (2021). Mapping mangrove zonation changes in Senegal with Landsat imagery using an OBIA
  approach combined with linear spectral unmixing. *Remote Sensing* 13(10), 1961.
- Oliver-Cabrera, E. & Wdowinski, S. (2016). InSAR-based mapping of tidal inundation extent and amplitude in
  Louisiana coastal wetlands. *Remote Sensing* 8(5), 393.
- Kellogg, K. et al. (2020). NASA-ISRO Synthetic Aperture Radar (NISAR) mission. *Proc. IEEE Aerospace Conference*, 1–21.
- Torres, R. et al. (2012). GMES Sentinel-1 mission. *Remote Sensing of Environment* 120, 9–24.
- Davidson, M. et al. (2023). ROSE-L – the Copernicus L-band synthetic aperture radar imaging mission.
  *Proc. IGARSS 2023*. First launch currently foreseen around 2029; the L-band results here are meant to carry
  over to it.
- ISCE3: InSAR Scientific Computing Environment, version 3 (NASA JPL/Caltech), the NISAR processor —
  <https://github.com/isce-framework/isce3>.
- Zanaga, D. et al. (2022). ESA WorldCover 10 m 2021 v200. Zenodo, doi:10.5281/zenodo.7254221.
- Caldwell, P. C., Merrifield, M. A. & Thompson, P. R. (2015). Sea level measured by tide gauges from global
  oceans — the Joint Archive for Sea Level holdings, version 5.5. NOAA NCEI.

## License

Apache-2.0 (see [LICENSE](LICENSE) and [NOTICE](NOTICE)). Data remain under the terms of their providers.
