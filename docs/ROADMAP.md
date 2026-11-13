# ROADMAP.md — PyWES term paper

> **For Claude:** read this file (together with `CLAUDE.md`) at the start of
> project work. It holds the plan, the current status, verified facts and the
> sanity checks. Rules: keep the MVP in focus, never invent numbers or results,
> report every filter with the number of samples it removes.
>
> **Keeping it up to date:** Mirza allows Claude to update this file **without
> asking** whenever it is useful — e.g. a step is finished, a decision is made,
> a new fact is verified, a new sanity check or risk shows up, or a plan
> changes. Rules for updates:
> - Change only what changed; don't rewrite or reorder sections for style.
> - Only verified facts and the team's own decisions — mark anything uncertain.
> - Never delete a decision silently: mark it "replaced by …" with the date.
> - Update the "Last updated" date and add one line to the change log (§11).
> - Tell Mirza in one line what was changed.
>
> Last updated: 2026-11-13.

---

## 1. Goal

**Big question:** How much do atmospheric stability and wind shear change wake
losses at Hill of Towie, and does a weather-aware wake model (one setting per
stability class) predict them better than one global setting?

Sub-goals:
1. Measure wake losses per stability class (and per shear class) in SCADA.
2. Compare FLORIS with one global tuning vs. one tuning per class.
3. (Extra) Simulate wake steering per class.
4. (Optional) Predict stability class / wake loss a few hours ahead.

Course term paper, not a publication: reuse published methods, but cite them.

---

## 2. Verified facts (from the dataset docs in `data/sample/`)

| Item | Value | Source |
|---|---|---|
| Turbines | 21 × Siemens SWT-2.3-VS-82, 2.3 MW, hub 59 m, rotor 82 m (rotor layer ~18–100 m) | `README_HillofTowie.md` |
| Timestamps | **All UTC**. SCADA `TimeStamp` = **end** of 10-min period | data dictionary |
| LiDAR timestamp | Start or end of the 10-min period is **not documented** (~5 min effect) | data dictionary, `analysis/stability_by_sector.py` |
| Directions | Meteorological (from north, clockwise) | data dictionary |
| Turbine directions | **Turbine north** → northing correction needed before any directional analysis; corrections can change in time | data dictionary |
| ZX300 (unit 2428) | Ground profiler ~280 m north of T11; heights 17…99, 181, 263 m above optical head (~1 m agl); 58 m ≈ hub | data dictionary |
| ZX300 period | 10-min 2025-08-21 – 2026-04-30 | README |
| ZXTM (unit 5060) | **Nacelle LiDAR on T07**; ranges 10, 48, 64, 85, 126, 208 m; rotor segment heights 26.2–91.8 m; FD fit gives shear exponent (assumes no veer) + fit residual for filtering | data dictionary |
| ZXTM period | 10-min 2025-12-12 – 2026-04-30 | README |
| LiDAR directions | **Instrument north** (ZX300 oriented to magnetic north per `lidar.py`) | data dictionary, ZXMS-RES029 |
| LiDAR missing values | 9987–9999 = bad data | data dictionary |
| Dynamic Yaw campaign | Jan–Apr 2026: Dynamic Yaw controller toggled on/off in cycles (fastlog) | `data/README.md`, README |
| Dynamic Yaw = what | RES wind farm controller: **wake steering + collective yaw control**, commands sent through the turbines' own controller (no turbine hardware added) | RES presentation (web) |
| Dynamic Yaw design | 13 turbines controlled with toggling, **50 min ON / 50 min OFF**; steering pairs with > 100 h data: **T11→T14, T03→T07, T02→T05**; T07 not controlled (benefiting turbine); reference turbines for collective yaw: T01, T13, T18 | RES presentation (web) |
| Dynamic Yaw results (RES) | Farm AEP +0.7 % (P50, ±0.4 %); wake steering +1.32 % (±0.54 %) in steerable directions only (~18 % of long-term weight); collective yaw +0.56 % (±0.47 %); site modelled wake loss ~7.6 % (low rotors, high TI → hard site for steering) | RES presentation (web) |
| Dynamic Yaw tags | **Only in fastlog** (not in 10-min SCADA): `computed_driver_post_processed_toggle_state` (0 = off, 1 = on), `..._core_wake_steering_offset_degrees` (signed, + = counter-clockwise), consensus wind direction (true north), yaw target (**turbine north**) | data dictionary |
| Campaign northing | Corrections used in the campaign: companion repo `scripts/wfc_analysis_2026/wind_up_config/northing_corrections/HOT_DY_Design_northing.yaml` | data dictionary |
| AeroUp | Retrofit works on all 21 turbines finished by **2023-05-26** (T13 first, 2021-09) → 2024–2026 data is all post-AeroUp | `Hill_of_Towie_AeroUp_install_dates.csv` |
| TuneUp | Control upgrade in **2024** on 9 test turbines (~+1.1 % P50, RES). **Turbines and exact dates not in our local docs** — check companion repo `tune_up.py` before using 2024 data | RES PEAK25 poster (web) |
| `ShutdownDuration` | Covers **only the v1 window (2016 – 2024-08)** → no shutdown file for 2025–2026; use alarm log + status/power filters there | README |
| `wtc_ActualWindDirection` | Present in 10-min SCADA from 2022-04 | CHANGELOG 2.1.0 |
| Dataset DOI | Concept DOI 10.5281/zenodo.14870021 (latest); v2.1.0 = 10.5281/zenodo.22662930 | CHANGELOG |
| Other filters available | `ShutdownDuration.csv` (to 2024-08 only), alarm logs, AeroUp retrofit dates | README |

Data on disk (`data/raw/`, 2026-11-10): SCADA 2025-01 – 2026-04, ERA5 2025-01 –
2025-07 (download script targets up to 2026-04), LiDAR 2428 + 5060, ShutdownDuration.

---

## 3. Status (2026-11-10)

**Done**
- Repo, `src/pywes` package, `pip install -e .`, Spyder `# %%` workflow.
- `scripts/download_era5.py`: monthly ERA5 single-level download (3×3 grid,
  fluxes, stress, u*, t2m, d2m, sp, u/v 10 m and 100 m, blh, sdfor).
- `src/pywes/stability.py`: 1/L (ECMWF method, incl. moisture flux), classes
  (neutral if |1/L| < 0.002 1/m, Tuna et al. 2026), 30° sectors, season,
  day/night, ERA5-hour matching, filters with report.
- `src/pywes/lidar.py`: ZX300 loader (sentinels → NaN), cleaning with report,
  shear exponent fit, veer, IEC-style disturbed sectors (**formula not yet
  verified against the standard**).
- `analysis/stability_by_sector.py`: stability per sector/season/day-night,
  sanity checks, LiDAR part (α, veer, TI per ERA5 class) — LiDAR part currently
  switched off (`RUN_LIDAR_PART = False`).
- `data/processed/era5_stability_hourly.csv` (from the ERA5 months on disk).

**In progress / next**
- ERA5 download for 2025-08 … 2026-04 (needed for the LiDAR overlap).
- Run the LiDAR part of `stability_by_sector.py` once ERA5 covers the LiDAR period.

**Small cleanups**
- `lidar.py` docstring says pyarrow is not in `requirements.txt` — it is now.
- `data/README.md` still has the "Claude generated" warning, placeholder
  citations/DOI and `scripts/<download_era5>.py` placeholder.

---

## 4. Roadmap

### Phase 1 — Data (MVP)
1. **ERA5:** finish download to 2026-04 (ask before starting large requests).
2. **SCADA loading and cleaning** → `data/processed/`
   - Apply northing corrections before using direction/yaw.
   - Remove: shutdowns (`ShutdownDuration` — **only to 2024-08**; for 2025–2026
     use stopping alarms from `tblAlarmLog` + power/status filters), stopping
     alarms, curtailment / derating, missing data, icing-like outliers; mark
     AeroUp periods (none after 2023-05-26).
   - **Mark or exclude Dynamic Yaw "on" periods (Jan–Apr 2026)** for the
     baseline wake analysis — yaw behaviour differs there. The ON/OFF flag is
     only in the fastlog (`toggle_state`, ~12 GB zip). Without the fastlog:
     exclude Jan–Apr 2026 for the 13 controlled turbines (simplest, safest).
   - Late 2025: RES says the trial ran "from late 2025" (web page) — **not
     verified** whether yaw offsets were active before 2026-01-01. Check yaw
     misalignment on the steering turbines (T02, T03, T11) in Oct–Dec 2025 vs
     earlier months before using that period as baseline.
   - Report each filter with removed samples.
3. **LiDAR loading and cleaning:** ZX300 (`lidar.py`), ZXTM loader still to
   write (use FD fit residual / fit flags to filter).

### Phase 2 — Stability classes (MVP)
1. ERA5 1/L → stable / neutral / unstable, per 30° sector (done for part of the period).
2. **ERA5 vs LiDAR stability check** (extra but small):
   - LiDAR cannot give L (no heat flux, no temperature gradient). Compare ERA5
     classes with LiDAR **proxies**: shear α, TI at hub (58 m), veer (28→99 m).
   - Join on nearest ERA5 hour (`nearest_era5_hour`), free-stream sectors only
     (`iec_disturbed_sectors`).
   - Plots: box plots per class, α per sector × class, z/L vs α and TI scatter,
     ERA5 vs LiDAR direction difference, diurnal cycle per month.
   - Keep it to comparison plots; do **not** replace ERA5 with LiDAR for the
     wake analysis (LiDAR covers only ~8 months).

### Phase 3 — Wind shear (MVP-light + extras)
1. **ERA5 shear:** α = ln(U100/U10) / ln(100/10) per hour; drop low-wind hours
   (report count). 10–100 m is close to the rotor layer (~18–100 m) here.
2. **Validate ERA5 α with LiDAR** (overlap period):
   - ZX300: fit α over the same layer (≈ 17–99 m) and over the rotor (28–99 m).
   - ZXTM on T07: `FD Vertical Wind Shear Exponent`, only when T07 is in
     free stream (otherwise it sees wake shear); filter by fit residual.
   - Bias, correlation, split by day/night and stability class.
3. **Shear classes** from data (e.g. terciles of α), not fixed literature thresholds.
4. Wake losses per shear class (Phase 4 method).
5. **Separate shear from stability:** 2D bins (stability × shear) if samples
   allow, or shear effect within one stability class. Always keep wind-speed bins.

### Phase 4 — Observed wake losses (MVP)
- Turbine pairs aligned with the wind; power ratio downstream / upstream
  (removes most of the shear/power-curve effect).
- Bin by wind direction (northing-corrected) and wind speed (upstream turbine).
- Per stability class, per sector; later per shear class.
- Option: OpenOA `WakeLosses` for a farm-level cross-check.
- Minimum samples per bin (e.g. `MIN_N_CELL`) — hide bins below it.

### Phase 5 — FLORIS (MVP)
- Model: Empirical Gaussian (`emgauss`); Gauss/GCH as second option.
- Tune the wake expansion (or TI) **once globally** and **once per class**.
- Output: comparison plots + error table (e.g. RMSE / bias of power ratio per class).
- `fmodel.set(wind_shear=...)` is **one float per model** (not per time step)
  → run once per shear class with the class-median α. `wind_veer` is one
  constant and only used in some models.
- Linking stability to TI / expansion is **our own contribution** (no FLORIS
  model has stability built in).

### Phase 6 — Paper
- Final figures → `results/figures/` (drafts in `draft/`, ask before overwriting).
- Methods with citations, filter tables, limitations.

### Extras (only after the MVP works — ask before starting)
- **Real wake-steering check (Dynamic Yaw, Jan–Apr 2026):** ON vs OFF power
  ratio for T11→T14, T03→T07, T02→T05 per stability class, and compare with
  FLORIS-predicted gains for the same pairs and offsets. Needs the fastlog
  (toggle state + commanded offsets). Expect large uncertainty: only ~100+ h
  per pair in steering directions, winter only (few unstable hours).
- Stability-dependent yaw lookup tables (FLORIS yaw optimization per class).
- Robustness to wind direction uncertainty (`UncertainFlorisModel`).
- Forecasting stability class / wake loss from weather data.
- Load proxies in FLORIS: `get_turbine_thrust_coefficients()` (→ thrust
  T = ½ρACtU²), `get_turbine_TIs()`, `load_optimization` module (`compute_lti`,
  VOC). These are **proxies, not structural loads**.

### Outlook only (write about, don't build)
- Dynamic wake models (OFF / FLORIDyn), helix, real loads (OpenFAST / FAST.Farm).

---

## 5. Sanity checks (run them, print OK / NOT OK)

**Stability (ERA5)**
- Unstable more frequent by day and in summer; stable more frequent at night
  and in winter. If not → check the sign of the ERA5 fluxes (positive downward).
- Our u* (from stress) vs ERA5 `zust`: ratio ≈ 1, high correlation.
- No jumps at month borders; no duplicate hours.
- Report `sdfor` (1/L less reliable over strong orography).

**LiDAR**
- Sentinels removed; Status Flags = "Fully Operational"; packet filter reported.
- Stable → higher α, more veer (clockwise with height), lower TI; unstable → opposite.
- α is unreliable at low wind speed and in unstable conditions.
- Direction difference ERA5 − LiDAR: expect an offset (instrument/magnetic north); don't mix frames.
- Use only free-stream sectors (ZX300 near T11; ZXTM only when T07 is unwaked).

**Shear**
- High shear mostly at night / stable, with low TI.
- Compare α only over the same height layer.

**SCADA / wakes**
- Wake loss ≥ 0 for the downstream turbine in an aligned pair; if not → check
  sign conventions and northing.
- Wake loss maximum at the pair's alignment direction; symmetric-ish around it.
- Wakes stronger and longer in stable conditions, faster recovery in unstable.
- Wake losses decrease above rated wind speed.
- Shear effect should get weaker once stability is held fixed.

**FLORIS / steering**
- Single turbine and zero yaw give the expected power/thrust curve.
- Per-class expansion: stable < neutral < unstable (slower recovery when stable).
- Steering gains small (a few %) and very sensitive to direction error
  (~4° error can turn a gain into a loss) — be suspicious of large gains.
- ~~Steering gains can only be shown with the model, not validated with real data.~~
  Replaced 2026-11-13: real ON/OFF data exists for Jan–Apr 2026 (3 pairs); RES
  measured +1.32 % in steering directions — our FLORIS gains for those pairs
  should be of a similar size (a few %, not tens of %).

---

## 6. Known risks
- ERA5 is coarse (hourly, ~30 km); 1/L less reliable over complex terrain.
- Some class × sector × speed bins may have too few samples.
- SCADA direction/yaw bias (northing) — biggest error source for wake binning.
- Dynamic Yaw periods (Jan–Apr 2026) can contaminate the baseline — and
  collective yaw control also changes yaw on the 13 controlled turbines, not
  only the steering ones.
- No `ShutdownDuration` for 2025–2026 → downtime filtering there depends on our
  own alarm/power filters.
- Possible controller testing in late 2025 (unverified, see Phase 1).
- Only ~16 months of SCADA on disk (dataset goes back to 2016) — download more
  years only if samples are too few (ask first, ~30 GB total).
- Scope creep.

## 7. Open questions for the team
- Shear classes: terciles or fixed thresholds?
- Which turbine pairs (spacing, free-stream upstream turbine)?
- How many SCADA years do we need for enough samples per class?
- Wind-speed source for binning: upstream SCADA, ERA5 100 m, or LiDAR?
- Use the Dynamic Yaw campaign as a real wake-steering check (extra)? Needs the
  ~12 GB fastlog download (ask first).
- Use the companion package `hot-open` (`pip install hot-open`, BSD-3) for
  loading SCADA/LiDAR/fastlog instead of our own loaders? New dependency → ask.
  Its ERA5 helper uses the Open-Meteo archive — **not checked** whether it
  includes the flux variables we need for 1/L, so keep our CDS download.

---

## 8. Decisions log
- Dataset: Hill of Towie open data (CC-BY-4.0).
- Stability: Obukhov length from ERA5 fluxes (ECMWF method); classes and
  filters after Tuna, Bingöl & Sofuoğlu (2026); also per 30° sector.
- Main library: FLORIS v4; PyWake only as optional cross-check.
- FLORIS model order: emgauss → Gauss/GCH → Cumulative Curl → TurbOPark
  (no deflection) → Jensen + Jiménez (baseline only).
- Steady-state FLORIS for 10-min data; dynamic models = outlook.
- IDE: Spyder, `.py` files with `# %%` cells in `analysis/`.
- LiDAR used to **check** ERA5 stability and shear, not to replace ERA5.

## 9. References to build on (cite in the paper)
- Tuna, Bingöl & Sofuoğlu (2026), J. Wind Eng. Ind. Aerodyn. 272, 106396 — sector-specific stability.
- Tuna & Bingöl (2020), Energies 13, 89 — shear exponent / stability background.
- ECMWF CKB, "ERA5: How to calculate Obukhov Length".
- WES preprint 2026-120 — Bayesian SCADA calibration of Jensen/TurbOPark by stability.
- BrightWind poster, WindEurope Tech Workshop 2025 (PO077) — wake losses per ERA5 stability class.
- WindEurope 2024 poster PO210 — condition-dependent wake decay coefficient.
- Howland et al. (2022), WES 7, 345 — wake steering over a diurnal cycle.
- NREL fy22osti/82125 — wind direction std. dev. as predictor of steering performance.
- TU Delft thesis (Lillgrund) — FLORIS calibrated on SCADA.
- WES 11, 2173 (2026) — direction error vs steering gain.
- WES 8, 1693 (2023) — preview wind direction and steering gains.
- Becker et al. (2025), WES 10, 1055 — OFF / FLORIDyn.
- FLORIS docs: https://natlabrockies.github.io/floris/
- Hill of Towie dataset: concept DOI 10.5281/zenodo.14870021; v2.1.0 DOI 10.5281/zenodo.22662930 (still to add in `data/README.md`).
- RES (Clerc, Torner et al.), Dynamic Yaw presentation, WindEurope Technology Workshop 2026 — campaign design and results: https://www.res-group.com/alex-clerc-dynamic-yaw-wind-europe-tech-presentation
- RES (Clerc), PEAK25 poster — AeroUp + TuneUp uplift with wind-up: https://www.res-group.com/wp-content/uploads/2026/05/PEAK25_Poster_Alex_Clerc.pdf

---

## 10. End-of-project checklist
- [ ] Decide whether to add the Spyder setup steps back to `README.md`.
- [ ] Fill attribution, citation and DOI in `data/README.md`; remove the "Claude generated" warning after checking.
- [ ] Verify the IEC disturbed-sector formula in `lidar.py` (or cite it as unverified).
- [ ] Final figures in `results/figures/`, no data or credentials committed.
- [ ] Every filter in the paper has a removed-samples count.

---

## 11. Change log
- 2026-11-10: Roadmap created (plan, verified data facts, status, sanity checks).
- 2026-11-10: Claude may now update this file without asking (rules at the top).
- 2026-11-13: Added Dynamic Yaw design/results/tags, AeroUp end date, TuneUp 2024, ShutdownDuration gap (2025–26), late-2025 check, DOIs; new extra "real wake-steering check"; replaced the "steering only in model" sanity check.
