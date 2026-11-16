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
> Last updated: 2026-11-16 (evening).

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
| Official northing files | Now in `data/sample/` (copied 2026-11-16 from RES companion repo, BSD-3, licence copy included): `optimized_northing_corrections.yaml` (from `northing.py` + `HoT_northing.yaml`, wind-up `optimize_northing_corrections`, ERA5T 100 m reference; data window **2016-01 – 2024-08**) and `HOT_DY_Design_northing.yaml` (entries from **2019-10-14**, incl. **2025** changes → **covers 2024–2026**). Each entry `[turbine, start_utc, correction_deg]` is valid from its start until the next entry for that turbine. **Apply as true = (turbine-north value + correction) mod 360** (RES `inspect_data.py`). Check: T11 −20.9°, T14 −11.6° → removes our +9.7° T11−T14 offset. 2025 changes: T02 2025-02-09 (+14.7°) and 2025-03-18 (+8.3°); T03 2025-08-07 (−177.5°, ~180° flip); T12 2025-08-27 (+140.6°); T07 2025-08-29 (+39.4°) | `data/sample/*northing*.yaml`, RES `inspect_data.py`, 2026-11-16 |
| AeroUp | Retrofit works on all 21 turbines finished by **2023-05-26** (T13 first, 2021-09) → 2024–2026 data is all post-AeroUp | `Hill_of_Towie_AeroUp_install_dates.csv` |
| TuneUp | **Pitch** TuneUp (not yaw), deployed **2024-03-14** (works 07:00–18:00 UTC; RES uses data from 17:00). Test turbines: **T02, T03, T06, T08, T09, T13, T16, T20, T21** (9). RES analysed it in 3 groups (east: T16/T20/T21 vs refs T17–T19; north: T08/T09/T13 vs T07/T11/T12/T14; south: T02/T03/T06/T08 vs T01/T04/T05/T07/T10). One-time change, ~+1.1 % P50 | RES `HoT_PitchTuneUp2024_*.yaml` (companion repo), PEAK25 poster |
| Yaw campaigns 2024-01 – 2025-12 | **Corrected 2026-11-16: Dynamic Yaw pre-campaign testing from 2025-10-28 13:00 UTC** (RES excludes ALL turbines to 2025-11-17 13:45). Then test/ref roles changed several times until **2026-02-03**; wrong steering commands on 2025-12-02 (T02, T03, T11); SAT test with **deliberately wrong north offset** on T15 + T18 (2025-11-19/20). RES analysis period starts 2026-01-07 13:00. Before 2025-10-28: no yaw campaign found. Lidars: ZX300 2025-08-21, ZXTM on T07 2025-12-12 | `data/sample/HOT_DY_exclusion_periods.yaml`, RES `HOT_dynamic_yaw.yaml`, lidar reports |
| Northing files agree | RES v1 file vs Dynamic Yaw file where both are valid (2019-10-15 – 2024-09-01, daily): median difference ≤ 0.8° for all 21 turbines; max 5–6° only on a few days around change dates (T10, T11, T14) | `pywes.northing.compare_tables`, 2026-11-16 |
| Curtailment periods (RES) | Farm-wide 'lots of curtailment' periods RES excludes in 2024–2025: 2024-01-28 – 02-01, 2024-04-09 – 04-13, 2024-06-28 – 06-29, **2024-08-09 – 09-28**, 2025-05-26 – 05-27, 2025-06-03 – 06-05 | `data/sample/HOT_DY_exclusion_periods.yaml` |
| `ShutdownDuration` | Covers **only the v1 window (2016 – 2024-08)** → no shutdown file for 2025–2026; use alarm log + status/power filters there | README |
| `ShutdownDuration.csv` format | Columns `TimeStamp_StartFormat`, `TurbineName` (T01…T21, not StationId), `ShutdownDuration` (s per 10 min, 0–600). Timestamp = **start** of period (SCADA = end → shift 10 min when joining). First row 2016-01-01 00:00, last row 2024-09-01 00:00 → covers our SCADA only for **2024-01 – 2024-08**; 2024-09 – 2026-04 has no shutdown file | file header/tail, checked 2026-11-13 |
| Downtime proxy | `600 s − wtc_ScReToOp_timeon` (time not ready to operate) equals `ShutdownDuration` in 99.1 % of 734,754 matched rows (Jan–Aug 2024); recall 0.994, precision 0.916 for downtime > 0 → used after 2024-08 | own check, 2026-11-15 |
| Raw file overlap | Monthly SCADA files hold TimeStamp month-start 00:00 **and** next month-start 00:00 → same keys in two files; last 2025 period (TimeStamp 2026-01-01 00:00) may sit in the 2026-01 file | file heads, 2026-11-15 |
| Sep 2025 gap | `tblSCTurbine_2025_09` holds only 50,506 of 90,720 expected (period, turbine) records (55.7 %) | test run, 2026-11-15 |
| Reduced power setpoint on T16/T17 | In the 2024-08 + 2025-09 test months, T16 runs mostly at `wtc_PowerRef_endvalue` 1993 or 2116 kW (only ~23 % at 2300 kW); T16 + T17 cause ~77 % of the `power_setpoint_reduced` removals → with the conservative rule T16 keeps < 20 % of its records. Cause not known; check on the full run | test run, 2026-11-16 |
| `wtc_ActualWindDirection` | Present in 10-min SCADA from 2022-04 | CHANGELOG 2.1.0 |
| Nacelle offsets differ between turbines | 2024–2025, no northing: median nacelle_T11 − nacelle_T14 = **+9.7°** (IQR +7.0…+11.6°, 69,028 10-min records); nacelle − ERA5 wd100 (hourly): T11 +19.8°, T14 +10.7° → northing correction needed before any pair analysis with a tight same-direction tolerance | own check (`turbine_comparison.py`), 2026-11-16 |
| Dataset DOI | Concept DOI 10.5281/zenodo.14870021 (latest); v2.1.0 = 10.5281/zenodo.22662930 | CHANGELOG |
| Other filters available | `ShutdownDuration.csv` (to 2024-08 only), alarm logs, AeroUp retrofit dates | README |

### Years before 2024 — are they good for the wake comparison? (checked 2026-11-16)

Sources: RES companion repo configs (AeroUp, TuneUp, northing, exclusions,
WeDoWind datapack), `Hill_of_Towie_AeroUp_install_dates.csv`, data dictionary.
Nothing pre-2024 is on disk yet.

| Period | Verdict | Why |
|---|---|---|
| 2016-01 – 2019-10-13 | usable with care | Only the v1 northing file covers it, with many large jumps (T01 ~+155° until 2016-06, T12 −143° / −163° in 2016–2020, T16 ~+80–89° 2017-05 – 2023-02, T05 2017–2018). RES used **2016–2019 as training data** in its WeDoWind challenge. Secondary-anemometer correction switched on 2017-08 / 2018-04 (our `ws` = active sensor `wtc_AcWindSp` — check). |
| **2019-10-14 – 2021-09-22** | **best extra period** | All turbines pre-AeroUp (T13 retrofit starts 2021-09-23), covered by both northing files (they agree), no RES farm-wide exclusion. T19 has a +107° correction until 2019-12-24 (corrected, but check). |
| 2021-09-23 – 2023-05-26 | avoid for power deficits | AeroUp rollout: T13 2021-09, the other 20 turbines 2022-07-14 – 2023-05-26 (power curves change turbine by turbine); RES farm-wide curtailment **2022-09-03 – 2023-02-12**. Speed deficits are less affected, but mixed states between A and B bias P_B/P_A. |
| 2023-05-27 – 2024-03-13 | usable | All post-AeroUp, pre-TuneUp; remove RES curtailment days in 2023 (04-09…04-13, 07-01…07-03, 08-18…08-19, 09-24…09-29, 10-10…10-12) and T11's northing change on 2023-12-06. |

Things that change for any pre-2024 use:
- `wtc_ActualWindDirection` only exists from 2022-04 → our cleaning lists `wd`
  as a required signal and would drop everything before. The wake comparison
  only needs `nacelle_pos`, so `wd` must be made optional for those years.
- `ShutdownDuration` (official downtime) covers all of 2016 – 2024-08 → better
  than the proxy there.
- ERA5 must be downloaded for the extra years (CDS request → ask first), and
  the yearly SCADA zips are ~1.4 GB each.
- Do not pool periods across AeroUp/TuneUp for **power** deficits; compare per
  period or use speed deficits.

Bonus found in the RES repo: a public FLORIS v4 model of Hill of Towie
(`scripts/wedowind_wake_model_2026/floris_inputs/hot_public.yaml` +
`swt_2_3_82_generic.yaml`, "wake model parameters are the FLORIS v4 defaults")
— useful starting point for Phase 5 (cite RES, BSD-3).

Data on disk (`data/raw/`, 2026-11-13): SCADA **2024-01 – 2026-04**, ERA5
**2025-01 – 2026-04** + 2024-01 (rest of 2024 not yet downloaded), LiDAR 2428 + 5060,
ShutdownDuration.

---

## 3. Status (2026-11-16)

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
  sanity checks, LiDAR part (α, veer, TI per ERA5 class) — LiDAR part now
  switched on (`RUN_LIDAR_PART = True`).
- ERA5 download 2025-01 … 2026-04 finished (covers the LiDAR period).
- `stability_by_sector.py` run on 2026-11-13: `data/processed/era5_stability_hourly.csv`,
  `stability_counts_by_sector*.csv` and draft figures 00–09 in
  `results/figures/draft/stability/` (incl. 07 ERA5 vs LiDAR direction, 08 LiDAR
  by class, 09 α per sector × class). **Results and sanity-check output not yet
  reviewed.**

- SCADA cleaning 2024–2025 written: `src/pywes/scada_cleaning.py` + `scripts/clean_scada_2024_2025.py` (single-turbine filters, no power-curve outliers, no wake filters; downtime = 0 s; outputs raw-schema CSVs, flags Parquet, CF NetCDF). Tested on 2024-08 (official downtime) and 2025-09 (proxy). **Full 2024–2025 run not done yet.** Before/after plots: `analysis/scada_cleaning_before_after.py` → `results/figures/draft/cleaning/` (tested on the 2 test months).
- Full 2024–2025 cleaning output found on disk (checked 2026-11-16): `scada_clean_core_2024_2025.nc` covers 2024-01-01 – 2026-01-01; report: 1,500,507 of 2,162,931 records kept (idling −21.8 %, reduced setpoint −5.1 %, downtime −3.5 %). T16 keeps only 27,712 records, T17 46,943 (see §2).
- Wake-pair script: `src/pywes/pairs.py` + `analysis/turbine_comparison.py` (upstream A vs downstream B, ERA5 direction window, nacelle same-direction check, ERA5 scaled to 59 m with ERA5 shear, speed + power deficit per stability class, sanity checks). Test run on T11 → T14 (238°, 4.9 D) only; **results not reviewed yet**.

- **Northing (2026-11-16):** `src/pywes/northing.py` reads RES's northing YAML
  files and applies true = (turbine north + correction) mod 360 (valid-from
  compared with the period START, like wind-up). Tested on real data, T11 → T14
  2024-01 – 2025-10-27: nacelle T11 − T14 **+9.8° → +0.5°**; nacelle − ERA5 wd100
  **T11 −1.3°, T14 −1.2°** after correction (before: about +20° / +11°).
- **Cleaning updated (2026-11-16):** step 16 = exclusion periods, one report
  line per group: `res_curtailment` (RES file), `tuneup_deployment`
  (2024-03-14 07:00–18:00, all turbines), `dynamic_yaw` (RES rows + all turbines
  from 2025-10-28 13:00). Northing applied → `wd_true`, `nacelle_pos_true`;
  marker column `post_tuneup`; NetCDF attributes `northing`, `exclusion_groups`.
  Tested on 2025-10 raw data (last kept period starts 2025-10-28 12:50).
  **Full 2024–2025 re-run still to do on Mirza's PC** (the NetCDF on disk is
  the old run; the comparison script warns about this).
- **Comparison updated (2026-11-16):** `analysis/turbine_comparison.py` ends on
  **27 Oct 2025** (`END = "2025-10-28 00:00"` = end of the last period of
  27 Oct, comment in the script), applies northing right after loading,
  warns about old NetCDF, Dynamic Yaw data and TuneUp pairs. Re-run T11 → T14:
  all sanity checks OK except the class order (unstable 0.185 > neutral 0.079,
  stable 0.231, power deficit in the core). Mean ws_A differs a lot per class
  (neutral 11.9 m/s, unstable/stable 7.7–7.8 m/s), so compare the classes in
  the same wind-speed bins (figure 07) before reading this as physics.

**In progress / next**
- Review the stability figures and the OK / NOT OK sanity-check printout.
- ERA5 2024: `era5_stability_hourly.csv` covers 2024-01-01 – 2026-04-30 (checked 2026-11-16).
- Re-run `scripts/clean_scada_2024_2025.py` (new exclusions + northing), then the comparison.
- Decide whether to add a pre-2024 period (see §2 "Years before 2024").

**Small cleanups**
- `lidar.py` docstring says pyarrow is not in `requirements.txt` — it is now.
- `data/README.md` still has the "Claude generated" warning, placeholder
  citations/DOI and `scripts/<download_era5>.py` placeholder.

---

## 4. Roadmap

### Phase 1 — Data (MVP)
1. **ERA5:** ~~finish download to 2026-04~~ done 2026-11-13; 2024 still open if
   2024 SCADA is used (ask before starting large requests).
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
   - **Late 2025 (verified 2026-11-16):** Dynamic Yaw testing started
     **2025-10-28 13:00 UTC** → baseline = **2024-01-01 – 2025-10-28 13:00** for all
     turbines (simplest; RES excludes ALL turbines until 2025-11-17 13:45, then
     turbine-specific periods). Also apply RES's curtailment periods (§2).
     For 2024: mark the pitch TuneUp step on 2024-03-14 for the 9 test turbines
     (power curve changes, yaw not).
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
- Dynamic Yaw testing from 2025-10-28 (verified) must be excluded from the baseline.
- Only ~16 months of SCADA on disk (dataset goes back to 2016) — download more
  years only if samples are too few (ask first, ~30 GB total).
- Scope creep.

## 7. Open questions for the team
- Shear classes: terciles or fixed thresholds?
- Which turbine pairs (spacing, free-stream upstream turbine)?
- How many SCADA years do we need for enough samples per class?
- Wind-speed source for binning: upstream SCADA, ERA5 100 m, or LiDAR?
- Add a pre-2024 period (best: 2019-10-14 – 2021-09-22)? Needs SCADA zips + ERA5 download
  and `wd` made optional in the cleaning (ask first).
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
- 2026-11-16: Added before/after cleaning plots script (status); T16/T17 reduced-setpoint finding (§2).
- 2026-11-16: Added wake-pair script (status), full 2024–2025 cleaning output (status), T11/T14 nacelle offset finding (§2).
- 2026-11-15: SCADA cleaning module + script added (status); verified downtime proxy, raw file overlap and Sep 2025 gap (§2).
- 2026-11-13: Added verified `ShutdownDuration.csv` format (start-of-period timestamps, TurbineName key, ends 2024-08).
- 2026-11-16: Added `pywes.northing`; cleaning step 16 (RES curtailment, TuneUp works, Dynamic Yaw from 2025-10-28 13:00) + northing columns; comparison END = 27 Oct 2025 + northing; pre-2024 assessment; RES public FLORIS model noted; status updated.
- 2026-11-16: Copied RES northing + exclusion files to `data/sample/`; corrected late 2025 (Dynamic Yaw testing from 2025-10-28); TuneUp = pitch, 2024-03-14, 9 turbines listed; northing sign convention; RES curtailment periods.
- 2026-11-16: Added official northing method/files from the companion repo (ERA5 − yaw, `optimized_northing_corrections.yaml`).
- 2026-11-16: TuneUp = one-time 2024 change (IEA Task 44 minutes); no documented yaw campaign 2024-01 – 2025-12; LiDAR commissioning times; yaw-step check also for 2024.
- 2026-11-13: Status updated — ERA5 to 2026-04 done, LiDAR part switched on and run (figures 00–09), SCADA 2024 + ERA5 2024-01 now on disk; next steps updated.
