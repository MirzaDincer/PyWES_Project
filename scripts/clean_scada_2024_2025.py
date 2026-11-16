"""Clean the 10-min Hill of Towie SCADA data for 2024 and 2025.

What it does
------------
Runs the single-turbine cleaning pipeline in ``pywes.scada_cleaning`` on all
records whose 10-min period STARTS between 2024-01-01 00:00 and
2025-12-31 23:50 (UTC). Nothing from 2023 or 2026 is kept. Only the first
lines of the 2026-01 files are read, to pick up the last 2025 period
(TimeStamp 2026-01-01 00:00 = period 2025-12-31 23:50 ... 2026-01-01 00:00).

Also removes the RES exclusion periods (curtailment, Dynamic Yaw tests), the
pitch TuneUp works day (2024-03-14) and everything from the start of Dynamic
Yaw testing (2025-10-28 13:00) - see "Exclusion periods" below - and adds
northing-corrected directions (wd_true, nacelle_pos_true, true north).

Not included here (on purpose): power-curve outlier filters and any wake
filter (free-stream checks, direction sectors, turbine pairs).

Inputs (read-only)
------------------
data/raw/2024/, data/raw/2025/, data/raw/2026/ (head only)
    tblSCTurbine, tblSCTurGrid, tblSCTurFlag, tblSCTurTemp, tblAlarmLog
data/raw/Hill_of_Towie_ShutdownDuration/ShutdownDuration.csv  (to 2024-08)
data/sample/Hill_of_Towie_turbine_metadata.csv, ..._alarms_description.csv
data/sample/HOT_DY_Design_northing.yaml      (RES northing corrections)
data/sample/HOT_DY_exclusion_periods.yaml    (RES exclusion periods)

Outputs (data/processed/)
-------------------------
scada_clean_2024_2025/<YYYY>/<table>_<YYYY>_<MM>.csv
    Kept records in the ORIGINAL raw schema (same header, lines copied
    byte-for-byte). Month = month of the period start.
scada_flags_2024_2025.parquet
    Every record with ERA5 keys (era5_time, latitude, longitude), core
    signals, downtime, one flag per filter, removed_by and keep.
scada_clean_core_2024_2025.nc
    Cleaned core signals on a (time, turbine) grid, CF/ERA5-style.
scada_cleaning_report_2024_2025.csv
    The filter table printed at the end.

Run from the repo root (takes several minutes; ~100 large CSV files):
    python scripts/clean_scada_2024_2025.py
"""

from __future__ import annotations

import time

from pywes import northing
from pywes.paths import PROCESSED_DIR, RAW_DIR, SAMPLE_DIR
from pywes.scada_cleaning import (
    DYNAMIC_YAW_FIRST,
    TUNEUP_WORKS,
    CleaningConfig,
    load_res_exclusion_periods,
    print_summary,
    run_pipeline,
    write_clean_tables,
    write_flags_parquet,
    write_netcdf,
)

# ----------------------------------------------------------------------------
# Settings - change thresholds here, not in the module
# ----------------------------------------------------------------------------
# Downtime: keep only periods with downtime <= this many seconds (of 600 s).
# 0 = zero tolerance. Try e.g. 60 later to see how many periods come back.
MAX_DOWNTIME_S = 0.0

# Where the downtime comes from:
#   "official_then_proxy" - ShutdownDuration (to 2024-08), then the proxy
#                           600 s - wtc_ScReToOp_timeon (validated on 2024-01..08)
#   "proxy_only"          - proxy for the whole period (one consistent rule)
#   "official_only"       - drops everything after 2024-08 (no official data)
DOWNTIME_SOURCE = "official_then_proxy"

# Curtailment: False = remove every period with a reduced power setpoint.
CURTAIL_ONLY_WHEN_LIMITING = False

# Optional filters, OFF by default (None):
FLAT_POWER_STD_KW = None     # e.g. 25.0 - removed normal low-wind running when tested
YAW_STD_MAX_DEG = None       # e.g. 5.0  - selection for wake work, not cleaning

START = "2024-01-01 00:00:00"        # first period START kept (UTC)
END = "2026-01-01 00:00:00"          # first period START NOT kept (UTC)

# Exclusion periods (step 16; times = period START, UTC, end excluded).
# Each group is its own line in the cleaning report.
#  1. RES exclusion file (HOT_DY_exclusion_periods.yaml, from RES's
#     scripts/wfc_analysis_2026/wind_up_config/exclusion_periods/HOT_DY.yaml):
#     farm-wide curtailment periods ("res_curtailment") and the Dynamic Yaw
#     test / role-change periods ("dynamic_yaw").
#  2. Pitch TuneUp works on 2024-03-14 07:00-18:00, all turbines (RES
#     HoT_PitchTuneUp2024_*.yaml). Not in the RES file above.
#  3. Dynamic Yaw from 2025-10-28 13:00 to the end of our window, ALL turbines:
#     RES only excludes all turbines until 2025-11-17 13:45 because RES wants
#     the toggle data. For our baseline (no wake control) we must remove it,
#     because the 10-min SCADA has no ON/OFF flag.
EXCLUDE_RES_PERIODS = True
EXCLUDE_TUNEUP_WORKS = True
EXCLUDE_DYNAMIC_YAW = True

EXCLUDE_PERIODS: list[dict] = []
if EXCLUDE_RES_PERIODS:
    EXCLUDE_PERIODS += load_res_exclusion_periods(SAMPLE_DIR / "HOT_DY_exclusion_periods.yaml",
                                                  START, END)
if EXCLUDE_TUNEUP_WORKS:
    EXCLUDE_PERIODS.append({"start": TUNEUP_WORKS[0], "end": TUNEUP_WORKS[1], "turbines": None,
                            "reason": "pitch TuneUp works (RES HoT_PitchTuneUp2024_*.yaml)",
                            "group": "tuneup_deployment"})
if EXCLUDE_DYNAMIC_YAW:
    EXCLUDE_PERIODS.append({"start": DYNAMIC_YAW_FIRST, "end": END, "turbines": None,
                            "reason": "Dynamic Yaw testing + campaign, ON/OFF unknown in 10-min data",
                            "group": "dynamic_yaw"})

# Northing (turbine north -> true north): RES corrections used for the Dynamic
# Yaw campaign. They start 2019-10-14 and include the 2025 changes (T02, T03,
# T07, T12), so they cover the whole window. Applied as
# true = (turbine north + correction) mod 360 -> columns wd_true, nacelle_pos_true.
NORTHING_FILE = northing.DEFAULT_NORTHING_FILE
NORTHING_CORRECTIONS_DEG = northing.as_config_dict(northing.load_northing_table(NORTHING_FILE))

CONFIG = CleaningConfig(
    raw_dir=RAW_DIR,
    sample_dir=SAMPLE_DIR,
    out_dir=PROCESSED_DIR,
    start=START,
    end=END,
    max_downtime_s=MAX_DOWNTIME_S,
    downtime_source=DOWNTIME_SOURCE,
    curtail_only_when_limiting=CURTAIL_ONLY_WHEN_LIMITING,
    flat_power_std_kw=FLAT_POWER_STD_KW,
    yaw_std_max_deg=YAW_STD_MAX_DEG,
    exclude_periods=EXCLUDE_PERIODS,
    northing_corrections_deg=NORTHING_CORRECTIONS_DEG,
    northing_source=NORTHING_FILE.name,
    tag="2024_2025",
)


def main() -> None:
    t0 = time.time()
    CONFIG.out_dir.mkdir(parents=True, exist_ok=True)

    df, report, info = run_pipeline(CONFIG)
    n_initial, n_final = len(df), int(df["keep"].sum())

    report_path = CONFIG.out_dir / f"scada_cleaning_report_{CONFIG.tag}.csv"
    report.to_csv(report_path, index=False)

    if CONFIG.write_flags_parquet:
        print(f"saved {write_flags_parquet(CONFIG, df)}")
    if CONFIG.write_netcdf:
        print(f"saved {write_netcdf(CONFIG, df)}")
    if CONFIG.write_clean_tables:
        print("Writing cleaned raw-schema tables ...")
        written = write_clean_tables(CONFIG, df)
        print(written.to_string(index=False))

    print_summary(report, info, n_initial, n_final)
    print(f"saved {report_path}")
    print(f"done in {(time.time() - t0) / 60:.1f} min")


if __name__ == "__main__":
    main()
