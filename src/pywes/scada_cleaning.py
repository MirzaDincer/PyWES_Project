"""Clean 10-minute Hill of Towie SCADA data (single-turbine filters only).

What this module does
---------------------
It runs a fixed, sequential cleaning pipeline on the 10-min SCADA tables and
returns one boolean flag per filter for every 10-min record. Nothing is
deleted from ``data/raw``. Each step only counts the rows it removes *after*
the previous steps, so the step counts add up to the total removed.

Pipeline (in this order)
    1  timestamp alignment   unparseable or off the 10-min grid
    2  temporal window       period START outside [start, end)
    3  duplicates            same (TimeStamp, StationId) twice -> keep first
    4  unknown turbine       StationId not in the turbine metadata
    5  missing signals       NaN in a required signal (or no row in a table)
    6  physical range        plausibility bounds, wind speed at logger floor
    7  frozen sensors        n consecutive values with ~zero std (FLASC)
    8  downtime              downtime [s] in the period > threshold
    9  turbine error         error timer > 0 s in the period
    10 stopping alarms       stopping alarm active in the period
    11 icing alarms          ice-detection alarm active in the period
    12 cable unwinding       yaw unwind alarm or unwind counter > 0
    13 idling / high pitch   mean power < idle limit, mean pitch > limit
    14 curtailment/derating  setpoint < rated, low power at high wind,
                             (optional) flat power
    15 optional filters      yaw-position std
    16 exclusion periods     one step per group, e.g. RES curtailment periods,
                             pitch TuneUp deployment, Dynamic Yaw testing

Exclusion periods and northing (since 2026-11-16)
-------------------------------------------------
* ``load_res_exclusion_periods`` reads RES's ``HOT_DY_exclusion_periods.yaml``
  (copied to ``data/sample/``; times = period START, like wind-up). Each
  period gets a group from its comment ("curtailment" -> res_curtailment,
  Dynamic Yaw tests -> dynamic_yaw), and every group is its own report step.
* Northing: ``pywes.northing`` (RES ``HOT_DY_Design_northing.yaml``);
  ``wd_true`` / ``nacelle_pos_true`` = (turbine north + correction) mod 360.

NOT included on purpose (see ROADMAP):
    * power-curve outlier filters (bin / MAD / IQR / clustering)
    * any wake filter: free-stream checks, direction sectors, turbine pairs

Downtime source
---------------
``ShutdownDuration.csv`` (RES, computed from SCADA and operator logs) only
covers the dataset v1 window, which ends 2024-08-31. For later periods we use
the proxy ``600 s - wtc_ScReToOp_timeon`` (time NOT "ready to operate").
Checked on Jan-Aug 2024 (734,754 matched rows): the proxy gives exactly the
same value as ShutdownDuration in 99.1 % of rows; recall of ShutdownDuration
> 0 is 0.994 and precision 0.916, so the proxy removes slightly MORE rows.

Time conventions (same as ERA5)
-------------------------------
All times are UTC and stored as timezone-naive datetimes, like the ``time`` /
``valid_time`` axis of ERA5 netCDF files. SCADA ``TimeStamp`` is the END of
the 10-min period. ERA5 hourly instantaneous fields are valid AT the hour, so
we match each 10-min record to the nearest hour of its period MIDDLE
(``pywes.stability.nearest_era5_hour``). Turbine positions are given as
``latitude`` / ``longitude`` in degrees (WGS84), the ERA5 coordinate names.

Sources for the filter rules
----------------------------
* Hill of Towie data dictionary + README (Zenodo 10.5281/zenodo.14870021):
  UTC end-labelled timestamps, 0.01 m/s logger floor, turbine-north frame.
* FLASC ``filter_by_sensor_stuck_faults`` defaults (n = 3, std < 0.001) and
  basic checks (NaN, power < 0, wind speed > 50 m/s):
  https://nrel.github.io/flasc/_autosummary/flasc.data_processing.filtering.html
* OpenOA WakeLosses derating rule: power < 1 % of rated above 4.5 m/s:
  https://openoa.readthedocs.io/en/latest/_modules/openoa/analysis/wake_losses.html
* Idling (mean power < 25 kW), high pitch (> 30 deg) and flat-power
  curtailment (mean < 2275 kW and std < 25 kW) from the FLASC-based baseline
  filters in arXiv:2607.13544 (Table 2; thresholds there are turbine-specific).
* Alarm codes and the stopping flag: ``Hill_of_Towie_alarms_description.csv``.
"""

from __future__ import annotations

import io
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from pywes import northing
from pywes.stability import nearest_era5_hour

# ----------------------------------------------------------------------------
# Constants of the Hill of Towie dataset
# ----------------------------------------------------------------------------
PERIOD = pd.Timedelta("10min")
PERIOD_S = 600.0
KEY = ["TimeStamp", "StationId"]
TIMESTAMP_FORMAT = "%Y-%m-%d %H:%M:%S"     # format of TimeStamp in the raw CSVs

# Columns we need for the filters, per raw table. Short names on the right
# are used inside the pipeline and in the compact (ERA5-ready) outputs.
CORE_COLUMNS: dict[str, dict[str, str]] = {
    "tblSCTurbine": {
        "wtc_AcWindSp_mean": "ws",                     # active sensor [m/s]
        "wtc_ActualWindDirection_mean": "wd",          # turbine north [deg]
        "wtc_NacelPos_mean": "nacelle_pos",            # turbine north [deg]
        "wtc_NacelPos_stddev": "nacelle_pos_std",      # [deg]
        "wtc_PitcPosA_mean": "pitch_a",                # [deg]
        "wtc_PitcPosB_mean": "pitch_b",
        "wtc_PitcPosC_mean": "pitch_c",
        "wtc_GenRpm_mean": "gen_rpm",                  # [rpm]
        "wtc_PowerRef_endvalue": "power_ref",          # setpoint [kW]
    },
    "tblSCTurGrid": {
        "wtc_ActPower_mean": "power",                  # [kW]
        "wtc_ActPower_min": "power_min",
        "wtc_ActPower_max": "power_max",
        "wtc_ActPower_stddev": "power_std",
    },
    "tblSCTurFlag": {
        "wtc_ScReToOp_timeon": "t_ready_s",            # ready to operate [s]
        "wtc_ScInOper_timeon": "t_running_s",          # running [s]
        "wtc_ScTurSto_timeon": "t_error_s",            # turbine error [s]
        "wtc_ScYawUnw_counts": "yaw_unwind_counts",    # [-]
    },
    "tblSCTurTemp": {
        "wtc_AmbieTmp_mean": "t_amb",                  # [deg C]
    },
}
BASE_TABLE = "tblSCTurbine"     # main table (README); defines the row universe

# Pitch TuneUp 2024 (RES companion repo, HoT_PitchTuneUp2024_{east,north,south}.yaml):
# works on 2024-03-14 07:00-18:00 UTC on all turbines; RES uses data from 17:00
# ("upgrade_first_dt_utc_start"). It is a pitch change, so the power curve of
# these 9 turbines changes on that day (yaw is not changed).
TUNEUP_TURBINES = ("T02", "T03", "T06", "T08", "T09", "T13", "T16", "T20", "T21")
TUNEUP_WORKS = ("2024-03-14 07:00:00", "2024-03-14 18:00:00")    # period START, UTC
TUNEUP_FIRST_AFTER = "2024-03-14 17:00:00"

# Dynamic Yaw: pre-campaign testing starts here (RES HOT_DY exclusion file).
# From then on yaw offsets / collective yaw can be active and the 10-min SCADA
# has no ON/OFF flag (only the fastlog has it).
DYNAMIC_YAW_FIRST = "2025-10-28 13:00:00"                         # period START, UTC


# ----------------------------------------------------------------------------
# Settings
# ----------------------------------------------------------------------------
@dataclass
class CleaningConfig:
    """All thresholds of the pipeline. Units are given in the comments.

    Bounds marked "plausibility" are our own generous choices to catch
    impossible values; they are NOT taken from a turbine datasheet.
    """

    raw_dir: Path
    sample_dir: Path
    out_dir: Path

    # --- Step 2: period START must be in [start, end) (UTC) ---
    start: str = "2024-01-01 00:00:00"
    end: str = "2026-01-01 00:00:00"

    # --- turbine ---
    rated_power_kw: float = 2300.0                     # SWT-2.3-VS-82 (metadata)

    # --- Step 5: signals that must not be NaN ---
    required: tuple[str, ...] = ("ws", "wd", "nacelle_pos", "power",
                                 "pitch_a", "pitch_b", "pitch_c", "power_ref")

    # --- Step 6: plausibility bounds (inclusive) ---
    ranges: dict[str, tuple[float, float]] = field(default_factory=lambda: {
        "ws": (0.0, 40.0),             # [m/s]
        "wd": (0.0, 360.0),            # [deg]
        "nacelle_pos": (0.0, 360.0),   # [deg]
        "power": (-50.0, 2530.0),      # [kW]  (-50 kW ... 1.1 x rated)
        "pitch_a": (-5.0, 95.0),       # [deg]
        "pitch_b": (-5.0, 95.0),
        "pitch_c": (-5.0, 95.0),
        "gen_rpm": (0.0, 2000.0),      # [rpm]
        "power_ref": (0.0, 2400.0),    # [kW]
        "t_amb": (-30.0, 40.0),        # [deg C]
    })
    ws_floor: float = 0.01            # [m/s] logger floor (data dictionary)

    # --- Step 7: frozen ("stuck") sensors, FLASC defaults ---
    frozen_columns: tuple[str, ...] = ("ws", "wd")
    frozen_n: int = 3                  # consecutive 10-min values
    frozen_std: float = 0.001          # std below this = frozen

    # --- Step 8: downtime ---
    max_downtime_s: float = 0.0        # keep only downtime <= this [s]
    # "official_then_proxy": ShutdownDuration where it exists, else the proxy
    # "official_only": periods without ShutdownDuration are removed
    # "proxy_only": proxy everywhere (consistent over the whole period)
    downtime_source: str = "official_then_proxy"

    # --- Step 9: turbine error timer ---
    max_error_s: float = 0.0           # [s]

    # --- Steps 10-12: alarm codes ---
    use_lookup_stopping_codes: bool = True       # Stopping == 1 in lookup
    extra_stopping_codes: tuple[int, ...] = ()
    ice_codes: tuple[int, ...] = (102, 8210, 8230, 8234, 8235, 8236)
    unwind_codes: tuple[int, ...] = (10105,)

    # --- Step 13: idling and high pitch (arXiv:2607.13544, Table 2) ---
    idle_power_kw: float = 25.0        # mean power below this = idling
    max_pitch_deg: float = 30.0        # mean pitch above this = not normal

    # --- Step 14: curtailment / derating ---
    power_ref_tol_kw: float = 0.5      # setpoint < rated - tol = derated
    # False: remove every period with a reduced setpoint (conservative).
    # True:  remove it only when the turbine actually hits the setpoint.
    curtail_only_when_limiting: bool = False
    limiting_fraction: float = 0.95    # "hits" = power >= 0.95 x setpoint
    low_power_ws: float = 4.5          # [m/s]  OpenOA derating rule
    low_power_fraction: float = 0.01   # power < 1 % of rated above low_power_ws
    # Flat power (arXiv:2607.13544): mean < 2275 kW and std < this. OFF by
    # default: on Sep 2025 it removed normal low-wind running (~100 kW).
    flat_power_std_kw: float | None = None
    flat_power_max_kw: float = 2275.0

    # --- Step 15: optional filters ---
    yaw_std_max_deg: float | None = None     # e.g. 5.0 (Mittelmeier et al. 2017)
    # --- Step 16: exclusion periods ---
    # dicts with keys start, end (period START, UTC, end excluded),
    # turbines (None = all), reason, group (report step name; default
    # "manual_exclusion"). One report step per group, in order of appearance.
    exclude_periods: list[dict] = field(default_factory=list)

    # --- northing (turbine north -> true north), applied to compact outputs ---
    # {"T01": [("2024-01-01", 3.2), ...]} = (valid_from UTC, offset [deg]),
    # e.g. northing.as_config_dict(northing.load_northing_table()).
    # Empty = no correction -> directions stay turbine north.
    northing_corrections_deg: dict[str, list[tuple[str, float]]] = field(default_factory=dict)
    northing_source: str = ""          # file name, written to the NetCDF attributes

    # --- outputs ---
    output_tables: tuple[str, ...] = ("tblSCTurbine", "tblSCTurGrid",
                                      "tblSCTurFlag", "tblSCTurTemp")
    write_clean_tables: bool = True
    write_flags_parquet: bool = True
    write_netcdf: bool = True
    tag: str = "2024_2025"


# ----------------------------------------------------------------------------
# Small helpers
# ----------------------------------------------------------------------------
def month_starts(start: pd.Timestamp, end: pd.Timestamp) -> list[pd.Timestamp]:
    """First day of every month that holds a period START in [start, end)."""
    last = end - PERIOD
    return list(pd.date_range(start.to_period("M").to_timestamp(),
                              last.to_period("M").to_timestamp(), freq="MS"))


def raw_file(raw_dir: Path, table: str, month: pd.Timestamp) -> Path:
    """Path of one monthly raw CSV, e.g. raw/2025/tblSCTurbine_2025_09.csv."""
    return Path(raw_dir) / f"{month.year}" / f"{table}_{month.year}_{month.month:02d}.csv"


def files_for_window(cfg: CleaningConfig, table: str) -> tuple[list[Path], Path | None]:
    """Monthly files that hold the window, plus the 'boundary' file.

    Raw files are labelled by TimeStamp (= period END) and also contain the
    first timestamp of the next month. The very last period of the window
    (start 2025-12-31 23:50) has TimeStamp 2026-01-01 00:00, which may sit
    only in the 2026-01 file. We read just the head of that file (rows with
    TimeStamp <= end); no other 2026 data is read.
    """
    start, end = pd.Timestamp(cfg.start), pd.Timestamp(cfg.end)
    files = [raw_file(cfg.raw_dir, table, m) for m in month_starts(start, end)]
    boundary = raw_file(cfg.raw_dir, table, end.to_period("M").to_timestamp())
    if boundary in files:
        boundary = None
    return files, boundary


def read_head_text(path: Path, max_timestamp: str, patience: int = 2000) -> str:
    """Header + all lines with TimeStamp <= max_timestamp from the top of a file.

    Assumes the file is sorted by time (true for the raw tables). Stops after
    ``patience`` lines in a row are later than ``max_timestamp``.
    """
    out, late = [], 0
    with open(path, newline="", encoding="utf-8", errors="surrogateescape") as fh:
        out.append(fh.readline())
        for line in fh:
            if line[:19] <= max_timestamp:
                out.append(line)
                late = 0
            else:
                late += 1
                if late >= patience:
                    break
    return "".join(out)


def load_core_table(cfg: CleaningConfig, table: str) -> tuple[pd.DataFrame, list[str]]:
    """Read TimeStamp, StationId and the core columns of one table.

    Returns (DataFrame with short column names + 'source_order', missing files).
    'source_order' keeps the read order, so "keep first" for duplicates is the
    same here and in the writer of the cleaned raw-schema files.
    """
    columns = CORE_COLUMNS[table]
    files, boundary = files_for_window(cfg, table)
    parts, missing = [], []

    def _read(src) -> pd.DataFrame:
        head = pd.read_csv(src, nrows=0)
        use = KEY + [c for c in columns if c in head.columns]
        if isinstance(src, io.StringIO):
            src.seek(0)
        df = pd.read_csv(src, usecols=use, dtype={c: "float32" for c in use[2:]})
        for c in columns:                     # schema change -> column of NaN
            if c not in df.columns:
                df[c] = np.float32(np.nan)
        return df

    for f in files:
        if not f.exists():
            missing.append(str(f))
            continue
        parts.append(_read(f))
    if boundary is not None and boundary.exists():
        text = read_head_text(boundary, cfg.end)
        parts.append(_read(io.StringIO(text)))

    if not parts:
        raise FileNotFoundError(f"No raw files found for {table} in {cfg.raw_dir}")
    df = pd.concat(parts, ignore_index=True)
    df["source_order"] = np.arange(len(df))
    df = df.rename(columns=columns)
    return df, missing


def load_turbine_metadata(sample_dir: Path) -> pd.DataFrame:
    """Turbine name, StationId, latitude, longitude [deg] from the metadata CSV."""
    meta = pd.read_csv(Path(sample_dir) / "Hill_of_Towie_turbine_metadata.csv",
                       encoding="utf-8-sig")
    meta = meta.rename(columns={"Turbine Name": "TurbineName", "Station ID": "StationId",
                                "Latitude": "latitude", "Longitude": "longitude"})
    return meta[["TurbineName", "StationId", "latitude", "longitude"]]


def load_shutdown_duration(raw_dir: Path, start: str, end: str) -> pd.DataFrame:
    """Official downtime [s] per 10-min period, keyed by period START.

    ``ShutdownDuration.csv`` uses ``TimeStamp_StartFormat`` (period START,
    "+00:00" suffix) and ``TurbineName``; it ends with August 2024.
    """
    path = Path(raw_dir) / "Hill_of_Towie_ShutdownDuration" / "ShutdownDuration.csv"
    if not path.exists():
        print(f"WARNING: {path} not found -> proxy downtime only")
        return pd.DataFrame(columns=["time_start", "TurbineName", "downtime_official_s"])
    parts = []
    reader = pd.read_csv(path, chunksize=2_000_000,
                         dtype={"TimeStamp_StartFormat": str, "TurbineName": str,
                                "ShutdownDuration": "float32"})
    for chunk in reader:
        # ISO strings compare correctly as text -> cheap filter before parsing
        ts = chunk["TimeStamp_StartFormat"].str[:19]
        parts.append(chunk[(ts >= start[:19]) & (ts < end[:19])])
    sd = pd.concat(parts, ignore_index=True)
    sd["time_start"] = pd.to_datetime(sd["TimeStamp_StartFormat"].str[:19],
                                      format=TIMESTAMP_FORMAT)
    sd = sd.rename(columns={"ShutdownDuration": "downtime_official_s"})
    sd = sd.drop_duplicates(["time_start", "TurbineName"], keep="first")
    return sd[["time_start", "TurbineName", "downtime_official_s"]]


def load_alarms(cfg: CleaningConfig) -> pd.DataFrame:
    """All alarm events (TimeOn, TimeOff, StationNr, Alarmcode) for the window.

    Alarms that started before the first loaded month are not seen.
    """
    files, boundary = files_for_window(cfg, "tblAlarmLog")
    parts = [pd.read_csv(f) for f in files + ([boundary] if boundary else []) if f.exists()]
    if not parts:
        return pd.DataFrame(columns=["TimeOn", "TimeOff", "StationNr", "Alarmcode"])
    al = pd.concat(parts, ignore_index=True).drop_duplicates()
    al["TimeOn"] = pd.to_datetime(al["TimeOn"], errors="coerce")
    al["TimeOff"] = pd.to_datetime(al["TimeOff"], errors="coerce")
    return al.dropna(subset=["TimeOn", "StationNr", "Alarmcode"])


def alarm_periods(alarms: pd.DataFrame, codes: set[int]) -> pd.DataFrame:
    """Expand alarms with the given codes to the 10-min periods they touch.

    A period [start, start + 10 min) is flagged if the alarm is active at any
    moment inside it. An alarm without TimeOff is treated as a single event
    at TimeOn (only its own period is flagged).
    Returns unique (time_start, StationId) pairs.
    """
    al = alarms[alarms["Alarmcode"].astype(int).isin(codes)]
    if al.empty:
        return pd.DataFrame(columns=["time_start", "StationId"])
    t_on = al["TimeOn"]
    t_off = al["TimeOff"].fillna(t_on)
    t_off = t_off.where(t_off >= t_on, t_on)               # broken TimeOff
    first = t_on.dt.floor(PERIOD)
    last = (t_off - pd.Timedelta("1ns")).where(t_off > t_on, t_on).dt.floor(PERIOD)
    n = ((last - first) // PERIOD).astype(int) + 1
    idx = np.repeat(np.arange(len(al)), n.to_numpy())
    offset = np.concatenate([np.arange(k) for k in n.to_numpy()])
    out = pd.DataFrame({
        "time_start": first.to_numpy()[idx] + offset * np.timedelta64(10, "m"),
        "StationId": al["StationNr"].astype(int).to_numpy()[idx],
    })
    return out.drop_duplicates()


def frozen_flag(df: pd.DataFrame, column: str, n: int, std_max: float) -> pd.Series:
    """True for every row inside a run of n consecutive near-constant values.

    Same idea as FLASC ``filter_by_sensor_stuck_faults``: a window of n
    consecutive 10-min values (same turbine, no time gap) whose std is below
    ``std_max`` marks all n rows as frozen. ``df`` must be sorted by
    (StationId, time_start).
    """
    x = df[column].astype("float64")
    step_ok = (df["StationId"].eq(df["StationId"].shift())
               & df["time_start"].diff().eq(PERIOD))
    window_ok = step_ok.astype(int).rolling(n - 1).sum().eq(n - 1)
    stuck_end = window_ok & x.rolling(n).std(ddof=0).lt(std_max)
    flag = stuck_end.copy()
    for k in range(1, n):                 # mark the earlier rows of the window
        flag |= stuck_end.shift(-k, fill_value=False)
    return flag


def apply_northing(df: pd.DataFrame, corrections: dict[str, list[tuple[str, float]]],
                   column: str) -> pd.Series:
    """Direction in true north: (turbine-north value + offset) mod 360.

    ``corrections`` maps TurbineName -> [(valid_from, offset_deg), ...];
    "valid from" is compared with the period START (``time_start``), like RES /
    wind-up. Turbines/times without an entry get NaN (unknown, not silently
    zero). The work is done by ``pywes.northing.to_true_north``.
    """
    table = northing.from_config_dict(corrections)
    return northing.to_true_north(df[column], df["TurbineName"], df["time_start"], table)


# One RES exclusion row, e.g. "- [ 'ALL', 2022-09-03 00:00:00, 2023-02-12 00:00:00 ]  # reason"
_EXCL_ROW = re.compile(
    r"^\s*-\s*\[\s*['\"]?(?P<who>[A-Za-z0-9_]+)['\"]?\s*,"
    r"\s*(?P<start>[0-9]{4}-[0-9]{2}-[0-9]{2}[ T][0-9:]+)\s*,"
    r"\s*(?P<end>[0-9]{4}-[0-9]{2}-[0-9]{2}[ T][0-9:]+)\s*\]\s*(?:#\s*(?P<reason>.*))?$"
)


def exclusion_group(reason: str) -> str:
    """Report group of a RES exclusion period, from its comment."""
    r = reason.lower()
    if "curtail" in r:
        return "res_curtailment"
    if any(k in r for k in ("dynamic yaw", "sat test", "wake steering", "test turbine",
                            " ref", "north offset")):
        return "dynamic_yaw"
    return "res_other"


def load_res_exclusion_periods(path: Path | str, start: str | None = None,
                               end: str | None = None) -> list[dict]:
    """Read RES's exclusion file (``HOT_DY_exclusion_periods.yaml``).

    Each row ``[turbine or 'ALL', start, end]  # reason`` becomes a dict for
    ``CleaningConfig.exclude_periods`` (turbines=None for 'ALL'). Times are the
    period START in UTC (wind-up convention), end excluded. With start/end
    given, only periods that overlap [start, end) are returned.
    Raises ValueError for a line that cannot be read (nothing skipped silently).
    """
    out = []
    lo = pd.Timestamp(start) if start else None
    hi = pd.Timestamp(end) if end else None
    for n, line in enumerate(Path(path).read_text(encoding="utf-8").splitlines(), start=1):
        if not line.strip() or line.strip().startswith("#"):
            continue
        m = _EXCL_ROW.match(line)
        if m is None:
            raise ValueError(f"{path}, line {n}: not an exclusion row: {line!r}")
        p_start, p_end = pd.Timestamp(m["start"]), pd.Timestamp(m["end"])
        if (hi is not None and p_start >= hi) or (lo is not None and p_end <= lo):
            continue
        reason = (m["reason"] or "").strip()
        out.append({"start": f"{p_start:%Y-%m-%d %H:%M:%S}", "end": f"{p_end:%Y-%m-%d %H:%M:%S}",
                    "turbines": None if m["who"].upper() == "ALL" else [m["who"]],
                    "reason": f"RES: {reason}", "group": exclusion_group(reason)})
    return out


# ----------------------------------------------------------------------------
# The pipeline
# ----------------------------------------------------------------------------
class StepLog:
    """Sequential filter log: each step removes rows from the 'kept' mask."""

    def __init__(self, n_initial: int):
        self.n_initial = n_initial
        self.rows: list[dict] = []

    def add(self, name: str, flag: pd.Series, keep: pd.Series, df: pd.DataFrame) -> pd.Series:
        """Remove rows where ``flag`` is True; store the flag and the count."""
        flag = flag.fillna(False).astype(bool)
        newly = flag & keep
        df[f"flag_{len(self.rows) + 1:02d}_{name}"] = flag
        df.loc[newly, "removed_by"] = name
        keep = keep & ~flag
        self.rows.append({"step": len(self.rows) + 1, "filter": name,
                          "flagged_total": int(flag.sum()),
                          "removed": int(newly.sum()),
                          "removed_pct": 100.0 * newly.sum() / self.n_initial,
                          "remaining": int(keep.sum())})
        return keep

    def table(self) -> pd.DataFrame:
        return pd.DataFrame(self.rows)


def run_pipeline(cfg: CleaningConfig) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Load the core signals, run all filters and return the results.

    Returns
    -------
    df     : one row per raw tblSCTurbine record read, with core signals,
             ERA5 keys, every filter flag, 'removed_by' and 'keep'
    report : one row per filter (removed rows, %, remaining)
    info   : extra numbers for the summary (missing files, downtime sources...)
    """
    info: dict = {}
    start, end = pd.Timestamp(cfg.start), pd.Timestamp(cfg.end)

    # --- load ---------------------------------------------------------------
    print("Loading core columns ...")
    base, missing = load_core_table(cfg, BASE_TABLE)
    info["missing_files"] = {BASE_TABLE: missing}
    df = base
    n0 = len(df)
    log = StepLog(n0)
    keep = pd.Series(True, index=df.index)
    df["removed_by"] = pd.Series(pd.NA, index=df.index, dtype="string")

    # --- Step 1: timestamp alignment -------------------------------------------
    # TimeStamp is the END of the period (UTC). Records off the 10-min grid
    # (e.g. 10:07:00) or unparseable are removed.
    df["time_end"] = pd.to_datetime(df["TimeStamp"], format=TIMESTAMP_FORMAT, errors="coerce")
    df["time_start"] = df["time_end"] - PERIOD
    on_grid = df["time_end"].notna() & df["time_end"].eq(df["time_end"].dt.floor(PERIOD))
    keep = log.add("timestamp_alignment", ~on_grid, keep, df)

    # --- Step 2: temporal window on the period START -------------------------
    # 2024-01-01 00:00 TimeStamp = period 2023-12-31 23:50 -> out.
    # 2026-01-01 00:00 TimeStamp = period 2025-12-31 23:50 -> in.
    in_window = (df["time_start"] >= start) & (df["time_start"] < end)
    keep = log.add("temporal_window", ~in_window, keep, df)

    # --- Step 3: duplicates ---------------------------------------------------
    # Raw files overlap by one timestamp at month borders; keep the first
    # record read (file order) of each (TimeStamp, StationId).
    sub = df[keep]
    dup = pd.Series(False, index=df.index)
    dup[sub.index[sub.duplicated(KEY, keep="first")]] = True
    keep = log.add("duplicates", dup, keep, df)
    info["slots_with_record"] = int(keep.sum())     # unique (period, turbine) records

    # --- Step 4: unknown turbine ----------------------------------------------
    meta = load_turbine_metadata(cfg.sample_dir)
    df = df.merge(meta, on="StationId", how="left")
    keep.index = df.index
    keep = log.add("unknown_turbine", df["TurbineName"].isna(), keep, df)

    # --- merge the other tables (dedup each the same way first) --------------
    for table in CORE_COLUMNS:
        if table == BASE_TABLE:
            continue
        other, miss = load_core_table(cfg, table)
        info["missing_files"][table] = miss
        other = other.drop_duplicates(KEY, keep="first").drop(columns="source_order")
        df = df.merge(other, on=KEY, how="left")
    keep.index = df.index

    # --- Step 5: missing required signals --------------------------------------
    # A missing row in tblSCTurGrid etc. also shows up here as NaN.
    missing_sig = df[list(cfg.required)].isna().any(axis=1)
    keep = log.add("missing_signals", missing_sig, keep, df)

    # --- Step 6: physical range -------------------------------------------------
    out_of_range = pd.Series(False, index=df.index)
    for col, (lo, hi) in cfg.ranges.items():
        out_of_range |= (df[col] < lo) | (df[col] > hi)      # NaN -> False
    out_of_range |= df["ws"] <= cfg.ws_floor                  # logger floor
    keep = log.add("physical_range", out_of_range, keep, df)

    # --- Step 7: frozen sensors ---------------------------------------------------
    df = df.sort_values(["StationId", "time_start", "source_order"]).reset_index(drop=True)
    keep = ~df.filter(like="flag_").any(axis=1)      # same mask, new row order
    frozen = pd.Series(False, index=df.index)
    for col in cfg.frozen_columns:
        # only check rows that are still kept, so stale duplicates don't count
        frozen |= frozen_flag_on_subset(df, keep, col, cfg.frozen_n, cfg.frozen_std)
    keep = log.add("frozen_sensor", frozen, keep, df)

    # --- Step 8: downtime ---------------------------------------------------------
    sd = load_shutdown_duration(cfg.raw_dir, cfg.start, cfg.end)
    df = df.merge(sd, on=["time_start", "TurbineName"], how="left")
    keep.index = df.index
    proxy = (PERIOD_S - df["t_ready_s"].astype("float64")).clip(0.0, PERIOD_S)
    if cfg.downtime_source == "official_then_proxy":
        df["downtime_s"] = df["downtime_official_s"].fillna(proxy)
        source = np.where(df["downtime_official_s"].notna(), "official",
                          np.where(proxy.notna(), "proxy", "missing"))
    elif cfg.downtime_source == "official_only":
        df["downtime_s"] = df["downtime_official_s"]
        source = np.where(df["downtime_official_s"].notna(), "official", "missing")
    elif cfg.downtime_source == "proxy_only":
        df["downtime_s"] = proxy
        source = np.where(proxy.notna(), "proxy", "missing")
    else:
        raise ValueError(f"unknown downtime_source: {cfg.downtime_source}")
    df["downtime_source"] = pd.Categorical(source, categories=["official", "proxy", "missing"])
    info["downtime_source_counts"] = df.loc[keep, "downtime_source"].value_counts().to_dict()
    # Missing downtime information counts as "not clean" (cannot be checked).
    down = df["downtime_s"].isna() | (df["downtime_s"] > cfg.max_downtime_s)
    keep = log.add("downtime", down, keep, df)

    # --- Step 9: turbine error timer -----------------------------------------------
    keep = log.add("turbine_error", df["t_error_s"] > cfg.max_error_s, keep, df)

    # --- Steps 10-12: alarms --------------------------------------------------------
    alarms = load_alarms(cfg)
    lookup = pd.read_csv(Path(cfg.sample_dir) / "Hill_of_Towie_alarms_description.csv")
    stop_codes = set(cfg.extra_stopping_codes)
    if cfg.use_lookup_stopping_codes:
        stop_codes |= set(lookup.loc[lookup["Stopping"] == 1, "Alarm Code"].astype(int))
    info["stopping_codes"] = sorted(stop_codes)
    for name, codes in [("stopping_alarm", stop_codes),
                        ("icing_alarm", set(cfg.ice_codes)),
                        ("cable_unwind", set(cfg.unwind_codes))]:
        hits = alarm_periods(alarms, codes).assign(_hit=True)
        flag = df[["time_start", "StationId"]].merge(hits, on=["time_start", "StationId"],
                                                     how="left")["_hit"].fillna(False)
        flag.index = df.index
        if name == "cable_unwind":                     # counter in tblSCTurFlag too
            flag |= df["yaw_unwind_counts"].fillna(0) > 0
        keep = log.add(name, flag, keep, df)

    # --- Step 13: idling and high pitch -----------------------------------------------
    df["pitch_mean"] = df[["pitch_a", "pitch_b", "pitch_c"]].mean(axis=1)
    keep = log.add("idling", df["power"] < cfg.idle_power_kw, keep, df)
    keep = log.add("high_pitch", df["pitch_mean"] > cfg.max_pitch_deg, keep, df)

    # --- Step 14: curtailment / derating ------------------------------------------------
    reduced_ref = df["power_ref"] < cfg.rated_power_kw - cfg.power_ref_tol_kw
    if cfg.curtail_only_when_limiting:
        reduced_ref &= df["power"] >= cfg.limiting_fraction * df["power_ref"]
    keep = log.add("power_setpoint_reduced", reduced_ref, keep, df)
    low_power = ((df["ws"] > cfg.low_power_ws)
                 & (df["power"] < cfg.low_power_fraction * cfg.rated_power_kw))
    keep = log.add("low_power_high_wind", low_power, keep, df)
    if cfg.flat_power_std_kw is not None:
        flat = (df["power"] < cfg.flat_power_max_kw) & (df["power_std"] < cfg.flat_power_std_kw)
        keep = log.add("flat_power", flat, keep, df)

    # --- Step 15: optional filters ---------------------------------------------------------
    if cfg.yaw_std_max_deg is not None:
        # NOTE: 10-min std of a wrapped angle may be wrong near 0/360 deg
        # (data dictionary) -> such periods are removed too (conservative).
        keep = log.add("yaw_unsteady", df["nacelle_pos_std"] > cfg.yaw_std_max_deg, keep, df)
    # --- Step 16: exclusion periods, one report step per group ------------------------
    groups: dict[str, pd.Series] = {}
    for p in cfg.exclude_periods:
        sel = (df["time_start"] >= pd.Timestamp(p["start"])) & \
              (df["time_start"] < pd.Timestamp(p["end"]))
        if p.get("turbines"):
            sel &= df["TurbineName"].isin(p["turbines"])
        name = p.get("group", "manual_exclusion")
        groups[name] = groups.get(name, pd.Series(False, index=df.index)) | sel
    for name, flag in groups.items():
        keep = log.add(name, flag, keep, df)
    info["exclusion_periods"] = len(cfg.exclude_periods)

    df["keep"] = keep

    # --- ERA5-compatible keys ------------------------------------------------------------
    df["era5_time"] = nearest_era5_hour(df["time_end"], period="10min", timestamp_is_end=True)
    df["cold_le_3C"] = df["t_amb"] <= 3.0      # marker only (IEA T19 reference uses > +3 C)
    df["post_tuneup"] = (df["TurbineName"].isin(TUNEUP_TURBINES)
                         & (df["time_start"] >= pd.Timestamp(TUNEUP_FIRST_AFTER)))   # marker only
    if cfg.northing_corrections_deg:
        df["wd_true"] = apply_northing(df, cfg.northing_corrections_deg, "wd")
        df["nacelle_pos_true"] = apply_northing(df, cfg.northing_corrections_deg, "nacelle_pos")
        n_nan = int((df["keep"] & df["nacelle_pos"].notna() & df["nacelle_pos_true"].isna()).sum())
        info["northing"] = (f"applied ({cfg.northing_source or 'corrections from config'}); "
                            f"kept records without a correction: {n_nan}")
    else:
        info["northing"] = "NOT applied; directions are turbine north"

    # completeness: how many (period, turbine) slots exist in the window
    n_periods = int((end - start) / PERIOD)
    info["expected_slots"] = n_periods * len(meta)
    info["window"] = f"period start {cfg.start} to {cfg.end} (end excluded), UTC"
    return df, log.table(), info


def frozen_flag_on_subset(df: pd.DataFrame, keep: pd.Series, column: str,
                          n: int, std_max: float) -> pd.Series:
    """Run ``frozen_flag`` only on kept rows; other rows get False."""
    sub = df.loc[keep, ["StationId", "time_start", column]]
    flag = pd.Series(False, index=df.index)
    flag.loc[sub.index] = frozen_flag(sub, column, n, std_max).to_numpy()
    return flag


# ----------------------------------------------------------------------------
# Outputs
# ----------------------------------------------------------------------------
def write_clean_tables(cfg: CleaningConfig, df: pd.DataFrame) -> pd.DataFrame:
    """Write the kept records in the ORIGINAL raw schema and text format.

    Lines are copied byte-for-byte from the raw files (same header, same
    number formatting). Output layout mirrors the raw layout:
        out_dir/scada_clean_<tag>/<YYYY>/<table>_<YYYY>_<MM>.csv
    with one difference: each record goes to the month of its period START,
    so no record appears in two files.

    Returns a table with the number of lines written per raw table.
    """
    kept = df.loc[df["keep"], ["TimeStamp", "StationId", "time_start"]]
    keys = (kept["TimeStamp"] + "," + kept["StationId"].astype(int).astype(str))
    month_of = dict(zip(keys, kept["time_start"].dt.strftime("%Y_%m")))
    root = Path(cfg.out_dir) / f"scada_clean_{cfg.tag}"
    counts = []

    for table in cfg.output_tables:
        files, boundary = files_for_window(cfg, table)
        handles: dict[str, object] = {}
        headers: dict[str, str] = {}
        written: set[str] = set()
        try:
            sources = [f for f in files if f.exists()]
            if boundary is not None and boundary.exists():
                sources.append(boundary)
            for src in sources:
                if src == boundary:
                    lines = io.StringIO(read_head_text(src, cfg.end), newline="")
                else:
                    lines = open(src, newline="", encoding="utf-8", errors="surrogateescape")
                with lines as fh:
                    header = fh.readline()
                    for line in fh:
                        i = line.find(",")
                        j = line.find(",", i + 1)
                        key = line[:j]
                        month = month_of.get(key)
                        if month is None or key in written:
                            continue                 # removed, or duplicate
                        if month not in handles:
                            path = root / month[:4] / f"{table}_{month}.csv"
                            path.parent.mkdir(parents=True, exist_ok=True)
                            handles[month] = open(path, "w", newline="", encoding="utf-8",
                                                  errors="surrogateescape")
                            handles[month].write(header)
                            headers[month] = header
                        elif headers[month] != header:
                            raise ValueError(
                                f"Schema change inside output month {month} of {table} "
                                f"(file {src.name}). Columns differ, cannot merge.")
                        handles[month].write(line)
                        written.add(key)
        finally:
            for h in handles.values():
                h.close()
        counts.append({"table": table, "rows_written": len(written),
                       "rows_kept": len(month_of)})
        if len(written) != len(month_of):
            print(f"WARNING: {table}: {len(month_of) - len(written)} kept records have "
                  f"no row in this table (written {len(written)} of {len(month_of)}).")
    return pd.DataFrame(counts)


def compact_columns(df: pd.DataFrame) -> list[str]:
    """Columns of the compact flags table (keys, ERA5 keys, signals, flags)."""
    front = ["TimeStamp", "time_start", "time_end", "era5_time", "StationId",
             "TurbineName", "latitude", "longitude"]
    signals = [s for t in CORE_COLUMNS.values() for s in t.values()]
    extra = ["pitch_mean", "downtime_s", "downtime_official_s", "downtime_source",
             "cold_le_3C", "post_tuneup"] + [c for c in ("wd_true", "nacelle_pos_true")
                                             if c in df.columns]
    flags = [c for c in df.columns if c.startswith("flag_")]
    return front + signals + extra + flags + ["removed_by", "keep"]


def write_flags_parquet(cfg: CleaningConfig, df: pd.DataFrame) -> Path:
    """All records (kept and removed) with signals and flags, as Parquet."""
    path = Path(cfg.out_dir) / f"scada_flags_{cfg.tag}.parquet"
    out = df[compact_columns(df)].copy()
    out["removed_by"] = out["removed_by"].astype("string")
    out.to_parquet(path, index=False)
    return path


def write_netcdf(cfg: CleaningConfig, df: pd.DataFrame) -> Path:
    """Cleaned core signals on a (time, turbine) grid, CF-style like ERA5.

    * ``time``: 10-min period END (UTC, = SCADA TimeStamp), full grid of the
      window; ``time_bnds`` gives [start, end] (CF cell bounds).
    * ``latitude`` / ``longitude`` per turbine (deg, WGS84).
    * Signals are NaN where the record was removed; ``keep`` = 1 kept,
      0 removed, missing record = 0.
    """
    import xarray as xr

    start, end = pd.Timestamp(cfg.start), pd.Timestamp(cfg.end)
    grid = pd.date_range(start + PERIOD, end, freq=PERIOD)      # period ends
    sub = df[df["keep"]].copy()
    var_units = {"power": "kW", "ws": "m s-1", "wd": "degree", "nacelle_pos": "degree",
                 "pitch_mean": "degree", "gen_rpm": "rpm", "power_ref": "kW",
                 "t_amb": "degC", "downtime_s": "s"}
    for c in ("wd_true", "nacelle_pos_true"):
        if c in sub.columns:
            var_units[c] = "degree"
    ds = (sub.set_index(["time_end", "TurbineName"])[list(var_units)]
          .astype("float32").to_xarray()
          .rename({"time_end": "time", "TurbineName": "turbine"}))
    meta = load_turbine_metadata(cfg.sample_dir).set_index("TurbineName")
    ds = ds.reindex(time=grid, turbine=meta.index.to_numpy())
    ds["keep"] = ds["power"].notnull().astype("int8")
    ds = ds.assign_coords(
        latitude=("turbine", meta["latitude"].to_numpy()),
        longitude=("turbine", meta["longitude"].to_numpy()),
        station_id=("turbine", meta["StationId"].to_numpy()),
    )
    ds["time_bnds"] = (("time", "bnds"), np.stack([grid - PERIOD, grid], axis=1))
    ds["time"].attrs.update(long_name="end of 10-min averaging period (UTC)",
                            bounds="time_bnds")
    ds["latitude"].attrs.update(units="degrees_north", standard_name="latitude")
    ds["longitude"].attrs.update(units="degrees_east", standard_name="longitude")
    names = {"wd": "wind direction, turbine north frame",
             "nacelle_pos": "nacelle position, turbine north frame",
             "wd_true": "wind direction, true north (northing corrected)",
             "nacelle_pos_true": "nacelle position, true north (northing corrected)"}
    for v, u in var_units.items():
        ds[v].attrs.update(units=u, long_name=names.get(v, v))
    ds["keep"].attrs.update(long_name="1 = record passed all cleaning steps")
    ds.attrs.update(
        title="Hill of Towie 10-min SCADA, cleaned (single-turbine filters)",
        source="Hill of Towie open dataset (RES/TRIG), CC-BY-4.0, doi:10.5281/zenodo.14870021",
        history=f"pywes.scada_cleaning, window {cfg.start} to {cfg.end} (period start)",
        Conventions="CF-1.8",
        northing=cfg.northing_source if cfg.northing_corrections_deg else "none",
        exclusion_groups=", ".join(dict.fromkeys(p.get("group", "manual_exclusion")
                                                 for p in cfg.exclude_periods)) or "none",
    )
    path = Path(cfg.out_dir) / f"scada_clean_core_{cfg.tag}.nc"
    enc = {v: {"zlib": True, "complevel": 4} for v in var_units}
    # same time units for time and its bounds (CF requirement)
    enc["time"] = {"units": "minutes since 1970-01-01 00:00:00", "calendar": "standard"}
    enc["time_bnds"] = {"units": "minutes since 1970-01-01 00:00:00", "calendar": "standard"}
    ds.to_netcdf(path, encoding=enc)
    return path


def print_summary(report: pd.DataFrame, info: dict, n_initial: int, n_final: int) -> None:
    """Print the filter log: initial rows, removed per step, final rows."""
    line = "-" * 78
    print(line)
    print(f"SCADA cleaning summary - tblSCTurbine records, {info.get('window', '')}")
    print(line)
    print(f"Initial rows read: {n_initial:>12,d}")
    print(f"{'step':>4}  {'filter':<26}{'removed':>12}{'% of initial':>14}{'remaining':>14}")
    for r in report.itertuples():
        print(f"{r.step:>4}  {r.filter:<26}{r.removed:>12,d}{r.removed_pct:>13.2f}%"
              f"{r.remaining:>14,d}")
    print(line)
    pct = 100.0 * n_final / n_initial if n_initial else 0.0
    print(f"Final retained rows: {n_final:>10,d}  ({pct:.2f} % of initial)")
    exp = info.get("expected_slots", 0)
    if exp:
        got = info.get("slots_with_record", 0)
        print(f"Completeness: {got:,d} of {exp:,d} expected (period, turbine) slots have a "
              f"record ({100.0 * got / exp:.2f} %); retained = {100.0 * n_final / exp:.2f} %")
    print(f"Downtime source (rows reaching step 8): {info.get('downtime_source_counts')}")
    print(f"Stopping alarm codes used: {info.get('stopping_codes')}")
    print(f"Northing: {info.get('northing')}")
    print(f"Exclusion periods used: {info.get('exclusion_periods', 0)}")
    for table, miss in info.get("missing_files", {}).items():
        if miss:
            print(f"Missing raw files for {table}: {len(miss)} (e.g. {Path(miss[0]).name})")
    print(line)
