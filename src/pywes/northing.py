"""Northing correction for Hill of Towie SCADA: turbine north -> true north.

Why
---
Yaw and direction signals in the SCADA (``wtc_NacelPos``, ``wtc_YawPos``,
``wtc_ActualWindDirection``) are in each turbine's own calibrated north
("turbine north"). The data dictionary warns that these calibration errors
"can be large and can change in time" and must be corrected before any
directional analysis. Example from our data (2024-2025, no correction):
nacelle T11 - nacelle T14 = +9.7 deg (median).

Where the corrections come from
-------------------------------
RES published the corrections in their companion repository
(github.com/resgroup/hill-of-towie-open-source-analysis, BSD-3, licence copy in
``data/sample/LICENSE_RES_hill-of-towie-open-source-analysis.txt``). We use
the files copied to ``data/sample/``:

* ``HOT_DY_Design_northing.yaml`` (DEFAULT) - used by RES for the Dynamic Yaw
  campaign. Entries start 2019-10-14 and include the changes in 2025, so it
  covers our whole 2024-2026 window.
* ``optimized_northing_corrections.yaml`` - result of RES ``northing.py``
  (wind-up ``optimize_northing_corrections`` against ERA5T 100 m wind
  direction) on dataset v1, so it only covers 2016-01 to 2024-08.

File format: one row per change, ``- ['T03', 2025-08-07 13:10:00, -177.46]``
= (turbine, valid from [UTC], correction [deg]). A row is valid from its time
until the next row of the same turbine.

How the correction is applied (checked in RES code)
---------------------------------------------------
    true_north = (turbine_north + correction) mod 360

This is how RES applies it in ``scripts/wfc_analysis_2026/inspect_data.py``.
It also matches our own check: T11 gets -20.9 deg and T14 -11.6 deg, which
removes the +9.7 deg T11 - T14 offset.

Time convention: RES / wind-up work with ``TimeStamp_StartFormat`` (period
START = SCADA TimeStamp - 10 min, see RES ``s02_generate_datapack.py``). So we
compare the "valid from" times with the period START. For data with period-END
timestamps (our NetCDF ``time``), pass ``timestamp_is_end=True``.

Limits
------
* The corrections are RES's estimates against ERA5, not survey values; RES
  says the turbine north calibrations "are apparently wrong quite often".
* Times before the first row of a turbine get NaN (unknown), never 0.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

from pywes.paths import SAMPLE_DIR

PERIOD = pd.Timedelta("10min")

DEFAULT_NORTHING_FILE = SAMPLE_DIR / "HOT_DY_Design_northing.yaml"   # 2019-10 -> 2026
V1_NORTHING_FILE = SAMPLE_DIR / "optimized_northing_corrections.yaml"  # 2016-01 -> 2024-08

# One correction row, e.g.  "    - ['T01', 2016-01-01 00:00:00, 156.77]"
_ROW = re.compile(
    r"^\s*-\s*\[\s*['\"]?(?P<turbine>[A-Za-z0-9_]+)['\"]?\s*,"
    r"\s*(?P<time>[0-9]{4}-[0-9]{2}-[0-9]{2}[ T][0-9:]+)\s*,"
    r"\s*(?P<corr>[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)\s*\]\s*(?:#.*)?$"
)


# ----------------------------------------------------------------------------
# Loading
# ----------------------------------------------------------------------------
def load_northing_table(path: Path | str = DEFAULT_NORTHING_FILE) -> pd.DataFrame:
    """Read a RES northing YAML file into a table.

    Returns columns turbine (str), valid_from (UTC, tz-naive), correction_deg
    (float), sorted by turbine and time. Raises ValueError for a line that is
    not empty, not a comment and not a correction row (nothing is skipped
    silently).
    """
    rows = []
    for n, line in enumerate(Path(path).read_text(encoding="utf-8").splitlines(), start=1):
        text = line.strip()
        if not text or text.startswith("#"):
            continue
        m = _ROW.match(line)
        if m is None:
            raise ValueError(f"{path}, line {n}: not a northing row: {line!r}")
        rows.append((m["turbine"], pd.Timestamp(m["time"]), float(m["corr"])))
    if not rows:
        raise ValueError(f"{path}: no northing rows found")
    table = pd.DataFrame(rows, columns=["turbine", "valid_from", "correction_deg"])
    table = table.sort_values(["turbine", "valid_from"]).reset_index(drop=True)
    if table.duplicated(["turbine", "valid_from"]).any():
        raise ValueError(f"{path}: two rows for the same turbine and time")
    table.attrs["source"] = Path(path).name
    return table


def as_config_dict(table: pd.DataFrame) -> dict[str, list[tuple[str, float]]]:
    """Convert the table to the format of ``CleaningConfig.northing_corrections_deg``:
    {"T01": [("2019-10-14 23:50:00", -23.57), ...], ...}."""
    out: dict[str, list[tuple[str, float]]] = {}
    for r in table.itertuples(index=False):
        out.setdefault(r.turbine, []).append((f"{r.valid_from:%Y-%m-%d %H:%M:%S}",
                                              float(r.correction_deg)))
    return out


def from_config_dict(corrections: dict[str, list[tuple[str, float]]]) -> pd.DataFrame:
    """Inverse of ``as_config_dict`` (for code that still uses the dict format)."""
    rows = [(t, pd.Timestamp(v), float(c)) for t, entries in corrections.items()
            for v, c in entries]
    table = pd.DataFrame(rows, columns=["turbine", "valid_from", "correction_deg"])
    return table.sort_values(["turbine", "valid_from"]).reset_index(drop=True)


# ----------------------------------------------------------------------------
# Applying
# ----------------------------------------------------------------------------
def correction_at(table: pd.DataFrame, turbines, times, timestamp_is_end: bool = False
                  ) -> pd.Series:
    """Correction [deg] valid for each (turbine, time).

    turbines : turbine names ("T01" ... "T21"), one per row
    times    : period START times (UTC); with timestamp_is_end=True they are
               period END times and 10 min are subtracted first
    Returns a float Series (same index as ``times`` if it is a Series). NaN
    where the turbine has no row in the table or the time is before its
    first row.
    """
    index = times.index if isinstance(times, pd.Series) else None
    t = pd.to_datetime(pd.Series(np.asarray(times), index=index))
    if timestamp_is_end:
        t = t - PERIOD
    names = pd.Series(np.asarray(turbines, dtype=object), index=t.index)
    out = pd.Series(np.nan, index=t.index, dtype="float64")
    for name, entries in table.groupby("turbine", sort=False):
        rows = names.eq(name).to_numpy()
        if not rows.any():
            continue
        starts = entries["valid_from"].to_numpy(dtype="datetime64[ns]")
        values = entries["correction_deg"].to_numpy(dtype="float64")
        # index of the last row with valid_from <= t (-1 = before the first row)
        k = np.searchsorted(starts, t[rows].to_numpy(dtype="datetime64[ns]"), side="right") - 1
        vals = np.where(k >= 0, values[np.clip(k, 0, None)], np.nan)
        vals[np.isnat(t[rows].to_numpy(dtype="datetime64[ns]"))] = np.nan
        out[rows] = vals
    return out


def to_true_north(directions, turbines, times, table: pd.DataFrame,
                  timestamp_is_end: bool = False) -> pd.Series:
    """Turbine-north directions [deg] -> true north [deg, 0-360).

    true = (turbine_north + correction) mod 360. NaN where the direction or
    the correction is unknown.
    """
    corr = correction_at(table, turbines, times, timestamp_is_end)
    d = pd.Series(np.asarray(directions, dtype="float64"), index=corr.index)
    return np.mod(d + corr, 360.0)


def apply_to_frame(df: pd.DataFrame, table: pd.DataFrame,
                   columns: tuple[str, ...] = ("wd", "nacelle_pos"),
                   turbine_col: str = "TurbineName", time_col: str = "time_start",
                   timestamp_is_end: bool = False, suffix: str = "_true") -> pd.DataFrame:
    """Add true-north copies ``<column><suffix>`` to a long table (one row per
    turbine and period). The original turbine-north columns stay unchanged."""
    out = df.copy()
    corr = correction_at(table, out[turbine_col], out[time_col], timestamp_is_end)
    for col in columns:
        out[f"{col}{suffix}"] = np.mod(out[col].astype("float64") + corr, 360.0)
    return out


def apply_to_pair(df: pd.DataFrame, table: pd.DataFrame, turb_a: str, turb_b: str,
                  columns: tuple[str, ...] = ("nacelle_pos",),
                  timestamp_is_end: bool = True) -> pd.DataFrame:
    """Correct a wide pair table from ``pywes.pairs.load_pair`` (returns a copy).

    ``<col>_A`` / ``<col>_B`` become TRUE north; the uncorrected values are kept
    as ``<col>_tn_A`` / ``<col>_tn_B`` (tn = turbine north). The index is the
    10-min period END (``timestamp_is_end=True``, NetCDF convention).
    Adds ``df.attrs["northing"]`` with the file name.
    """
    out = df.copy()
    for turb, tag in ((turb_a, "A"), (turb_b, "B")):
        corr = correction_at(table, np.repeat(turb, len(out)), pd.Series(out.index, index=out.index),
                             timestamp_is_end)
        for col in columns:
            out[f"{col}_tn_{tag}"] = out[f"{col}_{tag}"]
            out[f"{col}_{tag}"] = np.mod(out[f"{col}_{tag}"].astype("float64") + corr, 360.0)
    out.attrs = dict(df.attrs, northing=table.attrs.get("source", "table"))
    return out


# ----------------------------------------------------------------------------
# Reporting and checks
# ----------------------------------------------------------------------------
def segments(table: pd.DataFrame, start, end, turbines=None) -> pd.DataFrame:
    """Correction segments that overlap the window [start, end) (UTC).

    Returns turbine, valid_from, valid_to (next row, or NaT = open end),
    correction_deg - useful to print which corrections an analysis uses and to
    spot big jumps (e.g. T03 ~180 deg from 2025-08-07).
    """
    t = table.sort_values(["turbine", "valid_from"]).copy()
    t["valid_to"] = t.groupby("turbine")["valid_from"].shift(-1)
    start, end = pd.Timestamp(start), pd.Timestamp(end)
    ok = (t["valid_from"] < end) & (t["valid_to"].isna() | (t["valid_to"] > start))
    if turbines is not None:
        ok &= t["turbine"].isin(list(turbines))
    return t.loc[ok, ["turbine", "valid_from", "valid_to", "correction_deg"]].reset_index(drop=True)


def missing_coverage(table: pd.DataFrame, start, turbines) -> list[str]:
    """Turbines that have no correction at ``start`` (first row later or none)."""
    first = table.groupby("turbine")["valid_from"].min()
    return [t for t in turbines if t not in first.index or first[t] > pd.Timestamp(start)]


def compare_tables(a: pd.DataFrame, b: pd.DataFrame, start, end, freq: str = "1D"
                   ) -> pd.DataFrame:
    """Compare two correction tables on a time grid (e.g. the v1 file vs the
    Dynamic Yaw file where both are valid).

    Returns per turbine: n grid points with both defined, median and maximum
    absolute circular difference a - b [deg].
    """
    grid = pd.date_range(start, end, freq=freq, inclusive="left")
    rows = []
    for turb in sorted(set(a["turbine"]) | set(b["turbine"])):
        names = np.repeat(turb, len(grid))
        ca = correction_at(a, names, pd.Series(grid))
        cb = correction_at(b, names, pd.Series(grid))
        d = ((ca - cb + 180.0) % 360.0 - 180.0).abs().dropna()
        rows.append({"turbine": turb, "n": len(d),
                     "median_abs_diff_deg": d.median() if len(d) else np.nan,
                     "max_abs_diff_deg": d.max() if len(d) else np.nan})
    return pd.DataFrame(rows).set_index("turbine")


def monthly_offset(times, turbines, directions_true, reference, min_n: int = 50
                   ) -> pd.DataFrame:
    """Check after correcting: monthly circular median of (direction - reference).

    directions_true : corrected directions [deg, true north]
    reference       : reference directions [deg, true north], e.g. ERA5 wd100
    Returns a table month x turbine [deg]; months with < min_n values are NaN.
    After a good northing correction the values should be close to 0 and show
    no steps; a step means a re-calibration that the table misses.
    """
    d = pd.DataFrame({
        "month": pd.to_datetime(pd.Series(np.asarray(times))).dt.to_period("M").to_numpy(),
        "turbine": np.asarray(turbines, dtype=object),
        "diff": (np.asarray(directions_true, dtype="float64")
                 - np.asarray(reference, dtype="float64") + 180.0) % 360.0 - 180.0,
    }).dropna()
    g = d.groupby(["month", "turbine"])["diff"]
    med = g.median().where(g.count() >= min_n)
    return med.unstack("turbine")
