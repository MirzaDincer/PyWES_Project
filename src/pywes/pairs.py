"""Compare an upstream (free-stream) turbine A with a downstream turbine B.

Turbine A gets free-stream wind; turbine B stands behind A and sits in A's
wake when the wind blows from A towards B. The helpers below load the pair from
the cleaned SCADA NetCDF, join ERA5 stability, filter by direction and compute
wake statistics (speed ratio, power ratio, deficits) per stability class.

Conventions
-----------
* Directions: meteorological, deg, from north, clockwise.
* SCADA ``time`` = END of the 10-min period (UTC).
* ERA5 ``wd100`` is in TRUE north. SCADA ``nacelle_pos`` from the NetCDF is in
  TURBINE north; correct it with ``pywes.northing.apply_to_pair`` right after
  ``load_pair`` (RES corrections). The helpers here do not correct anything
  themselves; they work in whatever frame the columns are in.
* Column names after ``load_pair``: ``<signal>_A`` and ``<signal>_B``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from pywes.stability import CLASS_ORDER, nearest_era5_hour

# ERA5 columns used from data/processed/era5_stability_hourly.csv
ERA5_COLUMNS = ["wd100", "ws100", "ws10", "inv_L", "stability", "kept"]


# ----------------------------------------------------------------------------
# Angles
# ----------------------------------------------------------------------------
def angle_diff(a, b):
    """Circular difference a - b [deg], wrapped to [-180, 180).

    Example: angle_diff(5, 355) = 10, angle_diff(355, 5) = -10.
    Works with scalars, NumPy arrays and pandas Series.
    """
    return (a - b + 180.0) % 360.0 - 180.0


def circular_mean(deg) -> float:
    """Mean direction [deg, 0-360) of an array of directions; NaN if empty."""
    rad = np.radians(np.asarray(deg, dtype=float))
    rad = rad[~np.isnan(rad)]
    if rad.size == 0:
        return np.nan
    mean = float(np.mod(np.degrees(np.arctan2(np.sin(rad).mean(), np.cos(rad).mean())), 360.0))
    return 0.0 if mean >= 360.0 else mean          # rounding can give exactly 360.0


# ----------------------------------------------------------------------------
# Geometry
# ----------------------------------------------------------------------------
def pair_geometry(layout: pd.DataFrame, turb_a: str, turb_b: str, rotor_d: float) -> dict:
    """Alignment direction, spacing and elevation difference of a pair.

    layout  : table with columns "Turbine Name", x_utm, y_utm [m], hub_elev_m [m a.s.l.]
    rotor_d : rotor diameter [m]

    Returns a dict with
      align_dir_deg : wind direction [deg, from north] that puts B straight
                      behind A (= compass bearing from B to A)
      spacing_m, spacing_D : distance A-B [m] and [rotor diameters]
      dz_hub_m      : hub elevation of B minus A [m]

    The bearing is on the UTM grid. Grid north differs from true north by the
    meridian convergence, which is < 0.1 deg at Hill of Towie (close to the
    UTM zone 30 central meridian, 3 deg W).
    """
    lay = layout.set_index("Turbine Name")
    a, b = lay.loc[turb_a], lay.loc[turb_b]
    dx, dy = a["x_utm"] - b["x_utm"], a["y_utm"] - b["y_utm"]      # vector B -> A
    align = float(np.mod(np.degrees(np.arctan2(dx, dy)), 360.0))
    dist = float(np.hypot(dx, dy))
    return {
        "align_dir_deg": align,
        "spacing_m": dist,
        "spacing_D": dist / rotor_d,
        "dz_hub_m": float(b["hub_elev_m"] - a["hub_elev_m"]),
    }


# ----------------------------------------------------------------------------
# Loading and joining
# ----------------------------------------------------------------------------
def load_pair(nc_path: Path, turb_a: str, turb_b: str, start: str, end: str,
              variables: tuple[str, ...] = ("ws", "nacelle_pos", "power", "keep")) -> pd.DataFrame:
    """Load two turbines from the cleaned SCADA NetCDF on one time index.

    The index is the 10-min period END (UTC), selected with start <= time <= end.
    Signals are NaN where the cleaning removed the record; ``keep_A`` /
    ``keep_B`` are 1 = kept, 0 = removed or missing.
    Units: ws [m/s], nacelle_pos [deg, turbine north], power [kW].
    ``df.attrs["nc_attrs"]`` holds the NetCDF global attributes (e.g. which
    northing and exclusion periods the cleaning used).
    """
    import xarray as xr

    with xr.open_dataset(nc_path) as ds:
        sub = ds[list(variables)].sel(turbine=[turb_a, turb_b], time=slice(start, end)).load()
        nc_attrs = dict(ds.attrs)
    frames = []
    for turb, tag in ((turb_a, "A"), (turb_b, "B")):
        part = sub.sel(turbine=turb).drop_vars(
            [c for c in ("turbine", "latitude", "longitude", "station_id") if c in sub.coords]
        ).to_dataframe()
        frames.append(part[list(variables)].add_suffix(f"_{tag}"))
    df = pd.concat(frames, axis=1)
    df.index.name = "time"
    df.attrs.update(turbine_A=turb_a, turbine_B=turb_b, nc_attrs=nc_attrs)
    return df


def load_era5_table(csv_path: Path) -> pd.DataFrame:
    """Hourly ERA5 stability table (from analysis/stability_by_sector.py), index = UTC hour."""
    era5 = pd.read_csv(csv_path, parse_dates=["time_utc"], index_col="time_utc")
    era5["stability"] = pd.Categorical(era5["stability"], categories=CLASS_ORDER, ordered=True)
    era5["kept"] = era5["kept"].astype(bool)
    return era5


def era5_to_hub_height(ws10: pd.Series, ws100: pd.Series, hub_height: float
                       ) -> tuple[pd.Series, pd.Series]:
    """Scale ERA5 wind speed to hub height with the ERA5 power-law shear.

    alpha = ln(ws100 / ws10) / ln(100 / 10)
    ws_hub = ws100 * (hub_height / 100) ** alpha

    For 10 m < hub_height < 100 m this is an interpolation (log-linear in
    height between the two ERA5 levels), not an extrapolation.
    Returns (ws_hub [m/s], alpha [-]); NaN where ws10 or ws100 <= 0.
    """
    ok = (ws10 > 0) & (ws100 > 0)
    alpha = pd.Series(np.nan, index=ws10.index)
    alpha[ok] = np.log(ws100[ok] / ws10[ok]) / np.log(100.0 / 10.0)
    ws_hub = ws100 * (hub_height / 100.0) ** alpha
    return ws_hub, alpha


def join_era5(df: pd.DataFrame, era5: pd.DataFrame, hub_height: float
              ) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Add the nearest ERA5 hour and its values to every 10-min record.

    Adds columns era5_time, wd100, ws100, ws10, inv_L, stability, kept,
    alpha_era5 and ws_hub_era5 [m/s]. Then removes, step by step:
      1. no ERA5 hour in the table
      2. ERA5 hour not kept by the stability filters (ws100 < 3 m/s or
         |10 m / L| > 1, Tuna et al. 2026 - see stability.apply_filters)
      3. ERA5 shear exponent not defined (ws10 <= 0)
    Returns (joined table, report of removed 10-min timestamps).
    """
    out = df.copy()
    out["era5_time"] = nearest_era5_hour(out.index, "10min", timestamp_is_end=True)
    out = out.join(era5[ERA5_COLUMNS], on="era5_time")
    out["ws_hub_era5"], out["alpha_era5"] = era5_to_hub_height(out["ws10"], out["ws100"], hub_height)

    steps = [
        ("ERA5 hour available", out["inv_L"].notna() | out["wd100"].notna()),
        ("ERA5 hour kept (stability filters)", out["kept"].eq(True)),
        ("ERA5 shear exponent defined", out["alpha_era5"].notna()),
    ]
    keep, report = sequential_filter(steps, len(out))
    return out[keep].copy(), report


# ----------------------------------------------------------------------------
# Filters
# ----------------------------------------------------------------------------
def sequential_filter(steps: list[tuple[str, pd.Series]], n_initial: int
                      ) -> tuple[pd.Series, pd.DataFrame]:
    """Apply boolean masks one after another and count what each step removes.

    steps : list of (name, mask); mask True = record passes this step.
    Returns (final mask, report with columns step, filter, removed,
    removed_pct (of n_initial), remaining).
    """
    keep = None
    rows = []
    remaining = n_initial
    for i, (name, mask) in enumerate(steps, start=1):
        mask = mask.fillna(False).astype(bool)
        keep = mask if keep is None else keep & mask
        removed = remaining - int(keep.sum())
        remaining = int(keep.sum())
        rows.append((i, name, removed, 100.0 * removed / max(n_initial, 1), remaining))
    report = pd.DataFrame(rows, columns=["step", "filter", "removed", "removed_pct", "remaining"])
    return keep, report


def filter_pair(df: pd.DataFrame, wind_dir_deg: float, dir_tol_deg: float,
                same_dir_tol_deg: float, ws_min: float, ws_max: float,
                p_min_kw: float) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Direction and speed filters for the wake pair (each step is reported).

    1. both turbines have a kept record at this timestamp
    2. ERA5 wd100 within wind_dir_deg +/- dir_tol_deg        [deg, true north]
    3. |nacelle_A - nacelle_B| <= same_dir_tol_deg            [deg, true north
       after pywes.northing.apply_to_pair, otherwise turbine north]
    4. ws_A within [ws_min, ws_max]                           [m/s]
    5. power_A >= p_min_kw (avoids dividing by ~0 in P_B / P_A) [kW]

    Adds columns rel_dir (ERA5 wd100 - wind_dir_deg) and dnac (nacelle_A - nacelle_B).
    """
    out = df.copy()
    out["rel_dir"] = angle_diff(out["wd100"], wind_dir_deg)
    out["dnac"] = angle_diff(out["nacelle_pos_A"], out["nacelle_pos_B"])
    steps = [
        ("both turbines kept", (out["keep_A"] == 1) & (out["keep_B"] == 1)),
        (f"ERA5 wd100 in {wind_dir_deg:.1f} +/- {dir_tol_deg:g} deg", out["rel_dir"].abs() <= dir_tol_deg),
        (f"|nacelle_A - nacelle_B| <= {same_dir_tol_deg:g} deg", out["dnac"].abs() <= same_dir_tol_deg),
        (f"ws_A in [{ws_min:g}, {ws_max:g}] m/s", out["ws_A"].between(ws_min, ws_max)),
        (f"power_A >= {p_min_kw:g} kW", out["power_A"] >= p_min_kw),
    ]
    keep, report = sequential_filter(steps, len(out))
    return out[keep].copy(), report


# ----------------------------------------------------------------------------
# SCADA vs ERA5 check
# ----------------------------------------------------------------------------
def hourly_mean(df: pd.DataFrame, tag: str, min_records: int) -> pd.DataFrame:
    """Average the kept 10-min records of one turbine per ERA5 hour.

    tag         : "A" or "B"
    min_records : minimum number of kept 10-min records in the hour (max 6)
    Returns one row per hour with ws [m/s], nacelle_pos [deg, circular mean],
    n_records, and the ERA5 values of that hour (ws_hub_era5, wd100, stability).
    """
    sub = df[df[f"keep_{tag}"] == 1]
    g = sub.groupby("era5_time")
    out = pd.DataFrame({
        "ws": g[f"ws_{tag}"].mean(),
        "nacelle_pos": g[f"nacelle_pos_{tag}"].agg(circular_mean),
        "n_records": g[f"ws_{tag}"].count(),
        "ws_hub_era5": g["ws_hub_era5"].first(),
        "wd100": g["wd100"].first(),
        "stability": g["stability"].first(),
    })
    return out[out["n_records"] >= min_records]


def agreement_stats(x: pd.Series, y: pd.Series) -> dict:
    """Agreement of y (e.g. SCADA) with x (e.g. ERA5): N, bias = mean(y - x),
    RMSE, Pearson correlation, mean ratio y/x. Units as the inputs."""
    ok = x.notna() & y.notna()
    x, y = x[ok], y[ok]
    if len(x) < 2:
        return {"n": int(len(x)), "bias": np.nan, "rmse": np.nan, "corr": np.nan, "ratio_mean": np.nan}
    d = y - x
    return {
        "n": int(len(x)),
        "bias": float(d.mean()),
        "rmse": float(np.sqrt((d ** 2).mean())),
        "corr": float(np.corrcoef(x, y)[0, 1]),
        "ratio_mean": float((y / x).mean()),
    }


def direction_offset(scada_dir: pd.Series, ref_dir: pd.Series) -> dict:
    """Circular difference scada - ref [deg]: N, median, 25 % and 75 % quantiles.

    Before northing: a ROUGH northing-offset estimate (turbine north vs ERA5
    true north). After northing: the remaining offset, which should be small.
    ERA5 is coarse and at 100 m, and terrain turns the flow, so this is an
    indication, not a correction value.
    """
    d = angle_diff(scada_dir, ref_dir).dropna()
    return {"n": int(len(d)), "median": float(d.median()),
            "q25": float(d.quantile(0.25)), "q75": float(d.quantile(0.75))}


def ratio_by_sector(hourly: pd.DataFrame, width: float = 30.0, min_n: int = 1) -> pd.Series:
    """Median SCADA / ERA5-hub speed ratio per ERA5 direction sector (sector centres)."""
    sector = np.floor(np.mod(hourly["wd100"] + width / 2.0, 360.0) / width) * width % 360.0
    ratio = hourly["ws"] / hourly["ws_hub_era5"]
    g = ratio.groupby(sector)
    med = g.median()
    return med[g.count() >= min_n].reindex(np.arange(0, 360, width))


# ----------------------------------------------------------------------------
# Wake statistics
# ----------------------------------------------------------------------------
def add_wake_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Add ws_ratio = ws_B / ws_A, p_ratio = power_B / power_A and the deficits
    ws_deficit = 1 - ws_ratio, p_deficit = 1 - p_ratio (all [-])."""
    out = df.copy()
    out["ws_ratio"] = out["ws_B"] / out["ws_A"]
    out["p_ratio"] = out["power_B"] / out["power_A"]
    out["ws_deficit"] = 1.0 - out["ws_ratio"]
    out["p_deficit"] = 1.0 - out["p_ratio"]
    return out


def _group_stats(g: pd.DataFrame) -> pd.Series:
    """Statistics of one group of paired 10-min records."""
    n = len(g)
    if n == 0:
        return pd.Series({"n": 0})
    fit = np.polyfit(g["ws_A"], g["ws_B"], 1) if n >= 3 else (np.nan, np.nan)
    return pd.Series({
        "n": n,
        "ws_A_mean": g["ws_A"].mean(),                       # [m/s]
        "ws_B_mean": g["ws_B"].mean(),                       # [m/s]
        "ws_ratio_mean": g["ws_ratio"].mean(),               # mean of 10-min ratios [-]
        "ws_deficit": 1.0 - g["ws_B"].mean() / g["ws_A"].mean(),   # from mean speeds [-]
        "ws_deficit_std": g["ws_deficit"].std(),
        "ws_fit_slope": fit[0], "ws_fit_intercept": fit[1],  # ws_B = slope * ws_A + intercept
        "power_A_mean": g["power_A"].mean(),                 # [kW]
        "power_B_mean": g["power_B"].mean(),                 # [kW]
        "p_ratio_mean": g["p_ratio"].mean(),                 # mean of 10-min ratios [-]
        "p_deficit": 1.0 - g["power_B"].sum() / g["power_A"].sum(),  # energy-weighted [-]
        "p_deficit_std": g["p_deficit"].std(),
    })


def pair_stats(df: pd.DataFrame, by: str = "stability") -> pd.DataFrame:
    """Wake statistics overall ("all") and per group (default: stability class).

    ws_deficit = 1 - mean(ws_B) / mean(ws_A); p_deficit = 1 - sum(P_B) / sum(P_A)
    (energy-weighted, like a wake loss). The *_ratio_mean columns are plain
    means of the 10-min ratios; *_std are standard deviations of the 10-min deficits.
    """
    rows = {"all": _group_stats(df)}
    for key, g in df.groupby(by, observed=False):
        rows[str(key)] = _group_stats(g)
    out = pd.DataFrame(rows).T
    out.index.name = by
    out["n"] = out["n"].fillna(0).astype(int)
    return out


def binned_deficit(df: pd.DataFrame, x: str, edges: np.ndarray, min_n: int,
                   by: str = "stability") -> pd.DataFrame:
    """Deficits per bin of column x and per group.

    Returns a long table: group, bin centre, n, ws_deficit (from mean speeds),
    p_deficit (energy-weighted), and their standard errors (std / sqrt(n),
    of the 10-min deficits). Bins with n < min_n get NaN values (n is kept).
    """
    centres = 0.5 * (edges[:-1] + edges[1:])
    cut = pd.cut(df[x], edges, labels=centres, include_lowest=True)
    rows = []
    for (grp, c), g in df.groupby([by, cut], observed=False):
        n = len(g)
        row = {by: grp, "bin": float(c), "n": n,
               "ws_deficit": np.nan, "p_deficit": np.nan, "ws_sem": np.nan, "p_sem": np.nan}
        if n >= max(min_n, 1):
            row.update(
                ws_deficit=1.0 - g["ws_B"].mean() / g["ws_A"].mean(),
                p_deficit=1.0 - g["power_B"].sum() / g["power_A"].sum(),
                ws_sem=g["ws_deficit"].std() / np.sqrt(n) if n > 1 else np.nan,
                p_sem=g["p_deficit"].std() / np.sqrt(n) if n > 1 else np.nan,
            )
        rows.append(row)
    return pd.DataFrame(rows)
