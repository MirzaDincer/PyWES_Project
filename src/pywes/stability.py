"""Atmospheric stability from ERA5: Obukhov length, classes and sectors.

Method
------
* Inverse Obukhov length 1/L from ERA5 surface fluxes, following the ECMWF
  guide "ERA5: How to calculate Obukhov Length" (Stull 1988, eq. 5.7c),
  including the moisture (virtual heat) flux.
* Stability classes use the neutral limit of Tuna, Bingöl & Sofuoğlu (2026),
  J. Wind Eng. Ind. Aerodyn. 272, 106396: neutral if |1/L| < 0.002 1/m.
* Wind direction sectors of 30 deg (N = 345-15 deg, ...), as in the same paper.

Sign convention: ERA5 vertical fluxes are positive DOWNWARD. Inside this
module we flip them, so ``H_up`` > 0 means the surface heats the air
(daytime, unstable) and 1/L < 0 means unstable.
"""

from __future__ import annotations

import tempfile
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

# ----------------------------------------------------------------------------
# Physical constants (same values as the ECMWF Obukhov-length script)
# ----------------------------------------------------------------------------
KAPPA = 0.4          # von Karman constant [-]
G = 9.81             # gravity [m/s2]
RD = 287.06          # gas constant of dry air [J/(kg K)]
RETV = 0.6078        # Tv = T * (1 + RETV * q) [-]
CP = 1004.7          # specific heat of air at constant pressure [J/(kg K)]

# Neutral limit for |1/L| [1/m] (Tuna et al. 2026, Section 3.2)
NEUTRAL_INV_L = 0.002

CLASS_ORDER = ["unstable", "neutral", "stable"]
SEASON_ORDER = ["DJF", "MAM", "JJA", "SON"]


# ----------------------------------------------------------------------------
# 1. Load ERA5
# ----------------------------------------------------------------------------
def _open_one_file(path: Path) -> xr.Dataset:
    """Open one downloaded ERA5 file.

    The new CDS sometimes returns a zip with several netCDF files inside
    (for example one per "step type"), even if the file ends in .nc.
    We handle both cases and return one merged Dataset.
    """
    if zipfile.is_zipfile(path):
        parts = []
        with zipfile.ZipFile(path) as zf, tempfile.TemporaryDirectory() as tmp:
            for name in zf.namelist():
                if name.endswith(".nc"):
                    zf.extract(name, tmp)
                    parts.append(xr.load_dataset(Path(tmp) / name))
        ds = xr.merge(parts, compat="override")
    else:
        ds = xr.load_dataset(path)

    # New CDS netCDF: the time axis is called "valid_time"
    if "valid_time" in ds.dims or "valid_time" in ds.coords:
        ds = ds.rename({"valid_time": "time"})
    # Drop helper coordinates we don't need
    ds = ds.drop_vars([v for v in ("number", "expver") if v in ds.variables])
    return ds


def load_era5_at_site(era5_dir: Path, lat: float, lon: float) -> pd.DataFrame:
    """Load all monthly ERA5 files and interpolate them to one point.

    Parameters
    ----------
    era5_dir : folder with ``era5_hot_YYYY_MM.nc`` files
    lat, lon : site position [deg] (WGS84)

    Returns
    -------
    DataFrame indexed by UTC time (hourly). Columns use the ERA5 short names:
    sp [Pa], t2m [K], d2m [K], ishf [W/m2], ie [kg/(m2 s)], iews, inss [N/m2],
    zust [m/s], u10, v10, u100, v100 [m/s], blh [m], sdfor [m].
    """
    files = sorted(Path(era5_dir).glob("era5_hot_*.nc"))
    if not files:
        raise FileNotFoundError(
            f"No ERA5 files in {era5_dir}. Run scripts/download_era5.py first."
        )

    datasets = [_open_one_file(f) for f in files]
    ds = xr.concat(datasets, dim="time").sortby("time")
    # Remove duplicate hours (can happen at month borders)
    _, keep = np.unique(ds["time"].values, return_index=True)
    ds = ds.isel(time=np.sort(keep))

    # Bilinear interpolation from the 3 x 3 grid to the site
    point = ds.interp(latitude=lat, longitude=lon, method="linear")

    df = point.to_dataframe()
    df = df.drop(columns=[c for c in ("latitude", "longitude") if c in df.columns])
    df.index = pd.DatetimeIndex(df.index, name="time_utc")
    return df


# ----------------------------------------------------------------------------
# 2. Obukhov length (ECMWF method)
# ----------------------------------------------------------------------------
def saturation_specific_humidity(t: np.ndarray, p: np.ndarray) -> np.ndarray:
    """Saturation specific humidity over water [kg/kg].

    Same formula as ``qswat`` in the ECMWF script. With t = 2 m DEWPOINT, this
    gives the actual 2 m specific humidity.

    t : temperature [K], p : pressure [Pa]
    """
    rkbol = 1.380658e-23
    rnavo = 6.0221367e23
    r = rnavo * rkbol
    rd = 1000.0 * r / 28.9644      # gas constant of dry air
    rv = 1000.0 * r / 18.0153      # gas constant of water vapour
    restt = 611.21
    r2es = restt * rd / rv
    r3les = 17.502
    r4les = 32.19
    retv = rv / rd - 1.0
    rtt = 273.16

    foeew = r2es * np.exp((r3les * (t - rtt)) / (t - r4les))
    qs = foeew / p
    zcor = 1.0 / (1.0 - retv * qs)
    return qs * zcor


def add_obukhov_length(df: pd.DataFrame) -> pd.DataFrame:
    """Add friction velocity, heat fluxes and 1/L to an ERA5 DataFrame.

    New columns
    -----------
    q2      2 m specific humidity [kg/kg]
    tv2     2 m virtual temperature [K]
    rho     air density [kg/m3]
    ust     friction velocity from the surface stress [m/s]
    H_up    sensible heat flux, positive UPWARD [W/m2]
    wtv     kinematic virtual heat flux, positive upward [K m/s]
    inv_L   inverse Obukhov length [1/m]  (< 0 unstable, > 0 stable)
    L       Obukhov length [m] (inf when inv_L = 0)
    """
    out = df.copy()

    out["q2"] = saturation_specific_humidity(out["d2m"].to_numpy(), out["sp"].to_numpy())
    out["tv2"] = out["t2m"] * (1.0 + RETV * out["q2"])
    out["rho"] = out["sp"] / (RD * out["tv2"])

    # Friction velocity from the stress magnitude, with a small floor
    tau = np.hypot(out["iews"], out["inss"])
    out["ust"] = np.maximum(np.sqrt(tau / out["rho"]), 0.001)

    # ERA5 fluxes are positive DOWNWARD -> minus sign gives positive upward
    out["H_up"] = -out["ishf"]                     # [W/m2]
    wt = -out["ishf"] / (out["rho"] * CP)          # kinematic heat flux [K m/s]
    wq = -out["ie"] / out["rho"]                   # kinematic moisture flux
    out["wtv"] = wt + RETV * out["t2m"] * wq       # virtual heat flux [K m/s]

    # 1/L = -kappa * g * wtv / (Tv * u*^3)
    out["inv_L"] = -KAPPA * G * out["wtv"] / (out["tv2"] * out["ust"] ** 3)
    with np.errstate(divide="ignore"):
        out["L"] = 1.0 / out["inv_L"]
    return out


# ----------------------------------------------------------------------------
# 3. Wind, classes, sectors, season, day/night
# ----------------------------------------------------------------------------
def wind_speed_direction(u: pd.Series, v: pd.Series) -> tuple[pd.Series, pd.Series]:
    """Wind speed [m/s] and meteorological direction [deg, FROM north].

    u = eastward, v = northward component. A wind blowing from the west
    (u > 0, v = 0) gives 270 deg.
    """
    speed = np.hypot(u, v)
    direction = np.mod(np.degrees(np.arctan2(-u, -v)), 360.0)
    return speed, direction


def classify_stability(inv_L: pd.Series, neutral_limit: float = NEUTRAL_INV_L) -> pd.Series:
    """Three stability classes from 1/L [1/m].

    unstable: 1/L < -limit,  neutral: |1/L| <= limit,  stable: 1/L > limit.
    Default limit 0.002 1/m (|L| > 500 m is neutral), Tuna et al. (2026).
    """
    cls = np.select(
        [inv_L < -neutral_limit, inv_L > neutral_limit],
        ["unstable", "stable"],
        default="neutral",
    )
    cls = pd.Series(cls, index=inv_L.index)
    cls[inv_L.isna()] = np.nan
    return pd.Categorical(cls, categories=CLASS_ORDER, ordered=True)


def assign_sector(direction: pd.Series, width: float = 30.0) -> pd.Series:
    """Sector CENTRE [deg] for each direction; first sector centred on north.

    width = 30 -> centres 0, 30, ..., 330 and sector 0 covers 345-15 deg
    (same as Tuna et al. 2026, Table 2).
    """
    n_sectors = int(round(360.0 / width))
    idx = np.floor(np.mod(direction + width / 2.0, 360.0) / width).astype("Int64") % n_sectors
    return idx * width


def sector_labels(width: float = 30.0) -> dict[float, str]:
    """Readable labels for sector centres, e.g. {0: 'N', 30: 'NNE', ...} for 30 deg."""
    names16 = ["N", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE",
               "S", "SSW", "SW", "WSW", "W", "WNW", "NW", "NNW"]
    labels = {}
    for c in np.arange(0, 360, width):
        k = int(round(c / 22.5)) % 16
        labels[float(c)] = f"{names16[k]} ({int(c)}°)"
    return labels


def season_of(index: pd.DatetimeIndex) -> pd.Categorical:
    """Meteorological season: DJF, MAM, JJA, SON."""
    month_to_season = {12: "DJF", 1: "DJF", 2: "DJF", 3: "MAM", 4: "MAM", 5: "MAM",
                       6: "JJA", 7: "JJA", 8: "JJA", 9: "SON", 10: "SON", 11: "SON"}
    seasons = [month_to_season[m] for m in index.month]
    return pd.Categorical(seasons, categories=SEASON_ORDER, ordered=True)


def solar_elevation(index_utc: pd.DatetimeIndex, lat: float, lon: float) -> np.ndarray:
    """Approximate solar elevation angle [deg] (accuracy ~1 deg, enough for day/night).

    Uses the standard declination / equation-of-time approximation
    (NOAA "General Solar Position Calculations").
    """
    doy = index_utc.dayofyear.to_numpy()
    hour = (index_utc.hour + index_utc.minute / 60.0).to_numpy()
    gamma = 2.0 * np.pi / 365.0 * (doy - 1 + (hour - 12.0) / 24.0)   # fractional year [rad]

    eqtime = 229.18 * (0.000075 + 0.001868 * np.cos(gamma) - 0.032077 * np.sin(gamma)
                       - 0.014615 * np.cos(2 * gamma) - 0.040849 * np.sin(2 * gamma))  # [min]
    decl = (0.006918 - 0.399912 * np.cos(gamma) + 0.070257 * np.sin(gamma)
            - 0.006758 * np.cos(2 * gamma) + 0.000907 * np.sin(2 * gamma)
            - 0.002697 * np.cos(3 * gamma) + 0.00148 * np.sin(3 * gamma))         # [rad]

    true_solar_time = hour * 60.0 + eqtime + 4.0 * lon          # [min]
    hour_angle = np.radians(true_solar_time / 4.0 - 180.0)       # [rad]
    phi = np.radians(lat)
    cos_zenith = (np.sin(phi) * np.sin(decl)
                  + np.cos(phi) * np.cos(decl) * np.cos(hour_angle))
    return 90.0 - np.degrees(np.arccos(np.clip(cos_zenith, -1.0, 1.0)))


def nearest_era5_hour(timestamps: pd.Series | pd.DatetimeIndex,
                      period: str = "10min",
                      timestamp_is_end: bool = True) -> pd.DatetimeIndex:
    """Match 10-min records to the nearest ERA5 hour.

    The 10-min record is first moved to its MIDDLE (e.g. end-labelled 10:10
    -> 10:05) and then rounded to the nearest full hour. Use this later to join
    the stability table to SCADA (SCADA TimeStamp = end of period).
    """
    ts = pd.DatetimeIndex(timestamps)
    half = pd.Timedelta(period) / 2
    mid = ts - half if timestamp_is_end else ts + half
    return mid.round("h")


def build_stability_table(era5: pd.DataFrame, lat: float, lon: float,
                          hub_height: float, sector_width: float = 30.0) -> pd.DataFrame:
    """Full hourly stability table from the site ERA5 DataFrame.

    Adds: inv_L, L, zeta_hub (= hub_height / L), zeta_10 (= 10 m / L),
    ws100, wd100, ws10, wd10, stability class, sector, season, solar elevation
    and an is_day flag.
    """
    tab = add_obukhov_length(era5)
    tab["zeta_hub"] = hub_height * tab["inv_L"]
    tab["zeta_10"] = 10.0 * tab["inv_L"]
    tab["ws100"], tab["wd100"] = wind_speed_direction(tab["u100"], tab["v100"])
    tab["ws10"], tab["wd10"] = wind_speed_direction(tab["u10"], tab["v10"])
    tab["stability"] = classify_stability(tab["inv_L"])
    tab["sector"] = assign_sector(tab["wd100"], sector_width)
    tab["season"] = season_of(tab.index)
    tab["sun_elev"] = solar_elevation(tab.index, lat, lon)
    tab["is_day"] = tab["sun_elev"] > 0.0
    return tab


def apply_filters(tab: pd.DataFrame, min_ws: float = 3.0,
                  max_abs_zeta10: float = 1.0) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Filters from Tuna et al. (2026), Section 4, and a report of what was removed.

    1. remove missing 1/L or wind
    2. remove low wind: ws100 < min_ws [m/s]
    3. remove extreme stability: |zeta_10| = |10 m / L| > max_abs_zeta10

    Returns (filtered table, report table with the number of rows removed).
    """
    steps = []
    n0 = len(tab)
    ok_missing = tab[["inv_L", "ws100", "wd100"]].notna().all(axis=1)
    ok_wind = tab["ws100"] >= min_ws
    ok_extreme = tab["zeta_10"].abs() <= max_abs_zeta10

    mask = ok_missing.copy()
    steps.append(("missing 1/L or wind", n0 - int(mask.sum())))
    before = int(mask.sum())
    mask &= ok_wind
    steps.append((f"ws100 < {min_ws} m/s", before - int(mask.sum())))
    before = int(mask.sum())
    mask &= ok_extreme
    steps.append((f"|10 m / L| > {max_abs_zeta10}", before - int(mask.sum())))

    report = pd.DataFrame(steps, columns=["filter", "removed"])
    report["removed_pct"] = 100.0 * report["removed"] / n0
    report.loc[len(report)] = ["kept", int(mask.sum()), 100.0 * mask.sum() / n0]
    return tab[mask].copy(), report


def class_counts(tab: pd.DataFrame, by: str | list[str]) -> pd.DataFrame:
    """Counts and percentages of each stability class per group.

    Returns a table with columns n_unstable, n_neutral, n_stable, n_total,
    pct_unstable, pct_neutral, pct_stable.
    """
    counts = (tab.groupby(by, observed=False)["stability"]
              .value_counts().unstack(fill_value=0)
              .reindex(columns=CLASS_ORDER, fill_value=0))
    out = counts.add_prefix("n_")
    out["n_total"] = counts.sum(axis=1)
    pct = counts.div(counts.sum(axis=1).replace(0, np.nan), axis=0) * 100.0
    for c in CLASS_ORDER:
        out[f"pct_{c}"] = pct[c]
    return out
