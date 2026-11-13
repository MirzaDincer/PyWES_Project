"""ZX300 ground LiDAR (unit 2428) at Hill of Towie: load, clean, shear and veer.

Facts used here (from the dataset's data dictionary and the ZX300 deployment
report ZXMS-RES029-REP01-01, both in data/sample/):
* 10-min files: ``lidar_data/timeseries/2428/Wind10_2428@Y..._M.._D...parquet``
* time column: ``Time and Date``; every file has a large ``Checksum:...`` column
* values 9987 ... 9999 are bad-data sentinels
* heights 17, 28, 38, 48, 58, 68, 78, 88, 99, 181, 263 m above the optical
  head, which is ~1 m above ground
* deployment location 57.514311 N, -3.079795 E (WGS84), ~280 m north of T11
* oriented to MAGNETIC north (no offset) -> directions are in instrument north

Reading parquet needs ``pyarrow`` (not in requirements.txt yet).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

ZX300_LAT = 57.514311
ZX300_LON = -3.079795
HEAD_ABOVE_GROUND = 1.0                                   # [m]
ZX300_HEIGHTS = [17, 28, 38, 48, 58, 68, 78, 88, 99, 181, 263]   # [m] above optical head
SENTINEL_MIN, SENTINEL_MAX = 9987, 9999


def load_zx300_10min(folder: Path, start: str | None = None,
                     end: str | None = None) -> pd.DataFrame:
    """Load all ZX300 10-min files into one DataFrame indexed by time.

    Bad-data sentinels (9987-9999) are set to NaN. The checksum column is dropped.
    start / end: optional date strings, e.g. "2025-09-01".
    """
    files = sorted(Path(folder).glob("Wind10_*.parquet"))
    if not files:
        raise FileNotFoundError(f"No Wind10_*.parquet files in {folder}")

    frames = []
    for f in files:
        df = pd.read_parquet(f)
        df = df.drop(columns=[c for c in df.columns if str(c).startswith("Checksum")])
        frames.append(df)
    df = pd.concat(frames, ignore_index=True)

    df["Time and Date"] = pd.to_datetime(df["Time and Date"])
    df = (df.drop_duplicates(subset="Time and Date")
            .set_index("Time and Date").sort_index())
    df.index.name = "time"

    num = df.select_dtypes("number").columns
    df[num] = df[num].mask((df[num] >= SENTINEL_MIN) & (df[num] <= SENTINEL_MAX))

    if start is not None:
        df = df[df.index >= pd.Timestamp(start)]
    if end is not None:
        df = df[df.index < pd.Timestamp(end)]
    return df


def clean_zx300(df: pd.DataFrame, heights: list[int],
                min_packet_fraction: float = 0.5) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Remove records that are not fully operational or have too few packets.

    A record is kept only if
    * ``Status Flags`` is "Fully Operational", and
    * at every height in ``heights``, ``Packets in Average`` >= fraction x
      the median packet count at that height.

    Returns (cleaned DataFrame, report of removed rows).
    """
    n0 = len(df)
    steps = []
    mask = pd.Series(True, index=df.index)

    if "Status Flags" in df.columns:
        ok = df["Status Flags"].astype(str).str.strip() == "Fully Operational"
        steps.append(("Status Flags not 'Fully Operational'", int((mask & ~ok).sum())))
        mask &= ok

    for h in heights:
        col = f"Packets in Average at {h}m"
        limit = min_packet_fraction * df[col].median()
        ok = df[col] >= limit
        steps.append((f"packets at {h} m < {limit:.0f}", int((mask & ~ok).sum())))
        mask &= ok

    report = pd.DataFrame(steps, columns=["filter", "removed"])
    report["removed_pct"] = 100.0 * report["removed"] / n0
    report.loc[len(report)] = ["kept", int(mask.sum()), 100.0 * mask.sum() / n0]
    return df[mask].copy(), report


def shear_exponent(df: pd.DataFrame, heights: list[int]) -> pd.Series:
    """Power-law shear exponent alpha from a least-squares fit of
    ln(U) = alpha * ln(z) + c over the given heights (Tuna & Bingöl 2020, Eq. 6).

    Heights are converted to metres above ground (+1 m optical head).
    Rows with any missing or non-positive speed give NaN.
    """
    z = np.log(np.array(heights, dtype=float) + HEAD_ABOVE_GROUND)
    U = df[[f"Horizontal Wind Speed (m/s) at {h}m" for h in heights]].to_numpy(dtype=float)
    U = np.where(U > 0, U, np.nan)
    lnU = np.log(U)

    # Slope of a straight-line fit, done for all rows at once
    zc = z - z.mean()
    alpha = (lnU - lnU.mean(axis=1, keepdims=True)) @ zc / np.sum(zc ** 2)
    alpha[np.isnan(lnU).any(axis=1)] = np.nan
    return pd.Series(alpha, index=df.index, name="alpha")


def wind_veer(df: pd.DataFrame, h_low: int, h_high: int) -> pd.Series:
    """Direction change between two heights [deg], wrapped to [-180, 180).

    Positive = direction turns clockwise with height (veering), which is the
    usual Northern Hemisphere case, strongest in stable conditions.
    """
    d_low = df[f"Wind Direction (deg) at {h_low}m"]
    d_high = df[f"Wind Direction (deg) at {h_high}m"]
    return ((d_high - d_low + 180.0) % 360.0 - 180.0).rename("veer")


def iec_disturbed_sectors(lat0: float, lon0: float, turbines: pd.DataFrame,
                          rotor_d: float, max_dist_d: float = 20.0) -> pd.DataFrame:
    """Wind-direction sectors where a turbine wake can reach the instrument.

    Uses the disturbed-sector width as commonly quoted from IEC 61400-12-1
    (Annex A):
        width = 1.3 * arctan(2.5 * D / L + 0.15) + 10   [deg]
    centred on the bearing from the instrument to the turbine.
    NOTE: not checked against the standard text (paywalled) - verify before
    citing it in the report.
    Only turbines closer than max_dist_d rotor diameters are used.

    turbines: DataFrame with columns 'Turbine Name', 'Latitude', 'Longitude'.
    Returns one row per turbine: name, distance_m, bearing_deg, width_deg.
    """
    r_earth = 6_371_000.0
    dlat = np.radians(turbines["Latitude"] - lat0)
    dlon = np.radians(turbines["Longitude"] - lon0) * np.cos(np.radians(lat0))
    north = dlat * r_earth
    east = dlon * r_earth
    dist = np.hypot(north, east)
    bearing = np.mod(np.degrees(np.arctan2(east, north)), 360.0)
    width = 1.3 * np.degrees(np.arctan(2.5 * rotor_d / dist + 0.15)) + 10.0

    out = pd.DataFrame({"turbine": turbines["Turbine Name"].values,
                        "distance_m": dist.values,
                        "bearing_deg": bearing.values,
                        "width_deg": width.values})
    out = out[out["distance_m"] < max_dist_d * rotor_d]
    return out.sort_values("distance_m").reset_index(drop=True)


def in_disturbed_sector(direction: pd.Series, sectors: pd.DataFrame) -> pd.Series:
    """True where the wind direction [deg] falls inside any disturbed sector."""
    flag = pd.Series(False, index=direction.index)
    for _, s in sectors.iterrows():
        diff = (direction - s["bearing_deg"] + 180.0) % 360.0 - 180.0
        flag |= diff.abs() <= s["width_deg"] / 2.0
    return flag
