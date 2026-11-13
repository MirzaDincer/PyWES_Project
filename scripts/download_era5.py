"""Download ERA5 hourly single-level data for the Hill of Towie site.

What it does
------------
Sends one CDS request per month and saves one netCDF file per month to
``data/raw/era5/era5_hot_YYYY_MM.nc``. Months that already exist are skipped,
so you can stop and restart the script.

The variables are the ones needed for the Obukhov length (ECMWF guide
"ERA5: How to calculate Obukhov Length") plus winds and boundary-layer height:
    surface pressure, 2 m temperature, 2 m dewpoint temperature,
    instantaneous surface sensible heat flux, instantaneous moisture flux,
    instantaneous eastward / northward turbulent surface stress,
    friction velocity (cross-check only), 10 m and 100 m wind components,
    boundary layer height, std. dev. of filtered sub-grid orography.

Before you run it
-----------------
1. Make a free account at https://cds.climate.copernicus.eu and accept the
   ERA5 licence on the dataset page.
2. Put your API key in ``%USERPROFILE%\\.cdsapirc`` (never commit this file).
3. ``pip install cdsapi`` (already in requirements.txt).

Size: 3 x 3 grid points, 14 variables, ~744 hours per month -> a few MB per
month. Requests can wait in the CDS queue for several minutes each.

Run from the repo root:
    python scripts/download_era5.py
"""

from __future__ import annotations

import calendar

import cdsapi

from pywes.paths import RAW_DIR

# ----------------------------------------------------------------------------
# Settings
# ----------------------------------------------------------------------------
# Period: same as the SCADA files in data/raw (2025-01 ... 2026-04).
START = (2025, 1)   # (year, month)
END = (2026, 4)     # (year, month), inclusive

# Area [North, West, South, East] in degrees. The farm centre is about
# 57.506 N, -3.067 E; this box gives the 3 x 3 ERA5 grid points around it
# (0.25 deg grid), so we can interpolate to the site later.
AREA = [57.75, -3.25, 57.25, -2.75]

VARIABLES = [
    "surface_pressure",
    "2m_temperature",
    "2m_dewpoint_temperature",
    "instantaneous_surface_sensible_heat_flux",
    "instantaneous_moisture_flux",
    "instantaneous_eastward_turbulent_surface_stress",
    "instantaneous_northward_turbulent_surface_stress",
    "friction_velocity",
    "10m_u_component_of_wind",
    "10m_v_component_of_wind",
    "100m_u_component_of_wind",
    "100m_v_component_of_wind",
    "boundary_layer_height",
    "standard_deviation_of_filtered_subgrid_orography",
]

OUT_DIR = RAW_DIR / "era5"


def month_list(start: tuple[int, int], end: tuple[int, int]) -> list[tuple[int, int]]:
    """Return all (year, month) pairs from start to end, inclusive."""
    months = []
    year, month = start
    while (year, month) <= end:
        months.append((year, month))
        month += 1
        if month == 13:
            year, month = year + 1, 1
    return months


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    client = cdsapi.Client()

    for year, month in month_list(START, END):
        target = OUT_DIR / f"era5_hot_{year}_{month:02d}.nc"
        if target.exists():
            print(f"skip  {target.name} (already downloaded)")
            continue

        n_days = calendar.monthrange(year, month)[1]
        request = {
            "product_type": ["reanalysis"],
            "variable": VARIABLES,
            "year": [str(year)],
            "month": [f"{month:02d}"],
            "day": [f"{d:02d}" for d in range(1, n_days + 1)],
            "time": [f"{h:02d}:00" for h in range(24)],
            "area": AREA,
            "data_format": "netcdf",
            "download_format": "unarchived",
        }
        print(f"request {target.name} ...")
        client.retrieve("reanalysis-era5-single-levels", request).download(str(target))
        print(f"saved {target}")


if __name__ == "__main__":
    main()
