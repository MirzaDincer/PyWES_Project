# Data

**_Claude Generated Text, might contain wrong information. Will be updated later_**

Raw and processed data are **not tracked by git** (too large, and owned by third parties).
Download them into `data/raw/` and run the processing scripts to create `data/processed/`.

| Folder | Content | Tracked by git? |
|---|---|---|
| `raw/` | original downloads (SCADA CSVs, ERA5 NetCDF) | no |
| `processed/` | cleaned and merged data | no |
| `sample/` | small example files (e.g. turbine coordinates, a few days of data) | yes |

## Sources

### SCADA data: Hill of Towie wind farm open dataset

| Item | Details |
|---|---|
| Site | Hill of Towie, Scotland — 21 × Siemens SWT-2.3-VS-82 |
| Period | January 2016 – April 2026 |
| Resolution | 10-minute SCADA statistics, timestamps in UTC (end of period) |
| Released by | RES, on behalf of the owner TRIG |
| License | CC-BY-4.0 |
| Link | https://zenodo.org/records/14870023 (always redirects to the latest version) |
| Size | about 30 GB in total (yearly zip files) — download only the years we need |

Useful files in the record:
- Yearly SCADA zips (10-min statistics and alarm logs)
- Turbine metadata (coordinates, rated power, rotor diameter, hub height) and the data dictionary
- `ShutdownDuration.zip` — downtime per turbine (for filtering)
- AeroUp install dates — turbine upgrades that change power performance (filter or mark these periods)
- Jan–Apr 2026 only: turbine fastlog with Dynamic Yaw controller data (toggled on/off in cycles) and LiDAR data — possible extension, large files

Example download and loading code (Python): https://gh.hhj520.top/resgroup/hill-of-towie-open-source-analysis
— includes a yaw "northing" correction, which is useful because SCADA yaw/wind direction can be biased.

### Weather / stability data

| Dataset | Use | Link |
|---|---|---|
| ERA5 hourly data on single levels (Copernicus C3S) | heat flux, friction velocity, temperature, winds → Obukhov length / stability classes | https://cds.climate.copernicus.eu |

## Download steps

1. Download the turbine metadata and the SCADA years we need from Zenodo → `data/raw/hill_of_towie/`
2. Create a CDS account, save the API key in `~/.cdsapirc`.
3. [Run `scripts/<download_era5>.py` → `data/raw/era5/`]

## Attribution (required)

- **SCADA:** Hill of Towie wind farm open dataset, released by RES on behalf of TRIG, licensed under CC-BY-4.0. — add the dataset citation and DOI of the version we use from the Zenodo page: [citation + DOI]
- **ERA5:** "Generated using Copernicus Climate Change Service information [year]." — add the ERA5 dataset citation from the CDS page.
