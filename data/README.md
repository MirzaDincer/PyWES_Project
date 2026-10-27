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

Final choice of wind farm is still open — keep only the ones we actually use.

### SCADA data (10-minute)

| Dataset | Turbines | License | Link |
|---|---|---|---|
| Kelmarsh wind farm (Cubico) | 6 | CC-BY-4.0 | [Zenodo – link] |
| Penmanshiel wind farm (Cubico) | 14 Senvion MM82 | CC-BY-4.0 | https://zenodo.org/records/8253010 |
| Hill of Towie wind farm | 21 | CC-BY-4.0 | [link] |

Overview of open SCADA datasets: https://github.com/sltzgs/OpenWindSCADA

### Weather / stability data

| Dataset | Use | Link |
|---|---|---|
| ERA5 hourly data on single levels (Copernicus C3S) | heat flux, friction velocity, temperature, winds → Obukhov length / stability classes | https://cds.climate.copernicus.eu |

## Download steps

1. [Step to download SCADA data → `data/raw/<farm>/`]
2. Create a CDS account, save the API key in `~/.cdsapirc`.
3. [Run `scripts/<download_era5>.py` → `data/raw/era5/`]

## Attribution (required)

- **SCADA:** "[Kelmarsh / Penmanshiel] wind farm data, Cubico Sustainable Investments Ltd, licensed under CC-BY-4.0." — add the dataset citation from the Zenodo page.
- **ERA5:** "Generated using Copernicus Climate Change Service information [year]." — add the ERA5 dataset citation from the CDS page.
