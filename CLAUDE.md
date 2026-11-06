# CLAUDE.md

## Project Context
- **Name and purpose:** PyWES_Project — team term paper for *Python for Wind Energy Systems*, Hochschule Flensburg, WiSe 2026/27. Working title: "Stability-aware wake modelling and wake steering with Python". We measure how atmospheric stability and wind shear change wake losses in open SCADA data, compare one global wake-model setting against one setting per stability class, and simulate wake steering per stability class. It is a course project, not a publication: reusing published methods is fine when cited.
- **Tech stack:** Python ≥ 3.10, NumPy, Pandas, SciPy, xarray/netCDF4, Matplotlib, FLORIS v4 (main wake model and yaw optimization), cdsapi (ERA5 download). IDE: Spyder. Optional: PyWake, OpenOA, scikit-learn. No database. Full list in `requirements.txt`.
- **Repository structure:**
  - `src/pywes/` — reusable project code (installed as package `pywes`); paths in `src/pywes/paths.py`
  - `scripts/` — download and processing scripts
  - `analysis/` — analysis scripts with Spyder cells (`# %%`)
  - `data/raw/`, `data/processed/` — local only, git-ignored; `data/sample/` — small tracked files; sources and attribution in `data/README.md`
  - `results/figures/` — final figures for the paper (`draft/` is ignored)
  - `docs/` — notes, course profile, paper drafts

## Model chain (MVP — keep this in focus)
SCADA + ERA5 → cleaning + stability classes (Obukhov length: stable/neutral/unstable) → observed wake losses per class (turbine pairs, binned by direction and wind speed) → FLORIS with global vs. per-class tuning → comparison plots + error table.
Extras (only after the MVP works): stability-dependent yaw lookup tables, wind direction uncertainty, forecasting.

## Commands (Windows first; run from repo root)
```powershell
python -m venv .venv
.venv\Scripts\activate            # macOS/Linux: source .venv/bin/activate
pip install -r requirements.txt
pip install -e .                  # makes `import pywes` work everywhere
python scripts/<script_name>.py   # run a pipeline step
python -c "import floris, pywes; print(floris.__version__)"   # quick environment check
```
- **Tests:** no test suite yet. When `tests/` exists, run `pytest`. Check new functions on simple cases with known answers (e.g. a single turbine, or zero yaw).
- **Linting/formatting:** no linter is configured. Follow the style of the existing code; do not add tooling without asking.

## Coding conventions
- Use paths from `pywes.paths` (`RAW_DIR`, `PROCESSED_DIR`, `FIGURES_DIR`); never hard-code absolute paths like `D:\...` and don't rely on the current working directory (Spyder may change it).
- Use FLORIS v4 API (`FlorisModel`, `fmodel.set(...)`, `fmodel.run()`). Check the current docs before using an API you are unsure of: https://natlabrockies.github.io/floris/
- Put reusable logic in `src/pywes/`; keep analysis scripts thin (load → call functions → plot).
- Write analysis code as `.py` files split into `# %%` cells for Spyder, not `.ipynb` notebooks.
- Write short docstrings with units (m/s, deg, kW). Keep wind directions in meteorological convention (deg, from north).
- Raw data is read-only: write cleaned output to `data/processed/`.
- Comments and docs in simple, clear English.

## Boundaries
**Never:**
- Commit or stage anything from `data/raw/`, `data/processed/`, large data files (`.nc`, `.grib`, `.parquet`, `.zip`) or credentials (`.env`, `.cdsapirc`, API keys). See `.gitignore`.
- Invent references, dataset details, numbers or results. If something can't be verified, say so.
- Modify or delete raw data files.
- Hide failing steps or silently drop data; report filters and how many samples they removed.

**Ask before:**
- Adding, removing or upgrading dependencies (`requirements.txt`, `pyproject.toml`).
- Starting large downloads or ERA5 requests (time and CDS quota).
- Committing, pushing, creating branches, or any history rewrite (`rebase`, `reset --hard`, force push).
- Changing `LICENSE`, `README.md` structure, data attribution in `data/README.md`, or the folder layout.
- Expanding scope beyond the MVP (new models, ML forecasting, extra wind farms).
- Deleting or overwriting files in `results/figures/` or `docs/`.

## Physics sanity checks
- Wakes should be stronger and longer in stable conditions and recover faster in unstable conditions.
- Wake losses must be ≥ 0 for the downstream turbine in a pair aligned with the wind direction; check sign conventions if not.
- Wake steering gains are small (typically a few percent) and very sensitive to wind direction error; be suspicious of large gains.

## Planned later
- The Spyder installation/setup steps (choosing the `.venv` interpreter, installing `spyder-kernels`) were removed from `README.md` for now. When the project is close to its end (final cleanup before submission), remind Mirza to decide whether to add them back.

## Workflow
- Work on feature branches; merge to `main` through a pull request reviewed by a teammate.
- Don't commit Spyder project settings (`.spyproject/` is git-ignored).
- Cite sources for methods and data in code comments or `docs/`.
