# PyWES_Project

**_Claude Generated Text, might contain wrong information. Will be updated later_**

**Stability-aware wake modelling and wake steering with Python**
Term paper for *Python for Wind Energy Systems (PyWES)*, Hochschule Flensburg, WiSe 2026/27.

Team: Mirza Dincer, [Name 2], [Name 3], [Name 4]

## About

Wind turbine wakes reduce the power of downstream turbines. How strong and how long a wake is depends on the weather:
in **stable** conditions (calm nights, winter) wakes stay strong for a long distance, while in **unstable** conditions
(sunny, turbulent days) they recover quickly. Standard engineering wake models usually use one fixed tuning parameter.

In this project we:

1. Measure how wake losses change with atmospheric stability and wind shear in open SCADA data.
2. Test whether a wake model tuned per stability class predicts better than one global setting.
3. Simulate whether wake steering (yaw offsets) gives more benefit in some conditions than in others.
4. *(Optional)* Predict the stability class / expected wake loss a few hours ahead from weather data.

**Model chain (MVP):**
SCADA + ERA5 → cleaning + stability classes → observed wake losses per class → FLORIS (global vs. per-class tuning) → comparison plots + error table.

## Repository structure

```
PyWES_Project/
├── data/
│   ├── raw/          # downloaded SCADA / ERA5 files (not tracked by git)
│   ├── processed/    # cleaned, merged data (not tracked by git)
│   ├── sample/       # small example files that are tracked
│   └── README.md     # data sources, download steps, attribution
├── analysis/         # analysis scripts with Spyder cells (# %%)
├── scripts/          # download and processing scripts
├── src/pywes/        # reusable project code
├── results/figures/  # final figures for the paper
├── docs/             # notes, profile, term paper drafts
├── requirements.txt
└── LICENSE
```

## Setup

```bash
# create and activate a virtual environment
python -m venv .venv
# Windows:
.venv\Scripts\activate
# macOS / Linux:
source .venv/bin/activate

# install dependencies and this project (editable)
pip install -r requirements.txt
pip install -e .
```

### Using Spyder

1. In Spyder: **Tools → Preferences → Python interpreter → Use the following interpreter** and choose
   `.venv\Scripts\python.exe` from this repo. Restart the console.
2. If Spyder says `spyder-kernels` is missing or has the wrong version, install the version it shows
   (inside the activated `.venv`): `pip install "spyder-kernels==<version>"`.
3. Optional: **Projects → New Project → Existing directory** → select this repo folder
   (the `.spyproject/` folder it creates is git-ignored).
4. Split analysis scripts into cells with `# %%` and run them cell by cell (Ctrl+Enter).
   Use the paths from `pywes.paths` so scripts work no matter which working directory Spyder uses.

To download ERA5 data you need a free Copernicus Climate Data Store (CDS) account and an API key.
Keep the key in your home folder (`.cdsapirc`) — **never commit it**.

## Data

The raw data is **not** stored in this repository. See [`data/README.md`](data/README.md) for sources,
download steps and the required attribution.

## Main tools

- [FLORIS](https://github.com/NatLabRockies/floris) — engineering wake models and yaw optimization
- [PyWake](https://topfarm.pages.windenergy.dtu.dk/PyWake/) — optional cross-check of wake models
- [OpenOA](https://openoa.readthedocs.io) — optional, SCADA-based wake loss analysis
- NumPy, Pandas, xarray, Matplotlib

## Workflow notes

- Work on a branch, open a pull request, and let one teammate review before merging to `main`.
- Prefer `.py` scripts with `# %%` cells over Jupyter notebooks — they give clean diffs and run directly in Spyder.
- Save only final figures in `results/figures/`; put drafts in `results/figures/draft/` (ignored).

## License

Code: [MIT License](LICENSE).
Data: belongs to the original providers and is used under their licenses — see [`data/README.md`](data/README.md).
