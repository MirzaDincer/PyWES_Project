# -*- coding: utf-8 -*-
"""Atmospheric stability at Hill of Towie by wind-direction sector.

Steps (run cell by cell in Spyder, or the whole file):
  1. Settings and site data (turbine metadata)
  2. Load ERA5, compute 1/L and stability classes
  3. Filters (Tuna et al. 2026) and save the hourly table
  4. Physics sanity checks
  5. Plots: time series of 1/L, then by sector, season and day/night
  6. LiDAR check: shear, veer and TI per ERA5 stability class (2025-08 ... 2026-04)

Inputs : data/raw/era5/era5_hot_YYYY_MM.nc (from scripts/download_era5.py)
         data/raw/lidar_data/lidar_data/timeseries/2428/Wind10_*.parquet
         data/sample/Hill_of_Towie_turbine_metadata.csv
Outputs: data/processed/era5_stability_hourly.csv, *_counts_*.csv
         results/figures/draft/stability/*.png

References
  Method for 1/L: ECMWF, "ERA5: How to calculate Obukhov Length" (CKB).
  Classes, sectors and filters: Tuna, Bingöl & Sofuoğlu (2026), JWEIA 272, 106396.
  Shear exponent: Tuna & Bingöl (2020), Energies 13, 89, Eq. (6).
"""

# %% 1. Settings and site data ------------------------------------------------
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.colors import LinearSegmentedColormap

from pywes.paths import RAW_DIR, PROCESSED_DIR, SAMPLE_DIR, FIGURES_DIR
from pywes import stability as st

# --- settings you may change --------------------------------------------------
SECTOR_WIDTH = 30.0        # [deg] sector width (Tuna et al. 2026 use 30 deg)
MIN_WS = 3.0               # [m/s] low-wind filter on ERA5 100 m wind
MAX_ABS_ZETA10 = 1.0       # [-] remove |10 m / L| > 1 (Tuna et al. 2026)
MIN_N_CELL = 30            # [-] hide sector/season cells with fewer hours
RUN_LIDAR_PART = True      # set False if pyarrow is not installed yet

ERA5_DIR = RAW_DIR / "era5"
LIDAR_DIR = RAW_DIR / "lidar_data" / "lidar_data" / "timeseries" / "2428"
FIG_DIR = FIGURES_DIR / "draft" / "stability"
FIG_DIR.mkdir(parents=True, exist_ok=True)
PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

# --- plot style (recessive grid and axes, thin marks) -------------------------
CLASS_COLORS = {"unstable": "#e34948", "neutral": "#a3a19a", "stable": "#2a78d6"}
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#898781", "#e1e0d9"
plt.rcParams.update({
    "figure.dpi": 110, "savefig.dpi": 300, "savefig.bbox": "tight",
    "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
    "axes.edgecolor": "#c3c2b7", "axes.labelcolor": INK2, "text.color": INK,
    "xtick.color": MUTED, "ytick.color": MUTED,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.spines.top": False, "axes.spines.right": False,
    "legend.frameon": False,
})


def save(fig, name):
    """Save a figure as PNG in the draft folder and print the path."""
    path = FIG_DIR / f"{name}.png"
    fig.savefig(path)
    print(f"saved {path}")


# --- site data ------------------------------------------------------------------
turbines = pd.read_csv(SAMPLE_DIR / "Hill_of_Towie_turbine_metadata.csv", encoding="utf-8-sig")
SITE_LAT = turbines["Latitude"].mean()          # farm centre [deg]
SITE_LON = turbines["Longitude"].mean()
HUB_HEIGHT = float(turbines["Hub Height (m)"].iloc[0])     # 59 m
ROTOR_D = float(turbines["Rotor Diameter (m)"].iloc[0])    # 82 m
print(f"Farm centre: {SITE_LAT:.4f} N, {SITE_LON:.4f} E | hub {HUB_HEIGHT} m | D {ROTOR_D} m")

SECTOR_NAMES = st.sector_labels(SECTOR_WIDTH)
SECTORS = list(SECTOR_NAMES.keys())


# %% 2. Load ERA5 and compute stability ---------------------------------------
era5 = st.load_era5_at_site(ERA5_DIR, SITE_LAT, SITE_LON)
print(f"ERA5: {era5.index.min()} -> {era5.index.max()}, {len(era5)} hours")

tab_all = st.build_stability_table(era5, SITE_LAT, SITE_LON, HUB_HEIGHT, SECTOR_WIDTH)

# Sub-grid orography: ECMWF notes that 1/L is less reliable over strong
# orography. Report the value so we can mention it in the paper.
if "sdfor" in tab_all:
    print(f"Std. dev. of filtered sub-grid orography at site: {tab_all['sdfor'].mean():.1f} m")

# Cross-check: our u* (from the stress) vs ERA5's own friction velocity
if "zust" in tab_all:
    ratio = (tab_all["ust"] / tab_all["zust"]).median()
    corr = tab_all[["ust", "zust"]].corr().iloc[0, 1]
    print(f"u* check: median(u*_stress / zust) = {ratio:.2f}, correlation = {corr:.2f}")


# %% 3. Filters and save --------------------------------------------------------
tab, filter_report = st.apply_filters(tab_all, min_ws=MIN_WS, max_abs_zeta10=MAX_ABS_ZETA10)
print("\nFilter report (hours):")
print(filter_report.to_string(index=False, float_format="%.1f"))

# Also report how many hours have a very small heat flux (|H| < 10 W/m2).
# Tuna & Bingöl (2020) treat these fluxes as too uncertain; we keep them
# (they fall mostly in the neutral class) but report the share.
small_H = (tab["H_up"].abs() < 10.0).mean() * 100
print(f"Hours with |H| < 10 W/m2 (kept): {small_H:.1f} %")

cols_out = ["inv_L", "L", "zeta_hub", "zeta_10", "ust", "H_up", "wtv",
            "ws10", "wd10", "ws100", "wd100", "blh", "stability", "sector",
            "season", "sun_elev", "is_day"]
cols_out = [c for c in cols_out if c in tab_all.columns]
out = tab_all[cols_out].copy()
out["kept"] = out.index.isin(tab.index)
out.to_csv(PROCESSED_DIR / "era5_stability_hourly.csv")
print(f"saved {PROCESSED_DIR / 'era5_stability_hourly.csv'}")

counts_sector = st.class_counts(tab, "sector")
counts_sector.index = [SECTOR_NAMES[s] for s in counts_sector.index]
counts_sector.to_csv(PROCESSED_DIR / "stability_counts_by_sector.csv", float_format="%.2f")
counts_sector_season = st.class_counts(tab, ["sector", "season"])
counts_sector_season.to_csv(PROCESSED_DIR / "stability_counts_by_sector_season.csv",
                            float_format="%.2f")
print("\nClass shares per sector (%):")
print(counts_sector[["n_total", "pct_unstable", "pct_neutral", "pct_stable"]]
      .to_string(float_format="%.1f"))


# %% 4. Physics sanity checks ---------------------------------------------------
# Expectation: more unstable by day and in summer, more stable at night and in
# winter. If this is NOT true, check the sign of the ERA5 fluxes first.
dn = st.class_counts(tab, "is_day")
dn.index = dn.index.map({True: "day", False: "night"})
print("\nDay vs night (%):")
print(dn[["n_total", "pct_unstable", "pct_neutral", "pct_stable"]].to_string(float_format="%.1f"))

se = st.class_counts(tab, "season")
print("\nSeason (%):")
print(se[["n_total", "pct_unstable", "pct_neutral", "pct_stable"]].to_string(float_format="%.1f"))

ok_day = dn.loc["day", "pct_unstable"] > dn.loc["night", "pct_unstable"]
ok_night = dn.loc["night", "pct_stable"] > dn.loc["day", "pct_stable"]
print(f"\nCHECK unstable more frequent by day:  {'OK' if ok_day else 'NOT OK - check signs!'}")
print(f"CHECK stable more frequent at night:  {'OK' if ok_night else 'NOT OK - check signs!'}")


# %% 5-0. Time series of 1/L -------------------------------------------------------
# Hourly 1/L (dots, coloured by class) and the daily median (line).
# Only hours that passed the filters are shown, so filtered hours appear as gaps.
# What to look for: a daily cycle (unstable by day, stable at night), more
# unstable hours in summer, and no sudden jumps (e.g. at month borders).
Y_LIM = 0.05                                      # [1/m] plot range; larger |1/L| are clipped
fig, ax = plt.subplots(figsize=(10, 3.4))
ax.axhspan(-st.NEUTRAL_INV_L, st.NEUTRAL_INV_L, color=GRID, zorder=0, label="neutral band")
for cls in st.CLASS_ORDER:
    sub = tab.loc[tab["stability"] == cls, "inv_L"]
    ax.scatter(sub.index, sub.clip(-Y_LIM, Y_LIM), s=2, color=CLASS_COLORS[cls],
               linewidths=0, label=cls, zorder=2)
daily = tab["inv_L"].resample("D").median()
ax.plot(daily.index, daily.clip(-Y_LIM, Y_LIM), color=INK, linewidth=1.0,
        label="daily median", zorder=3)
ax.axhline(0, color=MUTED, linewidth=0.6, zorder=1)
ax.set_ylim(-Y_LIM * 1.05, Y_LIM * 1.05)
ax.set_ylabel("1/L [1/m]")
ax.set_title(f"ERA5 inverse Obukhov length, hourly (|1/L| > {Y_LIM} clipped; "
             "< 0 unstable, > 0 stable)")
locator = mdates.AutoDateLocator()
ax.xaxis.set_major_locator(locator)
ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(locator))   # short, unrotated dates
ax.legend(ncol=5, loc="upper center", bbox_to_anchor=(0.5, -0.16), markerscale=4)
save(fig, "00_invL_timeseries")


# %% 5a. Histogram of 1/L (all sectors) -------------------------------------------
fig, ax = plt.subplots(figsize=(6.5, 3.2))
bins = np.linspace(-0.05, 0.05, 101)
for cls in st.CLASS_ORDER:
    sub = tab.loc[tab["stability"] == cls, "inv_L"].clip(bins[0], bins[-1])
    ax.hist(sub, bins=bins, color=CLASS_COLORS[cls], label=f"{cls} (n={len(sub)})",
            edgecolor="white", linewidth=0.3)
for x in (-st.NEUTRAL_INV_L, st.NEUTRAL_INV_L):
    ax.axvline(x, color=INK2, linewidth=0.8, linestyle="--")
ax.set_xlabel("1/L [1/m]  (values outside ±0.05 are put in the end bins)")
ax.set_ylabel("Hours")
ax.set_title("ERA5 inverse Obukhov length at Hill of Towie (after filters)")
ax.legend()
save(fig, "01_invL_histogram")


# %% 5b. Class shares per sector (stacked bars) -----------------------------------
c = st.class_counts(tab, "sector").reindex(SECTORS)
x = np.arange(len(SECTORS))
fig, ax = plt.subplots(figsize=(7.5, 3.6))
bottom = np.zeros(len(SECTORS))
for cls in st.CLASS_ORDER:
    vals = c[f"pct_{cls}"].fillna(0).to_numpy()
    ax.bar(x, vals, bottom=bottom, width=0.7, color=CLASS_COLORS[cls],
           edgecolor="white", linewidth=1.0, label=cls)
    bottom += vals
for xi, n in zip(x, c["n_total"].fillna(0).astype(int)):
    ax.text(xi, 101.5, f"n={n}", ha="center", va="bottom", fontsize=7, color=INK2)
ax.set_xticks(x, [SECTOR_NAMES[s] for s in SECTORS], rotation=45, ha="right")
ax.set_ylim(0, 110)
ax.set_ylabel("Share of hours [%]")
ax.set_title(f"Stability class per {SECTOR_WIDTH:.0f}° sector (ERA5 100 m wind direction)")
ax.legend(ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.32))
ax.grid(axis="x", visible=False)
save(fig, "02_class_share_by_sector")


# %% 5c. Stability rose (hours per sector, stacked by class) -----------------------
fig = plt.figure(figsize=(5.2, 5.2))
ax = fig.add_subplot(projection="polar")
ax.set_theta_zero_location("N")
ax.set_theta_direction(-1)                       # clockwise, like a compass
theta = np.radians(SECTORS)
width = np.radians(SECTOR_WIDTH) * 0.9
bottom = np.zeros(len(SECTORS))
for cls in st.CLASS_ORDER:
    vals = c[f"n_{cls}"].fillna(0).to_numpy()
    ax.bar(theta, vals, width=width, bottom=bottom, color=CLASS_COLORS[cls],
           edgecolor="white", linewidth=1.0, label=cls)
    bottom += vals
ax.set_xticks(np.radians(np.arange(0, 360, 45)), ["N", "NE", "E", "SE", "S", "SW", "W", "NW"])
ax.set_title("Stability rose: hours per sector", pad=18)
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.06), ncol=3)
save(fig, "03_stability_rose")


# %% 5d. 1/L distribution in each sector (small multiples, like Tuna et al. Fig. S.3)
n_cols = 4
n_rows = int(np.ceil(len(SECTORS) / n_cols))
fig, axes = plt.subplots(n_rows, n_cols, figsize=(9, 2.0 * n_rows), sharex=True, sharey=True)
bins = np.linspace(-0.03, 0.03, 61)
for ax, s in zip(axes.flat, SECTORS):
    sub = tab[tab["sector"] == s]
    for cls in st.CLASS_ORDER:
        v = sub.loc[sub["stability"] == cls, "inv_L"].clip(bins[0], bins[-1])
        ax.hist(v, bins=bins, color=CLASS_COLORS[cls], density=False, linewidth=0)
    ax.set_title(f"{SECTOR_NAMES[s]}  n={len(sub)}", fontsize=8)
    for xv in (-st.NEUTRAL_INV_L, st.NEUTRAL_INV_L):
        ax.axvline(xv, color=INK2, linewidth=0.5, linestyle="--")
for ax in axes.flat[len(SECTORS):]:
    ax.set_visible(False)
fig.supxlabel("1/L [1/m]")
fig.supylabel("Hours")
handles = [plt.Rectangle((0, 0), 1, 1, color=CLASS_COLORS[k]) for k in st.CLASS_ORDER]
fig.suptitle("1/L per sector (values outside ±0.03 are put in the end bins)", y=1.0)
fig.tight_layout()
fig.legend(handles, st.CLASS_ORDER, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.01))
save(fig, "04_invL_hist_per_sector")


# %% 5e. Season and day/night (sanity check figure) ----------------------------------
g = st.class_counts(tab, ["season", "is_day"])
labels = [f"{s}\n{'day' if d else 'night'}" for s, d in g.index]
x = np.arange(len(g))
fig, ax = plt.subplots(figsize=(7.0, 3.4))
bottom = np.zeros(len(g))
for cls in st.CLASS_ORDER:
    vals = g[f"pct_{cls}"].fillna(0).to_numpy()
    ax.bar(x, vals, bottom=bottom, width=0.7, color=CLASS_COLORS[cls],
           edgecolor="white", linewidth=1.0, label=cls)
    bottom += vals
ax.set_xticks(x, labels)
ax.set_ylabel("Share of hours [%]")
ax.set_ylim(0, 100)
ax.set_title("Stability class by season and day/night (sun above horizon = day)")
ax.legend(ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.22))
ax.grid(axis="x", visible=False)
save(fig, "05_class_share_season_daynight")


# %% 5f. Sector x season heatmap: (stable % - unstable %) ------------------------------
# +100 = always stable, -100 = always unstable, 0 = balanced or mostly neutral.
cs = st.class_counts(tab, ["sector", "season"])
diff = (cs["pct_stable"] - cs["pct_unstable"]).unstack("season").reindex(SECTORS)
n_cell = cs["n_total"].unstack("season").reindex(SECTORS)
diff_masked = diff.where(n_cell >= MIN_N_CELL)

cmap = LinearSegmentedColormap.from_list("unstable_stable", ["#e34948", "#f0efec", "#2a78d6"])
cmap.set_bad("#ffffff")
fig, ax = plt.subplots(figsize=(5.0, 6.0))
im = ax.imshow(diff_masked.to_numpy(dtype=float), cmap=cmap, vmin=-100, vmax=100, aspect="auto")
ax.set_xticks(range(len(st.SEASON_ORDER)), st.SEASON_ORDER)
ax.set_yticks(range(len(SECTORS)), [SECTOR_NAMES[s] for s in SECTORS])
ax.grid(False)
for i in range(len(SECTORS)):
    for j in range(len(st.SEASON_ORDER)):
        v, n = diff.iat[i, j], n_cell.iat[i, j]
        if pd.isna(n) or n < MIN_N_CELL:
            txt, color = f"n={0 if pd.isna(n) else int(n)}", MUTED
        else:
            txt = f"{v:+.0f}\nn={int(n)}"
            color = "white" if abs(v) > 60 else INK
        ax.text(j, i, txt, ha="center", va="center", fontsize=7, color=color)
cb = fig.colorbar(im, ax=ax, shrink=0.8)
cb.set_label("stable % − unstable %")
ax.set_title(f"Sector × season (cells with n < {MIN_N_CELL} h left blank)")
save(fig, "06_sector_season_heatmap")

plt.show()


# %% 6. LiDAR check: does the measured profile agree with the ERA5 class? -----------
# Expectation (physics): stable -> high shear alpha, strong veer, low TI;
#                         unstable -> low alpha, little veer, high TI.
# Remember Tuna & Bingöl (2020): shear works for neutral/stable, not unstable.
if RUN_LIDAR_PART:
    from pywes import lidar as ld

    FIT_HEIGHTS = [28, 38, 48, 58, 68, 78, 88, 99]   # [m] ~ rotor range 18-100 m
    LIDAR_TIMESTAMP_IS_END = False   # NOT documented for the LiDAR - check! (only 5 min effect)

    lid_raw = ld.load_zx300_10min(LIDAR_DIR)
    print(f"\nZX300: {lid_raw.index.min()} -> {lid_raw.index.max()}, {len(lid_raw)} records")
    print("Status Flags found:", lid_raw["Status Flags"].value_counts().to_dict())

    lid, lid_report = ld.clean_zx300(lid_raw, FIT_HEIGHTS)
    print("LiDAR filter report:")
    print(lid_report.to_string(index=False, float_format="%.1f"))

    lid["alpha"] = ld.shear_exponent(lid, FIT_HEIGHTS)
    lid["veer"] = ld.wind_veer(lid, 28, 99)
    lid["TI58"] = lid["TI at 58m"]
    lid["wd_lidar"] = lid["Wind Direction (deg) at 58m"]

    # Remove directions where a turbine wake can reach the LiDAR
    disturbed = ld.iec_disturbed_sectors(ld.ZX300_LAT, ld.ZX300_LON, turbines, ROTOR_D)
    print("\nTurbines that can wake the LiDAR:")
    print(disturbed.to_string(index=False, float_format="%.1f"))
    waked = ld.in_disturbed_sector(lid["wd_lidar"], disturbed)
    print(f"LiDAR records in disturbed sectors (removed): {waked.mean() * 100:.1f} %")
    lid = lid[~waked]

    # Join each 10-min record to the nearest ERA5 hour
    lid["era5_hour"] = st.nearest_era5_hour(lid.index, "10min", LIDAR_TIMESTAMP_IS_END)
    m = lid.join(tab[["inv_L", "stability", "sector", "wd100", "ws100"]], on="era5_hour", how="inner")
    print(f"LiDAR records with a valid ERA5 class: {len(m)}")

    # --- 6a. ERA5 vs LiDAR wind direction (is the ERA5 direction rotated?) ---
    d_dir = (m["wd100"] - m["Wind Direction (deg) at 99m"] + 180) % 360 - 180
    fig, ax = plt.subplots(figsize=(5.5, 3.0))
    ax.hist(d_dir.dropna(), bins=np.arange(-60, 61, 2), color="#2a78d6", edgecolor="white", linewidth=0.3)
    ax.axvline(d_dir.median(), color=INK, linewidth=1)
    ax.set_xlabel("ERA5 wd100 − LiDAR wd99 [deg]")
    ax.set_ylabel("10-min records")
    ax.set_title(f"Direction difference, median = {d_dir.median():+.1f}° "
                 "(LiDAR in magnetic north)")
    save(fig, "07_era5_vs_lidar_direction")

    # --- 6b. alpha, veer and TI per ERA5 stability class (box plots) ---
    fig, axes = plt.subplots(1, 3, figsize=(9, 3.2))
    for ax, col, lab in zip(axes, ["alpha", "veer", "TI58"],
                            ["Shear exponent α [-]", "Veer 28→99 m [deg]", "TI at 58 m [-]"]):
        data = [m.loc[m["stability"] == k, col].dropna() for k in st.CLASS_ORDER]
        bp = ax.boxplot(data, showfliers=False, widths=0.5, patch_artist=True,
                        medianprops={"color": INK, "linewidth": 1.2})
        for patch, k in zip(bp["boxes"], st.CLASS_ORDER):
            patch.set_facecolor(CLASS_COLORS[k])
            patch.set_edgecolor("white")
        ax.set_xticks([1, 2, 3], [f"{k}\nn={len(d)}" for k, d in zip(st.CLASS_ORDER, data)])
        ax.set_ylabel(lab)
        ax.grid(axis="x", visible=False)
    fig.suptitle("ZX300 LiDAR profiles grouped by ERA5 stability class")
    fig.tight_layout()
    save(fig, "08_lidar_by_class")

    # --- 6c. median alpha per sector and class ---
    piv = m.pivot_table(index="sector", columns="stability", values="alpha",
                        aggfunc="median", observed=False).reindex(index=SECTORS, columns=st.CLASS_ORDER)
    npiv = m.pivot_table(index="sector", columns="stability", values="alpha",
                         aggfunc="count", observed=False).reindex(index=SECTORS, columns=st.CLASS_ORDER)
    piv_masked = piv.where(npiv >= MIN_N_CELL)
    cmap_seq = LinearSegmentedColormap.from_list("blue_seq", ["#cde2fb", "#2a78d6", "#0d366b"])
    cmap_seq.set_bad("#ffffff")
    fig, ax = plt.subplots(figsize=(4.6, 6.0))
    im = ax.imshow(piv_masked.to_numpy(dtype=float), cmap=cmap_seq, aspect="auto")
    ax.set_xticks(range(3), st.CLASS_ORDER)
    ax.set_yticks(range(len(SECTORS)), [SECTOR_NAMES[s] for s in SECTORS])
    ax.grid(False)
    vmax = np.nanmax(piv_masked.to_numpy(dtype=float)) if piv_masked.notna().any().any() else 1
    for i in range(len(SECTORS)):
        for j in range(3):
            v, n = piv.iat[i, j], npiv.iat[i, j]
            if pd.isna(n) or n < MIN_N_CELL:
                ax.text(j, i, f"n={0 if pd.isna(n) else int(n)}", ha="center", va="center",
                        fontsize=7, color=MUTED)
            else:
                color = "white" if v > 0.6 * vmax else INK
                ax.text(j, i, f"{v:.2f}\nn={int(n)}", ha="center", va="center", fontsize=7, color=color)
    fig.colorbar(im, ax=ax, shrink=0.8).set_label("median α [-]")
    ax.set_title("Median shear exponent α (LiDAR) per sector and ERA5 class")
    save(fig, "09_lidar_alpha_sector_class")

    # Agreement check in numbers
    med = m.groupby("stability", observed=False)[["alpha", "veer", "TI58"]].median()
    print("\nMedian LiDAR values per ERA5 class:")
    print(med.to_string(float_format="%.3f"))
    print("CHECK alpha stable > neutral > unstable:",
          "OK" if med.loc["stable", "alpha"] > med.loc["neutral", "alpha"] > med.loc["unstable", "alpha"]
          else "NOT OK - look at sectors / terrain")

    plt.show()
