# -*- coding: utf-8 -*-
"""Wake pair: upstream turbine A (free stream) vs downstream turbine B, per stability class.

Turbine A gets free-stream wind; turbine B stands behind A. For wind from the
alignment direction (A -> B), B is in A's wake. We compare wind speed and power
of B with A and group the result by ERA5 stability class.

    velocity deficit = 1 - mean(ws_B) / mean(ws_A)
    power deficit    = 1 - sum(P_B)  / sum(P_A)     (energy-weighted, like a wake loss)

Steps (run cell by cell in Spyder, or the whole file):
  0. Settings  - the only part you need to edit
  1. Load data and print the pair geometry
  2. Join ERA5 (nearest hour) and scale the ERA5 wind speed to hub height
  3. Check SCADA against ERA5 (wind speed, direction, time series)
  4. Direction and speed filters (each step reports the removed samples)
  5. Wake comparison per stability class (tables + figures)
  6. Physics sanity checks (OK / NOT OK) and save the tables

Inputs : data/processed/scada_clean_core_2024_2025.nc   (scripts/clean_scada_2024_2025.py)
         data/processed/era5_stability_hourly.csv      (analysis/stability_by_sector.py)
         data/processed/turbine_layout_utm.csv         (analysis/site_map.py)
         data/sample/Hill_of_Towie_turbine_metadata.csv
         data/sample/HOT_DY_Design_northing.yaml       (RES northing corrections)
Outputs: results/figures/draft/turbine_comparison/<A>_<B>_<dir>_*.png and *.csv

Period: 2024-01-01 to 2025-10-27 (inclusive). Dynamic Yaw testing starts on
2025-10-28 13:00 UTC; after that, yaw offsets can be active and the 10-min
SCADA has no ON/OFF flag, so those data are not a "no wake control" baseline.

Limitations (keep in mind when reading the results)
  * Northing: nacelle positions are corrected to true north with RES's
    corrections (pywes.northing, cell 1). The remaining offset to ERA5 is
    reported in cell 3; it should be small (ERA5 grid, terrain turning).
  * The direction window uses ERA5 (hourly, ~30 km grid). The measured wake
    profile is therefore wider and shallower than the real one - for all
    classes, so the comparison between classes stays fair.
  * The nacelle anemometer of B sits behind B's own rotor inside A's wake, so
    the velocity deficit is an indicator, not a free-stream wake deficit. The
    power deficit is the better wake-loss measure (ROADMAP Phase 4).
  * You choose a pair where A has no turbine upstream and B only has A
    upstream for the analysed direction; the script does not check this.

References
  Stability classes and ERA5 filters: Tuna, Bingöl & Sofuoğlu (2026), JWEIA 272, 106396.
  Shear exponent (power law): Tuna & Bingöl (2020), Energies 13, 89.
"""

# %% 0. Settings ------------------------------------------------------------------
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.colors import LinearSegmentedColormap

from pywes.paths import PROCESSED_DIR, SAMPLE_DIR, FIGURES_DIR
from pywes import northing
from pywes import pairs
from pywes import stability as st
from pywes.scada_cleaning import DYNAMIC_YAW_FIRST, TUNEUP_FIRST_AFTER, TUNEUP_TURBINES

# --- pair and direction --------------------------------------------------------------
TURBINE_A = "T11"           # UPSTREAM turbine (free stream)   - EXAMPLE, change it
TURBINE_B = "T14"           # DOWNSTREAM turbine (in A's wake) - EXAMPLE, change it
WIND_DIR_DEG = None         # centre of the analysed direction [deg, from north, ERA5 100 m]
                            # None = use the A -> B alignment direction (printed in cell 1)
DIR_TOL_DEG = 30.0          # ERA5 direction window: WIND_DIR_DEG +/- DIR_TOL_DEG [deg]
SAME_DIR_TOL_DEG = 20.0     # max |nacelle_A - nacelle_B| [deg] ("facing the same way")
CORE_HALF_WIDTH_DEG = 20.0  # "wake core" = |ERA5 dir - WIND_DIR_DEG| <= this [deg];
                            # used for the per-class table, the ws-binned plot and the checks

# --- wind speed, power, bins ----------------------------------------------------------
WS_MIN, WS_MAX = 4.0, 25.0  # free-stream speed range, from turbine A [m/s]
P_MIN_KW = 50.0             # remove records with power_A below this [kW] (power ratio)
WS_BIN = 1.0                # bin width for deficit vs wind speed [m/s]
DIR_BIN = 5.0               # bin width for the wake profile (deficit vs direction) [deg]
MIN_N_BIN = 30              # hide bins with fewer 10-min samples

# --- SCADA vs ERA5 check -----------------------------------------------------------------
MIN_REC_PER_HOUR = 4        # min. kept 10-min records per hour for the hourly mean (max 6)
MIN_CORR = 0.8              # our own rule of thumb for "SCADA agrees with ERA5"
TS_START = None             # start of the time-series plot, e.g. "2025-02-01"; None = auto
TS_DAYS = 31                 # length of the time-series plot [days]

# --- period (SCADA TimeStamp = end of the 10-min period, UTC) ------------------------------
# End date: 27 October 2025 (inclusive). Dynamic Yaw testing starts on
# 2025-10-28 13:00 UTC (RES exclusion file); after that, yaw offsets can be active
# and the 10-min SCADA has no ON/OFF flag. TimeStamp 2025-10-28 00:00 is the END
# of the last 10-min period of 27 Oct (23:50-00:00), so END = "2025-10-28 00:00".
START = "2024-01-01 00:10"   # first period 2024-01-01 00:00-00:10
END = "2025-10-28 00:00"     # last period 2025-10-27 23:50-00:00 (= end of 27 Oct 2025)

# --- northing (turbine north -> true north) -------------------------------------------------
APPLY_NORTHING = True                          # False only to compare with the old results
NORTHING_FILE = northing.DEFAULT_NORTHING_FILE  # RES HOT_DY_Design_northing.yaml

NC_PATH = PROCESSED_DIR / "scada_clean_core_2024_2025.nc"
ERA5_PATH = PROCESSED_DIR / "era5_stability_hourly.csv"
LAYOUT_PATH = PROCESSED_DIR / "turbine_layout_utm.csv"
FIG_DIR = FIGURES_DIR / "draft" / "turbine_comparison"
FIG_DIR.mkdir(parents=True, exist_ok=True)

# --- plot style (same as analysis/stability_by_sector.py) ---------------------------------
CLASS_COLORS = {"unstable": "#e34948", "neutral": "#a3a19a", "stable": "#2a78d6"}
TURB_COLORS = {"A": "#2a78d6", "B": "#eb6834"}     # categorical slots 1 and 2 (validated)
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
    """Save a figure as PNG in the draft folder (file name starts with the pair tag)."""
    path = FIG_DIR / f"{TAG}_{name}.png"
    fig.savefig(path)
    print(f"saved {path}")


def label_line_ends(ax, lines):
    """Direct labels at the last valid point of each line, spread so they don't overlap.

    lines: list of (x, y, text). Text is in ink, not in the series color.
    """
    ends = []
    for x, y, text in lines:
        y = np.asarray(y, dtype=float)
        ok = np.isfinite(y)
        if ok.any():
            ends.append([np.asarray(x, dtype=float)[ok][-1], y[ok][-1], text])
    if not ends:
        return
    y0, y1 = ax.get_ylim()
    gap = 0.06 * (y1 - y0)                      # min. vertical distance between labels
    ends.sort(key=lambda e: e[1])
    ys = [e[1] for e in ends]
    for i in range(1, len(ys)):
        ys[i] = max(ys[i], ys[i - 1] + gap)
    dx = 0.015 * (ax.get_xlim()[1] - ax.get_xlim()[0])
    for (x, _, text), y_lab in zip(ends, ys):
        ax.text(x + dx, y_lab, text, va="center", ha="left", fontsize=7, color=INK2, clip_on=False)


def break_wraps(series):
    """Insert NaN where a direction jumps over 0/360 deg, so plotted lines don't draw a vertical jump."""
    s = series.copy()
    jump = s.diff().abs() > 180.0
    s[jump] = np.nan
    return s


# %% 1. Load data and pair geometry ----------------------------------------------------
meta = pd.read_csv(SAMPLE_DIR / "Hill_of_Towie_turbine_metadata.csv", encoding="utf-8-sig")
HUB_HEIGHT = float(meta["Hub Height (m)"].iloc[0])        # 59 m
ROTOR_D = float(meta["Rotor Diameter (m)"].iloc[0])       # 82 m
P_RATED = float(meta["Rated power (kW)"].iloc[0])         # 2300 kW

layout = pd.read_csv(LAYOUT_PATH)
geo = pairs.pair_geometry(layout, TURBINE_A, TURBINE_B, ROTOR_D)
WIND_DIR = geo["align_dir_deg"] if WIND_DIR_DEG is None else float(WIND_DIR_DEG)
TAG = f"{TURBINE_A}_{TURBINE_B}_{WIND_DIR:03.0f}"

print(f"Pair {TURBINE_A} (upstream) -> {TURBINE_B} (downstream)")
print(f"  alignment direction (B behind A): {geo['align_dir_deg']:.1f} deg")
print(f"  spacing: {geo['spacing_m']:.0f} m = {geo['spacing_D']:.2f} D (D = {ROTOR_D:.0f} m)")
print(f"  hub elevation B - A: {geo['dz_hub_m']:+.1f} m")
print(f"  analysed direction: {WIND_DIR:.1f} +/- {DIR_TOL_DEG:g} deg (ERA5 100 m)")
off_align = abs(pairs.angle_diff(WIND_DIR, geo["align_dir_deg"]))
if off_align > DIR_TOL_DEG:
    print(f"  WARNING: analysed direction is {off_align:.1f} deg away from the alignment "
          f"direction (> DIR_TOL_DEG) - B is probably NOT in A's wake.")

df_raw = pairs.load_pair(NC_PATH, TURBINE_A, TURBINE_B, START, END)
print(f"\nSCADA: {df_raw.index.min()} -> {df_raw.index.max()}, {len(df_raw)} 10-min timestamps")

# Which cleaning made the NetCDF? (attributes written since 2026-11-16)
nc_attrs = df_raw.attrs.get("nc_attrs", {})
groups = nc_attrs.get("exclusion_groups", "")
if "res_curtailment" not in groups or "dynamic_yaw" not in groups:
    print("  WARNING: this NetCDF was made without the RES exclusion periods "
          f"(exclusion_groups = {groups or 'not recorded'}). Re-run "
          "scripts/clean_scada_2024_2025.py to remove the RES curtailment periods.")

# Dynamic Yaw: no record after the start of testing may be used as baseline
last_start = df_raw.index.max() - pd.Timedelta("10min")
if last_start >= pd.Timestamp(DYNAMIC_YAW_FIRST):
    print(f"  WARNING: data reach into Dynamic Yaw testing (from {DYNAMIC_YAW_FIRST} UTC) - "
          "set END to 2025-10-28 00:00 or earlier.")

# Pitch TuneUp 2024-03-14: power curve of these turbines changes on that day
tuned = [t for t in (TURBINE_A, TURBINE_B) if t in TUNEUP_TURBINES]
if len(tuned) == 1 and pd.Timestamp(START) < pd.Timestamp(TUNEUP_FIRST_AFTER):
    print(f"  NOTE: {tuned[0]} got the pitch TuneUp on 2024-03-14, the other turbine did not. "
          f"The power ratio P_B/P_A can change on that day; for power deficits consider "
          f"START = '{TUNEUP_FIRST_AFTER}' (period start) or later.")

# Northing: turbine north -> true north (RES corrections, see pywes.northing)
if APPLY_NORTHING:
    north_table = northing.load_northing_table(NORTHING_FILE)
    df_raw = northing.apply_to_pair(df_raw, north_table, TURBINE_A, TURBINE_B,
                                    columns=("nacelle_pos",), timestamp_is_end=True)
    seg = northing.segments(north_table, pd.Timestamp(START) - pd.Timedelta("10min"), END,
                            turbines=[TURBINE_A, TURBINE_B])
    print(f"Northing applied ({NORTHING_FILE.name}): true = (turbine north + correction) mod 360")
    print(seg.to_string(index=False, float_format="%.2f"))
    for t in (TURBINE_A, TURBINE_B):
        if (seg["turbine"] == t).sum() > 1:
            print(f"  NOTE: {t} has a northing change inside the period (see valid_from above).")
    n_nan = int((df_raw["nacelle_pos_tn_A"].notna() & df_raw["nacelle_pos_A"].isna()).sum()
                + (df_raw["nacelle_pos_tn_B"].notna() & df_raw["nacelle_pos_B"].isna()).sum())
    if n_nan:
        print(f"  WARNING: {n_nan} records have no northing correction (NaN direction).")
    FRAME = "true north"
else:
    print("Northing NOT applied: nacelle positions are in turbine north.")
    FRAME = "turbine north"
for t in ("A", "B"):
    n = int((df_raw[f"keep_{t}"] == 1).sum())
    print(f"  kept records {t} ({df_raw.attrs[f'turbine_{t}']}): {n} ({100 * n / len(df_raw):.1f} %)")
n_both = int(((df_raw["keep_A"] == 1) & (df_raw["keep_B"] == 1)).sum())
print(f"  both kept: {n_both} ({100 * n_both / len(df_raw):.1f} %)")

era5 = pairs.load_era5_table(ERA5_PATH)
print(f"ERA5 table: {era5.index.min()} -> {era5.index.max()}, {len(era5)} hours, "
      f"{era5['kept'].mean() * 100:.1f} % kept by the stability filters")


# %% 2. Join ERA5 and scale ERA5 wind speed to hub height --------------------------------
df, era5_report = pairs.join_era5(df_raw, era5, HUB_HEIGHT)
print("ERA5 join report (10-min timestamps):")
print(era5_report.to_string(index=False, float_format="%.2f"))
print(f"ERA5 shear exponent alpha (10-100 m): median {df['alpha_era5'].median():.3f}, "
      f"5-95 %: {df['alpha_era5'].quantile(0.05):.3f} ... {df['alpha_era5'].quantile(0.95):.3f}")
print(f"ERA5 speed scaled from 100 m to {HUB_HEIGHT:.0f} m: mean ws_hub / ws100 = "
      f"{(df['ws_hub_era5'] / df['ws100']).mean():.3f}")


# %% 3. Check SCADA against ERA5 (all directions, before the direction filters) ----------
# Hourly means of the kept 10-min records (ERA5 wind is one value per hour).
hourly = {t: pairs.hourly_mean(df, t, MIN_REC_PER_HOUR) for t in ("A", "B")}
# B is waked by A around the analysed direction - leave that sector out of B's speed check,
# otherwise the wake would look like a sensor bias.
in_sector_B = pairs.angle_diff(hourly["B"]["wd100"], WIND_DIR).abs() <= DIR_TOL_DEG
hourly_ws = {"A": hourly["A"], "B": hourly["B"][~in_sector_B]}

rows = []
for t in ("A", "B"):
    s = pairs.agreement_stats(hourly_ws[t]["ws_hub_era5"], hourly_ws[t]["ws"])
    d = pairs.direction_offset(hourly[t]["nacelle_pos"], hourly[t]["wd100"])
    rows.append({"turbine": df_raw.attrs[f"turbine_{t}"], "role": t,
                 "ws_n_hours": s["n"], "ws_bias_mps": s["bias"], "ws_rmse_mps": s["rmse"],
                 "ws_corr": s["corr"], "ws_ratio_mean": s["ratio_mean"],
                 "dir_n_hours": d["n"], "dir_offset_median_deg": d["median"],
                 "dir_offset_q25_deg": d["q25"], "dir_offset_q75_deg": d["q75"]})
verify = pd.DataFrame(rows).set_index("turbine")
print("SCADA vs ERA5 (hourly means; ws vs ERA5 scaled to hub height; "
      f"B without the sector {WIND_DIR:.0f} +/- {DIR_TOL_DEG:g} deg):")
print(verify.drop(columns="role").T.to_string(float_format="%.3f"))
if APPLY_NORTHING:
    print("dir_offset = circular median of (nacelle_pos - ERA5 wd100) AFTER northing: should be "
          "small; what is left comes from ERA5 (30 km grid, 100 m) and terrain turning the flow.")
else:
    print("dir_offset = circular median of (nacelle_pos - ERA5 wd100): a ROUGH northing-offset "
          "estimate (northing not applied).")

both = df[(df["keep_A"] == 1) & (df["keep_B"] == 1)]
rel_off = pairs.direction_offset(both["nacelle_pos_A"], both["nacelle_pos_B"])
print(f"Relative offset nacelle_A - nacelle_B (all directions, 10-min, {FRAME}): median "
      f"{rel_off['median']:+.1f} deg, IQR {rel_off['q25']:+.1f} ... {rel_off['q75']:+.1f} deg "
      f"(n = {rel_off['n']}).")
if APPLY_NORTHING:
    rel_raw = pairs.direction_offset(both["nacelle_pos_tn_A"], both["nacelle_pos_tn_B"])
    print(f"  before northing (turbine north): median {rel_raw['median']:+.1f} deg")

# --- 3a. scatter SCADA vs ERA5 hub-height wind speed -------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(7.6, 4.3), sharex=True, sharey=True)
lim = [0, max(25.0, np.nanpercentile(hourly["A"]["ws"], 99.5))]
for ax, t in zip(axes, ("A", "B")):
    h = hourly_ws[t]
    s = verify.iloc[0 if t == "A" else 1]
    ax.scatter(h["ws_hub_era5"], h["ws"], s=3, color=TURB_COLORS[t], alpha=0.25, linewidths=0)
    ax.plot(lim, lim, color=INK2, linewidth=0.8, linestyle="--", label="1:1")
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_aspect("equal")
    ax.set_xlabel(f"ERA5 wind speed scaled to {HUB_HEIGHT:.0f} m [m/s]")
    extra = "" if t == "A" else f"\n(sector {WIND_DIR:.0f} ± {DIR_TOL_DEG:g}° removed)"
    ax.set_title(f"{s.name} ({'upstream' if t == 'A' else 'downstream'}){extra}\n"
                 f"n = {int(s['ws_n_hours'])} h, bias {s['ws_bias_mps']:+.2f} m/s, "
                 f"RMSE {s['ws_rmse_mps']:.2f} m/s, r = {s['ws_corr']:.2f}", fontsize=8)
axes[0].set_ylabel("SCADA nacelle wind speed, hourly mean [m/s]")
axes[0].legend(loc="upper left")
fig.suptitle("SCADA vs ERA5 wind speed (all directions)")
fig.tight_layout()
save(fig, "01_scada_vs_era5_ws")

# --- 3b. SCADA / ERA5 speed ratio per direction sector ------------------------------------
# Shows where terrain or wakes of OTHER turbines change the ratio. B's wake sector is kept
# here on purpose (it shows up as a dip at the analysed direction).
fig, ax = plt.subplots(figsize=(7.0, 3.4))
ends = []
for t in ("A", "B"):
    r = pairs.ratio_by_sector(hourly[t], 30.0, MIN_N_BIN)
    ax.plot(r.index, r.values, color=TURB_COLORS[t], linewidth=2, marker="o", markersize=5,
            label=df_raw.attrs[f"turbine_{t}"])
    ends.append((r.index, r.values, df_raw.attrs[f"turbine_{t}"]))
label_line_ends(ax, ends)
ax.axvspan(WIND_DIR - DIR_TOL_DEG, WIND_DIR + DIR_TOL_DEG, color=GRID, alpha=0.6, zorder=0,
           label="analysed sector")
ax.axhline(1.0, color=MUTED, linewidth=0.6)
ax.set_xticks(np.arange(0, 360, 30))
ax.set_xlabel("ERA5 wind direction at 100 m, 30° sector centre [deg]")
ax.set_ylabel("median SCADA ws / ERA5 ws(hub) [-]")
ax.set_title(f"SCADA / ERA5 wind speed ratio by sector (sectors with n < {MIN_N_BIN} h hidden)")
ax.legend(ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.2))
save(fig, "02_scada_era5_ratio_by_sector")

# --- 3c. direction differences -------------------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(10, 3.0), sharey=False)
bins = np.arange(-60, 61, 2)
items = [(pairs.angle_diff(hourly["A"]["nacelle_pos"], hourly["A"]["wd100"]),
          f"{TURBINE_A} nacelle − ERA5 wd100", TURB_COLORS["A"]),
         (pairs.angle_diff(hourly["B"]["nacelle_pos"], hourly["B"]["wd100"]),
          f"{TURBINE_B} nacelle − ERA5 wd100", TURB_COLORS["B"]),
         (pairs.angle_diff(both["nacelle_pos_A"], both["nacelle_pos_B"]),
          f"{TURBINE_A} − {TURBINE_B} nacelle (10-min)", INK2)]
for ax, (d, lab, col) in zip(axes, items):
    d = d.dropna()
    ax.hist(d.clip(bins[0], bins[-1]), bins=bins, color=col, edgecolor="white", linewidth=0.3)
    ax.axvline(d.median(), color=INK, linewidth=1)
    ax.set_title(f"{lab}\nmedian {d.median():+.1f}°, n = {len(d)}", fontsize=8)
    ax.set_xlabel("difference [deg] (end bins hold values beyond ±60°)")
axes[0].set_ylabel("count")
fig.suptitle(f"Direction differences, all directions (nacelle in {FRAME}, ERA5 in true north)")
fig.tight_layout()
save(fig, "03_direction_differences")

# --- 3d. time series: SCADA (10-min) vs ERA5 (hourly) ---------------------------------------
if TS_START is None:
    # first day from which both turbines have kept data for TS_DAYS days in a row (or the start)
    daily_ok = both["ws_A"].resample("D").count() >= 100
    run = daily_ok.rolling(TS_DAYS).sum() >= TS_DAYS
    ts0 = run[run].index[0] - pd.Timedelta(days=TS_DAYS - 1) if run.any() else df.index[0].normalize()
else:
    ts0 = pd.Timestamp(TS_START)
ts1 = ts0 + pd.Timedelta(days=TS_DAYS)
w = df_raw.loc[ts0:ts1]
e = era5.loc[ts0:ts1]
e_hub, _ = pairs.era5_to_hub_height(e["ws10"], e["ws100"], HUB_HEIGHT)
fig, axes = plt.subplots(2, 1, figsize=(10, 5.0), sharex=True)
for t in ("A", "B"):
    axes[0].plot(w.index, w[f"ws_{t}"], color=TURB_COLORS[t], linewidth=1.2,
                 label=f"{df_raw.attrs[f'turbine_{t}']} SCADA (10-min)")
    axes[1].plot(w.index, break_wraps(w[f"nacelle_pos_{t}"]), color=TURB_COLORS[t], linewidth=1.2,
                 label=f"{df_raw.attrs[f'turbine_{t}']} nacelle ({FRAME})")
axes[0].plot(e.index, e_hub, color=INK, linewidth=1.5, marker="o", markersize=3,
             label=f"ERA5 scaled to {HUB_HEIGHT:.0f} m (hourly)")
axes[1].plot(e.index, break_wraps(e["wd100"]), color=INK, linewidth=1.5, marker="o", markersize=3,
             label="ERA5 wd100 (true north)")
axes[0].set_ylabel("wind speed [m/s]")
axes[1].set_ylabel("direction [deg]")
axes[1].set_ylim(0, 360)
axes[1].set_yticks(np.arange(0, 361, 90))
for ax in axes:
    ax.legend(ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.12), fontsize=7)
locator = mdates.AutoDateLocator()
axes[1].xaxis.set_major_locator(locator)
axes[1].xaxis.set_major_formatter(mdates.ConciseDateFormatter(locator))
axes[0].set_title(f"SCADA vs ERA5, {ts0:%Y-%m-%d} + {TS_DAYS} days (gaps = removed records)")
fig.tight_layout()
save(fig, "04_timeseries_scada_era5")
plt.show()


# %% 4. Direction and speed filters -----------------------------------------------------------
dfw, pair_report = pairs.filter_pair(df, WIND_DIR, DIR_TOL_DEG, SAME_DIR_TOL_DEG,
                                     WS_MIN, WS_MAX, P_MIN_KW)
dfw = pairs.add_wake_columns(dfw)
print("Pair filter report (10-min timestamps, after the ERA5 join):")
print(pair_report.to_string(index=False, float_format="%.2f"))
core = dfw[dfw["rel_dir"].abs() <= CORE_HALF_WIDTH_DEG]
print(f"Wake core |dir - {WIND_DIR:.0f}| <= {CORE_HALF_WIDTH_DEG:g} deg: {len(core)} records")
print("Samples per class (window / core):")
print(pd.DataFrame({"window": dfw["stability"].value_counts(sort=False),
                    "core": core["stability"].value_counts(sort=False)}).to_string())


# %% 5. Wake comparison per stability class ------------------------------------------------------
stats_window = pairs.pair_stats(dfw)
stats_core = pairs.pair_stats(core)
show = ["n", "ws_A_mean", "ws_B_mean", "ws_deficit", "ws_deficit_std", "ws_fit_slope",
        "ws_fit_intercept", "p_deficit", "p_deficit_std", "p_ratio_mean"]
print(f"Whole window ({WIND_DIR:.0f} +/- {DIR_TOL_DEG:g} deg):")
print(stats_window[show].to_string(float_format="%.3f"))
print(f"\nWake core ({WIND_DIR:.0f} +/- {CORE_HALF_WIDTH_DEG:g} deg):")
print(stats_core[show].to_string(float_format="%.3f"))

# Direction bins centred on 0 (the analysed direction)
half = DIR_BIN / 2.0
k = int(np.ceil((DIR_TOL_DEG - half) / DIR_BIN))
dir_edges = np.arange(-k, k + 2) * DIR_BIN - half
profile = pairs.binned_deficit(dfw, "rel_dir", dir_edges, MIN_N_BIN)
profile_all = pairs.binned_deficit(dfw.assign(group="all"), "rel_dir", dir_edges, MIN_N_BIN, by="group")
ws_edges = np.arange(WS_MIN, WS_MAX + WS_BIN / 2, WS_BIN)
by_ws = pairs.binned_deficit(core, "ws_A", ws_edges, MIN_N_BIN)

# --- 5a. scatter ws_A vs ws_B (wake core), coloured by class -----------------------------------
fig, ax = plt.subplots(figsize=(5.2, 5.0))
for cls in st.CLASS_ORDER:
    sub = core[core["stability"] == cls]
    ax.scatter(sub["ws_A"], sub["ws_B"], s=4, color=CLASS_COLORS[cls], alpha=0.35,
               linewidths=0, label=f"{cls} (n={len(sub)})")
lim = [WS_MIN - 1, WS_MAX + 1]
ax.plot(lim, lim, color=INK2, linewidth=0.8, linestyle="--", label="1:1 (no wake)")
ax.set_xlim(lim)
ax.set_ylim(lim)
ax.set_aspect("equal")
ax.set_xlabel(f"ws {TURBINE_A} (upstream) [m/s]")
ax.set_ylabel(f"ws {TURBINE_B} (downstream) [m/s]")
ax.set_title(f"Wake core {WIND_DIR:.0f} ± {CORE_HALF_WIDTH_DEG:g}°: points below 1:1 = wake")
ax.legend(markerscale=3, loc="upper left")
save(fig, "05_scatter_wsA_wsB")

# --- 5b. wake profile: deficit vs direction ------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), sharex=True)
for ax, col, sem, lab in zip(axes, ["ws_deficit", "p_deficit"], ["ws_sem", "p_sem"],
                             ["velocity deficit 1 − ws_B/ws_A [-]", "power deficit 1 − P_B/P_A [-]"]):
    ends = []
    for cls in st.CLASS_ORDER:
        p = profile[profile["stability"] == cls]
        ax.plot(p["bin"], p[col], color=CLASS_COLORS[cls], linewidth=2, marker="o", markersize=5,
                label=cls)
        ax.fill_between(p["bin"], p[col] - p[sem], p[col] + p[sem], color=CLASS_COLORS[cls],
                        alpha=0.15, linewidth=0)
        ends.append((p["bin"], p[col], cls))
    ax.axhline(0, color=MUTED, linewidth=0.6)
    ax.axvline(0, color=INK2, linewidth=0.6, linestyle="--")
    label_line_ends(ax, ends)
    ax.set_xlabel(f"ERA5 direction − {WIND_DIR:.0f}° [deg]")
    ax.set_ylabel(lab)
fig.suptitle(f"Wake profile {TURBINE_A} → {TURBINE_B} ({geo['spacing_D']:.1f} D), "
             f"band = ± 1 standard error, bins with n < {MIN_N_BIN} hidden")
fig.tight_layout()
fig.legend(*axes[0].get_legend_handles_labels(), ncol=3, loc="upper center", bbox_to_anchor=(0.5, 0.0))
save(fig, "06_wake_profile")

# --- 5c. deficit vs free-stream wind speed (wake core) -------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), sharex=True)
for ax, col, sem, lab in zip(axes, ["ws_deficit", "p_deficit"], ["ws_sem", "p_sem"],
                             ["velocity deficit [-]", "power deficit [-]"]):
    ends = []
    for cls in st.CLASS_ORDER:
        b = by_ws[by_ws["stability"] == cls]
        ax.plot(b["bin"], b[col], color=CLASS_COLORS[cls], linewidth=2, marker="o", markersize=5,
                label=cls)
        ax.fill_between(b["bin"], b[col] - b[sem], b[col] + b[sem], color=CLASS_COLORS[cls],
                        alpha=0.15, linewidth=0)
        ends.append((b["bin"], b[col], cls))
    ax.axhline(0, color=MUTED, linewidth=0.6)
    label_line_ends(ax, ends)
    ax.set_xlabel(f"ws {TURBINE_A} (free stream) [m/s]")
    ax.set_ylabel(lab)
fig.suptitle(f"Deficit vs wind speed, wake core {WIND_DIR:.0f} ± {CORE_HALF_WIDTH_DEG:g}° "
             f"(bins with n < {MIN_N_BIN} hidden)")
fig.tight_layout()
fig.legend(*axes[0].get_legend_handles_labels(), ncol=3, loc="upper center", bbox_to_anchor=(0.5, 0.0))
save(fig, "07_deficit_vs_ws")

# --- 5d. box plots of the 10-min deficits per class (wake core) -----------------------------------
fig, axes = plt.subplots(1, 2, figsize=(8, 3.4))
for ax, col, lab in zip(axes, ["ws_deficit", "p_deficit"], ["velocity deficit [-]", "power deficit [-]"]):
    data = [core.loc[core["stability"] == c, col].dropna() for c in st.CLASS_ORDER]
    bp = ax.boxplot(data, showfliers=False, widths=0.5, patch_artist=True,
                    medianprops={"color": INK, "linewidth": 1.2})
    for patch, c in zip(bp["boxes"], st.CLASS_ORDER):
        patch.set_facecolor(CLASS_COLORS[c])
        patch.set_edgecolor("white")
    ax.set_xticks([1, 2, 3], [f"{c}\nn={len(d)}" for c, d in zip(st.CLASS_ORDER, data)])
    ax.axhline(0, color=MUTED, linewidth=0.6)
    ax.set_ylabel(lab)
    ax.grid(axis="x", visible=False)
fig.suptitle(f"10-min deficits per class, wake core {WIND_DIR:.0f} ± {CORE_HALF_WIDTH_DEG:g}° "
             "(box = 25-75 %, whiskers = 1.5 IQR)")
fig.tight_layout()
save(fig, "08_deficit_boxplot")

# --- 5e. nacelle difference after the filters -------------------------------------------------------
fig, ax = plt.subplots(figsize=(5.5, 3.0))
bins = np.arange(-SAME_DIR_TOL_DEG - 1, SAME_DIR_TOL_DEG + 1.01, 0.5)
ax.hist(dfw["dnac"], bins=bins, color=INK2, edgecolor="white", linewidth=0.3)
for xv in (-SAME_DIR_TOL_DEG, SAME_DIR_TOL_DEG):
    ax.axvline(xv, color=MUTED, linewidth=0.8, linestyle="--")
ax.axvline(dfw["dnac"].median(), color=INK, linewidth=1)
ax.set_xlabel(f"nacelle {TURBINE_A} − nacelle {TURBINE_B} [deg] (dashed = tolerance)")
ax.set_ylabel("10-min records")
ax.set_title(f"After filters: median {dfw['dnac'].median():+.1f}°, n = {len(dfw)}")
save(fig, "09_nacelle_difference_filtered")

# --- 5f. number of samples per class and direction bin ------------------------------------------------
cnt = profile.pivot(index="stability", columns="bin", values="n").reindex(st.CLASS_ORDER)
cmap_seq = LinearSegmentedColormap.from_list("blue_seq", ["#f0efec", "#cde2fb", "#2a78d6", "#0d366b"])
fig, ax = plt.subplots(figsize=(max(6, 0.45 * cnt.shape[1] + 1.5), 2.4))
im = ax.imshow(cnt.to_numpy(dtype=float), cmap=cmap_seq, aspect="auto")
ax.set_xticks(range(cnt.shape[1]), [f"{c:+.0f}" for c in cnt.columns], fontsize=7)
ax.set_yticks(range(3), st.CLASS_ORDER)
ax.grid(False)
vmax = np.nanmax(cnt.to_numpy(dtype=float)) if cnt.size else 1
for i in range(cnt.shape[0]):
    for j in range(cnt.shape[1]):
        n = cnt.iat[i, j]
        ax.text(j, i, f"{int(n)}", ha="center", va="center", fontsize=6,
                color=MUTED if n < MIN_N_BIN else ("white" if n > 0.6 * vmax else INK))
ax.set_xlabel(f"ERA5 direction − {WIND_DIR:.0f}° [deg] (grey numbers: n < {MIN_N_BIN}, hidden in plots)")
ax.set_title("10-min samples per class and direction bin")
fig.colorbar(im, ax=ax, shrink=0.8).set_label("n")
save(fig, "10_samples_per_bin")
plt.show()


# %% 6. Physics sanity checks and save ----------------------------------------------------------
def check(name, ok, detail=""):
    """Print one sanity check: OK, NOT OK, or NOT ENOUGH DATA (ok = None)."""
    state = "NOT ENOUGH DATA" if ok is None else ("OK" if ok else "NOT OK")
    print(f"CHECK {name}: {state}  {detail}")


print("\nSanity checks")
c_all = stats_core.loc["all"]
# 1. B is slower than A in the wake core
check("deficit >= 0 in the wake core (B slower than A)",
      None if c_all["n"] < MIN_N_BIN else bool(c_all["ws_deficit"] >= 0 and c_all["p_deficit"] >= 0),
      f"ws {c_all.get('ws_deficit', np.nan):.3f}, power {c_all.get('p_deficit', np.nan):.3f} "
      "(if NOT OK: check pair order, direction, northing)")

# 2. wake profile peak near the analysed direction
pa = profile_all.dropna(subset=["p_deficit"])
if len(pa):
    peak = pa.loc[pa["p_deficit"].idxmax(), "bin"]
    limit = max(2 * DIR_BIN, CORE_HALF_WIDTH_DEG)
    left = pa.loc[pa["bin"] < 0, "p_deficit"].mean()
    right = pa.loc[pa["bin"] > 0, "p_deficit"].mean()
    check("wake profile peak near the analysed direction", bool(abs(peak) <= limit),
          f"power-deficit peak at {peak:+.0f} deg (limit ±{limit:g}); mean left {left:.3f} / "
          f"right {right:.3f} (roughly equal = symmetric)")
else:
    check("wake profile peak near the analysed direction", None)

# 3. deficit order stable > neutral > unstable (wake core)
enough = all(stats_core.loc[c, "n"] >= MIN_N_BIN for c in st.CLASS_ORDER)
for col in ("ws_deficit", "p_deficit"):
    v = {c: stats_core.loc[c, col] for c in st.CLASS_ORDER}
    check(f"{col}: stable > neutral > unstable",
          bool(v["stable"] > v["neutral"] > v["unstable"]) if enough else None,
          " / ".join(f"{c} {v[c]:.3f}" for c in st.CLASS_ORDER))

# 4. smaller deficit above rated power (wake core)
above = core[core["power_A"] >= 0.95 * P_RATED]
below = core[core["power_A"] < 0.70 * P_RATED]
if len(above) >= MIN_N_BIN and len(below) >= MIN_N_BIN:
    d_above = 1 - above["power_B"].sum() / above["power_A"].sum()
    d_below = 1 - below["power_B"].sum() / below["power_A"].sum()
    check("power deficit smaller near/above rated", bool(d_above < d_below),
          f"P_A >= 95 % rated: {d_above:.3f} (n={len(above)}) vs P_A < 70 % rated: "
          f"{d_below:.3f} (n={len(below)})")
else:
    check("power deficit smaller near/above rated", None,
          f"n above = {len(above)}, n below = {len(below)}")

# 5. SCADA vs ERA5 agreement
for t in ("A", "B"):
    r = verify.iloc[0 if t == "A" else 1]["ws_corr"]
    check(f"SCADA vs ERA5 wind speed correlation {verify.index[0 if t == 'A' else 1]}",
          bool(r >= MIN_CORR), f"r = {r:.2f} (rule of thumb >= {MIN_CORR})")

# 6. nacelle difference centred near 0
check("relative nacelle offset A - B small",
      bool(abs(rel_off["median"]) <= SAME_DIR_TOL_DEG / 2),
      f"median {rel_off['median']:+.1f} deg ({FRAME}, limit ±{SAME_DIR_TOL_DEG / 2:g}); if NOT OK, "
      "the same-direction filter removes many good records - check the northing of A and B")

# 7. enough samples per class
for c in st.CLASS_ORDER:
    n = int(stats_core.loc[c, "n"])
    check(f"samples in the wake core, {c}", bool(n >= MIN_N_BIN), f"n = {n} (min {MIN_N_BIN})")

# --- save tables -----------------------------------------------------------------------------------
settings = pd.Series({
    "turbine_A": TURBINE_A, "turbine_B": TURBINE_B, "wind_dir_deg": WIND_DIR,
    "align_dir_deg": geo["align_dir_deg"], "spacing_D": geo["spacing_D"],
    "dz_hub_m": geo["dz_hub_m"], "dir_tol_deg": DIR_TOL_DEG, "same_dir_tol_deg": SAME_DIR_TOL_DEG,
    "core_half_width_deg": CORE_HALF_WIDTH_DEG, "ws_min": WS_MIN, "ws_max": WS_MAX,
    "p_min_kw": P_MIN_KW, "start": START, "end": END,
    "northing": NORTHING_FILE.name if APPLY_NORTHING else "none",
    "nc_exclusion_groups": nc_attrs.get("exclusion_groups", "not recorded"),
}, name="value")
outputs = {
    "settings": settings.to_frame(),
    "stats_window": stats_window,
    "stats_core": stats_core,
    "profile": profile,
    "deficit_vs_ws": by_ws,
    "scada_vs_era5": verify,
    "filters_era5": era5_report,
    "filters_pair": pair_report,
}
for name, table in outputs.items():
    path = FIG_DIR / f"{TAG}_{name}.csv"
    table.to_csv(path, index=name not in ("profile", "deficit_vs_ws", "filters_era5", "filters_pair"),
                 float_format="%.5g")
    print(f"saved {path}")
