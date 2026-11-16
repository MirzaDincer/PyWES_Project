# -*- coding: utf-8 -*-
"""SCADA cleaning 2024-2025: before vs after comparison plots.

Steps (run cell by cell in Spyder, or the whole file):
  1. Settings, plot style, load the flags table
  2. Power curve before vs after (all turbines or one turbine)
  3. What each filter group removed (small multiples on the power curve)
  4. Rows removed per cleaning step (bar chart)
  5. Records per month before vs after (line chart)
  6. Share of records kept per turbine and month (heatmap)

Input : data/processed/scada_flags_2024_2025.parquet
        (written by scripts/clean_scada_2024_2025.py - run that first)
Output: results/figures/draft/cleaning/*.png

"Before" = all unique 10-min records of 2024-2025 that belong to a known
turbine, i.e. what is left after the structural steps 1-4 (timestamp
alignment, temporal window, duplicates, unknown turbine). Those steps only
remove records that are not part of the period, so they are not "cleaning"
in the physical sense. "After" = records with keep == True.
Records with missing wind speed or power cannot be drawn on the power curve;
they are still counted in the bar chart (cell 4).
"""

# %% 1. Settings, plot style, load data ------------------------------------------
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.colors import LinearSegmentedColormap

from pywes.paths import PROCESSED_DIR, FIGURES_DIR

# --- settings you may change ------------------------------------------------------
FLAGS_FILE = PROCESSED_DIR / "scada_flags_2024_2025.parquet"
TURBINE = None             # e.g. "T01" for one turbine; None = all 21 turbines
MAX_POINTS = 300_000       # random sample per scatter panel (keeps plots fast)
WS_BIN = 0.5               # [m/s] bin width of the median power curve
MIN_N_BIN = 20             # [-] hide median-curve bins with fewer records
RATED_POWER_KW = 2300.0    # SWT-2.3-VS-82
SEED = 42                  # fixed random sample -> same figure every run

FIG_DIR = FIGURES_DIR / "draft" / "cleaning"
FIG_DIR.mkdir(parents=True, exist_ok=True)

# --- plot style (same as analysis/stability_by_sector.py) ----------------------------
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#898781", "#e1e0d9"
BEFORE_COLOR = "#a3a19a"   # neutral gray = reference ("all records")
AFTER_COLOR = "#2a78d6"    # blue = cleaned data (categorical slot 1)
REMOVED_COLOR = "#eb6834"  # orange = removed records (categorical slot 2)
BACKGROUND_COLOR = "#e1e0d9"
plt.rcParams.update({
    "figure.dpi": 110, "savefig.dpi": 300, "savefig.bbox": "tight",
    "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
    "axes.edgecolor": "#c3c2b7", "axes.labelcolor": INK2, "text.color": INK,
    "xtick.color": MUTED, "ytick.color": MUTED,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.spines.top": False, "axes.spines.right": False,
    "legend.frameon": False,
})

# Filter steps grouped by what they catch (names = 'removed_by' values).
STRUCTURAL = ["timestamp_alignment", "temporal_window", "duplicates", "unknown_turbine"]
GROUPS = {
    "Data quality": ["missing_signals", "physical_range", "frozen_sensor"],
    "Downtime, errors & alarms": ["downtime", "turbine_error", "stopping_alarm",
                                  "icing_alarm", "cable_unwind"],
    "Idling & high pitch": ["idling", "high_pitch"],
    "Curtailment & derating": ["power_setpoint_reduced", "low_power_high_wind", "flat_power"],
    "Optional filters": ["yaw_unsteady"],
    "Exclusion periods": ["res_curtailment", "tuneup_deployment", "dynamic_yaw",
                          "res_other", "manual_exclusion"],
}


def save(fig, name):
    """Save a figure as PNG in the draft folder and print the path."""
    path = FIG_DIR / f"{name}.png"
    fig.savefig(path)
    print(f"saved {path}")


def sample(df, n=MAX_POINTS, seed=SEED):
    """Random subset of at most n rows (all rows if fewer)."""
    return df if len(df) <= n else df.sample(n, random_state=seed)


def median_curve(df, ws_bin=WS_BIN, min_n=MIN_N_BIN):
    """Median power [kW] per wind-speed bin [m/s]; bins with < min_n rows dropped."""
    edges = np.arange(0.0, df["ws"].max() + ws_bin, ws_bin)
    centre = pd.cut(df["ws"], edges, labels=edges[:-1] + ws_bin / 2)
    g = df.groupby(centre, observed=True)["power"].agg(["median", "size"])
    g = g[g["size"] >= min_n]
    return g.index.astype(float), g["median"].to_numpy()


def step_order(columns):
    """Filter names in pipeline order, from the 'flag_NN_<name>' columns."""
    flags = sorted(c for c in columns if c.startswith("flag_"))
    return [c.split("_", 2)[2] for c in flags]


# --- load ----------------------------------------------------------------------------
cols = ["time_start", "TurbineName", "ws", "power", "removed_by", "keep"]
flags = pd.read_parquet(FLAGS_FILE)
STEPS = step_order(flags.columns)
df = flags[cols].copy()
df["removed_by"] = df["removed_by"].astype("object")
if TURBINE is not None:
    df = df[df["TurbineName"] == TURBINE]
label = f"turbine {TURBINE}" if TURBINE else "all turbines"
tag = f"_{TURBINE}" if TURBINE else ""

before = df[~df["removed_by"].isin(STRUCTURAL)]
after = before[before["keep"]]
period = (f"{before['time_start'].min():%Y-%m-%d} to "
          f"{before['time_start'].max():%Y-%m-%d}")
print(f"{label}: before {len(before):,d} records, after {len(after):,d} "
      f"({100 * len(after) / len(before):.1f} % kept), {period}")


# %% 2. Power curve before vs after ------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.9), sharex=True, sharey=True)
panels = [("Before cleaning", before, BEFORE_COLOR), ("After cleaning", after, AFTER_COLOR)]
for ax, (title, data, color) in zip(axes, panels):
    pts = sample(data.dropna(subset=["ws", "power"]))
    ax.scatter(pts["ws"], pts["power"], s=1.0, color=color, alpha=0.25,
               linewidths=0, rasterized=True)
    x, y = median_curve(data.dropna(subset=["ws", "power"]))
    ax.plot(x, y, color=INK, linewidth=1.5, label="median per 0.5 m/s bin")
    ax.axhline(RATED_POWER_KW, color=MUTED, linewidth=0.8, linestyle="--")
    shown = f" (sample of {MAX_POINTS:,d})" if len(data) > MAX_POINTS else ""
    ax.set_title(f"{title}: {len(data):,d} records{shown}", loc="left")
    ax.set_xlabel("Nacelle wind speed [m/s]")
axes[0].set_ylabel("Active power [kW]")
axes[0].text(0.3, RATED_POWER_KW + 40, "rated 2300 kW", color=MUTED, fontsize=7, va="bottom")
axes[1].legend(loc="lower right")
axes[0].set_xlim(0, 25)
axes[0].set_ylim(-100, 2600)
fig.suptitle(f"Power curve before and after SCADA cleaning - {label}, {period}",
             x=0.01, ha="left", fontsize=10)
fig.tight_layout()
save(fig, f"00_power_curve_before_after{tag}")


# %% 3. What each filter group removed (small multiples) ---------------------------------
present = {g: s for g, s in GROUPS.items() if before["removed_by"].isin(s).any()}
ncol = 2
nrow = int(np.ceil(len(present) / ncol))
fig, axes = plt.subplots(nrow, ncol, figsize=(9.5, 3.4 * nrow), sharex=True, sharey=True,
                         squeeze=False)
bg = sample(after.dropna(subset=["ws", "power"]), n=MAX_POINTS // 3)
for ax, (group, steps) in zip(axes.flat, present.items()):
    removed = before[before["removed_by"].isin(steps)]
    pts = sample(removed.dropna(subset=["ws", "power"]))
    ax.scatter(bg["ws"], bg["power"], s=1.0, color=BACKGROUND_COLOR, alpha=0.5,
               linewidths=0, rasterized=True, label="kept")
    ax.scatter(pts["ws"], pts["power"], s=1.5, color=REMOVED_COLOR, alpha=0.35,
               linewidths=0, rasterized=True, label="removed by this group")
    pct = 100 * len(removed) / len(before)
    ax.set_title(f"{group}: {len(removed):,d} removed ({pct:.1f} %)", loc="left")
    names = ", ".join(s for s in steps if (removed["removed_by"] == s).any())
    ax.text(0.99, 0.03, names, transform=ax.transAxes, ha="right", va="bottom",
            fontsize=7, color=INK2)
for ax in axes.flat[len(present):]:
    ax.set_visible(False)
for ax in axes[-1]:
    ax.set_xlabel("Nacelle wind speed [m/s]")
for ax in axes[:, 0]:
    ax.set_ylabel("Active power [kW]")
axes[0, 0].set_xlim(0, 25)
axes[0, 0].set_ylim(-100, 2600)
handles = [plt.Line2D([], [], marker="o", linestyle="", color=c, markersize=5)
           for c in (BACKGROUND_COLOR, REMOVED_COLOR)]
fig.legend(handles, ["kept (background)", "removed by this group"], loc="upper right",
           ncol=2, bbox_to_anchor=(0.99, 1.0))
fig.suptitle(f"Records removed by each filter group - {label}, {period}",
             x=0.01, ha="left", fontsize=10)
fig.tight_layout(rect=(0, 0, 1, 0.97))
save(fig, f"01_removed_by_group{tag}")


# %% 4. Rows removed per cleaning step ------------------------------------------------------
counts = before["removed_by"].value_counts()
steps = [s for s in STEPS if s not in STRUCTURAL]
n_removed = np.array([counts.get(s, 0) for s in steps])
pct = 100 * n_removed / len(before)

fig, ax = plt.subplots(figsize=(7.5, 0.32 * len(steps) + 1.2))
y = np.arange(len(steps))[::-1]                 # pipeline order from top to bottom
ax.barh(y, pct, height=0.6, color=AFTER_COLOR, edgecolor="white", linewidth=1.0)
for yi, n, p in zip(y, n_removed, pct):
    ax.text(p + pct.max() * 0.01, yi, f"{n:,d} ({p:.1f} %)", va="center",
            fontsize=7, color=INK2)
ax.set_yticks(y, steps)
ax.set_xlim(0, pct.max() * 1.25 if pct.max() > 0 else 1)
ax.set_xlabel("Records removed [% of records before cleaning]")
ax.grid(axis="y", visible=False)
ax.set_title(f"Records removed per cleaning step (in pipeline order) - {label}\n"
             f"before {len(before):,d}  ->  after {len(after):,d} "
             f"({100 * len(after) / len(before):.1f} % kept)", loc="left")
fig.tight_layout()
save(fig, f"02_removed_per_step{tag}")


# %% 5. Records per month before vs after ---------------------------------------------------
month = before["time_start"].dt.to_period("M")
monthly = pd.DataFrame({
    "before": before.groupby(month).size(),
    "after": after.groupby(after["time_start"].dt.to_period("M")).size(),
})
# full monthly axis: months without any record show as gaps, not as straight lines
full = pd.period_range(monthly.index.min(), monthly.index.max(), freq="M")
monthly = monthly.reindex(full)
monthly.loc[monthly["before"].notna(), "after"] = monthly["after"].fillna(0)
n_turb = 1 if TURBINE else before["TurbineName"].nunique()
expected = pd.Series([p.days_in_month * 144 * n_turb for p in monthly.index], index=monthly.index)
x = monthly.index.to_timestamp()

fig, ax = plt.subplots(figsize=(9.5, 3.6))
ax.plot(x, expected, color=MUTED, linewidth=1.0, linestyle="--", label="expected (full month)")
ax.plot(x, monthly["before"], color=BEFORE_COLOR, linewidth=2.0, marker="o", markersize=4,
        label="before cleaning")
ax.plot(x, monthly["after"], color=AFTER_COLOR, linewidth=2.0, marker="o", markersize=4,
        label="after cleaning")
last = monthly["before"].last_valid_index()
for col, color in (("before", INK2), ("after", INK2)):         # direct labels at the end
    xl = last.to_timestamp()
    ax.annotate(f"{col}: {monthly.loc[last, col]:,.0f}", (xl, monthly.loc[last, col]),
                xytext=(6, 0), textcoords="offset points", va="center", fontsize=7, color=color)
ax.set_ylim(0, expected.max() * 1.1)
ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:,.0f}"))
ax.set_ylabel("10-min records per month")
ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))
ax.legend(loc="lower left", ncol=3)
ax.set_title(f"Records per month before and after cleaning - {label}", loc="left")
fig.autofmt_xdate()
fig.tight_layout()
save(fig, f"03_records_per_month{tag}")


# %% 6. Share kept per turbine and month (heatmap) ---------------------------------------------
if TURBINE is None:
    slot_month = before["time_start"].dt.to_period("M")
    kept_share = (before.assign(month=slot_month)
                  .groupby(["TurbineName", "month"])["keep"].mean().unstack())
    # share of the FULL month that is kept (missing records count as not kept)
    days = np.array([p.days_in_month for p in kept_share.columns])
    kept_full = (after.assign(month=after["time_start"].dt.to_period("M"))
                 .groupby(["TurbineName", "month"]).size().unstack()
                 .reindex(index=kept_share.index, columns=kept_share.columns).fillna(0)
                 / (days * 144) * 100)
    blues = LinearSegmentedColormap.from_list(
        "blues", ["#cde2fb", "#86b6ef", "#3987e5", "#256abf", "#184f95", "#0d366b"])

    fig, ax = plt.subplots(figsize=(max(6.0, 0.38 * kept_full.shape[1] + 2.5), 5.6))
    im = ax.imshow(kept_full.to_numpy(), aspect="auto", cmap=blues, vmin=0, vmax=100)
    ax.set_xticks(range(kept_full.shape[1]), [p.strftime("%Y-%m") for p in kept_full.columns],
                  rotation=90)
    ax.set_yticks(range(kept_full.shape[0]), kept_full.index)
    ax.grid(False)
    cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    cb.set_label("Records kept [% of a full month]")
    cb.outline.set_visible(False)
    ax.set_title("Share of each month kept after cleaning, per turbine\n"
                 "(missing raw records count as not kept)", loc="left")
    fig.tight_layout()
    save(fig, "04_kept_share_turbine_month")

    # table view for the paper / for checking single cells
    kept_full.round(1).to_csv(FIG_DIR / "04_kept_share_turbine_month.csv")
    print(f"saved {FIG_DIR / '04_kept_share_turbine_month.csv'}")

plt.show()
