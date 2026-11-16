# -*- coding: utf-8 -*-
"""Site map of Hill of Towie: turbine positions, coordinate grid and orography.

Steps (run cell by cell in Spyder, or the whole file):
  1. Settings and turbine coordinates (lat/lon -> UTM zone 30N)
  2. Download (once) and read the elevation tile(s)
  3. Cut the map area, compute hillshade, sample ground height at each turbine
  4. Sanity checks
  5. Plot the map and save the turbine table

Inputs : data/sample/Hill_of_Towie_turbine_metadata.csv
         data/raw/dem/N57W004.hgt.gz  (downloaded automatically, ~10-25 MB)
Outputs: results/figures/draft/site/hill_of_towie_site_map.png (+ .pdf)
         data/processed/turbine_layout_utm.csv

Elevation data
  "Terrain Tiles" (Mapzen / Linux Foundation), Skadi format = SRTM-style .hgt,
  1 arc-second (~31 m north-south, ~17 m east-west at 57.5 N), WGS84 lat/lon.
  Free, no account needed. How to cite: "Terrain Tiles was accessed on DATE from
  https://registry.opendata.aws/terrain-tiles". The tiles are merged from several
  open sources; for the exact source and attribution in this area see
  https://github.com/tilezen/joerd/blob/master/docs/attribution.md
  Note: satellite DEMs can include part of the forest canopy height.

Coordinates
  UTM zone 30N (EPSG:32630), WGS84. The central meridian of zone 30 is 3 deg W,
  which passes through the farm, so grid north = true north here (difference
  < 0.1 deg). UTM formula: Krueger series, see Karney (2011), J. Geodesy 85,
  475-485 (same method as PROJ). Accuracy far below 1 mm, so no pyproj needed.
"""

# %% 1. Settings and turbine coordinates ---------------------------------------
import gzip
import urllib.request

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import FuncFormatter, MultipleLocator

from pywes.paths import RAW_DIR, PROCESSED_DIR, SAMPLE_DIR, FIGURES_DIR

# --- settings you may change --------------------------------------------------
MARGIN = 1500.0          # [m] map border around the outer turbines
GRID_STEP = 500.0        # [m] UTM grid line spacing
CONTOUR_STEP = 10.0      # [m] contour line interval
LABEL_STEP = 50.0        # [m] labelled (thicker) contour interval
GRATICULE_LAT = 0.01     # [deg] latitude line spacing (~1.1 km)
GRATICULE_LON = 0.02     # [deg] longitude line spacing (~1.2 km at 57.5 N)
SHOW_ROTOR_CIRCLE = False  # True: draw rotor diameter to scale (D = 82 m)

DEM_DIR = RAW_DIR / "dem"
DEM_URL = "https://s3.amazonaws.com/elevation-tiles-prod/skadi/{lat_dir}/{name}.hgt.gz"
FIG_DIR = FIGURES_DIR / "draft" / "site"
FIG_DIR.mkdir(parents=True, exist_ok=True)
DEM_DIR.mkdir(parents=True, exist_ok=True)
PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

# --- plot style (same as the other analysis scripts) ---------------------------
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#898781", "#e1e0d9"
plt.rcParams.update({
    "figure.dpi": 110, "savefig.dpi": 300, "savefig.bbox": "tight",
    "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
    "axes.edgecolor": "#c3c2b7", "axes.labelcolor": INK2, "text.color": INK,
    "xtick.color": MUTED, "ytick.color": MUTED, "legend.frameon": False,
})
# one-hue sequential ramp for height (light = low, dark = high); no blue, so it
# is not mistaken for water
ELEV_CMAP = LinearSegmentedColormap.from_list(
    "elev", ["#f7f2e8", "#e8d6b8", "#d1b084", "#b0845a", "#865c3a", "#5c3b22"])


def latlon_to_utm(lat, lon, zone=30):
    """WGS84 lat/lon [deg] -> UTM easting, northing [m] (northern hemisphere).

    Krueger series to 6th order (Karney 2011). Works on scalars and arrays.
    """
    a, f = 6378137.0, 1 / 298.257223563          # WGS84 ellipsoid
    k0, e0 = 0.9996, 500000.0                    # UTM scale and false easting
    lon0 = np.radians(-183.0 + 6.0 * zone)       # central meridian (zone 30: -3 deg)
    n = f / (2 - f)
    A = a / (1 + n) * (1 + n**2 / 4 + n**4 / 64 + n**6 / 256)
    alpha = [n / 2 - 2 * n**2 / 3 + 5 * n**3 / 16 + 41 * n**4 / 180
             - 127 * n**5 / 288 + 7891 * n**6 / 37800,
             13 * n**2 / 48 - 3 * n**3 / 5 + 557 * n**4 / 1440
             + 281 * n**5 / 630 - 1983433 * n**6 / 1935360,
             61 * n**3 / 240 - 103 * n**4 / 140 + 15061 * n**5 / 26880
             + 167603 * n**6 / 181440,
             49561 * n**4 / 161280 - 179 * n**5 / 168 + 6601661 * n**6 / 7257600,
             34729 * n**5 / 80640 - 3418889 * n**6 / 1995840,
             212378941 * n**6 / 319334400]
    phi = np.radians(np.asarray(lat, dtype=float))
    lam = np.radians(np.asarray(lon, dtype=float)) - lon0
    e = np.sqrt(f * (2 - f))
    t = np.sinh(np.arctanh(np.sin(phi)) - e * np.arctanh(e * np.sin(phi)))
    xi_p = np.arctan2(t, np.cos(lam))
    eta_p = np.arctanh(np.sin(lam) / np.sqrt(1 + t**2))
    xi, eta = xi_p.copy(), eta_p.copy()
    for j, aj in enumerate(alpha, start=1):
        xi = xi + aj * np.sin(2 * j * xi_p) * np.cosh(2 * j * eta_p)
        eta = eta + aj * np.cos(2 * j * xi_p) * np.sinh(2 * j * eta_p)
    return e0 + k0 * A * eta, k0 * A * xi


turbines = pd.read_csv(SAMPLE_DIR / "Hill_of_Towie_turbine_metadata.csv", encoding="utf-8-sig")
turbines["x_utm"], turbines["y_utm"] = latlon_to_utm(turbines["Latitude"], turbines["Longitude"])
ROTOR_D = float(turbines["Rotor Diameter (m)"].iloc[0])     # 82 m
HUB_HEIGHT = float(turbines["Hub Height (m)"].iloc[0])      # 59 m
print(f"{len(turbines)} turbines, D = {ROTOR_D} m, hub height = {HUB_HEIGHT} m")

# map area in UTM [m], rounded out to whole grid steps
x_min = np.floor((turbines["x_utm"].min() - MARGIN) / GRID_STEP) * GRID_STEP
x_max = np.ceil((turbines["x_utm"].max() + MARGIN) / GRID_STEP) * GRID_STEP
y_min = np.floor((turbines["y_utm"].min() - MARGIN) / GRID_STEP) * GRID_STEP
y_max = np.ceil((turbines["y_utm"].max() + MARGIN) / GRID_STEP) * GRID_STEP

# lat/lon box that covers the UTM frame (local linear conversion around the farm
# centre + a small pad; checked by an assert in cell 4)
lat_c, lon_c = turbines["Latitude"].mean(), turbines["Longitude"].mean()
x_c, y_c = latlon_to_utm(lat_c, lon_c)
dlat_per_m = 1 / 111.2e3                                       # [deg/m]
dlon_per_m = 1 / (111.32e3 * np.cos(np.radians(lat_c)))        # [deg/m]
lat_lo, lat_hi = lat_c + (y_min - y_c) * dlat_per_m - 0.002, lat_c + (y_max - y_c) * dlat_per_m + 0.002
lon_lo, lon_hi = lon_c + (x_min - x_c) * dlon_per_m - 0.003, lon_c + (x_max - x_c) * dlon_per_m + 0.003


# %% 2. Download and read the elevation tile(s) --------------------------------
def tile_name(lat_floor, lon_floor):
    """SRTM tile name from its south-west corner, e.g. (57, -4) -> 'N57W004'."""
    ns = "N" if lat_floor >= 0 else "S"
    ew = "E" if lon_floor >= 0 else "W"
    return f"{ns}{abs(lat_floor):02d}{ew}{abs(lon_floor):03d}"


def read_hgt_tile(lat_floor, lon_floor):
    """Download (once) and read one 1x1 deg .hgt tile.

    Returns elevation [m] as float array (row 0 = north edge, NaN for voids)
    plus the latitude and longitude of the rows and columns [deg].
    """
    name = tile_name(lat_floor, lon_floor)
    path = DEM_DIR / f"{name}.hgt.gz"
    if not path.exists():
        url = DEM_URL.format(lat_dir=name[:3], name=name)
        print(f"downloading {url}")
        try:
            urllib.request.urlretrieve(url, path)
        except Exception as err:
            if path.exists():
                path.unlink()          # no half-downloaded file
            raise RuntimeError(f"Download failed ({err}). Download the file by "
                               f"hand from {url} and save it as {path}") from err
    with gzip.open(path, "rb") as fh:
        raw = np.frombuffer(fh.read(), dtype=">i2")       # big-endian int16
    size = int(round(np.sqrt(raw.size)))                  # 3601 for 1 arc-second
    if size * size != raw.size:
        raise ValueError(f"{path} is not a square .hgt grid ({raw.size} values)")
    z = raw.reshape(size, size).astype(float)
    z[z <= -32768] = np.nan                               # SRTM void value
    lats = lat_floor + 1 - np.arange(size) / (size - 1)   # north -> south
    lons = lon_floor + np.arange(size) / (size - 1)       # west -> east
    print(f"read {name}: {size} x {size} points")
    return z, lats, lons


def dem_for_box(lat_lo, lat_hi, lon_lo, lon_hi):
    """Elevation grid [m] for a lat/lon box; joins tiles if the box crosses a tile edge."""
    rows = []
    for lat_f in range(int(np.floor(lat_hi)), int(np.floor(lat_lo)) - 1, -1):   # north first
        cols = []
        for lon_f in range(int(np.floor(lon_lo)), int(np.floor(lon_hi)) + 1):
            z, lats, lons = read_hgt_tile(lat_f, lon_f)
            i = (lats >= lat_lo) & (lats <= lat_hi)
            j = (lons >= lon_lo) & (lons <= lon_hi)
            cols.append((z[np.ix_(i, j)], lats[i], lons[j]))
        # neighbouring tiles share their edge column -> drop the duplicate
        z_row = np.hstack([c[0] if k == 0 else c[0][:, 1:] for k, c in enumerate(cols)])
        lon_row = np.concatenate([c[2] if k == 0 else c[2][1:] for k, c in enumerate(cols)])
        rows.append((z_row, cols[0][1], lon_row))
    z = np.vstack([r[0] if k == 0 else r[0][1:] for k, r in enumerate(rows)])
    lats = np.concatenate([r[1] if k == 0 else r[1][1:] for k, r in enumerate(rows)])
    return z, lats, rows[0][2]


elev, lat_vec, lon_vec = dem_for_box(lat_lo, lat_hi, lon_lo, lon_hi)
LON2D, LAT2D = np.meshgrid(lon_vec, lat_vec)
X2D, Y2D = latlon_to_utm(LAT2D, LON2D)     # UTM of every DEM point [m]


# %% 3. Hillshade and ground height at the turbines ----------------------------
def hillshade(z, dx, dy, azimuth=315.0, altitude=45.0):
    """Hillshade (0 = dark, 1 = bright) for a grid with spacing dx (east), dy (north) [m].

    Light from the north-west (315 deg), the usual cartographic choice.
    z must have row 0 at the north edge.
    """
    dz_dy, dz_dx = np.gradient(z, -dy, dx)          # -dy: rows go north -> south
    slope = np.arctan(np.hypot(dz_dx, dz_dy))
    aspect = np.arctan2(-dz_dx, -dz_dy)              # direction the slope faces
    az, alt = np.radians(azimuth), np.radians(altitude)
    hs = np.sin(alt) * np.cos(slope) + np.cos(alt) * np.sin(slope) * np.cos(az - aspect)
    return np.clip(hs, 0, 1)


dy_m = np.radians(abs(lat_vec[1] - lat_vec[0])) * 6371e3                 # ~31 m
dx_m = np.radians(abs(lon_vec[1] - lon_vec[0])) * 6371e3 * np.cos(np.radians(lat_vec.mean()))
shade = hillshade(np.nan_to_num(elev, nan=np.nanmean(elev)), dx_m, dy_m)

# ground height at each turbine: bilinear interpolation on the lat/lon grid
from scipy.interpolate import RegularGridInterpolator
interp = RegularGridInterpolator((lat_vec[::-1], lon_vec), elev[::-1, :])
turbines["ground_elev_m"] = interp(np.c_[turbines["Latitude"], turbines["Longitude"]])
turbines["hub_elev_m"] = turbines["ground_elev_m"] + HUB_HEIGHT


# %% 4. Sanity checks ------------------------------------------------------------
inside = (X2D >= x_min) & (X2D <= x_max) & (Y2D >= y_min) & (Y2D <= y_max)
n_void = int(np.isnan(elev[inside]).sum())
z_in = elev[inside]
print(f"map area: {(x_max - x_min) / 1e3:.1f} km x {(y_max - y_min) / 1e3:.1f} km, "
      f"{inside.sum()} DEM points, {n_void} voids")
print(f"terrain in map area: {np.nanmin(z_in):.0f} ... {np.nanmax(z_in):.0f} m a.s.l.")
print(f"ground at turbines : {turbines['ground_elev_m'].min():.0f} ... "
      f"{turbines['ground_elev_m'].max():.0f} m a.s.l. "
      f"(spread {np.ptp(turbines['ground_elev_m']):.0f} m)")
# the DEM box must really cover the whole map frame
assert X2D.min() <= x_min and X2D.max() >= x_max and Y2D.min() <= y_min and Y2D.max() >= y_max, \
    "DEM box does not cover the map frame - increase the lat/lon box"
assert turbines["ground_elev_m"].notna().all(), "DEM void at a turbine position"
assert 0 < np.nanmin(z_in) and np.nanmax(z_in) < 1500, "elevation out of plausible range"
# spacing between nearest neighbours, in rotor diameters (layout check)
dxy = np.hypot(turbines["x_utm"].values[:, None] - turbines["x_utm"].values[None, :],
               turbines["y_utm"].values[:, None] - turbines["y_utm"].values[None, :])
np.fill_diagonal(dxy, np.inf)
nn = dxy.min(axis=1) / ROTOR_D
print(f"nearest-neighbour spacing: {nn.min():.1f} ... {nn.max():.1f} D (median {np.median(nn):.1f} D)")


# %% 5. Plot the map ---------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7.2, 6.0 * (y_max - y_min) / (x_max - x_min) + 0.8))

# orography: filled height bands, hillshade on top, then contour lines
z_lo = np.floor(np.nanmin(z_in) / CONTOUR_STEP) * CONTOUR_STEP
z_hi = np.ceil(np.nanmax(z_in) / CONTOUR_STEP) * CONTOUR_STEP
levels = np.arange(z_lo, z_hi + CONTOUR_STEP, CONTOUR_STEP)
cf = ax.contourf(X2D, Y2D, elev, levels=levels, cmap=ELEV_CMAP, zorder=1)
ax.pcolormesh(X2D, Y2D, shade, cmap="gray", vmin=0, vmax=1, alpha=0.28,
              shading="nearest", antialiased=False, linewidth=0, zorder=2, rasterized=True)
minor = [lv for lv in levels if lv % LABEL_STEP != 0]
major = [lv for lv in levels if lv % LABEL_STEP == 0]
ax.contour(X2D, Y2D, elev, levels=minor, colors="#6b5a48", linewidths=0.25,
           alpha=0.6, zorder=3)
cs = ax.contour(X2D, Y2D, elev, levels=major, colors="#4a3b2c", linewidths=0.6, zorder=3)
ax.clabel(cs, fmt="%.0f", fontsize=6.5, inline_spacing=2)

# lat/lon graticule (thin dashed lines, labelled at the top and right edges)
for lat_g in np.arange(np.ceil(lat_lo / GRATICULE_LAT) * GRATICULE_LAT, lat_hi, GRATICULE_LAT):
    lon_line = np.linspace(lon_lo, lon_hi, 100)
    xg, yg = latlon_to_utm(np.full_like(lon_line, lat_g), lon_line)
    ax.plot(xg, yg, color=INK2, lw=0.4, ls=(0, (4, 3)), alpha=0.7, zorder=4)
    if y_min < yg[-1] < y_max:
        ax.annotate(f"{lat_g:.2f}°N", xy=(1, (np.interp(x_max, xg, yg) - y_min) / (y_max - y_min)),
                    xycoords="axes fraction", xytext=(3, 0), textcoords="offset points",
                    va="center", ha="left", fontsize=7, color=INK2)
for lon_g in np.arange(np.ceil(lon_lo / GRATICULE_LON) * GRATICULE_LON, lon_hi, GRATICULE_LON):
    lat_line = np.linspace(lat_lo, lat_hi, 100)
    xg, yg = latlon_to_utm(lat_line, np.full_like(lat_line, lon_g))
    ax.plot(xg, yg, color=INK2, lw=0.4, ls=(0, (4, 3)), alpha=0.7, zorder=4)
    x_top = np.interp(y_max, yg, xg)
    if x_min < x_top < x_max:
        ax.annotate(f"{abs(lon_g):.2f}°W", xy=((x_top - x_min) / (x_max - x_min), 1),
                    xycoords="axes fraction", xytext=(0, 3), textcoords="offset points",
                    va="bottom", ha="center", fontsize=7, color=INK2)

# turbines
ax.scatter(turbines["x_utm"], turbines["y_utm"], s=34, marker="^", color=INK,
           edgecolor="white", linewidth=0.8, zorder=6, label="Turbine (Siemens SWT-2.3-82)")
if SHOW_ROTOR_CIRCLE:
    for xt, yt in zip(turbines["x_utm"], turbines["y_utm"]):
        ax.add_patch(plt.Circle((xt, yt), ROTOR_D / 2, fill=False, lw=0.6, color=INK, zorder=6))
# label position per turbine (default: centred above the marker). Add an entry
# here (e.g. "T05": "below") if a label covers a neighbour.
LABEL_POS = {}
OFFSETS = {"right": ((4, 3), "left", "bottom"), "left": ((-4, 3), "right", "bottom"),
           "below": ((0, -6), "center", "top"), "above": ((0, 6), "center", "bottom")}
halo = [pe.withStroke(linewidth=2.2, foreground="white")]
for _, row in turbines.iterrows():
    off, ha, va = OFFSETS[LABEL_POS.get(row["Turbine Name"], "above")]
    ax.annotate(row["Turbine Name"], (row["x_utm"], row["y_utm"]), xytext=off,
                textcoords="offset points", ha=ha, va=va, fontsize=6.5, color=INK,
                zorder=7, path_effects=halo)

# UTM grid, axes in km
ax.set_xlim(x_min, x_max)
ax.set_ylim(y_min, y_max)
ax.set_aspect("equal")
ax.xaxis.set_major_locator(MultipleLocator(GRID_STEP))
ax.yaxis.set_major_locator(MultipleLocator(GRID_STEP))
km = FuncFormatter(lambda v, _: f"{v / 1e3:.1f}")
ax.xaxis.set_major_formatter(km)
ax.yaxis.set_major_formatter(km)
ax.grid(True, color="white", lw=0.6, alpha=0.8, zorder=5)
ax.set_axisbelow(False)
ax.tick_params(labelsize=7)
plt.setp(ax.get_xticklabels(), rotation=90)
ax.set_xlabel("UTM 30N easting [km]")
ax.set_ylabel("UTM 30N northing [km]")

# north arrow and 1 km scale bar (lower left)
ax.annotate("N", xy=(0.05, 0.95), xytext=(0.05, 0.87), xycoords="axes fraction",
            ha="center", va="center", fontsize=9, fontweight="bold",
            arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.2), zorder=8)
sb_x, sb_y = x_min + 0.04 * (x_max - x_min), y_min + 0.04 * (y_max - y_min)
ax.plot([sb_x, sb_x + 1000], [sb_y, sb_y], color=INK, lw=2.5, solid_capstyle="butt", zorder=8)
ax.text(sb_x + 500, sb_y + 0.012 * (y_max - y_min), "1 km", ha="center", va="bottom",
        fontsize=7.5, zorder=8)

cbar = fig.colorbar(cf, ax=ax, shrink=0.7, pad=0.1)
cbar.set_label("Terrain height [m a.s.l.]")
cbar.ax.tick_params(labelsize=7)
cbar.outline.set_visible(False)
ax.legend(loc="lower right", fontsize=7.5, frameon=True, facecolor="white",
          edgecolor="none", framealpha=0.8)
ax.set_title("Hill of Towie wind farm: layout and terrain", loc="left", pad=16)
ax.annotate(f"Contours every {CONTOUR_STEP:.0f} m (labelled every {LABEL_STEP:.0f} m). "
         "Grid: UTM zone 30N, WGS84; dashed: lat/lon.\n"
         "Elevation: Terrain Tiles (registry.opendata.aws/terrain-tiles), 1 arc-second.",
            xy=(0, 0), xycoords="axes fraction", xytext=(0, -48), textcoords="offset points",
            fontsize=6.5, color=MUTED, ha="left", va="top")

for ext in ("png", "pdf"):
    path = FIG_DIR / f"hill_of_towie_site_map.{ext}"
    fig.savefig(path)
    print(f"saved {path}")
plt.show()

# turbine table for later use (e.g. FLORIS layout: x = x_utm - x_utm.mean(), y likewise)
cols = ["Turbine Name", "Latitude", "Longitude", "x_utm", "y_utm", "ground_elev_m", "hub_elev_m"]
out = turbines[cols].round({"x_utm": 1, "y_utm": 1, "ground_elev_m": 1, "hub_elev_m": 1})
out_path = PROCESSED_DIR / "turbine_layout_utm.csv"
out.to_csv(out_path, index=False)
print(f"saved {out_path}")
print(out.to_string(index=False))
