#!/usr/bin/env python3
"""ENA sky maps in ecliptic J2000 coordinates.

    python3 enacpp/tools/ena_plot.py OUT/flux_band13.h5                      # nose-centred (UPWIND, V1, V2)
    python3 enacpp/tools/ena_plot.py OUT/flux_band13.h5 --view tail          # tail-centred (DOWNWIND)
    python3 enacpp/tools/ena_plot.py OUT/flux_band*.h5 --config enacpp/configs/plot_settings.toml
    python3 enacpp/tools/ena_plot.py OUT/flux_band13.h5 --scale 1.8 --vmin 0 --vmax 12
    python3 enacpp/tools/ena_plot.py OUT/flux_band03_2010-2011.h5 --component heliosheath_pui --log
    python3 enacpp/tools/ena_plot.py run_b/flux_band13.h5 --ratio run_a/flux_band13.h5
    python3 enacpp/tools/ena_plot.py OUT/flux_band13.h5 --write-ecliptic     # all components on a 6-deg ecliptic grid

Every option can be set in a plot settings file (--config; configs/plot_settings.toml
lists them all); command-line options override the file. Input: the flux files of ena_run
(flux_bandNN.h5) or of the time-dependent method (flux_bandNN_<window>.h5); one figure per
file. The legacy count-rate files (...<TAG>CountRateKappa<k>EbandNN.dat) are read as well,
turned back into flux (count rate / E_central / G).

Coordinates
-----------
The MHD solution is in the simulation frame: z along the solar rotation axis, x 5 deg
above the interstellar flow (in the solar equatorial plane), y completing the right-handed
system. The maps are drawn in ecliptic J2000. The rotation is fixed by two directions known
in both frames (TRIAD): the solar rotation axis (kept exactly; inclination 7.25 deg,
ascending node 75.76 deg, J2000) and the interstellar flow ((26.3, 0.3, -2.3) km/s in the
simulation; ecliptic longitude 75.4 deg, latitude -5.2 deg). The two agree to 0.1 deg; the
upwind direction lands at (255.4, +5.2), the downwind direction at (75.4, -5.2).

Cells
-----
Grid values are cell centred (theta = i dtheta +- dtheta/2, phi = j dphi +- dphi/2). Every
cell is rotated and drawn as its own polygon, without resampling: the boxes of the regions
with missing flux outline exactly the flagged cells. A 3-deg cap around the simulation
south pole has no cell and stays blank. --idl uses the placement of the original IDL plotting
code instead (no rotation, a longitude shift only), for comparison with old figures.

Map
---
Mollweide projection (--projection hammer: Hammer), longitude increasing to the left. The
map is centred on the upwind direction (--view nose) or the downwind direction (--view
tail), including its latitude (+5.2 / -5.2 deg) - so the ecliptic equator appears slightly
curved; --view ecliptic centres the map on ecliptic (0, 0); --center LON LAT sets any centre
(e.g. --center 255 0). Linear colour scale (--log for a log scale) with 4 divisions; limits
from the map, or --vmin/--vmax, or per band in the settings file, or one range for all maps
of a band (--common-range). Values outside the limits take the end colours. --scale
multiplies the flux (the title then says "(Flux x S)"). Regions with missing flux are grey
with a black box per cell. Markers: UPWIND, V1, V2 on nose maps, DOWNWIND on tail maps.
--limits LONMIN LONMAX LATMIN LATMAX shows only a part of the sky (degrees relative to the
map centre, longitude positive to the left).

Needs numpy, h5py and matplotlib (pip install numpy h5py matplotlib).
"""

import argparse
import ast
import glob
import os
import re
import sys

import numpy as np

try:
    import h5py
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection, LineCollection
    from matplotlib.colors import LogNorm, Normalize
    from matplotlib.path import Path
    from matplotlib.patches import Ellipse, PathPatch
    from matplotlib.ticker import FuncFormatter
except ImportError:
    sys.exit("ena_plot.py needs numpy, h5py and matplotlib (pip install numpy h5py matplotlib)")

DEG = np.pi / 180.0

# ---------------------------------------------------------------------------------
# Settings: defaults < settings file (--config) < command line
# ---------------------------------------------------------------------------------
DEFAULTS = {
    # [frame]
    "flow_sim": [26.3, 0.3, -2.3],        # km/s, downstream ISM flow in the simulation frame
    "flow_ecliptic": [75.4, -5.2],        # deg, the same direction in ecliptic J2000
    "solar_axis": [7.25, 75.76],          # deg, inclination of the solar equator and its ascending node (J2000)
    # [plot]
    "component": "total",
    "view": "nose",                       # nose (centre on upwind), tail (on downwind), ecliptic (on lon 0, lat 0)
    "center": None,                       # [lon, lat] of the map centre (overrides view)
    "projection": "mollweide",            # mollweide or hammer
    "idl": False,                         # placement of the original IDL plots (no rotation)
    "scale": 1.0,                         # flux multiplied by this
    "log": False,                         # log colour scale
    "vmin": None, "vmax": None,           # colour limits for every band
    "common_range": False,                # one range per band over all input files
    "legend_style": "auto",               # colour-bar numbers: "integer" (e.g. 190), "float" (to the
                                          #   nearest hundredth, 0.85) or "exponent" (mantissa to the nearest hundredth,
                                          #   9.00e-03); "auto" picks integer for maxima >= 10, float for >= 0.1,
                                          #   exponent below
    "cmap": "jet",                        # matplotlib colour map
    "limits": None,                       # [lonmin, lonmax, latmin, latmax] relative to the centre
    "boxes": True,                        # boxes over regions with missing flux
    "grid": True,
    "grid_dlon": 30.0, "grid_dlat": 30.0,
    "edge_meridians": 4.0,                # extra meridians this far either side of the map edge (0: none)
    "lon_labels": "signed",               # signed (-180..180) or 0-360
    "title": None,
    "units_label": "ENA flux (cm$^{-2}$ s$^{-1}$ sr$^{-1}$ keV$^{-1}$)",
    "markers_nose": ["UPWIND", "V1", "V2"],
    "markers_tail": ["DOWNWIND"],
    "markers_ecliptic": ["UPWIND", "DOWNWIND", "V1", "V2"],
    "marker_positions": {"V1": [254.0, 35.0], "V2": [289.0, -32.0]},   # ecliptic deg
    # band -> [vmin, vmax]: maxima for bands 9-17, minima 0
    "band_limits": {9: [0.0, 265.0], 10: [0.0, 140.0], 11: [0.0, 60.0], 12: [0.0, 30.0], 13: [0.0, 12.0],
                    14: [0.0, 0.90], 15: [0.0, 0.08], 16: [0.0, 0.025], 17: [0.0, 0.009]},
    "nsub": 4,                            # points per cell edge
    # [output]
    "format": "png", "dpi": 200, "out": None, "figsize": [9.0, 6.0], "fontsize": 12,
    "write_ecliptic": False, "ecliptic_res": 6.0,
}


def parse_settings(path):
    """Plot settings file: TOML-like ([section], key = value, strings, numbers, booleans, arrays).
    Read without a TOML library (works with any Python 3)."""
    out, section = {}, ""
    for n, raw in enumerate(open(path), 1):
        line = re.sub(r'#(?=(?:[^"]*"[^"]*")*[^"]*$).*', "", raw).strip()  # comments outside quotes
        if not line:
            continue
        if line.startswith("["):
            section = line.strip("[] ")
            continue
        if "=" not in line:
            sys.exit(f"{path}:{n}: expected key = value")
        key, val = (x.strip() for x in line.split("=", 1))
        key = key.strip('"')
        val = re.sub(r"\btrue\b", "True", re.sub(r"\bfalse\b", "False", val))
        try:
            value = ast.literal_eval(val)
        except (ValueError, SyntaxError):
            sys.exit(f"{path}:{n}: cannot read the value '{val}'")
        if section in ("markers",):
            out.setdefault("marker_positions", {})[key] = list(value)
        elif section == "band_limits":
            out.setdefault("band_limits", {})[int(key)] = list(value)
        else:
            out[key] = value
    return out


# ---------------------------------------------------------------------------------
# Frames
# ---------------------------------------------------------------------------------
def unit(lon_deg, lat_deg):
    lon, lat = np.asarray(lon_deg) * DEG, np.asarray(lat_deg) * DEG
    return np.stack([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)], axis=-1)


def lonlat(v):
    """Longitude in [0, 360) and latitude (deg) of vector(s) v [..., 3]."""
    v = np.asarray(v, dtype=float)
    r = np.linalg.norm(v, axis=-1)
    lon = np.degrees(np.arctan2(v[..., 1], v[..., 0])) % 360.0
    lat = np.degrees(np.arcsin(np.clip(v[..., 2] / r, -1.0, 1.0)))
    return lon, lat


def solar_pole_ecliptic(inclination, node):
    """Pole of the solar equator in ecliptic coordinates: 90 deg behind the ascending node."""
    return unit(node - 90.0, 90.0 - inclination)


def sim_to_ecliptic(flow_sim=DEFAULTS["flow_sim"], flow_ecliptic=DEFAULTS["flow_ecliptic"],
                    solar_axis=DEFAULTS["solar_axis"]):
    """Rotation matrix R with v_ecliptic = R @ v_simulation (TRIAD: solar axis exact, flow fixes the rest)."""
    def triad(z, f):
        z = z / np.linalg.norm(z)
        x = f - (f @ z) * z
        x /= np.linalg.norm(x)
        return np.vstack([x, np.cross(z, x), z])
    f_sim = np.asarray(flow_sim, float) / np.linalg.norm(flow_sim)
    t_sim = triad(np.array([0.0, 0.0, 1.0]), f_sim)
    t_ecl = triad(solar_pole_ecliptic(*solar_axis), unit(*flow_ecliptic))
    return t_ecl.T @ t_sim


def frame_check(flow_sim=DEFAULTS["flow_sim"], flow_ecliptic=DEFAULTS["flow_ecliptic"], solar_axis=DEFAULTS["solar_axis"]):
    """Angle between the flow and the solar axis in both frames (deg): they should agree."""
    f_sim = np.asarray(flow_sim, float) / np.linalg.norm(flow_sim)
    a_sim = np.degrees(np.arccos(f_sim[2]))
    a_ecl = np.degrees(np.arccos(np.clip(unit(*flow_ecliptic) @ solar_pole_ecliptic(*solar_axis), -1, 1)))
    return a_sim, a_ecl


def view_matrix(center_lon, center_lat):
    """Rotation of ecliptic vectors that brings the map centre to (0, 0) (oblique map)."""
    a, b = center_lon * DEG, center_lat * DEG
    rz = np.array([[np.cos(a), np.sin(a), 0], [-np.sin(a), np.cos(a), 0], [0, 0, 1]])
    ry = np.array([[np.cos(b), 0, np.sin(b)], [0, 1, 0], [-np.sin(b), 0, np.cos(b)]])
    return ry @ rz


def to_display(v_ecl, V):
    """Display longitude/latitude (deg) of ecliptic vectors: longitude relative to the map centre in
    (-180, 180], positive to the east; the projections draw east to the LEFT (sky-map convention)."""
    lon, lat = lonlat(v_ecl @ V.T)
    return ((lon + 180.0) % 360.0) - 180.0, lat


# ---------------------------------------------------------------------------------
# Projections (x, y in units of the radius)
# ---------------------------------------------------------------------------------
def mollweide(lon, lat):
    """Mollweide projection; lon/lat in radians, lon = display longitude (copies beyond +-pi are clipped)."""
    lat = np.asarray(lat, float)
    t = lat.copy()
    target = np.pi * np.sin(lat)
    for _ in range(50):  # Newton for 2t + sin 2t = pi sin(lat)
        f = 2 * t + np.sin(2 * t) - target
        df = 2 + 2 * np.cos(2 * t)
        step = np.where(df > 1e-12, f / np.where(df > 1e-12, df, 1.0), 0.0)
        t = t - step
        if np.all(np.abs(step) < 1e-12):
            break
    t = np.where(np.abs(lat) >= np.pi / 2 - 1e-12, np.sign(lat) * np.pi / 2, t)
    return -2 * np.sqrt(2) / np.pi * np.asarray(lon) * np.cos(t), np.sqrt(2) * np.sin(t)


def hammer(lon, lat):
    """Hammer-Aitoff projection; lon/lat in radians (copies beyond +-pi are pushed outside the map and clipped)."""
    lon, lat = np.asarray(lon, float), np.asarray(lat, float)
    lc = np.clip(lon, -np.pi, np.pi)
    d = np.sqrt(1 + np.cos(lat) * np.cos(lc / 2))
    x = -2 * np.sqrt(2) * np.cos(lat) * np.sin(lc / 2) / d
    return np.where(np.abs(lon) > np.pi + 1e-9, -np.sign(lon) * 3.0, x), np.sqrt(2) * np.sin(lat) / d


PROJECTIONS = {"mollweide": mollweide, "hammer": hammer}
# Note: display longitude is positive to the left, so the projections return x = -(...).


# ---------------------------------------------------------------------------------
# Cells
# ---------------------------------------------------------------------------------
def cell_outline(theta0, theta1, phi0, phi1, nsub):
    """Boundary of a theta/phi cell (deg), nsub points per edge, counter-clockwise."""
    t = np.linspace(theta0, theta1, nsub + 1)
    p = np.linspace(phi0, phi1, nsub + 1)
    th = np.concatenate([t, np.full(nsub, theta1), t[::-1][1:], np.full(nsub - 1, theta0)])
    ph = np.concatenate([np.full(nsub + 1, phi0), p[1:], np.full(nsub, phi1), p[::-1][1:-1]])
    return th, ph


def cell_vectors(ntheta, nphi, nsub=4, idl=False):
    """Outline of every cell: simulation-frame unit vectors (cell centred), or with idl=True already
    ecliptic unit vectors of the placement of the original IDL plots. Returns [((i, j), vectors)]."""
    dth, dph = 180.0 / ntheta, 360.0 / nphi
    out = []
    for i in range(ntheta):
        for j in range(nphi):
            if idl:
                # original IDL plots: column x = j + nphi/2 - round(281/dphi) (mod nphi), drawn from lon -180 + x dphi;
                # row y = ntheta-1-i, drawn from lat -90 + y dtheta (pixels start at the grid point)
                x = (j + nphi // 2 - int(round(281.0 / dph))) % nphi
                lon0, lat0 = -180.0 + x * dph, -90.0 + (ntheta - 1 - i) * dth
                th, ph = cell_outline(90.0 - (lat0 + dth), 90.0 - lat0, lon0, lon0 + dph, nsub)
            else:
                th0, th1 = max(0.0, i * dth - dth / 2), min(180.0, i * dth + dth / 2)
                th, ph = cell_outline(th0, th1, j * dph - dph / 2, j * dph + dph / 2, nsub)
            out.append(((i, j), unit(ph, 90.0 - th)))
    return out


def display_polygons(x, y):
    """Closed outline in display lon/lat (deg) -> list of (x, y): unwrapped along the outline,
    closed through the pole if it turns around one, with copies shifted by +-360 deg where it
    crosses the map edge (removed by the clipping to the map outline)."""
    loop = np.degrees(np.unwrap(np.radians(np.append(x, x[0]))))
    x, total, y = loop[:-1], loop[-1] - loop[0], np.asarray(y, float)
    if abs(total) > 180.0:  # around a pole
        pole = 90.0 if np.mean(y) > 0 else -90.0
        x = np.concatenate([x, [x[0] + total, x[0] + total, x[0]]])
        y = np.concatenate([y, [y[0], pole, pole]])
        return [(x + s, y) for s in (-360.0, 0.0, 360.0)]
    polys = [(x, y)]
    if x.min() < -180.0:
        polys.append((x + 360.0, y))
    if x.max() > 180.0:
        polys.append((x - 360.0, y))
    return polys


# ---------------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------------
# IBEX energy bands (bands 0-7 IBEX-Lo, 8-13 IBEX-Hi): central energy (keV) and
# geometric factor (cm^2 sr keV/keV, triples), used to turn the legacy count-rate files back into flux
LEGACY_ECENTRAL = [0.014, 0.027, 0.052, 0.102, 0.197, 0.451, 0.908, 1.903, 0.45, 0.71, 1.11, 1.74, 2.73, 4.29]
LEGACY_G = [2.6e-5, 5.3e-5, 8.1e-5, 9.1e-5, 9.0e-5, 1.0e-4, 1.9e-4, 2.7e-4, 1.3e-4, 4.1e-4, 7.5e-4, 1.3e-3, 2.4e-3, 4.5e-3]


def read_legacy(path):
    """Legacy count-rate file (...<TAG>CountRateKappa<k>EbandNN.dat, one value per line, theta-major)
    -> flux map: count rate / E_central / G.
    The grid is 2:1 (nphi = 2 ntheta), e.g. 30 x 60 for 1800 values. The tag in the name (PUI1, PUI2,
    PUI3, SW, CombPUIs) is the component; there is no coverage information in these files."""
    v = np.loadtxt(path, dtype=float).ravel()
    nt = int(round(np.sqrt(v.size / 2)))
    if 2 * nt * nt != v.size:
        sys.exit(f"{path}: {v.size} values do not make a 2:1 theta-phi grid")
    m = re.search(r"Eband(\d+)", os.path.basename(path))
    if not m or int(m.group(1)) >= len(LEGACY_G):
        sys.exit(f"{path}: no IBEX band number (EbandNN) in the file name")
    band = int(m.group(1))
    # the tag between the model name and "CountRate" (the model name may itself contain the tag letters,
    # so take the known tags first, then whatever follows the last underscore)
    stem = os.path.basename(path).split("CountRate")[0]
    tag = next((t for t in ("CombPUIs", "PUI1", "PUI2", "PUI3", "SW") if stem.endswith(t)), stem.split("_")[-1])
    info = {"band": band, "ecentral_kev": LEGACY_ECENTRAL[band], "coverage": None, "legacy_component": tag}
    return v.reshape(nt, 2 * nt) / LEGACY_ECENTRAL[band] / LEGACY_G[band], info


def read_flux(path, component):
    if path.endswith(".dat"):  # legacy count-rate file: one component per file
        return read_legacy(path)
    with h5py.File(path, "r") as f:
        if component not in f:
            names = [k for k in f.keys() if f[k].ndim == 2 and not k.startswith("coverage")]
            sys.exit(f"{path}: no component '{component}' (available: {', '.join(names)})")
        m = np.array(f[component][()], dtype=float)
        info = {k: f[k][()] for k in ("band", "emin_kev", "emax_kev", "ecentral_kev", "window_start_yr",
                                       "window_end_yr", "include_incomplete_lines") if k in f}
        info["coverage"] = np.array(f["coverage_complete"][()]) if "coverage_complete" in f else None
    return m, info


def map_values(path, s):
    m, info = read_flux(path, s["component"])
    if s.get("ratio"):
        ref, _ = read_flux(s["ratio"], s["component"])
        with np.errstate(divide="ignore", invalid="ignore"):
            m = np.where(ref != 0, m / ref, np.nan)
    else:
        m = m * float(s["scale"])
    return m, info


def title_of(info, s):
    if s["title"]:
        return s["title"]
    t = f"Energy = {float(info['ecentral_kev']):.3g} keV" if "ecentral_kev" in info else os.path.basename(s["_path"])
    if not s.get("ratio") and float(s["scale"]) != 1.0:
        t += f" (Flux x {float(s['scale']):g})"
    comp = info.get("legacy_component", s["component"])  # legacy files: the tag in the file name
    if comp not in ("total", "CombPUIs"):
        t += f", {comp}"
    if s.get("ratio"):
        t += f", ratio to {os.path.basename(s['ratio'])}"
    if "window_start_yr" in info:
        t += f", {float(info['window_start_yr']):g}-{float(info['window_end_yr']):g}"
    return t


# ---------------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------------
def map_outline(proj):
    """Outline of the whole map (closed path) in projected coordinates."""
    lat = np.linspace(-90, 90, 361)
    x1, y1 = proj(np.full_like(lat, -180.0) * DEG, lat * DEG)
    x2, y2 = proj(np.full_like(lat, 180.0) * DEG, lat[::-1] * DEG)
    return np.column_stack([np.concatenate([x1, x2]), np.concatenate([y1, y2])])


def split_line(x, y):
    """Split a display polyline where it jumps across the map edge."""
    cut = np.where(np.abs(np.diff(x)) > 90.0)[0] + 1
    return [(a, b) for a, b in zip(np.split(x, cut), np.split(y, cut)) if len(a) > 1]


def draw_grid(ax, proj, V, s, clip, center_lon):
    segs = []
    lat = np.linspace(-90, 90, 721)
    lon_lines = list(np.arange(0.0, 360.0, s["grid_dlon"]))
    if s["edge_meridians"]:
        anti = (center_lon + 180.0) % 360.0
        lon_lines += [anti - s["edge_meridians"], anti + s["edge_meridians"]]
    for lon in lon_lines:
        x, y = to_display(unit(np.full_like(lat, lon), lat), V)
        for a, b in split_line(x, y):
            segs.append(np.column_stack(proj(a * DEG, b * DEG)))
    lons = np.linspace(0, 360, 1441)
    for la in np.arange(-90 + s["grid_dlat"], 90.0, s["grid_dlat"]):
        x, y = to_display(unit(lons, np.full_like(lons, la)), V)
        for a, b in split_line(x, y):
            segs.append(np.column_stack(proj(a * DEG, b * DEG)))
    lc = LineCollection(segs, colors="black", linewidths=0.8, zorder=3)
    lc.set_clip_path(clip)
    ax.add_collection(lc)
    fs = s["fontsize"] * 0.75
    box = dict(boxstyle="round,pad=0.08", fc="white", ec="none", alpha=0.55)
    for lon in np.arange(0.0, 360.0, s["grid_dlon"]):  # longitude labels along the ecliptic equator
        x, y = to_display(unit(lon, 0.0), V)
        if abs(x) > 165.0:
            continue
        px, py = proj(np.array(x) * DEG, np.array(y) * DEG)
        name = f"{((lon + 180) % 360) - 180:.0f}" if s["lon_labels"] == "signed" else f"{lon:.0f}"
        ax.text(px, py, name, ha="center", va="center", fontsize=fs, zorder=4, bbox=box)
    lab_lon = (center_lon + 60.0) % 360.0  # latitude labels along a meridian left of the centre
    for la in np.arange(-90 + s["grid_dlat"], 90.0, s["grid_dlat"]):
        if abs(la) < 1e-9:
            continue
        x, y = to_display(unit(lab_lon, la), V)
        px, py = proj(np.array(x) * DEG, np.array(y) * DEG)
        ax.text(px, py, f"{la:.0f}", ha="center", va="center", fontsize=fs, zorder=4, bbox=box)


def marker_list(s, R):
    flow = np.asarray(s["flow_sim"], float)
    up = lonlat(-(R @ (flow / np.linalg.norm(flow))))
    pos = {"UPWIND": up, "DOWNWIND": ((up[0] + 180.0) % 360.0, -up[1])}
    for k, v in s["marker_positions"].items():
        pos[k] = tuple(v)
    names = s["markers_" + s["view"]]
    out = []
    for n in names:
        if n not in pos:
            sys.exit(f"marker '{n}' has no position (set it under [markers] in the settings file)")
        out.append((n, float(pos[n][0]), float(pos[n][1])))
    return out, up


def colour_limits(vals, s, band, common):
    """Colour-bar range: vmin/vmax > band_limits > common range > the map's own min/max. Returns (lo, hi,
    lo_auto, hi_auto); an end is 'auto' when it comes from the data (map or common range), so it may be
    rounded outwards to the label precision (nice_limits)."""
    good = np.isfinite(vals)
    if s["log"]:
        good &= vals > 0
    lo = hi = None
    lo_auto = hi_auto = False
    if band is not None and band in s["band_limits"]:
        lo, hi = s["band_limits"][band]
    if common is not None and lo is None:
        lo, hi = common
        lo_auto = hi_auto = True
    if s["vmin"] is not None:
        lo = s["vmin"]
    if s["vmax"] is not None:
        hi = s["vmax"]
    if lo is None:
        lo = float(np.min(vals[good])) if good.any() else (1.0 if s["log"] else 0.0)
        lo_auto = True
    if hi is None:
        hi = float(np.max(vals[good])) if good.any() else lo * 10 if s["log"] else 1.0
        hi_auto = True
    if hi <= lo:
        hi = lo * 10 if s["log"] else lo + 1.0
    return float(lo), float(hi), lo_auto, hi_auto


def legend_style(s, lo, hi):
    """The colour-bar number style: the one asked for, or for "auto" by the size of the numbers."""
    if s["legend_style"] != "auto":
        return s["legend_style"]
    big = max(abs(lo), abs(hi))
    return "integer" if big >= 10 else "float" if big >= 0.1 else "exponent"


def label_step(style, lo, hi):
    """Smallest step the labels show: 1 (integer), 0.01 (float) or 0.01 of the leading power of ten (exponent)."""
    if style == "integer":
        return 1.0
    if style == "float":
        return 0.01
    big = max(abs(lo), abs(hi))
    return 0.01 * 10.0 ** np.floor(np.log10(big)) if big > 0 else 0.01


def nice_limits(lo, hi, lo_auto, hi_auto, style, ndiv=4):
    """Round the ends that came from the data outwards to the label step, and extend the top so the ndiv
    divisions fall on that step too: every colour-bar label is then exact (e.g. 62.8-186.3 -> 62-190 by 32,
    integer). Limits that were set by hand are kept as they are."""
    q = label_step(style, lo, hi)
    if lo_auto:
        lo = np.floor(lo / q + 1e-9) * q
    if hi_auto:
        hi = lo + ndiv * q * max(1.0, np.ceil((hi - lo) / (ndiv * q) - 1e-9))
    return float(lo), float(hi)


def format_label(v, style):
    """One colour-bar number: integers to the nearest whole number, floats and exponent mantissas to the
    nearest hundredth, halves rounded up (132.5 -> 133, 0.225 -> 0.23; the 1e-9 absorbs binary round-off)."""
    if style == "integer":
        return "%d" % np.floor(v + 0.5 + 1e-9)
    if style == "float":
        return "%.2f" % (np.floor(v * 100 + 0.5 + 1e-9) / 100)
    if v == 0:
        return "0.00e+00"
    e = int(np.floor(np.log10(abs(v))))
    m = np.sign(v) * np.floor(abs(v) / 10.0 ** e * 100 + 0.5 + 1e-9) / 100
    if abs(m) >= 10:  # 9.995 -> 10.00: carry into the exponent
        m, e = m / 10, e + 1
    return "%.2fe%+03d" % (m, e)


def legend_formatter(style):
    """Colour-bar label formatter for the style (see format_label)."""
    return FuncFormatter(lambda v, pos=None: format_label(v, style))


def plot_map(path, s, R, common=None):
    s = dict(s, _path=path)
    m, info = map_values(path, s)
    ntheta, nphi = m.shape
    proj = PROJECTIONS[s["projection"]]
    markers, up = marker_list(s, R)
    if s["center"] is not None:
        c_lon, c_lat = float(s["center"][0]), float(s["center"][1])
    elif s["view"] == "nose":
        c_lon, c_lat = float(up[0]), float(up[1])
    elif s["view"] == "ecliptic":
        c_lon, c_lat = 0.0, 0.0
    else:
        c_lon, c_lat = (float(up[0]) + 180.0) % 360.0, -float(up[1])
    V = view_matrix(c_lon, c_lat)

    cov = info["coverage"] if s["boxes"] else None
    polys, vals, flags = [], [], []
    for (i, j), pts in cell_vectors(ntheta, nphi, nsub=int(s["nsub"]), idl=s["idl"]):
        x, y = to_display(pts if s["idl"] else pts @ R.T, V)
        for px, py in display_polygons(x, y):
            polys.append(np.column_stack(proj(px * DEG, py * DEG)))
            vals.append(m[i, j])
            flags.append(cov is not None and cov[i, j] == 0)
    vals, flags = np.array(vals), np.array(flags, bool)

    fig, ax = plt.subplots(figsize=tuple(s["figsize"]))
    ax.set_aspect("equal")
    ax.axis("off")
    outline = map_outline(proj)
    clip_path = Path(outline)
    region = None
    if s["limits"]:  # part of the sky only (relative to the centre; longitude positive to the left)
        lo1, lo2, la1, la2 = map(float, s["limits"])
        ll = np.linspace(lo1, lo2, 181)
        bb = np.linspace(la1, la2, 91)
        rx = np.concatenate([ll, np.full_like(bb, lo2), ll[::-1], np.full_like(bb, lo1)])
        ry = np.concatenate([np.full_like(ll, la1), bb, np.full_like(ll, la2), bb[::-1]])
        region = np.column_stack(proj(rx * DEG, ry * DEG))
        clip_path = Path(region)
    clip = PathPatch(clip_path, transform=ax.transData, facecolor="none", edgecolor="none")
    ax.add_patch(clip)
    ext = region if region is not None else outline
    pad = 0.03 * (ext[:, 0].max() - ext[:, 0].min())
    ax.set_xlim(ext[:, 0].min() - pad, ext[:, 0].max() + pad)
    ax.set_ylim(ext[:, 1].min() - pad, ext[:, 1].max() + pad)

    lo, hi, lo_auto, hi_auto = colour_limits(vals[~flags], s, int(info["band"]) if "band" in info else None, common)
    style = legend_style(s, lo, hi)
    if not s["log"]:
        lo, hi = nice_limits(lo, hi, lo_auto, hi_auto, style)
    norm = LogNorm(lo, hi) if s["log"] else Normalize(lo, hi)
    cmap = plt.get_cmap(s["cmap"])
    colors = np.tile(np.array([0.75, 0.75, 0.75, 1.0]), (len(vals), 1))  # grey: missing flux or no value
    ok = np.isfinite(vals) & ~flags
    v = np.clip(vals[ok], lo, hi)  # outside the limits: end colours
    colors[ok] = cmap(norm(v))
    pc = PolyCollection(polys, facecolors=colors, edgecolors=colors, linewidths=0.3, zorder=1)
    pc.set_clip_path(clip)
    ax.add_collection(pc)
    if flags.any():
        boxes = PolyCollection([p for p, f in zip(polys, flags) if f], facecolors="none", edgecolors="black",
                               linewidths=0.9, zorder=2)
        boxes.set_clip_path(clip)
        ax.add_collection(boxes)
    if s["grid"]:
        draw_grid(ax, proj, V, s, clip, c_lon)
    ax.plot(*(region if region is not None else outline).T, color="black", lw=1.2, zorder=5)
    if region is None:
        ax.plot(*outline[[0, -1]].T, color="black", lw=1.2, zorder=5)

    for name, lo_, la_ in markers:
        x, y = to_display(unit(lo_, la_), V)
        if region is None and abs(x) > 178:
            continue
        px, py = proj(np.array(x) * DEG, np.array(y) * DEG)
        ax.plot(px, py, "o", color="black", ms=7, zorder=6)
        below = name == "DOWNWIND"
        ax.text(px + 0.05, py - 0.13 if below else py + 0.06, name, fontsize=s["fontsize"], fontweight="bold",
                ha="left", va="top" if below else "bottom", zorder=6)

    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    cb = fig.colorbar(sm, ax=ax, orientation="horizontal", fraction=0.05, pad=0.03, aspect=40)
    if not s["log"]:
        cb.set_ticks(np.linspace(lo, hi, 5))  # 4 divisions
    if not s["log"] or s["legend_style"] != "auto":  # log scale: decade labels unless a style is chosen
        cb.ax.xaxis.set_major_formatter(legend_formatter(style))
    if not s["log"]:
        labels = [legend_formatter(style)(t) for t in np.linspace(lo, hi, 5)]
        if len(set(labels)) < len(labels):  # e.g. integer labels on a 0-0.9 range: 0 0 0 1 1
            print(f"WARNING: {os.path.basename(path)}: legend_style '{style}' is too coarse for the range "
                  f"{lo:g}-{hi:g} (labels {' '.join(labels)}); use float or exponent")
    cb.ax.tick_params(labelsize=s["fontsize"] * 0.9)
    cb.set_label("flux ratio" if s.get("ratio") else s["units_label"], fontsize=s["fontsize"] * 0.85)
    ax.set_title(title_of(info, s), loc="left", fontsize=s["fontsize"] * 1.25)
    if flags.any():
        kept = int(info.get("include_incomplete_lines", 0)) == 1
        fig.text(0.5, 0.0, f"boxes: regions with missing flux ({int((cov == 0).sum())} cells"
                 f"{'' if kept else ', flux set to 0'})", ha="center", fontsize=s["fontsize"] * 0.7)
    note = "placement of the original IDL plots" if s["idl"] else "ecliptic J2000"
    fig.text(0.99, 0.0, f"{note}; {s['projection']}; centre ({c_lon:.1f}, {c_lat:+.1f})", ha="right",
             fontsize=s["fontsize"] * 0.55, color="0.4")

    stem = os.path.splitext(os.path.basename(path))[0] + ("" if path.endswith(".dat") else f"_{s['component']}") + ("_ratio" if s.get("ratio") else "") + \
        ("_idl" if s["idl"] else "") + f"_{s['view'] if s['center'] is None else 'custom'}"
    out = s["out"] or os.path.dirname(path) or "."
    os.makedirs(out, exist_ok=True)
    name = os.path.join(out, f"{stem}.{s['format']}")
    fig.savefig(name, dpi=int(s["dpi"]), bbox_inches="tight")
    plt.close(fig)
    return name


def write_ecliptic(path, R, res):
    """All components on a regular ecliptic grid of `res` deg (bin centres): each bin takes the
    simulation cell that contains its centre; NaN where there is no cell."""
    with h5py.File(path, "r") as f:
        data = {k: np.array(f[k][()]) for k in f.keys() if f[k].ndim == 2}
    ntheta, nphi = next(iter(data.values())).shape
    dth, dph = 180.0 / ntheta, 360.0 / nphi
    lon = np.arange(res / 2, 360.0, res)
    lat = np.arange(-90.0 + res / 2, 90.0, res)
    LON, LAT = np.meshgrid(lon, lat, indexing="xy")
    v_sim = unit(LON, LAT) @ R
    ph, la = lonlat(v_sim)
    i = np.rint((90.0 - la) / dth).astype(int)
    j = np.rint(ph / dph).astype(int) % nphi
    inside = i < ntheta
    i = np.clip(i, 0, ntheta - 1)
    out = os.path.splitext(path)[0] + "_ecliptic.h5"
    with h5py.File(out, "w") as g:
        g["lon_deg"], g["lat_deg"] = lon, lat
        for k, a in data.items():
            g[k] = np.where(inside, a[i, j], np.nan)
        g.attrs["note"] = (f"regular ecliptic J2000 grid, {res} deg bins (centres in lon_deg, lat_deg); arrays [lat][lon]; "
                           "each bin holds the simulation cell containing its centre; NaN: no cell")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="+", help="flux_bandNN*.h5 files or legacy ...CountRate...EbandNN.dat files (wildcards allowed)")
    ap.add_argument("--config", help="plot settings file (configs/plot_settings.toml lists every setting)")
    ap.add_argument("--component")
    ap.add_argument("--view", choices=["nose", "tail", "ecliptic"],
                    help="map centre: nose (upwind), tail (downwind) or ecliptic (lon 0, lat 0)")
    ap.add_argument("--center", type=float, nargs=2, metavar=("LON", "LAT"), help="map centre, ecliptic deg")
    ap.add_argument("--projection", choices=sorted(PROJECTIONS))
    ap.add_argument("--idl", action="store_const", const=True, help="placement of the original IDL plots (no rotation)")
    ap.add_argument("--scale", type=float, help="multiply the flux by this")
    ap.add_argument("--log", action="store_const", const=True, help="log colour scale")
    ap.add_argument("--linear", dest="log", action="store_const", const=False, help="linear colour scale (default)")
    ap.add_argument("--vmin", type=float)
    ap.add_argument("--vmax", type=float)
    ap.add_argument("--common-range", action="store_const", const=True, help="one range per band over all files")
    ap.add_argument("--legend-style", choices=["auto", "integer", "float", "exponent"])
    ap.add_argument("--cmap")
    ap.add_argument("--limits", type=float, nargs=4, metavar=("LONMIN", "LONMAX", "LATMIN", "LATMAX"))
    ap.add_argument("--ratio", help="divide by the same component of this file")
    ap.add_argument("--no-boxes", dest="boxes", action="store_const", const=False)
    ap.add_argument("--no-grid", dest="grid", action="store_const", const=False)
    ap.add_argument("--no-markers", action="store_true")
    ap.add_argument("--title")
    ap.add_argument("--format")
    ap.add_argument("--dpi", type=int)
    ap.add_argument("--out", help="output directory (default: next to each input file)")
    ap.add_argument("--flow-sim", type=float, nargs=3, metavar=("VX", "VY", "VZ"))
    ap.add_argument("--flow-ecliptic", type=float, nargs=2, metavar=("LON", "LAT"))
    ap.add_argument("--solar-axis", type=float, nargs=2, metavar=("INCLINATION", "NODE"))
    ap.add_argument("--write-ecliptic", action="store_const", const=True)
    ap.add_argument("--ecliptic-res", type=float)
    a = ap.parse_args()

    s = {k: (dict(v) if isinstance(v, dict) else v) for k, v in DEFAULTS.items()}
    if a.config:
        for k, v in parse_settings(a.config).items():
            if k not in s:
                sys.exit(f"{a.config}: unknown setting '{k}'")
            if isinstance(s[k], dict):
                s[k].update(v)
            else:
                s[k] = v
    for k, v in vars(a).items():
        if k in ("files", "config", "no_markers", "ratio") or v is None:
            continue
        s[k] = v
    s["ratio"] = a.ratio
    if a.no_markers:
        s["markers_nose"], s["markers_tail"], s["markers_ecliptic"] = [], [], []

    a_sim, a_ecl = frame_check(s["flow_sim"], s["flow_ecliptic"], s["solar_axis"])
    if abs(a_sim - a_ecl) > 1.0:
        print(f"WARNING: the flow makes {a_sim:.2f} deg with the solar axis in the simulation frame but {a_ecl:.2f} deg "
              f"in ecliptic J2000: check flow_sim, flow_ecliptic and solar_axis")
    R = sim_to_ecliptic(s["flow_sim"], s["flow_ecliptic"], s["solar_axis"])
    files = [g for f in a.files for g in (sorted(glob.glob(f)) or [f])]

    common = {}
    if s["common_range"]:  # one colour range per band over all files
        for f in files:
            m, info = map_values(f, s)
            band = int(info["band"]) if "band" in info else None
            cov = info["coverage"]
            vals = m[(cov != 0)] if cov is not None else m.ravel()
            vals = vals[np.isfinite(vals) & ((vals > 0) if s["log"] else True)]
            if vals.size:
                lo, hi = common.get(band, (np.inf, -np.inf))
                common[band] = (min(lo, vals.min()), max(hi, vals.max()))
    for f in files:
        band = None
        if s["common_range"]:
            with h5py.File(f, "r") as h:
                band = int(h["band"][()]) if "band" in h else None
        print(plot_map(f, s, R, common.get(band) if s["common_range"] else None))
        if s["write_ecliptic"]:
            print(write_ecliptic(f, R, float(s["ecliptic_res"])))


if __name__ == "__main__":
    main()
