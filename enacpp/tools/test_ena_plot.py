#!/usr/bin/env python3
"""Checks of the coordinate handling of ena_plot.py. Run: python3 enacpp/tools/test_ena_plot.py"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ena_plot as ep  # noqa: E402

fails = 0


def check(name, ok, detail=""):
    global fails
    print(f"{'ok  ' if ok else 'FAIL'} {name}{(': ' + detail) if detail else ''}")
    fails += 0 if ok else 1


R = ep.sim_to_ecliptic()
FLOW_SIM, FLOW_ECLIPTIC, SOLAR_AXIS = ep.DEFAULTS["flow_sim"], ep.DEFAULTS["flow_ecliptic"], ep.DEFAULTS["solar_axis"]
check("rotation is orthonormal", np.allclose(R @ R.T, np.eye(3), atol=1e-12) and abs(np.linalg.det(R) - 1) < 1e-12)

# the solar rotation axis (simulation z) lands exactly on the pole of the solar equator
pole = R @ np.array([0, 0, 1.0])
check("simulation z = solar rotation axis", np.allclose(pole, ep.solar_pole_ecliptic(*SOLAR_AXIS), atol=1e-12),
      "ecliptic lon %.2f, lat %.2f" % ep.lonlat(pole))

# the flow: simulation (26.3, 0.3, -2.3) -> close to ecliptic (75.4, -5.2); exact up to the 0.1-deg mismatch
f = R @ (np.array(FLOW_SIM) / np.linalg.norm(FLOW_SIM))
sep = np.degrees(np.arccos(np.clip(f @ ep.unit(*FLOW_ECLIPTIC), -1, 1)))
a_sim, a_ecl = ep.frame_check()
check("flow direction within the frame mismatch", sep <= abs(a_sim - a_ecl) + 1e-9,
      f"separation {sep:.3f} deg, mismatch {a_sim - a_ecl:+.3f} deg")
nose = ep.lonlat(-f)
check("nose (upstream) at about (255.4, +5.2)", abs(nose[0] - 255.4) < 0.15 and abs(nose[1] - 5.2) < 0.15,
      "lon %.2f, lat %.2f" % nose)

# simulation equator: tilted by the 7.25-deg inclination; flanks reach +-7.2 deg
lon_f1, lat_f1 = ep.lonlat(R @ ep.unit(90.0, 0.0))
lon_f2, lat_f2 = ep.lonlat(R @ ep.unit(270.0, 0.0))
check("flanks at +-7.2 deg latitude", abs(lat_f1 - 7.2) < 0.1 and abs(lat_f2 + 7.2) < 0.1,
      f"phi 90 -> ({lon_f1:.1f}, {lat_f1:.1f}), phi 270 -> ({lon_f2:.1f}, {lat_f2:.1f})")

# projections: Mollweide equal area, edges, pole
x, y = ep.mollweide(np.array([np.pi, -np.pi, 0.0, 0.0]), np.array([0.0, 0.0, np.pi / 2, -np.pi / 2]))
check("Mollweide edges and poles (longitude positive to the left)", np.allclose(x, [-2 * np.sqrt(2), 2 * np.sqrt(2), 0, 0]) and
      np.allclose(y, [0, 0, np.sqrt(2), -np.sqrt(2)]))
# equal area: band between latitudes b1, b2 has area proportional to sin b2 - sin b1
lat = np.linspace(-np.pi / 2, np.pi / 2, 20001)
_, yy = ep.mollweide(np.zeros_like(lat), lat)
xx, _ = ep.mollweide(np.full_like(lat, np.pi), lat)
area = np.trapezoid(2 * np.abs(xx), yy) if hasattr(np, "trapezoid") else np.trapz(2 * np.abs(xx), yy)
check("Mollweide total area = 4 pi (unit sphere)", abs(area - 4 * np.pi) < 1e-3, f"{area:.5f}")

# cells cover the sphere once (cell-centred grid: only the 3-deg cap around the simulation south pole is empty)
ntheta, nphi = 30, 60
solid = 0.0
for (i, j), pts in ep.cell_vectors(ntheta, nphi, nsub=1):
    th0, th1 = max(0, i * 6 - 3), min(180, i * 6 + 3)
    solid += (np.cos(np.radians(th0)) - np.cos(np.radians(th1))) * np.radians(6)
cap = 2 * np.pi * (1 - np.cos(np.radians(3)))
check("cells cover the sphere except the south cap", abs(solid + cap - 4 * np.pi) < 1e-9,
      f"cells {solid:.6f} sr + cap {cap:.6f} sr")

# placement of the original IDL plots (--idl): the nose cell (theta 84, phi 180) drawn at lon 258-264, lat 0-6 (pixel starting at its grid point)
cells = dict(ep.cell_vectors(ntheta, nphi, nsub=1, idl=True))
lon, lat = ep.lonlat(cells[(14, 30)])
lon = np.where(lon < 200, lon + 360, lon)
check("--idl placement of the nose cell", np.isclose(lon.min() % 360, 258.0) and np.isclose(lon.max() % 360, 264.0)
      and np.isclose(lat.min(), 0.0) and np.isclose(lat.max(), 6.0),
      f"lon {lon.min() % 360:.0f}-{lon.max() % 360:.0f}, lat {lat.min():.0f}-{lat.max():.0f}")

# the two cells that contain the ecliptic poles: their outlines turn by 360 deg in longitude and are
# closed through the pole; every other cell is an ordinary polygon
enclosing = []
for (i, j), pts in ep.cell_vectors(ntheta, nphi, nsub=4):
    x, y = ep.to_display(pts @ R.T, ep.view_matrix(255.0, 0.0))
    polys = ep.display_polygons(x, y)
    if len(polys) == 3 and abs(polys[1][1][-1]) == 90.0:
        enclosing.append(((i, j), polys[1][1][-1]))
check("one cell around each ecliptic pole", sorted(p for _, p in enclosing) == [-90.0, 90.0],
      ", ".join(f"cell {c} -> pole {p:+.0f}" for c, p in enclosing))

# view rotation: the chosen centre goes to display (0, 0)
up = ep.lonlat(-f)
for (clon, clat) in [(up[0], up[1]), (0.0, 0.0), ((up[0] + 180) % 360, -up[1])]:
    x, y = ep.to_display(ep.unit(clon, clat), ep.view_matrix(clon, clat))
    check(f"centre ({clon:.1f}, {clat:+.1f}) at display (0, 0)", abs(x) < 1e-9 and abs(y) < 1e-9)
# longitude increases to the left: a point 10 deg east of the centre has display x > 0 and projected x < 0
x, y = ep.to_display(ep.unit(10.0, 0.0), ep.view_matrix(0.0, 0.0))
px, _ = ep.mollweide(np.array(x) * ep.DEG, np.array(y) * ep.DEG)
check("east is to the left", np.isclose(x, 10.0) and px < 0)

# legacy count-rate files: theta-major, 2:1 grid, flux = count rate / E_central / G
import tempfile  # noqa: E402
with tempfile.TemporaryDirectory() as d:
    cr = np.arange(1800, dtype=float).reshape(30, 60)
    path = os.path.join(d, "Opher_multiion_withneutralsPUI2CountRateKappa1.00Eband13.dat")
    np.savetxt(path, cr.ravel(), fmt="%16.8f")
    m, info = ep.read_flux(path, "total")
    check("legacy count-rate file -> flux", m.shape == (30, 60) and np.allclose(m, cr / 4.29 / 4.5e-3)
          and info["legacy_component"] == "PUI2" and info["band"] == 13)

# colour-bar labels: automatic ends rounded outwards so the 4 divisions fall on the label step
for style, lo, hi, want in [("integer", 62.8, 186.3, (62.0, 190.0)), ("float", 0.0988, 0.8306, (0.09, 0.85)),
                            ("exponent", 0.00123, 0.00869, (0.00123, 0.00871))]:
    got = ep.nice_limits(lo, hi, True, True, style)
    ticks = np.linspace(*got, 5) / ep.label_step(style, *got)
    check(f"{style} labels exact ({lo}-{hi} -> {got[0]:g}-{got[1]:g})", np.allclose(got, want)
          and np.allclose(ticks, np.round(ticks)))
check("limits set by hand are kept", ep.nice_limits(0.0, 0.9, False, False, "float") == (0.0, 0.9))
fm = [ep.format_label(v, st) for st, v in (("integer", 66.25), ("integer", 132.5), ("float", 0.225), ("float", 0.675),
                                          ("exponent", 0.00225), ("exponent", 0.009995), ("exponent", 0.0))]
check("label formats (halves round up)", fm == ["66", "133", "0.23", "0.68", "2.25e-03", "1.00e-02", "0.00e+00"],
      ", ".join(fm))
s_ = dict(ep.DEFAULTS)
check("default band limits (bands 9-17)", s_["band_limits"][13] == [0.0, 12.0]
      and s_["band_limits"][9] == [0.0, 265.0] and s_["band_limits"][17] == [0.0, 0.009])

print(f"{fails} failures")
sys.exit(1 if fails else 0)
