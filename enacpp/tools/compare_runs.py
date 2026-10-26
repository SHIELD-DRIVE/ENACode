#!/usr/bin/env python3
"""Compare the maps of two runs made on different machines (e.g. your workstation and Athena).

    python3 enacpp/tools/compare_runs.py REF_DIR OTHER_DIR [--rtol 1e-9] [--floor 1e-12]

same_maps.py asks for bit-identical maps, which holds for any thread count or MPI layout on
one machine. Between machines it need not: another compiler, another maths library (exp, log,
pow can differ in the last bit between glibc versions) or another CPU can change the last
digits. This tool reports how far apart the runs are:

- every flux_bandNN*.h5 of REF_DIR is matched with the file of the same name in OTHER_DIR;
- the run description must agree: band, energies, model, version, settings, number of
  streamlines and of streamlines with a TS, r_TS along V2;
- every map (2D dataset) is compared cell by cell. The relative difference of a cell is
  |a - b| / max(|a|, floor * F), F the largest total flux of the band: values that are
  negligible in the band's total (e.g. a solar-wind component of 1e-276) are measured against
  floor * F, so round-off in them does not count; any change that matters for the total does;
- the legacy count-rate files (*.dat, written to 17 digits) are compared the same way, with F
  the largest count rate over all files of the band.

Verdict per map: IDENTICAL (bit for bit), AGREE (largest relative difference <= rtol) or
DIFFER (the worst cells are listed with theta, phi and both values). Exit status 0 when nothing
DIFFERs and no file is missing, 1 otherwise.
"""
import argparse
import glob
import os
import sys

import numpy as np

try:
    import h5py
except ImportError:
    sys.exit("compare_runs.py needs numpy and h5py (pip install numpy h5py)")

# scalar entries that describe the run and must agree (when present in both files)
DESCRIPTION = ("band", "emin_kev", "emax_kev", "ecentral_kev", "geometric_factor", "model", "version", "settings",
               "streamlines", "streamlines_with_ts", "r_ts_v2_au", "window_start_yr", "window_end_yr")


def rel_diff(a, b, floor_abs):
    """Cell-wise relative difference, measured against floor_abs for negligible cells."""
    scale = np.maximum(np.abs(a), floor_abs)
    with np.errstate(divide="ignore", invalid="ignore"):
        d = np.where(scale > 0, np.abs(a - b) / scale, np.where(a == b, 0.0, np.inf))
    both_nan = np.isnan(a) & np.isnan(b)
    d[both_nan] = 0.0
    d[np.isnan(a) ^ np.isnan(b)] = np.inf
    return d


def compare_map(name, a, b, args, band_max, dtheta=None, dphi=None):
    """Print one line for the map (and its worst cells if it differs); return 'identical', 'agree' or 'differ'."""
    if a.shape != b.shape:
        print(f"  DIFFER    {name}: shape {a.shape} vs {b.shape}")
        return "differ"
    if np.array_equal(a, b, equal_nan=True):
        return "identical"
    d = rel_diff(a.astype(float), b.astype(float), args.floor * band_max)
    worst = float(np.max(d))
    if worst <= args.rtol:
        print(f"  AGREE     {name}: max relative difference {worst:.1e}")
        return "agree"
    over = int(np.sum(d > args.rtol))
    print(f"  DIFFER    {name}: max relative difference {worst:.1e}, {over} of {d.size} cells over {args.rtol:g}")
    order = [k for k in np.argsort(d, axis=None)[::-1][:args.show] if d.flat[k] > args.rtol]
    for k in order:
        idx = tuple(int(i) for i in np.unravel_index(k, d.shape))
        where = f"cell {idx if len(idx) > 1 else idx[0]}"
        if d.ndim == 2 and dtheta:
            where += f" (theta {idx[0] * dtheta:g}, phi {idx[1] * dphi:g})"
        print(f"              {where}: {a[idx]:.10g} vs {b[idx]:.10g} (relative {d[idx]:.1e})")
    return "differ"


def as_text(v):
    v = v.decode() if isinstance(v, bytes) else v
    return v.item() if isinstance(v, np.generic) else v


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("ref")
    ap.add_argument("other")
    ap.add_argument("--rtol", type=float, default=1e-9, help="largest relative difference accepted (default 1e-9)")
    ap.add_argument("--floor", type=float, default=1e-12,
                    help="cells below floor x the map's maximum are compared against that level (default 1e-12)")
    ap.add_argument("--show", type=int, default=5, help="worst cells listed per differing map (default 5)")
    args = ap.parse_args()

    counts = {"identical": 0, "agree": 0, "differ": 0}
    missing = []
    files = sorted(glob.glob(os.path.join(args.ref, "flux_band*.h5")))
    if not files:
        sys.exit(f"no flux_band*.h5 files in {args.ref}")
    for f in files:
        g = os.path.join(args.other, os.path.basename(f))
        print(os.path.basename(f))
        if not os.path.exists(g):
            print("  MISSING in", args.other)
            missing.append(g)
            continue
        with h5py.File(f, "r") as a, h5py.File(g, "r") as b:
            for k in DESCRIPTION:  # the two runs must be the same run
                if k in a and k in b and as_text(a[k][()]) != as_text(b[k][()]):
                    print(f"  DIFFER    description '{k}': {as_text(a[k][()])!r} vs {as_text(b[k][()])!r}")
                    counts["differ"] += 1
            maps = [k for k in a if isinstance(a[k], h5py.Dataset) and a[k].ndim == 2]
            # the band's scale: its total flux, else the largest of its maps
            band_max = float(np.nanmax(np.abs(a["total"][()]))) if "total" in maps else \
                max(float(np.nanmax(np.abs(a[k][()]))) for k in maps)
            for k in maps:
                if k not in b:
                    print(f"  MISSING   map {k} in {g}")
                    missing.append(f"{g}:{k}")
                    continue
                counts[compare_map(k, a[k][()], b[k][()], args, band_max, 6.0, 6.0)] += 1
    dats = sorted(glob.glob(os.path.join(args.ref, "*.dat")))
    if dats:
        print("legacy count-rate files")
    ref_vals = {f: np.loadtxt(f, ndmin=1) for f in dats}
    band_of = {f: (os.path.basename(f).split("Eband")[1] if "Eband" in os.path.basename(f) else "") for f in dats}
    band_max = {}
    for f, v in ref_vals.items():  # the band's scale: its largest count rate over all its files
        band_max[band_of[f]] = max(band_max.get(band_of[f], 0.0), float(np.nanmax(np.abs(v))) if v.size else 0.0)
    for f in dats:
        g = os.path.join(args.other, os.path.basename(f))
        if not os.path.exists(g):
            print("  MISSING", os.path.basename(f))
            missing.append(g)
            continue
        counts[compare_map(os.path.basename(f), ref_vals[f], np.loadtxt(g, ndmin=1), args, band_max[band_of[f]])] += 1

    print(f"\n{counts['identical']} identical, {counts['agree']} agree within {args.rtol:g}, "
          f"{counts['differ']} differ, {len(missing)} missing")
    if counts["differ"] or missing:
        print("RESULT: DIFFERENT - see the maps marked DIFFER or MISSING above")
        return 1
    print("RESULT: SAME" + (" (bit-identical)" if counts["agree"] == 0 else f" (within relative {args.rtol:g})"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
