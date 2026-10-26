#!/usr/bin/env python3
"""Steady-state limit test of the time-dependent method (after td_repeat_snapshot.py):
compares the time-dependent maps made from one snapshot repeated N times with the
steady-state maps of that snapshot.

With identical snapshots the only difference is the flux made before the first
snapshot, which is missing from the early windows. Three checks, per band and window:

 1. Missing share (all windows, all lines of sight whose flux is written): 1 -
    total_td/total_steady must equal coverage_missing_flux_fraction, the code's own
    estimate of the missing share. Expect about 1e-15. Lines of sight flagged as
    regions with missing flux are 0 by default (include_incomplete_lines = false) and
    are then skipped; run with --set time_dependent.include_incomplete_lines=true to
    check them too.
 2. Fully covered lines of sight (missing share exactly 0): every component must
    equal the steady-state map. Expect about 1e-15 for every value that is not
    negligible (at least 1e-12 of its line's total flux). Negligible values, e.g. a SW
    component of 1e-276, can differ by up to about 1e-9: in the steady-state code the
    survival factor enters the cell flux before the other factors, and the intermediate
    product of such a value drops below the smallest normal double (about 1e-308), where
    it has only a few significant digits. The time-dependent code applies the survival
    factor last and keeps the precision.
 3. Lines of sight marked complete but not fully covered (missing share > 0, below
    the coverage thresholds): not an error. Listed with the missing share of the
    total and of the component that misses the most (a single component, e.g. the
    heliosheath-created PUIs, can miss a larger share than the total).

usage: python3 td_vs_steady.py STEADY_DIR TD_DIR"""
import sys, glob, os, re
import h5py, numpy as np

steady, td = sys.argv[1], sys.argv[2]
worst_share = worst_full = 0.0
full_cases = []  # (relative difference, band, window, component, theta index, phi index, td, steady, line total)
n_files = n_full_lines = 0
partial = []
for f in sorted(glob.glob(f"{td}/flux_band*_*.h5")):
    m = re.match(r"flux_band(\d+)_(.+)\.h5", os.path.basename(f))
    ref = f"{steady}/flux_band{m.group(1)}.h5"
    if not os.path.exists(ref):
        continue
    n_files += 1
    with h5py.File(f) as a, h5py.File(ref) as b:
        comps = [c for c in b.keys() if b[c].ndim == 2 and c in a]
        share = a["coverage_missing_flux_fraction"][()]
        complete = a["coverage_complete"][()] == 1
        written = complete | (a["include_incomplete_lines"][()] == 1 if "include_incomplete_lines" in a else True)
        # 1. missing share of the total
        x, y = a["total"][()], b["total"][()]
        nz = (y != 0) & written
        worst_share = max(worst_share, float(np.max(np.abs((1 - x[nz] / y[nz]) - share[nz]))) if nz.any() else 0.0)
        # 2. fully covered lines: every component equal
        full = share == 0
        n_full_lines += int(full.sum())
        for c in comps:
            x, y = a[c][()], b[c][()]
            sel = full & (y != 0)
            if sel.any():
                r = np.where(sel, np.abs(x / np.where(y == 0, 1, y) - 1), 0.0)
                worst_full = max(worst_full, float(r.max()))
                for i, j in zip(*np.nonzero(r > 1e-13)):
                    full_cases.append((float(r[i, j]), m.group(1), m.group(2), c, int(i), int(j), float(x[i, j]),
                                       float(y[i, j]), float(b["total"][i, j])))
            if np.any(x[full & (y == 0)] != 0):
                worst_full = np.inf
        # 3. complete but partly covered
        part = complete & ~full
        if part.any():
            comp_worst, comp_name = 0.0, ""
            for c in comps:
                x, y = a[c][()], b[c][()]
                sel = part & (y != 0)
                if sel.any():
                    r = float(np.max(1 - x[sel] / y[sel]))
                    if r > comp_worst:
                        comp_worst, comp_name = r, c
            partial.append((m.group(1), m.group(2), int(part.sum()), float(share[part].max()), comp_name, comp_worst))

print(f"{n_files} band/window files compared")
print(f"1. missing share, 1 - total_td/total_steady vs coverage_missing_flux_fraction (all lines): max difference {worst_share:.3g}")
sig = [c for c in full_cases if c[8] and c[7] / c[8] >= 1e-12]
worst_sig = max((c[0] for c in sig), default=0.0)
print(f"2. fully covered lines of sight ({n_full_lines} line-windows): max relative difference, all components {worst_full:.3g}")
print(f"   values of at least 1e-12 of their line's total flux: max relative difference "
      f"{worst_sig if worst_sig > 1e-13 else 0.0:.3g} (differences up to 1e-13 counted as 0); "
      f"negligible values differing by more than 1e-13: {len(full_cases) - len(sig)}")
if full_cases:
    full_cases.sort(reverse=True)
    print(f"   {len(full_cases)} values differ by more than 1e-13; the largest (theta/phi as grid indices; "
          f"'share' = the component's share of that line's total flux):")
    for r, band, win, c, i, j, x, y, tot in full_cases[:12]:
        print(f"     {r:.3g}  band {band}, window {win}, {c}, line ({i},{j}) = theta {i * 6} phi {j * 6} deg: "
              f"td {x:.10g}, steady {y:.10g}, share {y / tot if tot else 0:.3g}")
print("3. lines marked complete but not fully covered (expected; missing share below the thresholds):")
if not partial:
    print("   none")
for band, win, n, s, c, r in partial:
    print(f"   band {band}, window {win}: {n} lines; missing share of the total up to {s:.3f}; "
          f"most missing component: {c}, up to {r:.3f}")
