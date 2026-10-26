#!/usr/bin/env python3
"""Compares the per-cell fields that ena_run builds from the streamlines (checkpoint C3, written
with --dump) between two runs, to see which quantity a change of the streamline file changes.

    python3 enacpp/tools/compare_dumps.py DUMP_A DUMP_B [--file C3_stream_final.h5] [--rtol 0.01]

DUMP_A and DUMP_B are the --dump folders of two ena_run runs. C3_stream_final.h5 holds the fields
after the gap filling; C3_stream_after_pass00.h5 (if written) holds them as binned from the
streamline points, before any filling. For each field ([theta][phi][r] cells):
  - cells with a value in one run only (NaN or empty in the other);
  - cells with a value in both that differ by more than rtol, by theta row and by r range;
  - the cells that differ the most.
The grid spacing for the printed angles and radii is read from C0_config.h5 if it holds them,
otherwise 6 deg, 6 deg, 2 AU from r = 2 AU (change with --dtheta, --dphi, --dr, --r-inner).
  --theta-max T     count only cells with theta <= T deg (e.g. the polar cap)
  --los THETA PHI   also print, along that line of sight, the main fields of both runs at every
                    radius where either run has a value (A / B)
Needs numpy and h5py (e.g. ~/pytecplot-env/bin/python3).
"""
import argparse
import os
import sys

import h5py
import numpy as np

R_RANGES = ((0, 100), (100, 200), (200, 400), (400, 800), (800, 1e9))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("a")
    ap.add_argument("b")
    ap.add_argument("--file", default="C3_stream_final.h5")
    ap.add_argument("--rtol", type=float, default=0.01)
    ap.add_argument("--show", type=int, default=5)
    ap.add_argument("--theta-max", type=float, default=None)
    ap.add_argument("--los", type=float, nargs=2, metavar=("THETA", "PHI"))
    ap.add_argument("--max-rows", type=int, default=120)
    ap.add_argument("--dtheta", type=float, default=6.0)
    ap.add_argument("--dphi", type=float, default=6.0)
    ap.add_argument("--dr", type=float, default=2.0)
    ap.add_argument("--r-inner", type=float, default=2.0)
    a = ap.parse_args()
    fa, fb = os.path.join(a.a, a.file), os.path.join(a.b, a.file)
    for f in (fa, fb):
        if not os.path.exists(f):
            print(f"missing: {f} (run ena_run with --dump on a steady-state settings file)")
            return 1
    with h5py.File(fa, "r") as ha, h5py.File(fb, "r") as hb:
        names = [n for n in ha.keys() if n in hb]
        print(f"{a.file}: {len(names)} fields in both runs; rtol {a.rtol:g}")
        for n in names:
            x, y = np.asarray(ha[n], dtype=float), np.asarray(hb[n], dtype=float)
            if x.shape != y.shape:
                print(f"\n{n}: shapes differ {x.shape} vs {y.shape}")
                continue
            has_x = np.isfinite(x) & (x != 0)
            has_y = np.isfinite(y) & (y != 0)
            both = has_x & has_y
            only_x, only_y = has_x & ~has_y, has_y & ~has_x
            rel = np.zeros_like(x)
            m = np.maximum(np.abs(x), np.abs(y))
            rel[both] = np.abs(x[both] - y[both]) / m[both]
            bad = both & (rel > a.rtol)
            if a.theta_max is not None and x.ndim == 3:
                keep = (np.arange(x.shape[0]) * a.dtheta <= a.theta_max + 1e-9)[:, None, None]
                both, only_x, only_y, bad = both & keep, only_x & keep, only_y & keep, bad & keep
            if not bad.any() and not only_x.any() and not only_y.any():
                print(f"\n{n}: same within {a.rtol:g} ({int(both.sum())} cells with a value)")
                continue
            print(f"\n{n}: {int(both.sum())} cells with a value in both, {int(bad.sum())} differ by more than "
                  f"{a.rtol:g} (median difference of those {np.median(rel[bad]) if bad.any() else 0:.2g}); "
                  f"value only in A: {int(only_x.sum())}, only in B: {int(only_y.sum())}")
            if x.ndim == 3:
                nt = x.shape[0]
                rows = bad.reshape(nt, -1).sum(axis=1)
                print("   differing cells by theta row (deg: count): " +
                      ", ".join(f"{t * a.dtheta:.0f}: {int(c)}" for t, c in enumerate(rows) if c))
                r = a.r_inner + a.dr * np.arange(x.shape[2])
                parts = []
                for lo, hi in R_RANGES:
                    sel = (r >= lo) & (r < hi)
                    if sel.any():
                        nb, nboth = int(bad[:, :, sel].sum()), int(both[:, :, sel].sum())
                        parts.append(f"{lo:.0f}-{'' if hi > 1e8 else f'{hi:.0f}'} AU: {nb} of {nboth}")
                print("   differing cells by r: " + ", ".join(parts))
                idx = np.argsort(-np.where(bad, rel, 0).ravel())[:a.show]
                for flat in idx:
                    if not bad.ravel()[flat]:
                        break
                    i, j, k = np.unravel_index(flat, x.shape)
                    print(f"   (theta {i * a.dtheta:.0f}, phi {j * a.dphi:.0f}, r {a.r_inner + k * a.dr:.0f} AU): "
                          f"{x[i, j, k]:.6g} vs {y[i, j, k]:.6g}")
        if a.los:
            i, j = int(round(a.los[0] / a.dtheta)), int(round(a.los[1] / a.dphi)) % ha[names[0]].shape[1]
            cols = [n for n in ("ptr_ct", "cool", "vshock", "ntr", "nref", "nen", "esw", "em", "kappa_ref") if n in names]
            ax = {n: np.asarray(ha[n][i, j, :], dtype=float) for n in cols}
            bx = {n: np.asarray(hb[n][i, j, :], dtype=float) for n in cols}
            print(f"\nline of sight theta {i * a.dtheta:.0f}, phi {j * a.dphi:.0f}: values A / B ('-' = none); * = differ by > {a.rtol:g}")
            print("   r (AU) " + "".join(f"{n:>26s}" for n in cols))
            rows = 0
            for k in range(len(ax[cols[0]])):
                vals = [(ax[n][k], bx[n][k]) for n in cols]
                if not any((np.isfinite(u) and u != 0) or (np.isfinite(v) and v != 0) for u, v in vals):
                    continue
                cells = []
                for u, v in vals:
                    m = max(abs(u), abs(v)) if np.isfinite(u) and np.isfinite(v) else 0
                    star = "*" if m > 0 and abs(u - v) / m > a.rtol else " "
                    f = lambda w: "-" if not np.isfinite(w) or w == 0 else f"{w:.4g}"
                    cells.append(f"{f(u)}/{f(v)}{star}")
                print(f"   {a.r_inner + k * a.dr:6.0f} " + "".join(f"{c:>26s}" for c in cells))
                rows += 1
                if rows >= a.max_rows:
                    print("   ... (more: --max-rows)")
                    break
    return 0


if __name__ == "__main__":
    sys.exit(main())
