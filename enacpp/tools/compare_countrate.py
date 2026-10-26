#!/usr/bin/env python3
"""Compare the legacy count-rate files (*CountRateKappa*EbandNN.dat) of two runs.

    python3 compare_countrate.py RUN_A_DIR RUN_B_DIR [--nphi 60]

Every file name found in RUN_A_DIR is looked up in RUN_B_DIR. For each file
(component tag and band) it reports whether the maps are identical and, if
not, the ratio B/A: median, range, number of pixels changed by more than 1e-6
and 1%, and where the largest change is (theta, phi indices).
Exit code 0 if all files are identical, 1 otherwise.
"""
import argparse
import glob
import os
import re
import sys

import numpy as np


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_a")
    ap.add_argument("run_b")
    ap.add_argument("--nphi", type=int, default=60, help="phi points per theta row (default 60)")
    args = ap.parse_args()

    files = sorted(glob.glob(os.path.join(args.run_a, "*CountRate*Eband*.dat")))
    if not files:
        sys.exit(f"no *CountRate*Eband*.dat files in {args.run_a}")
    all_same = True
    print(f"{'file':60s} {'result':>10s}  details")
    for fa in files:
        name = os.path.basename(fa)
        fb = os.path.join(args.run_b, name)
        if not os.path.exists(fb):
            print(f"{name:60s} {'MISSING':>10s}  not in {args.run_b}")
            all_same = False
            continue
        a = np.loadtxt(fa).ravel()
        b = np.loadtxt(fb).ravel()
        if a.shape != b.shape:
            print(f"{name:60s} {'SHAPE':>10s}  {a.size} vs {b.size} values")
            all_same = False
            continue
        if np.array_equal(a, b):
            print(f"{name:60s} {'identical':>10s}")
            continue
        all_same = False
        den = np.maximum(np.abs(a), np.abs(b))
        rel = np.where(den > 0, np.abs(b - a) / np.where(den > 0, den, 1), 0.0)
        m = a > 0
        ratio = b[m] / a[m]
        i = int(np.argmax(rel))
        t, p = divmod(i, args.nphi)
        # B/A only where A > 0; a map that is zero in A (e.g. a component switched off there) has no ratio
        rtxt = (f"B/A median {np.median(ratio):.4g}, range {ratio.min():.4g}-{ratio.max():.4g}" if ratio.size
                else f"A is zero everywhere, B max {np.abs(b).max():.4e} (median {np.median(b):.4e})")
        print(f"{name:60s} {'DIFFERS':>10s}  {rtxt}; "
              f"{int((rel > 1e-6).sum())} pixels > 1e-6, {int((rel > 0.01).sum())} > 1%; "
              f"largest at theta index {t}, phi index {p} ({a[i]:.4e} -> {b[i]:.4e})")
    print("\nAll files identical." if all_same else "\nSome files differ.")
    return 0 if all_same else 1


if __name__ == "__main__":
    sys.exit(main())
