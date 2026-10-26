#!/usr/bin/env python3
"""Bit-identity check of the flux_band*.h5 maps and legacy .dat files of run directories.

    python3 same_maps.py REF_DIR OTHER_DIR [OTHER_DIR ...]

Prints, for each other directory, the maps and files that differ from REF_DIR (exit 1 if any).
Used to check that the maps do not depend on the thread count.
"""
import sys, glob, os, filecmp
import h5py, numpy as np
ref = sys.argv[1]; ok = True
for other in sys.argv[2:]:
    n = 0; bad = []
    for f in sorted(glob.glob(f"{ref}/flux_band*.h5")):
        g = os.path.join(other, os.path.basename(f))
        with h5py.File(f) as a, h5py.File(g) as b:
            for k in a:
                if a[k].ndim == 2:
                    n += 1
                    if not np.array_equal(a[k][()], b[k][()]): bad.append((os.path.basename(f), k))
    dats = sorted(glob.glob(f"{ref}/*.dat"))
    dbad = [os.path.basename(f) for f in dats if not filecmp.cmp(f, os.path.join(other, os.path.basename(f)), shallow=False)]
    print(f"{other}: {n} maps, {len(dats)} .dat files; differing: {bad + dbad if bad or dbad else 'none'}")
    ok &= not bad and not dbad
sys.exit(0 if ok else 1)
