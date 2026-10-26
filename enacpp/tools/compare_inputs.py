#!/usr/bin/env python3
"""Compare ENA input files made two ways (e.g. by hand in Tecplot and by tools/tecplot_prep.py).

    python3 enacpp/tools/compare_inputs.py MINE.dat DRIVER.dat [MINE2 DRIVER2 ...] [--rtol 1e-6]
                                           [--atol-frac 1e-6]

Pairs are given one after the other. For each pair:
  1. IDENTICAL if the files are byte for byte the same (the expected result: the driver keeps every
     row Tecplot writes as it is, apart from the header and the theta > 180 deg rows it drops);
  2. otherwise the rows of numbers are compared in order. Lines that are not rows of numbers
     (TITLE, VARIABLES, ZONE, ...) are compared only by count, so a header kept in one file and
     removed in the other shows up as a note, not as a difference in the data. Reported: the
     number of rows, the first rows that differ, and the largest difference in each column.
  3. streamline files (several ZONEs in both files) are compared streamline by streamline, so one
     streamline with a different number of points does not shift all the rows after it. Reported:
     how many streamlines match, and for each one that differs its seed point and lengths.
A value counts as different if it differs by more than rtol (relative) AND by more than
atol-frac times the largest absolute value of its column (an absolute floor, so that values near
zero, e.g. a velocity component of 0.1 km/s next to 400 km/s, do not fail on their last digit).
Needs only the Python standard library. Exit status 0 if every pair matches.
"""
import argparse
import re
import sys

# a row of numbers: only numbers (with optional exponent) separated by blanks
NUMBER_ROW = re.compile(rb"^\s*[-+]?(\d+\.?\d*|\.\d+)([eEdD][-+]?\d+)?(\s+[-+]?(\d+\.?\d*|\.\d+)([eEdD][-+]?\d+)?)*\s*$")


def rows(path):
    """Yields ('row', [floats]) for each row of numbers and ('other', line) for every other line."""
    with open(path, "rb") as f:
        for line in f:
            if NUMBER_ROW.match(line):
                yield "row", [float(x.replace(b"d", b"e").replace(b"D", b"e")) for x in line.split()]
            elif line.strip():
                yield "other", line.strip()


def scan(path):
    """One pass over a file: the largest absolute value of each column and the number of ZONEs."""
    scale, zones = [], 0
    for kind, v in rows(path):
        if kind == "row":
            if len(v) > len(scale):
                scale += [0.0] * (len(v) - len(scale))
            for c, x in enumerate(v):
                if abs(x) > scale[c]:
                    scale[c] = abs(x)
        elif v.startswith(b"ZONE"):
            zones += 1
    return scale, zones


def differs(u, w, c, rtol, floor):
    """True if u and w differ by more than rtol (relative) and more than the column's floor."""
    d = abs(u - w)
    return d > rtol * max(abs(u), abs(w)) and d > floor[c]


def zones_of(path):
    """The rows of numbers of a streamline file, grouped by ZONE (one list per streamline)."""
    out = []
    for kind, v in rows(path):
        if kind == "row":
            if not out:
                out.append([])
            out[-1].append(v)
        elif v.startswith(b"ZONE"):
            out.append([])
    return [z for z in out if z]


def compare_zones(a, b, rtol, floor, show=5):
    """Compares two streamline files streamline by streamline. Returns True if all match."""
    za, zb = zones_of(a), zones_of(b)
    if len(za) != len(zb):
        print(f"   {len(za)} streamlines in the first file, {len(zb)} in the second; compared in order")
    same = within = 0
    bad = []
    for i, (p, q) in enumerate(zip(za, zb)):
        if p == q:
            same += 1
            continue
        if len(p) == len(q) and not any(differs(u, w, c, rtol, floor)
                                        for x, y in zip(p, q) for c, (u, w) in enumerate(zip(x, y))):
            within += 1
            continue
        bad.append((i, p, q))
    print(f"   {min(len(za), len(zb))} streamlines compared: {same} identical, {within} within tolerance, "
          f"{len(bad)} different")
    for i, p, q in bad[:show]:
        n = min(len(p), len(q))
        first = next((k for k in range(n) if any(differs(u, w, c, rtol, floor)
                                                  for c, (u, w) in enumerate(zip(p[k], q[k])))), n)
        print(f"   streamline {i + 1}: seed ({p[0][0]:.4g}, {p[0][1]:.3g}, {p[0][2]:.4g}); "
              f"{len(p)} points in the first file, {len(q)} in the second; they part at point {first + 1}")
        for label, z in (("first ", p), ("second", q)):
            end = z[-1]
            print(f"      {label} file ends at ({end[0]:.4g}, {end[1]:.4g}, {end[2]:.4g})")
    if len(bad) > show:
        print(f"   ... and {len(bad) - show} more")
    if bad:
        seam_summary(za, zb, bad)
    return not bad and len(za) == len(zb)


def on_seam(z):
    """True if a streamline starts on the y = 0 plane (phi = 0 or 180 deg): the seam of the
    spherical grid, where the direction of the first step can depend on the last digit."""
    x, y, zz = z[0][0], z[0][1], z[0][2]
    return abs(y) <= 1e-6 * max(1e-30, (x * x + zz * zz) ** 0.5)


def stuck(z, tol=0.01):
    """True if a streamline never gets farther than tol AU from its seed."""
    x0, y0, z0 = z[0][0], z[0][1], z[0][2]
    return all(((p[0] - x0) ** 2 + (p[1] - y0) ** 2 + (p[2] - z0) ** 2) ** 0.5 < tol for p in z)


def seed_angles(z):
    """theta, phi (deg) of a streamline's seed."""
    import math
    x, y, zz = z[0][0], z[0][1], z[0][2]
    r = (x * x + y * y + zz * zz) ** 0.5
    return math.degrees(math.acos(max(-1.0, min(1.0, zz / r)))), math.degrees(math.atan2(y, x)) % 360.0


def seam_summary(za, zb, bad, show=5):
    """Where the differing streamlines start, and which ones never leave their seed."""
    seam_all = sum(on_seam(z) for z in za)
    seam_bad = [i for i, p, q in bad if on_seam(p)]
    print(f"   seeds on the y = 0 plane (phi 0 or 180 deg, the seam of the grid): {seam_all} of {len(za)}; "
          f"of the {len(bad)} different streamlines, {len(seam_bad)} start there")
    others = [(i, p, q) for i, p, q in bad if not on_seam(p)]
    if others:
        print(f"   different streamlines starting elsewhere: {len(others)}, e.g. seeds at (theta, phi) "
              + ", ".join("(%.0f, %.0f)" % seed_angles(p) for _, p, _ in others[:show]))
    stuck_a = {i for i, z in enumerate(za) if stuck(z)}
    stuck_b = {i for i, z in enumerate(zb) if stuck(z)}
    for label, zs, st in (("first", za, stuck_a), ("second", zb, stuck_b)):
        # how many stuck seeds per seed theta (the seeds lie on rings of constant theta)
        by_theta = {}
        for i in st:
            t = round(seed_angles(zs[i])[0])
            by_theta[t] = by_theta.get(t, 0) + 1
        print(f"   streamlines that never leave their seed (within 0.01 AU) in the {label} file: {len(st)}"
              + ("; by seed theta: " + ", ".join(f"{t} deg: {n}" for t, n in sorted(by_theta.items())) if st else ""))
    if stuck_a != stuck_b:
        print(f"   stuck in one file only: {len(stuck_a ^ stuck_b)} streamlines "
              f"({len(stuck_a - stuck_b)} only in the first, {len(stuck_b - stuck_a)} only in the second)")


def same_bytes(a, b):
    with open(a, "rb") as f, open(b, "rb") as g:
        while True:
            x, y = f.read(1 << 20), g.read(1 << 20)
            if x != y:
                return False
            if not x:
                return True


def compare(a, b, rtol, floor, show=5):
    """Compares the numeric rows of a and b in order. Returns True if they agree within tolerance."""
    ra, rb = rows(a), rows(b)
    n = nother_a = nother_b = ndiff = 0
    zones_a = zones_b = 0
    maxabs, maxrel = {}, {}
    ok = True
    while True:
        # next numeric row of each file, counting the other lines on the way
        xa = xb = None
        for kind, v in ra:
            if kind == "row":
                xa = v
                break
            nother_a += 1
            zones_a += v.startswith(b"ZONE")
        for kind, v in rb:
            if kind == "row":
                xb = v
                break
            nother_b += 1
            zones_b += v.startswith(b"ZONE")
        if xa is None or xb is None:
            break
        n += 1
        if len(xa) != len(xb):
            print(f"   row {n}: {len(xa)} columns in the first file, {len(xb)} in the second")
            return False
        bad = False
        for c, (u, w) in enumerate(zip(xa, xb)):
            d = abs(u - w)
            r = d / max(abs(u), abs(w)) if d else 0.0
            if d > maxabs.get(c, 0.0):
                maxabs[c] = d
            if r > maxrel.get(c, 0.0):
                maxrel[c] = r
            bad |= differs(u, w, c, rtol, floor)
        if bad:
            ndiff += 1
            if ndiff <= show:
                print(f"   row {n} differs:\n      {' '.join(f'{v:.9g}' for v in xa)}\n      {' '.join(f'{v:.9g}' for v in xb)}")
    rest_a = xa is not None or any(k == "row" for k, _ in ra)
    rest_b = xb is not None or any(k == "row" for k, _ in rb)
    if rest_a or rest_b:
        print(f"   the {'first' if rest_a else 'second'} file has more rows (both agree on the first {n})")
        ok = False
    print(f"   {n} rows compared; {ndiff} differ beyond the tolerance")
    if nother_a != nother_b:
        print(f"   note: {nother_a} non-number lines (header/ZONE) in the first file, {nother_b} in the second")
    if zones_a != zones_b:
        # a grid file with the Tecplot header in one copy only has 1 ZONE vs 0: just a note.
        # Two streamline files with different numbers of zones (streamlines) are a difference.
        print(f"   {'ZONES differ' if zones_a and zones_b else 'note: ZONE lines'}: "
              f"{zones_a} in the first file, {zones_b} in the second")
        ok &= not (zones_a and zones_b)
    if maxrel:
        print("   largest relative difference per column: "
              + "  ".join(f"{c + 1}:{maxrel[c]:.2g}" for c in sorted(maxrel)))
    return ok and ndiff == 0


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("files", nargs="+", help="pairs: MINE DRIVER [MINE DRIVER ...]")
    ap.add_argument("--rtol", type=float, default=1e-6, help="relative tolerance per value (default 1e-6)")
    ap.add_argument("--atol-frac", type=float, default=1e-6,
                    help="absolute floor, as a fraction of the largest |value| of the column (default 1e-6)")
    a = ap.parse_args()
    if len(a.files) % 2:
        sys.exit("give the files in pairs: MINE DRIVER [MINE DRIVER ...]")
    all_ok = True
    for mine, drv in zip(a.files[::2], a.files[1::2]):
        print(f"{mine}  vs  {drv}")
        if same_bytes(mine, drv):
            print("   IDENTICAL")
            continue
        scale, zones_a = scan(mine)
        _, zones_b = scan(drv)
        floor = [a.atol_frac * x for x in scale] + [0.0] * 100
        if zones_a > 1 and zones_b > 1:  # streamline files: one zone per streamline
            ok = compare_zones(mine, drv, a.rtol, floor)
        else:
            ok = compare(mine, drv, a.rtol, floor)
        print(f"   SAME within tolerance (rtol {a.rtol:g}, atol-frac {a.atol_frac:g})" if ok else "   DIFFERENT")
        all_ok &= ok
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
