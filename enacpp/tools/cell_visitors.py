#!/usr/bin/env python3
"""Which streamlines put points into a range of grid cells, as the ENA code bins them.

    python3 enacpp/tools/cell_visitors.py STREAMLINES.dat --theta 12 --phi 18 --r 126 180 [--other OTHER.dat]

The ENA code bins each streamline point to the nearest grid point: theta from the streamline file's
theta column (theta0), r from its r column (r0) and phi from x and y. This lists the streamlines
with points in the cells at the given theta and phi (deg) between the two radii (AU): their seed
(theta, phi), how many points, and the radii they span. It also bins the same points by the
theta computed from x, y, z, to show whether the theta column agrees with the position.
With --other, it shows for each listed seed what the same seed's streamline in the other file
does in those cells (pairs by seed, within 0.01 AU).
Grid: --dtheta 6 --dphi 6 --dr 2 --r-inner 2 (defaults). Needs only the Python standard library.
"""
import argparse
import math
import sys


def zones(path):
    """Yields (names, rows) per zone; rows are lists of floats."""
    names, cur, in_vars = [], None, False
    with open(path) as f:
        for line in f:
            s = line.lstrip()
            if s.upper().startswith("VARIABLES"):
                in_vars = True
            if line.startswith("ZONE"):
                in_vars = False
                if cur:
                    yield names, cur
                cur = []
            elif cur is not None and s and (s[0].isdigit() or s[0] in "+-."):
                cur.append([float(v) for v in line.split()])
            elif in_vars:
                names += line.split('"')[1::2]
    if cur:
        yield names, cur


def col(names, name, fallback):
    return names.index(name) if name in names else fallback


def visits(path, a):
    """seed (rounded x, y, z) -> (seed theta, phi, hits by column theta, hits by xyz theta, r range)"""
    out = {}
    it, ip = round(a.theta / a.dtheta), round(a.phi / a.dphi) % round(360 / a.dphi)
    k0, k1 = round((a.r[0] - a.r_inner) / a.dr), round((a.r[1] - a.r_inner) / a.dr)
    nphi = round(360 / a.dphi)
    for names, rows in zones(path):
        cth, cr = col(names, "theta0", 16), col(names, "r0", 18)
        hit_col = hit_xyz = 0
        rr = []
        for p in rows:
            x, y, z = p[0], p[1], p[2]
            ph = math.degrees(math.atan2(y, x))
            if ph < 0:
                ph = 360 - abs(ph)
            pi = round(ph / a.dphi) % nphi
            if pi != ip:
                continue
            k = round((p[cr] - a.r_inner) / a.dr)
            if not (k0 <= k <= k1):
                continue
            r = math.sqrt(x * x + y * y + z * z)
            if round(p[cth] / a.dtheta) == it:
                hit_col += 1
                rr.append(p[cr])
            if round(math.degrees(math.acos(max(-1, min(1, z / r)))) / a.dtheta) == it:
                hit_xyz += 1
        s = rows[0]
        rs = math.sqrt(s[0] ** 2 + s[1] ** 2 + s[2] ** 2)
        key = (s[0], s[1], s[2])
        st = math.degrees(math.acos(max(-1, min(1, s[2] / rs))))
        sp = math.degrees(math.atan2(s[1], s[0])) % 360
        out[key] = (st, sp, hit_col, hit_xyz, (min(rr), max(rr)) if rr else None)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("streamlines")
    ap.add_argument("--theta", type=float, required=True)
    ap.add_argument("--phi", type=float, required=True)
    ap.add_argument("--r", type=float, nargs=2, required=True, metavar=("RMIN", "RMAX"))
    ap.add_argument("--other")
    ap.add_argument("--dtheta", type=float, default=6.0)
    ap.add_argument("--dphi", type=float, default=6.0)
    ap.add_argument("--dr", type=float, default=2.0)
    ap.add_argument("--r-inner", type=float, default=2.0)
    a = ap.parse_args()
    va = visits(a.streamlines, a)
    vb = visits(a.other, a) if a.other else {}
    print(f"cells theta {a.theta:g}, phi {a.phi:g}, r {a.r[0]:g}-{a.r[1]:g} AU")
    hits = {k: v for k, v in va.items() if v[2] or v[3]}
    print(f"{a.streamlines}: {len(hits)} streamlines with points there "
          f"({sum(v[2] for v in va.values())} points binned by the theta column, "
          f"{sum(v[3] for v in va.values())} by the theta of x, y, z)")
    if a.other:
        hb = {k: v for k, v in vb.items() if v[2] or v[3]}
        print(f"{a.other}: {len(hb)} streamlines with points there "
              f"({sum(v[2] for v in vb.values())} points by the theta column, {sum(v[3] for v in vb.values())} by x, y, z)")
    print("   seed (theta, phi)   points (theta column / xyz theta)   r range (AU)" + ("   | other file: points, r range" if a.other else ""))
    for k, (st, sp, hc, hx, rr) in sorted(hits.items(), key=lambda kv: -kv[1][2]):
        line = f"   ({st:5.1f}, {sp:5.1f})          {hc:4d} / {hx:4d}                    " + (f"{rr[0]:.0f}-{rr[1]:.0f}" if rr else "-")
        if a.other:
            o = None
            for kb, vv in vb.items():  # the same seed in the other file: nearest, within 0.01 AU
                if math.dist(k, kb) < 0.01:
                    o = vv
                    break
            line += "   | " + ("no streamline" if o is None else f"{o[2]} / {o[3]}, " + (f"{o[4][0]:.0f}-{o[4][1]:.0f}" if o[4] else "-"))
        print(line)
    return 0


if __name__ == "__main__":
    sys.exit(main())
