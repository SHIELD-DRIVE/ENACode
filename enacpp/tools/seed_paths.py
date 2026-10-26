#!/usr/bin/env python3
"""Where the streamlines from chosen seeds go: number of points, the point closest to the north
and south poles, and the last point, for each streamline whose seed lies on the given meridian.

    python3 enacpp/tools/seed_paths.py STREAMLINES.dat --seed-phi 198 [--seed-theta-max 90]

A streamline that passes over a pole shows a small minimum (or large maximum) theta and a last
point on the far side (phi changed by about 180 deg). Needs only the Python standard library.
"""
import argparse
import math
import sys


def angles(p):
    r = math.sqrt(p[0] ** 2 + p[1] ** 2 + p[2] ** 2)
    return r, math.degrees(math.acos(max(-1, min(1, p[2] / r)))), math.degrees(math.atan2(p[1], p[0])) % 360


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("streamlines")
    ap.add_argument("--seed-phi", type=float, required=True)
    ap.add_argument("--seed-theta-max", type=float, default=180.0)
    a = ap.parse_args()
    print(f"{a.streamlines}: streamlines seeded at phi {a.seed_phi:g} deg")
    print("   seed theta   points   closest to the north pole (r, theta, phi)   last point (r, theta, phi)")
    zone = []

    def report(z):
        if not z:
            return
        r0, t0, p0 = angles(z[0])
        dp = abs(p0 - a.seed_phi) % 360
        if min(dp, 360 - dp) > 0.5 or t0 > a.seed_theta_max + 0.5:
            return
        best = min(z, key=lambda p: angles(p)[1])
        rb, tb, pb = angles(best)
        rl, tl, pl = angles(z[-1])
        print(f"   {t0:8.1f}   {len(z):6d}   ({rb:6.1f}, {tb:5.1f}, {pb:5.1f})                       ({rl:6.1f}, {tl:5.1f}, {pl:5.1f})")

    with open(a.streamlines) as f:
        for line in f:
            if line.startswith("ZONE"):
                report(zone)
                zone = []
                continue
            s = line.lstrip()
            if s and (s[0].isdigit() or s[0] in "+-."):
                v = s.split(None, 3)
                zone.append((float(v[0]), float(v[1]), float(v[2])))
    report(zone)
    return 0


if __name__ == "__main__":
    sys.exit(main())
