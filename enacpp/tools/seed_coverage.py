#!/usr/bin/env python3
"""Which seeds of a seed file have a streamline in a streamline file (e.g. Tecplot's), per ring of
seed theta, and optionally the seed files to trace the same seeds again with ena_trace.

    python3 enacpp/tools/seed_coverage.py STREAMLINES.dat SEEDS.txt [--write-matched OUT.txt]
                                          [--write-unique OUT.txt]

A streamline belongs to the seed nearest to its first point (within 0.01 AU). Reported per seed
theta: the seeds, how many have a streamline, and how many seeds sit at the same point as an
earlier seed (the copies at the poles).
  --write-matched OUT.txt  the seeds that have a streamline, in the order of the seed file
                           (to trace with ena_trace exactly the seeds the file was made from)
  --write-unique OUT.txt   the seed file without repeated seeds (each pole once)
  --list-missing           list theta and phi of every seed without a streamline
Needs only the Python standard library.
"""
import argparse
import math
import sys


def read_seeds(path):
    seeds = []
    for line in open(path):
        line = line.split("#")[0].split()
        if len(line) >= 3:
            seeds.append(tuple(float(v) for v in line[:3]))
    return seeds


def first_points(path):
    """The first point (x, y, z) of each zone of a Tecplot streamline file."""
    out, want = [], False
    for line in open(path):
        s = line.strip()
        if s.startswith("ZONE"):
            want = True
            continue
        if want and s and (s[0].isdigit() or s[0] in "+-."):
            v = s.split()
            out.append((float(v[0]), float(v[1]), float(v[2])))
            want = False
    return out


def theta_deg(p):
    r = math.sqrt(p[0] ** 2 + p[1] ** 2 + p[2] ** 2)
    return math.degrees(math.acos(max(-1.0, min(1.0, p[2] / r))))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("streamlines")
    ap.add_argument("seeds")
    ap.add_argument("--write-matched")
    ap.add_argument("--write-unique")
    ap.add_argument("--list-missing", action="store_true")
    a = ap.parse_args()
    seeds = read_seeds(a.seeds)
    firsts = first_points(a.streamlines)
    # repeated seeds: same point as an earlier seed
    key = lambda p: (round(p[0], 6), round(p[1], 6), round(p[2], 6))
    seen, repeat = {}, [False] * len(seeds)
    for i, s in enumerate(seeds):
        k = key(s)
        if k in seen:
            repeat[i] = True
        else:
            seen[k] = i
    # each streamline to its nearest seed (repeated seeds share a point: count each streamline once,
    # to the first free copy)
    has = [0] * len(seeds)
    unmatched = 0
    by_point = {}
    for i, s in enumerate(seeds):
        by_point.setdefault(key(s), []).append(i)
    for f in firsts:
        best, bd = None, 1e30
        for k, idx in by_point.items():
            d = math.dist(f, seeds[idx[0]])
            if d < bd:
                best, bd = k, d
        if bd > 0.01:
            unmatched += 1
            continue
        free = [i for i in by_point[best] if not has[i]]
        has[(free or by_point[best])[0]] += 1
    rings = {}
    for i, s in enumerate(seeds):
        t = round(theta_deg(s))
        r = rings.setdefault(t, [0, 0, 0])
        r[0] += 1
        r[1] += has[i] > 0
        r[2] += repeat[i]
    print(f"{a.streamlines}: {len(firsts)} streamlines; {a.seeds}: {len(seeds)} seeds "
          f"({sum(repeat)} at the same point as an earlier seed)")
    print("   seed theta   seeds   with a streamline   repeated seeds")
    for t in sorted(rings):
        n, w, rep = rings[t]
        flag = "" if w == n else f"   <- {n - w} without"
        print(f"   {t:6d} deg   {n:5d}   {w:17d}   {rep:14d}{flag}")
    print(f"   total: {sum(1 for h in has if h)} seeds with a streamline, {len(seeds) - sum(1 for h in has if h)} without"
          + (f"; {unmatched} streamlines start at no seed" if unmatched else ""))
    if a.list_missing:
        # theta and phi (degrees, phi 0-360) of the seeds without a streamline, grouped by theta
        miss = {}
        for i, s in enumerate(seeds):
            if not has[i]:
                ph = math.degrees(math.atan2(s[1], s[0])) % 360.0
                miss.setdefault(round(theta_deg(s)), []).append(round(ph) % 360)
        print("   seeds without a streamline (phi in deg, per seed theta):")
        for t in sorted(miss):
            print(f"   {t:6d} deg: " + " ".join(str(v) for v in sorted(miss[t])))
    if a.write_matched:
        with open(a.write_matched, "w") as g:
            g.write(f"# the seeds of {a.seeds} that have a streamline in {a.streamlines}\n")
            for s, h in zip(seeds, has):
                if h:
                    g.write(f"{s[0]!r} {s[1]!r} {s[2]!r}\n")
        print(f"wrote {a.write_matched}")
    if a.write_unique:
        with open(a.write_unique, "w") as g:
            g.write(f"# {a.seeds} without repeated seeds (each pole once)\n")
            for s, rep in zip(seeds, repeat):
                if not rep:
                    g.write(f"{s[0]!r} {s[1]!r} {s[2]!r}\n")
        print(f"wrote {a.write_unique}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
