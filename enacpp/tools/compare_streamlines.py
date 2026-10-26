#!/usr/bin/env python3
"""Compares two streamline files seed by seed: which streamlines are stuck (never leave their
seed), how long they are, and how far apart the paths of the same seed run.

    python3 enacpp/tools/compare_streamlines.py A.dat B.dat [--write-a-like-b OUT.dat]

Streamlines are paired by their first point (the seed, within 0.01 AU), so the two files may hold
the streamlines in any order and need not have the same seeds. Reported:
  - the streamlines of each file that are stuck, by seed theta;
  - the number of points per streamline, how far out they reach, and the spacing of the points
    (the ENA code averages the points in each grid cell, so the spacing changes the maps);
  - for seeds whose streamline moves in both files: the angle (deg, seen from the Sun) between the
    two paths where they first reach 20, 50, 100, 200 and 400 AU, and the seeds where the paths
    part the most;
  - the variable names of both files, and how much each variable differs between the two files
    at the seed and where the streamlines first reach 20, 50 and 100 AU (where the paths still
    agree, so the values should too).
  --write-a-like-b OUT.dat  writes A's streamlines for only the seeds whose streamline in B moves
                            (drops the seeds B has no streamline for or a stuck one), to test with
                            ena_run how much of the map difference those seeds explain.
Needs only the Python standard library. Reads files in Tecplot point format (x, y, z first).
"""
import argparse
import math
import sys

RADII = (20.0, 50.0, 100.0, 200.0, 400.0)
BINS = (2.0, 20.0, 50.0, 100.0, 200.0, 400.0, 800.0, 1e9)   # r bins of the spacing table (AU)
STUCK_AU = 0.01      # a streamline that never gets farther than this from its seed is stuck
MATCH_AU = 0.01      # seeds closer than this are the same seed


def is_number_row(line):
    s = line.lstrip()
    return bool(s) and (s[0].isdigit() or s[0] in "+-.")


SAMPLE_R = (20.0, 50.0, 100.0)   # radii where all variables are compared (besides the seed)


def read_paths(path):
    """The (x, y, z) points of each zone of a streamline file, in file order, the variable names,
    and for each zone all values at the seed and at the first point at or beyond each SAMPLE_R
    (a dict: "seed" or radius -> list of values)."""
    zones, samples, names = [], [], []
    prev_r, prev_line = 0.0, ""
    in_vars = False
    with open(path) as f:
        for line in f:
            if line.lstrip().upper().startswith("VARIABLES"):
                in_vars = True
            if line.startswith("ZONE"):
                in_vars = False
                zones.append([])
                samples.append({})
            elif zones and is_number_row(line):
                v = line.split(None, 3)
                p = (float(v[0]), float(v[1]), float(v[2]))
                z, smp = zones[-1], samples[-1]
                z.append(p)
                r = math.sqrt(p[0] ** 2 + p[1] ** 2 + p[2] ** 2)
                if len(z) == 1:
                    smp["seed"] = [float(x) for x in line.split()]
                for R in SAMPLE_R:
                    if R not in smp and r >= R:
                        # the values at exactly r = R: linear between this point and the previous one
                        cur = [float(x) for x in line.split()]
                        if len(z) > 1 and prev_r < R:
                            w = (R - prev_r) / (r - prev_r)
                            pv = [float(x) for x in prev_line.split()]
                            cur = [(1 - w) * a + w * b for a, b in zip(pv, cur)]
                        smp[R] = cur
                prev_r, prev_line = r, line
            elif in_vars:
                # VARIABLES = "X AU" "Y AU" ... (one or several lines, before the first ZONE)
                parts = line.split('"')
                names += parts[1::2]
    keep = [i for i, z in enumerate(zones) if z]
    return [zones[i] for i in keep], names, [samples[i] for i in keep]


def angles(p):
    """theta, phi (deg) of a point."""
    r = math.sqrt(p[0] ** 2 + p[1] ** 2 + p[2] ** 2)
    return math.degrees(math.acos(max(-1.0, min(1.0, p[2] / r)))), math.degrees(math.atan2(p[1], p[0])) % 360.0


def is_stuck(z):
    x0, y0, z0 = z[0]
    return all(math.dist(p, (x0, y0, z0)) < STUCK_AU for p in z)


def first_at_radius(z, R):
    """The first point of a path at distance >= R from the Sun, or None."""
    for p in z:
        if p[0] ** 2 + p[1] ** 2 + p[2] ** 2 >= R * R:
            return p
    return None


def angle_between(p, q):
    """The angle (deg) between the directions of two points, seen from the Sun."""
    dot = sum(a * b for a, b in zip(p, q))
    n = math.sqrt(sum(a * a for a in p)) * math.sqrt(sum(b * b for b in q))
    return math.degrees(math.acos(max(-1.0, min(1.0, dot / n))))


def spacing_table(label, zs, stuck):
    """Median distance between consecutive points, per bin of r, and how far the streamlines reach."""
    steps = [[] for _ in BINS[:-1]]
    rmax = []
    for i, z in enumerate(zs):
        if i in stuck:
            continue
        rs = [math.sqrt(p[0] ** 2 + p[1] ** 2 + p[2] ** 2) for p in z]
        rmax.append(max(rs))
        for k in range(1, len(z)):
            rm = 0.5 * (rs[k] + rs[k - 1])
            for b in range(len(BINS) - 1):
                if BINS[b] <= rm < BINS[b + 1]:
                    steps[b].append(math.dist(z[k], z[k - 1]))
                    break
    cells = "  ".join(f"{percentile(v, 0.5):7.2f}" if v else "      -" for v in steps)
    counts = "  ".join(f"{len(v):7d}" for v in steps)
    print(f"   {label}  median spacing (AU): {cells}")
    print(f"   {label}  points:              {counts}")
    if rmax:
        print(f"   {label}  farthest point of a streamline (AU): min {min(rmax):.0f}, median {percentile(rmax, 0.5):.0f}, "
              f"max {max(rmax):.0f}")


def pair_by_seed(za, zb):
    """For each zone of A, the index of the zone of B with the same seed (or None). Uses a grid of
    cells of MATCH_AU, so it is fast for thousands of streamlines."""
    cell = lambda p: tuple(int(math.floor(c / MATCH_AU)) for c in p)
    grid = {}
    for j, z in enumerate(zb):
        grid.setdefault(cell(z[0]), []).append(j)
    used, out = set(), []
    for z in za:
        c, best, bd = cell(z[0]), None, MATCH_AU
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for j in grid.get((c[0] + dx, c[1] + dy, c[2] + dz), ()):
                        d = math.dist(z[0], zb[j][0])
                        if j not in used and d <= bd:
                            best, bd = j, d
        if best is not None:
            used.add(best)
        out.append(best)
    return out


def by_theta(zs, idx):
    counts = {}
    for i in idx:
        t = round(angles(zs[i][0])[0])
        counts[t] = counts.get(t, 0) + 1
    return ", ".join(f"{t} deg: {n}" for t, n in sorted(counts.items())) or "none"


def percentile(v, q):
    v = sorted(v)
    return v[min(len(v) - 1, int(q * (len(v) - 1) + 0.5))] if v else float("nan")


def write_subset(src, out, keep):
    """Copies src to out with only the zones whose index (0-based) is in keep."""
    n, write = -1, True
    with open(src) as f, open(out, "w") as g:
        for line in f:
            if line.startswith("ZONE"):
                n += 1
                write = n in keep
            if write:
                g.write(line)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("a")
    ap.add_argument("b")
    ap.add_argument("--write-a-like-b", metavar="OUT")
    ap.add_argument("--show", type=int, default=10, help="how many of the most different seeds to list")
    a = ap.parse_args()
    za, names_a, sa = read_paths(a.a)
    zb, names_b, sb = read_paths(a.b)
    pair = pair_by_seed(za, zb)
    stuck_a = {i for i, z in enumerate(za) if is_stuck(z)}
    stuck_b = {j for j, z in enumerate(zb) if is_stuck(z)}
    print(f"A = {a.a}: {len(za)} streamlines, {len(stuck_a)} stuck")
    print(f"B = {a.b}: {len(zb)} streamlines, {len(stuck_b)} stuck")
    print(f"same seed in both: {sum(j is not None for j in pair)}; "
          f"only in A: {sum(j is None for j in pair)}; only in B: {len(zb) - sum(j is not None for j in pair)}")
    print(f"stuck in A, by seed theta: {by_theta(za, stuck_a)}")
    print(f"stuck in B, by seed theta: {by_theta(zb, stuck_b)}")
    sa_only = [i for i, j in enumerate(pair) if j is not None and i in stuck_a and j not in stuck_b]
    sb_only = [i for i, j in enumerate(pair) if j is not None and j in stuck_b and i not in stuck_a]
    print(f"same seed, stuck only in A: {len(sa_only)}; stuck only in B: {len(sb_only)}"
          + (f" (by seed theta: {by_theta(za, sb_only)})" if sb_only else ""))
    for label, zs in (("A", za), ("B", zb)):
        lens = [len(z) for z in zs]
        if lens:
            print(f"points per streamline in {label}: min {min(lens)}, median {percentile(lens, 0.5)}, max {max(lens)}")
    print("\npoint spacing along the streamlines that move, by r of the points:")
    print("              r (AU):            " + "  ".join(f"{int(BINS[b]):>3d}-{('' if BINS[b + 1] > 1e8 else int(BINS[b + 1])):<3}"
                                                       for b in range(len(BINS) - 1)))
    spacing_table("A", za, stuck_a)
    spacing_table("B", zb, stuck_b)
    # the paths of the seeds that move in both files
    both = [(i, j) for i, j in enumerate(pair) if j is not None and i not in stuck_a and j not in stuck_b]
    print(f"\nseeds whose streamline moves in both files: {len(both)}")
    print("   angle between the two paths (deg) where they first reach r:")
    print("     r (AU)   compared   median    90%     max    over 1 deg   over 5 deg")
    worst = {}
    for R in RADII:
        ang = []
        for i, j in both:
            p, q = first_at_radius(za[i], R), first_at_radius(zb[j], R)
            if p is not None and q is not None:
                d = angle_between(p, q)
                ang.append(d)
                if R == 100.0:
                    worst[i] = d
        if ang:
            print(f"     {R:6.0f}   {len(ang):8d}   {percentile(ang, 0.5):6.3f}  {percentile(ang, 0.9):6.3f}  "
                  f"{max(ang):6.2f}   {sum(d > 1 for d in ang):10d}   {sum(d > 5 for d in ang):10d}")
        else:
            print(f"     {R:6.0f}   no streamline reaches it in both files")
    if worst:
        print(f"   the {a.show} seeds whose paths part the most at 100 AU (seed theta, phi: angle):")
        for i in sorted(worst, key=worst.get, reverse=True)[:a.show]:
            t, ph = angles(za[i][0])
            print(f"     ({t:5.1f}, {ph:5.1f}): {worst[i]:.2f} deg")
    # ---- the variables ----
    print("\nvariables (column: A | B):")
    for c in range(max(len(names_a), len(names_b))):
        na = names_a[c] if c < len(names_a) else "-"
        nb = names_b[c] if c < len(names_b) else "-"
        print(f"   {c + 1:2d}: {na:<22s} | {nb:<22s}{'' if na == nb else '   <- different name'}")
    ncol = min(len(names_a), len(names_b)) or 3
    print("   relative difference |A - B| / max(|A|, |B|) of each column, over the seeds that move in both")
    print("   files: median / 90% / max, at the seed and where the streamlines first reach r (values")
    print("   interpolated to exactly that r)")
    keys = ["seed"] + list(SAMPLE_R)
    print("   " + " " * 24 + "".join(f"{('seed' if k == 'seed' else f'r = {k:.0f} AU'):>26s}" for k in keys))
    for c in range(ncol):
        cells = []
        for k in keys:
            d = []
            for i, j in both:
                ra, rb = sa[i].get(k), sb[j].get(k)
                if ra is None or rb is None or c >= len(ra) or c >= len(rb):
                    continue
                m = max(abs(ra[c]), abs(rb[c]))
                d.append(abs(ra[c] - rb[c]) / m if m > 0 else 0.0)
            cells.append(f"{percentile(d, 0.5):8.1e}/{percentile(d, 0.9):7.1e}/{max(d):7.1e}" if d else f"{'-':>24s}")
        name = names_a[c] if c < len(names_a) else str(c + 1)
        print(f"   {name[:22]:<22s}  " + "".join(f"{x:>26s}" for x in cells))
    if a.write_a_like_b:
        keep = {i for i, j in enumerate(pair) if j is not None and j not in stuck_b}
        write_subset(a.a, a.write_a_like_b, keep)
        print(f"\nwrote {a.write_a_like_b}: the {len(keep)} streamlines of A whose seed moves in B")
    return 0


if __name__ == "__main__":
    sys.exit(main())
