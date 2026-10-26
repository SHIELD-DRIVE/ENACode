#!/usr/bin/env python3
"""Writes a time-dependent configuration that repeats the first snapshot of a
steady-state configuration N times (the steady-state limit test of the
time-dependent method). Every [[snapshot]] entry is replaced by a
[time_dependent] section whose files are those of the first snapshot.

usage: python3 td_repeat_snapshot.py STEADY.toml OUT.toml [--n 30] [--dt 1.0]
           [--start 2000] [--cadence 2.0]
OUT.toml is written next to STEADY.toml (relative paths stay valid)."""
import argparse, os, re

ap = argparse.ArgumentParser()
ap.add_argument("steady")
ap.add_argument("out")
ap.add_argument("--n", type=int, default=30, help="number of snapshots")
ap.add_argument("--dt", type=float, default=1.0, help="snapshot step (years)")
ap.add_argument("--start", type=float, default=2000.0)
ap.add_argument("--cadence", type=float, default=2.0, help="output window length and cadence (years)")
a = ap.parse_args()

text = open(a.steady).read()
lines = text.splitlines()
# split off the [[snapshot]] blocks; keep the first one's file entries
head, snap, in_snap, seen = [], {}, False, 0
for ln in lines:
    s = ln.strip()
    if s.startswith("[[snapshot]]"):
        in_snap, seen = True, seen + 1
        continue
    if s.startswith("[") and not s.startswith("[["):
        in_snap = False
    if s.startswith("[time_dependent"):
        raise SystemExit("the configuration already has a [time_dependent] section")
    if in_snap:
        m = re.match(r'\s*(plasma|neutral|boundaries|streamlines|file)\s*=\s*(".*?")', ln)
        if m and seen == 1:
            snap[m.group(1)] = m.group(2)
        continue
    head.append(ln)
if not snap:
    raise SystemExit("no [[snapshot]] entry found")
end = a.start + a.n * a.dt
out = "\n".join(head).rstrip() + f"""

# Steady-state limit test: snapshot 1 of {os.path.basename(a.steady)} repeated {a.n} times.
[time_dependent]
start = {a.start!r}
end = {end!r}
dt = {a.dt!r}
output_start = {a.start!r}
output_end = {end!r}
cadence = {a.cadence!r}

[time_dependent.files]
""" + "".join(f"{k} = {v}\n" for k, v in snap.items())
path = os.path.join(os.path.dirname(os.path.abspath(a.steady)), os.path.basename(a.out))
open(path, "w").write(out)
print(f"wrote {path}: {a.n} snapshots, {a.start:g} to {end:g}, windows of {a.cadence:g} yr")
