#!/usr/bin/env python3
"""Writes a copy of a run configuration that reads the HDF5 snapshots made by
ena_convert instead of the text files: [inputs] format = "hdf5" and dir = H5DIR,
every [[snapshot]] entry gets file = "snapshot_NNN.h5" (in order, from 000), and a
[time_dependent.files] section gets file = "snapshot_{index:03}.h5" with
first_index = 0 (or file = "snapshot_000.h5" if its file names have no {index},
i.e. one snapshot repeated, as written by td_repeat_snapshot.py). Everything else
is copied unchanged.

usage: python3 use_hdf5.py CONFIG.toml OUT.toml H5DIR
OUT.toml is written next to CONFIG.toml; a relative H5DIR is taken relative to the
current directory and written as an absolute path."""
import os, re, sys

cfg, out, h5dir = sys.argv[1], sys.argv[2], os.path.abspath(sys.argv[3])
lines = open(cfg).read().splitlines()
# does [time_dependent.files] use {index}? (if not, every snapshot is the same file)
indexed, sec = False, ""
for ln in lines:
    if ln.strip().startswith("["):
        sec = ln.strip()
    elif sec == "[time_dependent.files]" and "{index" in ln:
        indexed = True
res, section, nsnap = [], "", 0
file_keys = re.compile(r'\s*(plasma|neutral|boundaries|streamlines|file|first_index)\s*=')
for ln in lines:
    s = ln.strip()
    if s.startswith("["):
        section = s
        res.append(ln)
        if s == "[[snapshot]]":
            res.append(f'file = "snapshot_{nsnap:03d}.h5"')
            nsnap += 1
        elif s == "[time_dependent.files]":
            if indexed:
                res.append('file = "snapshot_{index:03}.h5"')
                res.append("first_index = 0")
            else:
                res.append('file = "snapshot_000.h5"')
        continue
    if section == "[inputs]" and re.match(r'\s*(format|dir)\s*=', ln):
        continue  # replaced below
    if (section == "[[snapshot]]" or section == "[time_dependent.files]") and file_keys.match(ln):
        continue  # text file names (and the old pattern) are replaced
    res.append(ln)
# [inputs] settings right after its header
i = next((k for k, l in enumerate(res) if l.strip() == "[inputs]"), None)
if i is None:
    res += ["", "[inputs]"]
    i = len(res) - 1
res[i + 1:i + 1] = ['format = "hdf5"', f'dir = "{h5dir}"']
path = os.path.join(os.path.dirname(os.path.abspath(cfg)), os.path.basename(out))
open(path, "w").write("\n".join(res) + "\n")
print(f"wrote {path}: HDF5 snapshots from {h5dir}")
