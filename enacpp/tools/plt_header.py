#!/usr/bin/env python3
"""Prints the header of a Tecplot binary data file (.plt, format version 112): the variables, every
zone (type, size, variable locations, sharing) and the data formats of the first zones, and checks
that the file size matches. Reads only the header and the start of each zone's data, so it is quick
even on files of many GB. Plain Python 3, no packages needed.

    python3 enacpp/tools/plt_header.py FILE.plt
"""
import os
import struct
import sys

ZONE_TYPES = ["ORDERED", "FELINESEG", "FETRIANGLE", "FEQUADRILATERAL", "FETETRAHEDRON", "FEBRICK",
              "FEPOLYGON", "FEPOLYHEDRON"]
NODES_PER_ELEMENT = {1: 2, 2: 3, 3: 4, 4: 4, 5: 8}
FORMATS = {1: ("float", 4), 2: ("double", 8), 3: ("int32", 4), 4: ("int16", 2), 5: ("byte", 1), 6: ("bit", 1)}


class Reader:
    def __init__(self, f, endian):
        self.f, self.e = f, endian

    def i32(self):
        return struct.unpack(self.e + "i", self.f.read(4))[0]

    def f32(self):
        return struct.unpack(self.e + "f", self.f.read(4))[0]

    def f64(self):
        return struct.unpack(self.e + "d", self.f.read(8))[0]

    def string(self):
        out = []
        while True:
            c = self.i32()
            if c == 0:
                return "".join(out)
            out.append(chr(c) if 0 < c < 0x110000 else "?")


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    path = sys.argv[1]
    size = os.path.getsize(path)
    with open(path, "rb") as fh:
        magic = fh.read(8)
        print(f"file: {path} ({size / 2**30:.2f} GB)")
        print(f"magic: {magic!r}")
        if not magic.startswith(b"#!TDV"):
            sys.exit("not a classic Tecplot binary file (a .szplt file starts with #!SZPLT)")
        endian = "<"
        r = Reader(fh, endian)
        if r.i32() != 1:
            endian = ">"
            r = Reader(fh, endian)
        filetype = r.i32()
        title = r.string()
        nv = r.i32()
        names = [r.string() for _ in range(nv)]
        print(f"byte order: {'little' if endian == '<' else 'big'}-endian; file type {filetype} "
              f"({['full', 'grid', 'solution'][filetype] if 0 <= filetype <= 2 else '?'}); title {title!r}")
        print(f"{nv} variables:")
        for i, n in enumerate(names, 1):
            print(f"  {i:3d} {n}")

        zones = []
        while True:
            marker = r.f32()
            if marker == 299.0:
                z = {"name": r.string(), "parent": r.i32(), "strand": r.i32(), "time": r.f64()}
                r.i32()  # not used
                z["type"] = r.i32()
                z["loc"] = [r.i32() for _ in range(nv)] if r.i32() == 1 else [0] * nv
                z["raw_faces"] = r.i32()
                nmisc = r.i32()
                if nmisc != 0:
                    r.i32()  # face neighbour mode
                    if z["type"] != 0:
                        r.i32()  # fully connected
                if z["type"] == 0:
                    z["ijk"] = (r.i32(), r.i32(), r.i32())
                else:
                    z["npts"] = r.i32()
                    if z["type"] in (6, 7):
                        z["poly"] = (r.i32(), r.i32(), r.i32(), r.i32())
                    z["nelem"] = r.i32()
                    r.i32(); r.i32(); r.i32()  # cell dims (future use)
                aux = []
                while r.i32() == 1:
                    name = r.string(); r.i32(); value = r.string()
                    aux.append((name, value))
                z["aux"] = aux
                zones.append(z)
            elif marker == 357.0:
                break
            elif marker == 799.0:  # data-set auxiliary data
                name = r.string(); r.i32(); value = r.string()
                print(f"data-set aux: {name} = {value!r}")
            elif marker == 899.0:  # variable auxiliary data
                var = r.i32(); name = r.string(); r.i32(); value = r.string()
                print(f"variable {var + 1} aux: {name} = {value!r}")
            else:
                sys.exit(f"header item {marker} (geometry, text or custom labels) is not handled by this script; "
                         f"the zones read so far: {len(zones)}")
        print(f"{len(zones)} zone(s); header ends at byte {fh.tell()}")

        # the data section: per zone, formats, sharing, min/max, the data, the connectivity
        for zi, z in enumerate(zones):
            if r.f32() != 299.0:
                sys.exit(f"zone {zi + 1}: data section marker not found at byte {fh.tell() - 4}: header misread?")
            fmt = [r.i32() for _ in range(nv)]
            passive = [r.i32() for _ in range(nv)] if r.i32() == 1 else [0] * nv
            share = [r.i32() for _ in range(nv)] if r.i32() == 1 else [-1] * nv
            share_conn = r.i32()
            nown = [v for v in range(nv) if not passive[v] and share[v] < 0]
            fh.seek(16 * len(nown), 1)  # min and max of each variable that is stored here
            if z["type"] == 0:
                ni, nj, nk = z["ijk"]
                nnode, ncell = ni * nj * nk, max(ni - 1, 1) * max(nj - 1, 1) * max(nk - 1, 1)
            else:
                nnode, ncell = z["npts"], z["nelem"]
            nbytes = 0
            for v in nown:
                count = ncell if z["loc"][v] == 1 else nnode
                nbytes += count * FORMATS.get(fmt[v], ("?", 4))[1]
            if z["type"] in NODES_PER_ELEMENT and share_conn < 0:
                nbytes += 4 * NODES_PER_ELEMENT[z["type"]] * z["nelem"]
            elif z["type"] in (6, 7):
                sys.exit("polygon/polyhedron zones: connectivity size not computed by this script")
            if zi < 3 or zi == len(zones) - 1:
                kind = ZONE_TYPES[z["type"]] if 0 <= z["type"] < len(ZONE_TYPES) else str(z["type"])
                size_txt = f"I x J x K = {z['ijk']}" if z["type"] == 0 else f"{z['npts']} nodes, {z['nelem']} elements"
                print(f"zone {zi + 1} {z['name']!r}: {kind}, {size_txt}, time {z['time']:g}, strand {z['strand']}")
                cells = [names[v] for v in range(nv) if z["loc"][v] == 1]
                print(f"  cell-centred variables: {', '.join(cells) if cells else 'none (all at the nodes)'}")
                print(f"  data formats: {', '.join(sorted(set(FORMATS.get(x, ('?', 0))[0] for x in fmt)))}; "
                      f"shared variables: {sum(1 for s in share if s >= 0)}; connectivity shared from zone "
                      f"{share_conn + 1 if share_conn >= 0 else '-'}; passive: {sum(passive)}")
                if z["aux"]:
                    print("  aux: " + "; ".join(f"{a}={b}" for a, b in z["aux"]))
            elif zi == 3:
                print(f"  ... ({len(zones) - 4} more zones) ...")
            fh.seek(nbytes, 1)
        end = fh.tell()
        print(f"data ends at byte {end}; file size {size}: "
              + ("match, the header was read correctly" if end == size else f"MISMATCH by {size - end} bytes"))


if __name__ == "__main__":
    main()
