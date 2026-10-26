#!/usr/bin/env python3
"""Times steps 2a + 2b (read the fixed neutrals, interpolate them onto the grid) on their own, to
find out why they take hours in PyTecplot batch mode but a minute in the Tecplot window.

    tec360-env -- python3 enacpp/tools/prep_batch_test.py GRID FIXED_NEUTRALS [options]

GRID is the spherical grid (e.g. ENA_sphere_2AU.szplt), FIXED_NEUTRALS the fixed-neutral file of
step 1 (e.g. prep_work/neutral_fixed_000.plt). Options choose how Tecplot is driven:
    (none)            batch mode, grid loaded with the PyTecplot load functions (as the driver)
    --load-by-macro   batch mode, grid loaded with a $!ReadDataSet macro (as the Tecplot window)
    --connect         in an open Tecplot window (Scripting > PyTecplot Connections, accept)
    --macro-only      no Python: writes test_batch.mcr (load + 2a + 2b) for Tecplot's own batch
                      mode:  tec360 -b -p test_batch.mcr
                      With --neutral-source SRC it also writes test_step1.mcr (load SRC + step 1,
                      which writes FIXED), and with --mhd MHD test_mhd.mcr (load the grid + steps
                      2c + 2d), so that step 1 and the MHD interpolation can be timed as well.
Prints the time of each part. If a part runs far beyond a few minutes, stop it (kill -9) and
note which variant hangs.
"""
import argparse
import os
import sys
import time

MACROS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "configs", "macros", "single_ion")


def body(name):
    """The commands of one step macro, without its '#!MC' line."""
    with open(os.path.join(MACROS, name)) as f:
        return "".join(ln for ln in f if not ln.startswith("#!MC"))


def load_commands(path):
    """The $!ReadDataSet commands the Tecplot window records for File > Load Data."""
    path = os.path.abspath(path)
    if path.endswith(".szplt"):
        return ("$!NewLayout\n"
                f"$!ReadDataSet  '\"STANDARDSYNTAX\" \"1.0\" \"FILENAME_FILE\" \"{path}\"'\n"
                "  DataSetReader = 'Tecplot Subzone Data Loader'\n  ReadDataOption = New\n  ResetStyle = Yes\n")
    return f"$!NewLayout\n$!ReadDataSet  '\"{path}\"'\n  ReadDataOption = New\n  ResetStyle = Yes\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("grid")
    ap.add_argument("fixed")
    ap.add_argument("--load-by-macro", action="store_true")
    ap.add_argument("--connect", nargs="?", const=7600, type=int, metavar="PORT")
    ap.add_argument("--macro-only", action="store_true")
    ap.add_argument("--neutral-source", metavar="SRC", help="--macro-only: also write test_step1.mcr")
    ap.add_argument("--mhd", metavar="MHD", help="--macro-only: also write test_mhd.mcr")
    a = ap.parse_args()
    # the fixed-neutral file may not exist yet when step 1 (test_step1.mcr) is to write it
    needed = [a.grid] + ([] if a.neutral_source else [a.fixed]) + [p for p in (a.neutral_source, a.mhd) if p]
    for p in needed:
        if not os.path.exists(p):
            sys.exit(f"not found: {p}")
    varset = f'$!VARSET |NEUTRALFIXED| = "{os.path.abspath(a.fixed)}"\n'

    if a.macro_only:
        with open("test_batch.mcr", "w") as g:
            g.write("#!MC 1410\n" + varset + load_commands(a.grid)
                    + body("2a_read_neutrals.mcr") + body("2b_interpolate_neutrals.mcr"))
        print("wrote test_batch.mcr (grid + steps 2a, 2b: the neutral interpolation)")
        if a.neutral_source:
            with open("test_step1.mcr", "w") as g:
                g.write("#!MC 1410\n" + varset + load_commands(a.neutral_source) + body("1_fix_neutrals.mcr"))
            print(f"wrote test_step1.mcr (load {a.neutral_source} + step 1; writes {a.fixed})")
        if a.mhd:
            with open("test_mhd.mcr", "w") as g:
                g.write("#!MC 1410\n" + f'$!VARSET |MHDFILE| = "{os.path.abspath(a.mhd)}"\n'
                        + load_commands(a.grid) + body("2c_read_mhd.mcr") + body("2d_interpolate_mhd.mcr"))
            print("wrote test_mhd.mcr (grid + steps 2c, 2d: the MHD interpolation)")
        print("run each with Tecplot's own batch mode, e.g.  time tec360 -b -p test_batch.mcr")
        return 0

    import tecplot as tp
    if a.connect:
        tp.session.connect(port=a.connect)
    mode = ("connected to the Tecplot window" if a.connect else "batch mode") + \
           (", grid loaded by macro" if a.load_by_macro else ", grid loaded by PyTecplot")
    print(f"PyTecplot {tp.version_info}, Tecplot {tp.sdk_version_info}; {mode}", flush=True)

    def run(name, text):
        """Runs macro text from a temporary file (PyTecplot runs $!VARSET only from files)."""
        path = os.path.abspath(f"test_{os.getpid()}.mcr")
        with open(path, "w") as g:
            g.write("#!MC 1410\n" + varset + text)
        t = time.time()
        tp.macro.execute_file(path)
        os.remove(path)
        print(f"  {name}: {time.time() - t:.0f} s", flush=True)

    t = time.time()
    if a.load_by_macro:
        run("load grid (macro)", load_commands(a.grid))
    else:
        tp.new_layout()
        (tp.data.load_tecplot_szl if a.grid.endswith(".szplt") else tp.data.load_tecplot)(a.grid)
        print(f"  load grid (PyTecplot): {time.time() - t:.0f} s", flush=True)
    run("2a read neutrals", body("2a_read_neutrals.mcr"))
    ds = tp.active_frame().dataset
    for z in ds.zones():
        print(f"     zone {z.index + 1}: {z.name!r}, {z.zone_type.name}, {z.num_points} points", flush=True)
    run("2b interpolate neutrals", body("2b_interpolate_neutrals.mcr"))
    print("done", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
