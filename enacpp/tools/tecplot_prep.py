#!/usr/bin/env python3
"""Prepare ENA input files from a series of MHD solutions with Tecplot macros.

    tec360-env -- python3 enacpp/tools/tecplot_prep.py enacpp/configs/prep_single_ion.toml
    python3 enacpp/tools/tecplot_prep.py CONFIG --batch         # Tecplot's own batch mode, no PyTecplot
    python3 enacpp/tools/tecplot_prep.py CONFIG --dry-run       # what would run; checks the macros (no Tecplot)
    python3 enacpp/tools/tecplot_prep.py CONFIG --check-only    # check existing output files (no Tecplot)
    options: --first N --last M (snapshot range), --redo (redo snapshots already done),
             --work-dir DIR (intermediate files), --stage-dir DIR (copy the inputs there first),
             --batch [--tecplot CMD] (run the macros with "tec360 -b -p" instead of PyTecplot)

The configuration (e.g. configs/prep_single_ion.toml) names:
  - the inputs of each snapshot ([input.NAME]): a numbered pattern ("mhd/run_{index:03}.szplt"),
    a sorted file listing (glob = "mhd/3d__var_2_*.szplt"), a list file (one path per line, in
    snapshot order), or one fixed file used by every snapshot (the spherical grid);
  - intermediate files ([work_file.NAME]) written by one macro and read by the next;
  - the steps, in order ([[step]]): each runs a Tecplot macro (macro = "...") or a program
    (command = ["program", "{NEUTRALSRC}", ...], e.g. bin/ena_fixneutrals; "{VAR}" stands for the
    file of macro variable VAR). A macro step can first load an input or an intermediate file
    into a new layout (load = "NAME"), as done by hand before running that macro;
  - the output files ([output.NAME]) with the layout the ENA code reads.

Every input, intermediate and output file is handed to the macros as a Tecplot macro variable
(|MHDFILE|, |PLASMAFILE|, ...: the macro_variable of each entry), as an absolute path built from
the directory the driver is started in, so the macros contain no paths of their own and work in
any directory on any machine. |INDEX| (e.g. 007) is set too. Before running, every macro is
scanned for absolute paths left in it; the run stops if one is found.

Tecplot is run in one of two ways:
  - with PyTecplot (the default): Tecplot runs inside this Python process, started through
    tec360-env, and each step's macro is run on its own;
  - with --batch: no PyTecplot is needed. The macro steps between two program steps are written
    into one macro (the macro variables, then for each step the commands that load its input and
    the step's macro), which is run by Tecplot's own batch mode, "tec360 -b -p MACRO", in the work
    folder. For the single-ion configuration that is one Tecplot run per snapshot (steps 1b to 3,
    between ena_fixneutrals and ena_trace). Tecplot writes its messages to batch.log in the work folder; a failed run shows
    the end of it. The resulting files are the same: the same macros run in the same order.

For each snapshot, the outputs are written under temporary names (<name>.part), checked
(grid files: header removed, the rows of theta index above theta_index_max dropped where set, then
one row per grid cell, the expected number of columns, no NaN/Inf; the streamline file: the
expected number of variables and at least one zone), then renamed.

With --stage-dir DIR (or [work] stage_dir), each snapshot's input files are first copied to DIR,
and Tecplot reads the copies. On a cluster, DIR should be the compute node's local disk: Tecplot
reads .szplt files piece by piece as it needs them, which is very slow from a network file system
(measured on NASA's Aitken Rome nodes: the MHD interpolation took 45 minutes reading from
/nobackup, 2 to 4 minutes from the node's disk; copying the 13 GB file took 12 s). The copies of a
snapshot are deleted when it is done; an input given as one fixed file (the grid) is copied once. A snapshot whose
outputs all exist is skipped, so a stopped run can simply be started again. Each finished
snapshot is recorded in prep_index.txt in the output folder, with the input files it came from.
"""
import argparse
import ast
import glob as globmod
import os
import re
import shutil
import subprocess
import sys
import time


# ---------------------------------------------------------------------------------------------
# configuration (a small TOML subset: [section], [[step]], key = value, # comments)

def strip_comment(line):
    """The line without its # comment (a # inside double quotes is kept)."""
    inside = False
    for i, ch in enumerate(line):
        if ch == '"':
            inside = not inside
        elif ch == "#" and not inside:
            return line[:i].strip()
    return line.strip()


def parse_config(path):
    cfg, steps, cur = {}, [], None
    for raw in open(path):
        line = strip_comment(raw)
        if not line:
            continue
        if line == "[[step]]":
            cur = {}
            steps.append(cur)
            continue
        m = re.match(r"^\[([A-Za-z0-9_.]+)\]$", line)
        if m:
            cur = cfg.setdefault(m.group(1), {})
            continue
        m = re.match(r"^([A-Za-z0-9_]+)\s*=\s*(.+)$", line)
        if not m or cur is None:
            sys.exit(f"{path}: cannot read line: {raw.rstrip()}")
        value = {"true": "True", "false": "False"}.get(m.group(2).strip(), m.group(2).strip())
        try:
            cur[m.group(1)] = ast.literal_eval(value)
        except (ValueError, SyntaxError):
            sys.exit(f"{path}: cannot read the value of '{m.group(1)}': {value}")
    cfg["steps"] = steps
    return cfg


def sections(cfg, prefix):
    return {k.split(".", 1)[1]: v for k, v in cfg.items() if k.startswith(prefix + ".")}


def input_resolver(name, spec):
    """(function index -> path, number of snapshots it provides or None for unlimited)."""
    if "pattern" in spec:
        return (lambda i, p=spec["pattern"]: p.format(index=i)), None
    if "file" in spec:
        return (lambda i, p=spec["file"]: p), None
    if "glob" in spec:
        files = sorted(globmod.glob(spec["glob"]))
        return (lambda i, f=files: f[i]), len(files)
    if "list" in spec:
        base = os.path.dirname(os.path.abspath(spec["list"]))
        files = [ln.strip() for ln in open(spec["list"]) if ln.strip() and not ln.strip().startswith("#")]
        files = [f if os.path.isabs(f) else os.path.join(base, f) for f in files]
        return (lambda i, f=files: f[i]), len(files)
    sys.exit(f"[input.{name}] needs one of: pattern, glob, list, file")


# ---------------------------------------------------------------------------------------------
# macro scan: no absolute paths may be left in the macros

PATH_IN_QUOTES = re.compile(r"""["']\s*((?:/|[A-Za-z]:[\\/]|~/)[^"'|]*)["']""")


def scan_macros(steps, variables):
    problems, used = [], set()
    for st in steps:
        if "command" in st:  # a program step: its "{VAR}" arguments count as uses
            used |= {m for arg in st["command"] for m in re.findall(r"\{([A-Za-z0-9_]+)\}", arg)}
            continue
        text = open(st["macro"]).read()
        used |= set(re.findall(r"\|([A-Za-z0-9_]+)\|", text))
        for n, line in enumerate(text.splitlines(), 1):
            for m in PATH_IN_QUOTES.finditer(line):
                problems.append(f"{st['macro']}, line {n}: absolute path {m.group(1)!r}; use a macro variable")
    return problems, used


# ---------------------------------------------------------------------------------------------
# checks of the output files

NUMBER = re.compile(rb"^\s*[-+]?(\d|\.\d)")


def check_grid_file(path, ncols, nrows):
    """One row per grid cell, ncols numbers per row, no NaN or Inf. Returns an error text or ''."""
    rows, bad_cols, first_bad = 0, 0, None
    with open(path, "rb") as f:
        for n, line in enumerate(f, 1):
            if not line.strip():
                continue
            if not NUMBER.match(line):
                return f"line {n} is not a row of numbers (a header? set strip_header = true): {line[:60]!r}"
            low = line.lower()
            if b"nan" in low or b"inf" in low:
                return f"line {n} has NaN or Inf: {line[:80]!r}"
            rows += 1
            if len(line.split()) != ncols:
                bad_cols += 1
                first_bad = first_bad or n
    if bad_cols:
        return f"{bad_cols} rows do not have {ncols} columns (first at line {first_bad})"
    if rows != nrows:
        return f"{rows} rows; expected {nrows} (one per grid cell)"
    return ""


def clean_grid_file(path, strip, column, maximum):
    """Rewrites a grid file in one pass:
    - strip: keep only the rows of numbers (removes the TITLE/VARIABLES/ZONE header Tecplot writes);
    - column, maximum: drop the rows whose value in that column (1-based) is above maximum. Used for
      the theta index: the spherical grid of the macros runs over theta indices 0-59, and indices
      30-59 (theta > 180 deg) wrap around and repeat the sky, so only 0-29 are kept, as done by hand.
    Returns the number of rows dropped."""
    tmp, dropped = path + ".clean", 0
    with open(path, "rb") as f, open(tmp, "wb") as g:
        for line in f:
            is_row = NUMBER.match(line)
            if not is_row:
                if not strip:
                    g.write(line)
                continue
            if column:
                fields = line.split()
                if len(fields) >= column and float(fields[column - 1]) > maximum + 0.5:
                    dropped += 1
                    continue
            g.write(line)
    os.replace(tmp, path)
    return dropped


def check_streamline_file(path, nvars):
    """Tecplot ASCII file: nvars variables, at least one zone, no NaN/Inf in the data."""
    with open(path, "rb") as f:
        text = f.read()
    m = re.search(rb"VARIABLES\s*=(.*?)\n\s*ZONE", text, re.S | re.I)
    if not m:
        return "no VARIABLES ... ZONE header (not a Tecplot ASCII file?)"
    nv = len(re.findall(rb'"[^"]*"', m.group(1))) or len(m.group(1).replace(b",", b" ").split())
    if nv != nvars:
        return f"{nv} variables; expected {nvars}"
    zones = len(re.findall(rb"^\s*ZONE\b", text, re.M | re.I))
    if zones == 0:
        return "no ZONE (no streamlines)"
    data = re.sub(rb'"[^"]*"', b"", text)  # names may contain anything; the numbers may not
    if re.search(rb"(?<![A-Za-z])(nan|inf)(?![A-Za-z])", data, re.I):
        return "NaN or Inf in the data"
    return f"ok:{zones}"


def check_output(o, path, grid_rows):
    if o["kind"] == "grid":
        column = int(o.get("theta_index_column", 0))
        dropped = 0
        if o.get("strip_header", False) or column:
            dropped = clean_grid_file(path, o.get("strip_header", False), column, float(o.get("theta_index_max", 0)))
        err = check_grid_file(path, int(o["columns"]), grid_rows)
        note = f"{grid_rows} rows x {o['columns']} columns" + (f" ({dropped} rows with theta index > "
                                                                f"{o.get('theta_index_max')} dropped)" if column else "")
        return err, "" if err else note
    res = check_streamline_file(path, int(o["variables"]))
    if res.startswith("ok:"):
        return "", f"{res[3:]} streamlines, {o['variables']} variables"
    return res, ""


# ---------------------------------------------------------------------------------------------
# macro text

def load_commands(path):
    """The commands the Tecplot window records for File > Load Data, into a new layout."""
    path = os.path.abspath(path)
    if path.endswith(".szplt"):
        return ("$!NewLayout\n"
                f"$!ReadDataSet  '\"STANDARDSYNTAX\" \"1.0\" \"FILENAME_FILE\" \"{path}\"'\n"
                "  DataSetReader = 'Tecplot Subzone Data Loader'\n  ReadDataOption = New\n  ResetStyle = Yes\n")
    return f"$!NewLayout\n$!ReadDataSet  '\"{path}\"'\n  ReadDataOption = New\n  ResetStyle = Yes\n"


def macro_body(path):
    """The commands of a step macro, without its '#!MC' first line."""
    with open(path) as f:
        return [ln for ln in f if not ln.startswith("#!MC")]


def varset_lines(variables):
    """$!VARSET lines that give every macro variable its file."""
    return [f'$!VARSET |{k}| = "{v}"\n' for k, v in variables.items()]


# ---------------------------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("config")
    ap.add_argument("--first", type=int)
    ap.add_argument("--last", type=int)
    ap.add_argument("--redo", action="store_true", help="redo snapshots whose files already exist")
    ap.add_argument("--dry-run", action="store_true", help="print what would run and check the macros; no Tecplot")
    ap.add_argument("--check-only", action="store_true", help="only check existing output files; no Tecplot")
    ap.add_argument("--connect", nargs="?", const=7600, type=int, metavar="PORT",
                    help="run in an open Tecplot window instead of batch mode (in Tecplot: Scripting > "
                         "PyTecplot Connections..., accept connections; default port 7600)")
    ap.add_argument("--load-by-macro", action="store_true",
                    help="load the inputs with $!ReadDataSet macro commands, as the Tecplot window does, "
                         "instead of the PyTecplot load functions")
    ap.add_argument("--batch", action="store_true",
                    help="run the macros with Tecplot's own batch mode (tec360 -b -p), without PyTecplot")
    ap.add_argument("--tecplot", metavar="CMD", default="tec360",
                    help="--batch: the Tecplot program (default tec360, as set up by the Tecplot module)")
    ap.add_argument("--work-dir", metavar="DIR", help="folder of the intermediate files (overrides [work] dir)")
    ap.add_argument("--stage-dir", metavar="DIR",
                    help="copy each snapshot's input files to DIR first and let Tecplot read the copies "
                         "(e.g. the compute node's local disk; overrides [work] stage_dir)")
    a = ap.parse_args()

    cfg = parse_config(a.config)
    run, grid, steps = cfg.get("run", {}), cfg["grid"], cfg["steps"]
    inputs, works, outputs = sections(cfg, "input"), sections(cfg, "work_file"), sections(cfg, "output_file")
    if not steps or not outputs:
        sys.exit("the configuration needs [[step]] entries and [output_file.NAME] sections")
    grid_rows = int(grid["ntheta"]) * int(grid["nphi"]) * int(grid["nr"])
    outdir = cfg.get("output", {}).get("dir", "prep_out")
    # environment variables such as $TMPDIR are expanded in these two paths
    workdir = os.path.expandvars(a.work_dir or cfg.get("work", {}).get("dir", "prep_work"))
    stagedir = a.stage_dir or cfg.get("work", {}).get("stage_dir")
    stagedir = os.path.expandvars(stagedir) if stagedir else None
    keep_work = cfg.get("work", {}).get("keep", False)

    resolvers = {n: input_resolver(n, s) for n, s in inputs.items()}
    available = [cnt for _, cnt in resolvers.values() if cnt is not None]
    # an empty listing would silently give zero snapshots: say which input found nothing, and where
    empty = [n for n, (_, cnt) in resolvers.items() if cnt == 0]
    if empty:
        sys.exit("no files found for: " + ", ".join(
            f"[input.{n}] {'glob' if 'glob' in inputs[n] else 'list'} = \"{inputs[n].get('glob', inputs[n].get('list'))}\""
            for n in empty) + f"\n(relative paths are taken from the directory the driver runs in: {os.getcwd()})")
    counts = {n: cnt for n, (_, cnt) in resolvers.items() if cnt is not None}
    if len(set(counts.values())) > 1:  # e.g. a neutral solution missing for one MHD solution
        print("WARNING: the input listings have different lengths: "
              + ", ".join(f"{n} {c}" for n, c in counts.items())
              + ". Snapshots pair the files by sorted position, so check the pairs below.")
    first = a.first if a.first is not None else int(run.get("first_index", 0))
    count = run.get("count", min(available) - first if available else None)
    if count is None:
        sys.exit("set [run] count (no input is a glob or list that fixes the number of snapshots)")
    last = a.last if a.last is not None else first + int(count) - 1
    if available and last >= min(available):
        sys.exit(f"snapshots {first}-{last} asked, but the input listings give only {min(available)} files")

    # every macro must exist and contain no absolute paths; the variables must all be used
    names = {s["macro_variable"] for s in list(inputs.values()) + list(works.values()) + list(outputs.values())}
    if not a.check_only:
        for st in steps:
            if ("macro" in st) == ("command" in st):
                sys.exit(f"step '{st.get('name', '?')}': give either macro = \"...\" or command = [...]")
            if "command" in st and st.get("load"):
                sys.exit(f"step '{st.get('name', '?')}': a command step cannot load a file")
            if st.get("load") and st["load"] not in inputs and st["load"] not in works:
                sys.exit(f"step '{st.get('name', '?')}': load = \"{st['load']}\" is neither an [input.] nor a [work_file.]")
        missing = [st["macro"] for st in steps if "macro" in st and not os.path.exists(st["macro"])]
        missing += [st["command"][0] for st in steps if "command" in st and not os.path.exists(st["command"][0])]
        if missing:
            sys.exit("not found: " + ", ".join(missing) + " (a program path is relative to where the driver starts)")
        problems, used = scan_macros(steps, names)
        loaded = {(inputs.get(st["load"]) or works[st["load"]])["macro_variable"] for st in steps if st.get("load")}
        for v in sorted(names - used - loaded):
            problems.append(f"no macro uses |{v}| (and no step loads it)")
        for v in sorted(used - names - {"INDEX", "OUTDIR", "WORKDIR", "NUMZONES", "NUMVARS"}):
            problems.append(f"a macro uses |{v}|, which the configuration does not define")
        if problems:
            print("Problems in the macros:\n  " + "\n  ".join(problems))
            if not a.dry_run:
                return 1
        elif a.dry_run:
            print("macros: ok (no absolute paths; every variable defined and used)")
    if a.dry_run:
        print(f"snapshots: {first}-{last} ({last - first + 1}); "
              + ", ".join(f"{n}: {c} files" for n, c in counts.items()))

    tp = None
    if a.batch and a.connect:
        sys.exit("--batch and --connect cannot be used together")
    if a.batch and not (a.dry_run or a.check_only) and not shutil.which(a.tecplot):
        sys.exit(f"--batch: the Tecplot program '{a.tecplot}' is not found (load the Tecplot module, or give --tecplot)")
    if not (a.dry_run or a.check_only):
        if not a.batch:
            try:
                import tecplot as tp  # PyTecplot (pip install pytecplot); run through tec360-env
            except ImportError:
                sys.exit("PyTecplot is not available: pip install --user pytecplot, and run through tec360-env,"
                         " or use --batch (Tecplot's own batch mode, no PyTecplot)")
        if a.connect:
            # connected mode: the commands run in an open Tecplot window (with its own settings),
            # not in a separate batch-mode Tecplot inside this Python process
            tp.session.connect(port=a.connect)
            print(f"connected to the Tecplot window on port {a.connect}", flush=True)
        os.makedirs(outdir, exist_ok=True)
        os.makedirs(workdir, exist_ok=True)
        if stagedir:
            os.makedirs(stagedir, exist_ok=True)
            print(f"inputs are copied to {stagedir} first ({shutil.disk_usage(stagedir).free / 2**30:.0f} GB free there)",
                  flush=True)

    def load(path):
        """Loads one input file into a new layout, as done by hand before running a macro."""
        if a.load_by_macro:
            loader = os.path.join(workdir, "load.mcr")
            with open(loader, "w") as g:
                g.write("#!MC 1410\n" + load_commands(path))
            tp.macro.execute_file(os.path.abspath(loader))
            os.remove(loader)
            return
        tp.new_layout()
        if path.endswith(".szplt"):
            tp.data.load_tecplot_szl(path)
        else:
            tp.data.load_tecplot(path)

    def run_batch(group, tag, variables, read, work):
        """--batch: runs a sequence of macro steps in one Tecplot process ("tec360 -b -p MACRO"), as
        one macro: the macro variables, then for each step the commands that load its input (as
        done by hand before running that macro) and the step's macro."""
        first_step = steps.index(group[0]) + 1
        path = os.path.abspath(os.path.join(
            workdir, f"tecplot_{tag}_steps{first_step}-{first_step + len(group) - 1}.mcr"))
        with open(path, "w") as g:
            g.write("#!MC 1410\n")
            g.write(f"# written by tecplot_prep.py for snapshot {tag}: the macro variables, then the steps\n")
            g.writelines(varset_lines(variables))
            for st in group:
                g.write(f"# ---- step '{st.get('name', st['macro'])}' ({st['macro']})\n")
                if st.get("load"):
                    g.write(load_commands(read[st["load"]] if st["load"] in read else work[st["load"]]))
                g.writelines(macro_body(st["macro"]))
        # Tecplot's batch mode writes its messages (and errors) to batch.log in the folder it runs in
        log = os.path.join(workdir, "batch.log")
        if os.path.exists(log):
            os.remove(log)
        names = ", ".join(st.get("name", st["macro"]) for st in group)
        ts = time.time()
        proc = subprocess.run([a.tecplot, "-b", "-p", path], cwd=workdir, stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT, universal_newlines=True)
        logtext = open(log, errors="replace").read().splitlines() if os.path.exists(log) else []
        errors = [ln for ln in logtext if re.search(r"\b(err|error|fatal)\b", ln, re.I)]
        if proc.returncode != 0 or errors:
            for ln in proc.stdout.splitlines()[-15:]:
                print(f"[{tag}]     tecplot: {ln}")
            for ln in (errors or logtext)[-15:]:
                print(f"[{tag}]     batch.log: {ln}")
            raise RuntimeError(f"Tecplot (batch) failed in the steps {names} (exit status {proc.returncode}); "
                               f"the macro is {path}")
        if not keep_work:
            os.remove(path)
            if os.path.exists(log):
                os.remove(log)
        print(f"[{tag}]   Tecplot (batch): {names}: {time.time() - ts:.0f} s", flush=True)

    done = failed = skipped = 0
    t_all = time.time()
    for idx in range(first, last + 1):
        tag = f"{idx:03d}"
        final = {n: os.path.join(outdir, o["file"].format(index=idx)) for n, o in outputs.items()}
        if a.check_only:
            for n, o in outputs.items():
                if not os.path.exists(final[n]):
                    print(f"[{tag}] {n}: MISSING {final[n]}")
                    failed += 1
                    continue
                err, note = check_output(o, final[n], grid_rows)
                print(f"[{tag}] {n}: {'FAIL: ' + err if err else 'ok, ' + note}  {final[n]}")
                failed += bool(err)
            continue
        if not a.redo and all(os.path.exists(p) for p in final.values()):
            print(f"[{tag}] skipped: all output files exist")
            skipped += 1
            continue
        src = {n: resolvers[n][0](idx) for n in inputs}
        work = {n: os.path.join(workdir, w["file"].format(index=idx)) for n, w in works.items()}
        variables = {"INDEX": tag, "OUTDIR": os.path.abspath(outdir), "WORKDIR": os.path.abspath(workdir)}
        variables.update({inputs[n]["macro_variable"]: os.path.abspath(p) for n, p in src.items()})
        variables.update({works[n]["macro_variable"]: os.path.abspath(p) for n, p in work.items()})
        variables.update({o["macro_variable"]: os.path.abspath(final[n] + ".part") for n, o in outputs.items()})
        if a.dry_run:
            print(f"[{tag}]")
            for n, p in src.items():
                print(f"   input {n}: {p}" + ("" if os.path.exists(p) else "   (FILE NOT FOUND)"))
            if stagedir:
                print(f"   (the inputs are copied to {stagedir} first, and Tecplot reads the copies)")
            for k, v in variables.items():
                print(f"   |{k}| = {v}")
            for st in steps:
                if "command" in st:
                    print("   run " + " ".join(arg.format(**variables) for arg in st["command"]))
                else:
                    print(f"   run {st['macro']}" + (f"  (after loading '{st['load']}' into a new layout)" if st.get("load") else ""))
            if a.batch:
                print(f"   (--batch: the macros between two programs run as one macro, with {a.tecplot} -b -p)")
            continue
        lost = [p for p in src.values() if not os.path.exists(p)]
        if lost:
            print(f"[{tag}] FAILED: input not found: {', '.join(lost)}")
            failed += 1
            continue
        t0 = time.time()
        current = None
        read = dict(src)  # the files Tecplot reads: the inputs, or their copies in the stage folder
        staged = []       # copies of this snapshot's inputs, deleted when it is done
        try:
            if stagedir:
                for n, p in src.items():
                    dst = os.path.join(stagedir, f"{n}_{os.path.basename(p)}")
                    if "file" not in inputs[n]:  # a fixed input (the grid) is kept for the next snapshots
                        staged.append(dst)
                    if not os.path.exists(dst):
                        ts = time.time()
                        shutil.copyfile(p, dst + ".copying")  # a stopped copy never looks complete
                        os.replace(dst + ".copying", dst)
                        print(f"[{tag}]   copied {n} ({os.path.getsize(dst) / 2**30:.1f} GB) to the stage folder: "
                              f"{time.time() - ts:.0f} s", flush=True)
                    read[n] = dst
                    variables[inputs[n]["macro_variable"]] = os.path.abspath(dst)

            group = []  # --batch: the macro steps since the last program step, run in one Tecplot
            for st in steps + [None]:
                if group and (st is None or "command" in st):
                    current = None
                    run_batch(group, tag, variables, read, work)
                    group = []
                if st is None:
                    break
                current = st  # for the failure report
                ts = time.time()
                if "command" in st:
                    # a program step, e.g. ena_fixneutrals: its output goes to the log as it comes
                    args = [arg.format(**variables) for arg in st["command"]]
                    proc = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
                    for ln in proc.stdout.splitlines():
                        print(f"[{tag}]     {ln}", flush=True)
                    if proc.returncode != 0:
                        raise RuntimeError(f"{os.path.basename(args[0])} failed (exit status {proc.returncode})")
                    print(f"[{tag}]   {st.get('name', args[0])}: {time.time() - ts:.0f} s", flush=True)
                    continue
                if a.batch:
                    group.append(st)
                    continue
                if st.get("load"):
                    load(read[st["load"]] if st["load"] in read else work[st["load"]])
                    print(f"[{tag}]   loaded {st['load']}: {time.time() - ts:.0f} s", flush=True)
                # PyTecplot accepts $!VARSET only inside a macro file, not through execute_command.
                # So each step runs a small wrapper macro in the work folder: the $!VARSET lines,
                # then the step's macro as it is (without its own "#!MC" first line).
                wrapper = os.path.join(workdir, f"step{steps.index(st) + 1}_{tag}.mcr")
                with open(wrapper, "w") as g:
                    g.write("#!MC 1410\n")
                    g.write(f"# written by tecplot_prep.py: the macro variables, then {st['macro']}\n")
                    g.writelines(varset_lines(variables))
                    g.writelines(macro_body(st["macro"]))
                tp.macro.execute_file(os.path.abspath(wrapper))
                if not keep_work:
                    os.remove(wrapper)
                print(f"[{tag}]   {st.get('name', st['macro'])}: {time.time() - ts:.0f} s", flush=True)
            errors = []
            for n, o in outputs.items():
                part = final[n] + ".part"
                if not os.path.exists(part):
                    errors.append(f"{n}: not written (does a macro write |{o['macro_variable']}|?)")
                    continue
                err, note = check_output(o, part, grid_rows)
                if err:
                    errors.append(f"{n}: {err}")
                else:
                    print(f"[{tag}]   {n}: ok, {note}")
            if errors:
                raise RuntimeError("; ".join(errors))
            for n in outputs:
                os.replace(final[n] + ".part", final[n])
            if not keep_work:
                for p in work.values():
                    if os.path.exists(p):
                        os.remove(p)
            with open(os.path.join(outdir, "prep_index.txt"), "a") as f:
                f.write(f"{tag}  " + "  ".join(f"{n}={os.path.abspath(p)}" for n, p in src.items()) + "\n")
            done += 1
            print(f"[{tag}] done in {time.time() - t0:.0f} s", flush=True)
        except Exception as e:  # report, keep the .part and work files for inspection, go on
            failed += 1
            print(f"[{tag}] FAILED: {e}", flush=True)
            if current is not None and tp is not None and "macro" in current:
                # PyTecplot does not say which macro line failed; the loaded variables usually
                # show why (e.g. an input that is not the file the macro was recorded on)
                print(f"[{tag}]   in step '{current.get('name', current['macro'])}'"
                      + (f" (after loading '{current['load']}')" if current.get("load") else ""))
                try:
                    names = [v.name for v in tp.active_frame().dataset.variables()]
                    print(f"[{tag}]   Tecplot had {len(names)} variables loaded: " + ", ".join(names))
                except Exception:
                    print(f"[{tag}]   (no dataset was loaded)")
                print(f"[{tag}]   To see the failing line: load the same input in the Tecplot window and play "
                      f"the macro {os.path.join(os.path.abspath(workdir), 'step' + str(steps.index(current) + 1) + '_' + tag + '.mcr')}",
                      flush=True)
        finally:
            for p in staged:  # this snapshot's input copies (a stopped copy is the .copying file)
                for q in (p, p + ".copying"):
                    if os.path.exists(q):
                        os.remove(q)
    if a.dry_run:
        return 0
    if a.check_only:
        print(f"\n{last - first + 1} snapshots checked: {failed} files missing or failing the checks")
        return 1 if failed else 0
    print(f"\n{done} snapshots done, {skipped} skipped (already done), {failed} failed; {time.time() - t_all:.0f} s")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
