#!/usr/bin/env python3
"""SUPERSEDED by ena_plot.py (ecliptic J2000 maps). Kept for simulation-frame quick looks.

Quick-look sky map of an ena_run flux file (flux_bandNN.h5), Hammer projection.

    python3 enacpp/tools/plot_flux.py cpp_output/flux_band13.h5
    python3 enacpp/tools/plot_flux.py run_b/flux_band13.h5 --ratio run_a/flux_band13.h5
    python3 enacpp/tools/plot_flux.py cpp_output/flux_band13.h5 --component accelerated --center 180

Components: transmitted, reflected, accelerated, solar_wind, total (default).
The map is drawn in simulation coordinates: longitude = phi, latitude =
90 - theta. --center sets the phi at the map centre (default 0).
"""

import argparse
import sys

import numpy as np

try:
    import h5py
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
except ImportError:
    sys.exit("plot_flux.py needs numpy, h5py and matplotlib (pip install numpy h5py matplotlib)")

COMPONENTS = ["transmitted", "reflected", "accelerated", "solar_wind", "total", "pu2"]  # pu2: multi-ion only


def load(path, component):
    with h5py.File(path, "r") as f:
        if component not in f:
            sys.exit(f"{path} has no dataset '{component}'")
        m = f[component][()]
        band = int(f["band"][()]) if "band" in f else None
        settings = f["settings"][()].decode() if "settings" in f else ""
    return m, band, settings


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("file")
    ap.add_argument("--component", default="total", choices=COMPONENTS)
    ap.add_argument("--ratio", metavar="REF.h5", help="plot file / REF instead of the flux")
    ap.add_argument("--center", type=float, default=0.0, help="phi (deg) at the map centre")
    ap.add_argument("--log", action="store_true", help="logarithmic colour scale for the flux")
    ap.add_argument("--out", help="output image (default: <file>_<component>.png)")
    args = ap.parse_args()

    m, band, settings = load(args.file, args.component)
    nt, npp = m.shape
    dth, dph = 180.0 / nt, 360.0 / npp
    # cell edges around the grid points theta = i*dth (i = 0: north pole), phi = j*dph
    lat_e = np.clip(90.0 - (np.arange(nt + 1) - 0.5) * dth, -90, 90)
    lat_e[-1] = -90.0  # the last row (theta = 180 - dth) also covers the south pole
    lon_c = np.arange(npp) * dph
    lon_c = (lon_c - args.center + 180.0) % 360.0 - 180.0
    order = np.argsort(lon_c)
    lon_c = lon_c[order]
    lon_e = np.concatenate([lon_c - dph / 2, [lon_c[-1] + dph / 2]])
    if lon_e[0] < -180.0:  # a cell straddles +/-180: draw its other half on the right
        order = np.concatenate([order, order[:1]])
        lon_e = np.concatenate([lon_e, [180.0]])
    lon_e = np.clip(lon_e, -180, 180)
    data = m[:, order]

    fig = plt.figure(figsize=(10, 5.6))
    fig.patch.set_facecolor("#fcfcfb")
    ax = fig.add_subplot(111, projection="hammer")
    LON, LAT = np.meshgrid(np.radians(lon_e), np.radians(lat_e))
    if args.ratio:
        ref, _, _ = load(args.ratio, args.component)
        ref = ref[:, order]  # same column order as data
        with np.errstate(divide="ignore", invalid="ignore"):
            r = np.where(ref > 0, data / ref, np.nan)
        lr = np.log2(r)
        good = np.abs(lr[np.isfinite(lr)])
        lim = max(float(np.percentile(good, 99)) if good.size else 0.0, 0.01)
        cmap = LinearSegmentedColormap.from_list("div", ["#104281", "#3987e5", "#f0efec", "#e66767", "#a31f1f"])
        im = ax.pcolormesh(LON, LAT, lr, cmap=cmap, vmin=-lim, vmax=lim, shading="flat")
        cb = fig.colorbar(im, ax=ax, orientation="horizontal", shrink=0.6, pad=0.06)
        ticks = np.linspace(-lim, lim, 5)
        cb.set_ticks(ticks)
        cb.set_ticklabels([f"{2 ** t:.3g}" for t in ticks])
        cb.set_label(f"{args.component}: ratio to {args.ratio}")
        med = np.nanmedian(r)
        title = f"{args.component} flux ratio, median {med:.4g}"
    else:
        if args.log:
            with np.errstate(divide="ignore"):
                shown = np.log10(np.where(data > 0, data, np.nan))
            label = f"log10 {args.component} flux (cm^-2 s^-1 sr^-1 keV^-1)"
        else:
            shown = data
            label = f"{args.component} flux (cm^-2 s^-1 sr^-1 keV^-1)"
        im = ax.pcolormesh(LON, LAT, shown, cmap="viridis", shading="flat")
        cb = fig.colorbar(im, ax=ax, orientation="horizontal", shrink=0.6, pad=0.06)
        cb.set_label(label)
        title = f"{args.component} flux"
    if band is not None:
        title += f", band {band:02d}"
    ax.set_title(title, pad=14)
    ax.grid(True, color="#999999", linewidth=0.4, alpha=0.6)
    ax.tick_params(labelsize=7, colors="#555555")
    if settings:
        fig.text(0.5, 0.01, settings, ha="center", fontsize=7, color="#555555")
    out = args.out or args.file.rsplit(".", 1)[0] + f"_{args.component}{'_ratio' if args.ratio else ''}.png"
    fig.savefig(out, dpi=130, facecolor=fig.get_facecolor(), bbox_inches="tight")
    print(f"written {out}")


if __name__ == "__main__":
    main()
