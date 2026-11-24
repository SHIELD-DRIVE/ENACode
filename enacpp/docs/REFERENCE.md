# ENA code (C++ version 2.8): full reference

Every program, setting and check in detail. Start with the user guide (`docs/ENA_USER_GUIDE.md`).

C++ ENA line-of-sight code: heliospheric ENA sky maps from MHD solutions, for the single-ion four-population model (Kornbleuth+2023), the extended single-ion model and the multi-ion model (Kornbleuth+2026), as steady-state or time-dependent maps, on many cores (OpenMP) and nodes (MPI).

The code started as a port of an IDL code. `CHANGES_FROM_IDL.md` lists every difference from it and its effect on the maps; to reproduce the IDL results, use the version 1.x archives (`enacpp_v1.0_faithful.zip`, or `enacpp_v1.2_esi.zip` with the multi-ion and extended single-ion models), whose faithful mode matches the IDL code.

Versions:
- **2.2:** v_shock with the regions TS is the flow speed at the first point past the region crossing with M^2 < 0.64 (`physics.vshock = "first_below"`); `ena_run` prints a diagnostic of the flow speeds around the crossings.
- **2.3:** the multi-ion model (below).
- **2.4:** the extended single-ion model (below).
- **2.5:** survival losses in region 6 as well as region 3, in all three models; the `[region6]` settings of the extended single-ion model.
- **2.6:** OpenMP threads; the results are bit-identical for any thread count (see Threads below).
- **2.7:** time-dependent maps from a sequence of snapshots (see Time-dependent maps below).
- **2.8:** MPI ranks for time-dependent maps, and HDF5 snapshots of every model (see MPI below); the maps are bit-identical for any layout.

What it contains:
- **Run configuration** (TOML) with a snapshot list for time-dependent runs.
- **Grid and field storage**, with r varying fastest.
- **Readers** for the plasma, neutral, `boundaries.txt` and Tecplot streamline files, plus an ASCII-to-HDF5 converter (`ena_convert`).
- **Streamline processing:** termination-shock partition fits, extinction, binning to the grid and gap filling.
- **Flux** and band integration, for all usable bands: IBEX-Lo 0-7, IBEX-Hi 9-13, INCA 14-17.
- **Checkpoint dumps** C0-C5 (HDF5) of the intermediate arrays (`ena_run --dump`), for debugging.

## Multi-ion model (version 2.3)

`configs/run_mi.toml` runs the multi-ion model (Kornbleuth+2026). The model works as follows:
- **Inputs:** three ion fluids (SW; pu3, PUIs picked up in the supersonic wind; pu2, PUIs created in the heliosheath), total H plus four neutral populations, and streamlines with pu3. These are the MI file layouts; the temperatures are not halved.
- **pu3 split:** at the region-based termination shock, pu3 is split into transmitted, reflected and DSA-accelerated parts using the hybrid fits. The values are carried along the streamline. The scaling is `multi_ion.ts_scaling = "compressed"` (default, tailored to the hybrid fits) or `"ratio"`.
- **Shock frame and injection:** v_shock is the pu3 speed at the first point past the crossing below 215 km/s. The DSA injection energy comes from the SW speed at 10 AU on the line of sight (`multi_ion.injection = "streamline"` carries the speed just upstream of the TS instead).
- **Reconnection region:** the pu2 and SW slots are replaced by two kappa populations. The settings are in `[multi_ion]`; the radius and longitude limits can be removed.
- **Outputs:** `flux_bandNN.h5` with `transmitted`, `reflected`, `accelerated`, `pu2`, `solar_wind`, `total`, and the legacy count-rate files with the tags PUI1, PUI2, PUI3, PUI4, SW, CombPUIs.
- **Test:** `./build/test_physics_mi PATH/U_physics_mi.h5` checks the MI physics functions against the unit-test dump of the original IDL code (agreement about 1e-7).

## Extended single-ion model (version 2.4)

`model.populations = "extended_si"` runs the extended single-ion model; `configs/run_esi.toml` is the configuration for the reference inputs. It is the single-ion model with these additions:
- **Inputs.** The single-ion plasma file (12 columns) and streamlines (19 variables), with the 41-column neutral file (total H and populations 1-4). Sensors: the multi-ion table (IBEX-Lo at the CG-corrected central energies, Ultra/INCA).
- **Components.** `transmitted`, `reflected`, `accelerated`, `heliosheath_pui`, `solar_wind`, `ne4`, `neu`, `total`; legacy tags PUI1, PUI2, PUI3, PUI4, SW, Ne4, Neu, CombPUIs.
- **Physics.**
  - Cross section at the ion-neutral relative speed.
  - Heliosheath PUIs: n_p − n1 − n2 − n3, T = ½ m_p |u_p − u_H|²/k.
  - Neutral populations 1 and 2 at the inner boundary.
  - E_l = ½ m_p v² with v the flow speed inside 12 AU on the streamline.
- **Other choices.** Region-6 cells get the region-3 populations (see Region 6 below), no survival losses in region 2, neutral velocities in km/s for all populations, no reconnection-region test populations.
- **Region 6 (version 2.5).** `[region6]` sets the kappas of the region-6 transmitted PUIs, reflected PUIs, heliosheath PUIs and SW, and an accelerated fraction of these ions (`accelerated_fraction` f, `accelerated_energy_share` e, `accelerated_kappa`, `accelerated_from`), with density and pressure conserved; output component `region6_accelerated`. `physics.kappa_heliosheath_pui` sets the region-3 heliosheath PUIs (default `kappa_oh`). All defaults give the version 2.4 populations. See `configs/run_esi.toml`.
- **Original IDL code.** The version 1.3 archive reproduces it.

## Threads (OpenMP, version 2.6)

The code runs on all cores of a node with OpenMP. The work is split into independent items:
- lines of sight in the flux, the survival exponent and the gap fill;
- grid cells for the extinction;
- streamline points;
- blocks of text lines when the inputs are read.

Each item writes only its own results, and sums along a line of sight stay in order. The maps and checkpoint dumps are therefore **bit-identical for any thread count**, including a build without OpenMP, and an input or physics error is reported for the same line or cell as in a serial run. Binning the streamline points into cells stays serial, because its sums run over the points in file order.

- **Thread count:** `ena_run ... --threads N`, or the `OMP_NUM_THREADS` environment variable. The default is all cores; the log prints the count.
- **MPI ranks × threads:** see MPI below.
- **Build without OpenMP:** `cmake -DENA_WITH_OPENMP=OFF`.

Measured on 4 cores with full-size synthetic inputs (30 × 60 × 736 grid, 1800 streamlines, about 3× the heliosheath cells of the reference model, 23 energies):
| Stage | Version 2.5 | Version 2.6, 1 thread | Version 2.6, 4 threads |
|---|---|---|---|
| Total | 21 s | about 16 s | about 7 s |
| Flux, 23 energies | 6.8 s | 6.6 s | 1.8 s |

Most of the single-thread gain comes from keeping freed memory for reuse (fewer page faults). The flux scales almost linearly with threads.

## Time-dependent maps (version 2.7)

With a `[time_dependent]` section in the configuration, `ena_run` makes ENA maps per energy band **and per time window** from a sequence of snapshots (plasma, neutral, boundaries and streamline files at a constant step dt). The method is guided by the time-dependent flux of Zirnstein et al. (2015, Sec. 2.2.4, Eq. 2); `include/ena/time_dependent.hpp` describes it in full and `configs/run_esi_td.toml` is a commented example. All times are in years.

**Time cells.** Snapshot s holds, unchanged, for [start + s·dt, start + (s+1)·dt), with no interpolation in time. N = (end − start)/dt snapshots must be given.

**What is computed.**
- **Production:** for every snapshot, energy, cell and component, the ENA production of the steady-state model, without the survival factor.
- **Travel:** ENAs move radially inward at v = (2E/m_p)^½ to the observer at `observer_r_au` (1 AU). From radius r they arrive after τ = (r − r_obs)/v. Gravity and radiation pressure are neglected: for H they nearly balance, and the deflection matters only below about 0.1 keV.
- **Survival along the path:** an ENA made in cell k at the centre of time cell s passes a closer cell k′ at the time t + (r_k − r_k′)/v. The loss term of k′ (regions 3 and 6, from `physics.survival_observer_r` outward, as in the steady state) is taken from the snapshot of that time; after the last snapshot, from the last one. The exponent is summed with per-snapshot prefix sums, so a path costs one difference per snapshot it crosses.
- **Arrival:** the ENAs of cell k made during time cell s arrive during [start + s·dt + τ_k, start + (s+1)·dt + τ_k), at that cell's steady flux rate. Each output window receives the part that overlaps it, divided by the window length: **the flux averaged over the window**.

**Output windows.** Windows run from `output_start` to `output_end` every `cadence`, each `window` long (default: the cadence). For example, 1-year maps (`cadence = 1`), 6-month maps (`cadence_months = 6`), or 1-year running averages every 6 months (`window = 1`, `cadence_months = 6`). The output range must lie inside [start, end].

**Coverage.** ENAs made before the first snapshot are missing from the first windows. A line of sight is marked **incomplete** (a region with missing flux) in a window by either of two tests:
- **Extinction test** (`coverage_extinction`, default 0.3, about the cooling length): ENAs arriving at the window start a from beyond r_cov = r_obs + v (a − start) are missing. The line is incomplete if the streamline extinction e^−τ of those heliosheath cells still reaches the limit.
- **Missing-flux test** (`coverage_missing_flux`, default 0 = off; the share is always written): the extinction does not reduce the heliosheath-created PUIs (extended single-ion model) or the SW share, so the missing flux itself is estimated from the first snapshot.
  - Each emitting cell's flux at each energy (all components, with its survival) is multiplied by the fraction of the window during which its ENAs were made before the first snapshot.
  - These are summed and divided by the window-averaged flux of all cells. The energies are combined with the band's count-rate weights.
  - The low-energy heliosheath PUIs therefore count for as much as they contribute at each energy.
  - With identical snapshots this share equals 1 − (time-dependent map / steady map), checked to 1e-15.
- **Regions with missing flux:** by default their flux is **set to 0** in every component map of that band and window (`include_incomplete_lines = false`). With `include_incomplete_lines = true` their incomplete flux is kept. Either way they are flagged: `coverage_complete = 0` in the HDF5 file, and a legacy flag file `Opher_multiion_withneutralsMissingFluxFlagKappa<k>EbandNN_<window>.dat` (1 = missing flux, same layout as the count-rate files), so the plots can box these cells as regions to be evaluated further (milestone 7).
- **Warnings:** each incomplete window is printed with a WARNING giving the number of lines flagged by each test, and is listed in `time_dependent_summary.txt`.

**Outputs** (in `output_dir` or `--out`):
- `flux_bandNN_<window>.h5`, e.g. `flux_band03_2020-2021.h5`: the window-averaged flux of every component, plus the coverage of that window:
  - `coverage_complete` (1 or 0 per line of sight);
  - `coverage_max_extinction_missing` and `coverage_missing_flux_fraction`;
  - `r_coverage_au` (at the lowest energy of the band), `incomplete_lines` (with `incomplete_lines_extinction` and `incomplete_lines_missing_flux`), the window and the input range.
- The legacy count-rate files, with `_<window>` before `.dat`.
- `time_dependent_summary.txt`: one line per band and window, with the coverage.

**Order of work and memory.**
- **Order:** the snapshots are read once, from the last to the first, so the loss terms a path needs (the emission snapshot and later ones) are always at hand. Each snapshot is read and its streamlines are processed as in a steady-state run, then the production of every energy is binned.
- **Memory:** the loss sums of all snapshots are kept, N × lines of sight × cells (survival_observer_r to r_upper) × 8 bytes per energy (about 0.5 GB per energy for 100 snapshots of a 30 × 60 grid out to 800 AU). If all energies do not fit in `memory_gb` (default 8), the energies are split into groups and the snapshots are read once per group.
- **Determinism:** every sum has a fixed order, so the maps are **bit-identical for any thread count and any energy grouping**.

**Checks.**
- **Steady-state limit:** `tools/td_repeat_snapshot.py` writes a configuration that repeats one snapshot N times. `tools/td_vs_steady.py` then compares the result with the steady-state maps:
  - the missing share of the total flux equals `coverage_missing_flux_fraction` on every line of sight (expect about 1e-15);
  - fully covered lines of sight equal the steady-state maps in every component (about 1e-15);
  - lines marked complete but not fully covered are listed with their missing shares. These are below the thresholds, so they are not errors.
- **Unit tests:** `test_m2` checks the configuration, the windows and the path sum against a cell-by-cell sum.

**Snapshot files.** List them as `[[snapshot]]` entries (`time`, if given, must be the start of the snapshot's time cell), or give name patterns in `[time_dependent.files]`: `{index}` is replaced by the snapshot number (from `first_index`), `{index:03}` pads it to 3 digits. HDF5 snapshots (`ena_convert`) work for every model; with them, use `file = "snapshot_{index:03}.h5"`.

## MPI (version 2.8)

Time-dependent runs can be spread over several MPI ranks, each with its own OpenMP threads:
```
mpiexec -n 8 enacpp/build/ena_run run_td.toml --threads 16
```

**How the work is split.**
1. **Production, over snapshots:** each rank reads a contiguous block of snapshots, processes their streamlines, and computes their loss sums and per-cell ENA production.
2. **Exchange, by line of sight:** after each round, each rank sends every line of sight's data to the rank that integrates that line of sight, using point-to-point messages.
3. **Integration, over lines of sight:** each rank integrates its block of lines of sight over all snapshots (path survival, arrival binning, coverage). Rank 0 collects the maps and writes the outputs.

**Same maps for any layout.** Every line of sight is summed in the same order on one rank, so the maps are **bit-identical for any number of ranks, threads and energy groups**, and to version 2.7. This was checked with 1–7 ranks, uneven blocks and energy groups, for the multi-ion and extended single-ion models, and from text and HDF5 inputs.

**Errors.** If a rank fails, for example on a bad input file, every rank stops together. Rank 0 reports the message with the snapshot and the rank, instead of the job hanging.

**Memory.** `memory_gb` is per rank. Each rank stores its own lines of sight for all snapshots: the loss sums and the per-cell production. The number of emitting cells is estimated from the first snapshot, which rank 0 reads once at the start. More ranks therefore also means less memory per rank.

**Choosing ranks and threads.** Reading and processing the streamlines of one snapshot are mostly serial, so several ranks with fewer threads each usually beat one rank with many threads. On 4 cores (full-size synthetic grid, 6 snapshots, 23 energies, HDF5 inputs):
- 1 rank × 4 threads: 40 s;
- 4 ranks × 1 thread: 32 s;
- the exchange took 0.6 s of it.

With 100 snapshots, use as many ranks per node as memory allows. `tools/run_athena.pbs` is the job-script template for a production run on Athena. `tools/run_athena_td_test.pbs` is the time-dependent test: one steady-state snapshot repeated 24 times, run from text and HDF5 inputs with 16 × 16, 4 × 64 and 1 × 256 ranks × threads. All four runs must give bit-identical files, and they are compared with the workstation (`compare_runs.py`) and with the steady-state maps (`td_vs_steady.py`).

**Steady-state test on Athena** (passed 2026-10-16: the extended single-ion maps of bands 3, 6, 9, 13 and 14 are bit-identical to the workstation, with 256 and with 16 threads). `tools/run_athena_steady.pbs` runs a steady-state configuration on one whole node, one thread per core, with the build without MPI (`cmake -S enacpp -B enacpp/build_nompi -DENA_WITH_MPI=OFF`). It then repeats the run with fewer threads, to check that the maps do not depend on the thread count. `tools/compare_runs.py REF OTHER` compares the maps of runs from two machines. It reports each map as bit-identical, agreeing within a relative tolerance (default 1e-9), or different, listing the worst cells. Values negligible in the band's total flux are measured against that total, so round-off in them does not count. Another compiler or maths library can change the last digits; on one machine, the maps are bit-identical for any thread count (`tools/same_maps.py`).

**Compilers.** Use GCC, or Clang or another compiler that keeps standard floating-point rules by default. The build turns off fused multiply-add (`-ffp-contract=off`); with Intel compilers it also sets `-fp-model=precise`, because their default fast model reorders sums.

**Steady-state runs** stay single-rank: under `mpiexec` only rank 0 works.

## Build

Requirements: a C++17 compiler (GCC 9 or newer, with OpenMP, as standard; tested with GCC 9, 10 and 13, which give byte-identical results), CMake 3.16 or newer, and the HDF5 C library with headers. HDF5 is optional, but without it there are no dumps, no `flux_bandNN.h5` output and no converter. MPI is optional: when CMake finds an MPI library (e.g. `module load` of the system MPI first), `ena_run` and `ena_convert` are built with it; `-DENA_WITH_MPI=OFF` builds without.

On Ubuntu:
```
sudo apt install g++ cmake libhdf5-dev
```

Build and test, from the `enacpp` directory:
```
cmake -S . -B build
cmake --build build -j
./build/test_m2 build/test_tmp      # expect "69 checks, 0 failures"
./build/test_physics PATH/TO/U_physics.h5   # optional: physics functions vs the IDL unit-test dump (agree to ~1e-7)
```

## Run

**Master configuration: `configs/all_settings.toml`.** It lists every setting, with its default and the models that use it. Copy it, edit the copy, and run it.
- **Model:** switch it with `model.populations`. Settings of the other models can stay in the file: they are skipped, with a `note:` in the log.
- **Time-dependent maps:** switch them on with `time_dependent.enabled = true`.
- **Defaults:** its active values reproduce `run_esi.toml`. The other files in `configs/` are shorter examples of the same format.

From the top directory of the ENA code, with `enacpp/` placed there:
```
enacpp/build/ena_run enacpp/configs/run_reference.toml
enacpp/build/ena_run enacpp/configs/run_reference.toml --bands 9,10,11,12,13 --out cpp_output_hi
```
The reference TS distance r_TS,V2 of the partition fits is computed from the model: the first TS crossing along the Voyager 2 trajectory in `configs/v2_trajectory_sim.txt` (r, theta, phi in the simulation frame of the BU model), set by `physics.v2_trajectory_file`. A fixed direction (`physics.v2_direction = [THETA, PHI]`) also works. Without either, `ena_run` uses the constant 78.5 AU and warns. For another model or frame, give its own trajectory file.

## Outputs of `ena_run`
Written to `output_dir` (or `--out DIR`):
- **`flux_bandNN.h5`**: the main output. Band-averaged ENA flux maps (cm^-2 s^-1 sr^-1 keV^-1), one dataset per component: `transmitted`, `reflected`, `accelerated`, `solar_wind`, `total`. Arrays are [theta][phi], with theta = i*dtheta (i = 0 is the north pole). The file also records the code version, the settings, the band limits and the geometric factor.
- **`Opher_multiion_withneutrals<TAG>CountRateKappa<k>EbandNN.dat`**: legacy count-rate text files (TAG = PUI1, PUI2, PUI3, SW, CombPUIs), one value per cell, for older plotting tools. Turn them off with `--no-legacy`.
- **`tec/denstemptest_EbandNN.dat`** (with `--tecplot` or `[run] tecplot_cells = true`; single-ion model, steady state): the per-cell values of each band at its centre energy E_c, for Tecplot. E_c is an extra `make_flux` evaluation, not one of the sub-band energies of the band sum, so every column belongs to E_c and the maps do not change. The 43 columns of the original code's file, in the same order and units: x, y, z (AU), region code, the indices i, j, k, θ, φ, r; density and temperature of the transmitted, reflected, accelerated and solar-wind populations; plasma density and temperature, v_r, v_t (km/s); f_test1, lc_test, the extinction; the streamline flow v_x, v_y, v_z (cm/s); the flux of each population; B_x, B_y, B_z (G); the neutral density; v_plas; the accelerated fraction, E_min and E_max of the accelerated ions (keV), kappa_ref and v_shock (cm/s); and a last column that is always 0, kept for the layout. The file has a Tecplot header (one ordered zone: I = r over the integration window, J = φ, K = θ, point format) and loads directly. Rows run from `r_lower` to `r_upper` (excluded), and the last θ row is repeated at θ = 180° so the grid closes at the south pole. About 420 MB per band; writing it takes several seconds per band.

## Maps: `tools/ena_plot.py`

Sky maps of the flux files in **ecliptic J2000** coordinates, in a Mollweide projection. It needs numpy, h5py and matplotlib.
```
python3 enacpp/tools/ena_plot.py OUT/flux_band13.h5                          # nose-centred map, total flux
python3 enacpp/tools/ena_plot.py OUT/flux_band13.h5 --view tail               # tail-centred (DOWNWIND marked)
python3 enacpp/tools/ena_plot.py OUT/flux_band13.h5 --view ecliptic           # centred on ecliptic (0, 0)
python3 enacpp/tools/ena_plot.py OUT/flux_band13.h5 --scale 1.8 --vmin 0 --vmax 12   # flux x 1.8, fixed colour range
python3 enacpp/tools/ena_plot.py OUT/flux_band*.h5 --config enacpp/configs/plot_settings.toml --out maps
python3 enacpp/tools/ena_plot.py run_b/flux_band13.h5 --ratio run_a/flux_band13.h5 --cmap RdBu_r
python3 enacpp/tools/ena_plot.py OUT/flux_band13.h5 --write-ecliptic          # all components on a regular 6-deg ecliptic grid
python3 enacpp/tools/ena_plot.py OUT/flux_band*_*.h5 --movie --frame-seconds 0.5   # time-dependent maps: one movie per band
python3 enacpp/tools/ena_plot.py OUT/flux_band13_*.h5 --movie --normalize      # each map divided by its own maximum
python3 enacpp/tools/ena_plot.py old/Opher_multiion_withneutralsCombPUIsCountRateKappa1.00Eband13.dat --idl --projection hammer
python3 enacpp/tools/test_ena_plot.py                                         # checks of the coordinate handling
```
The legacy count-rate files (`...<TAG>CountRateKappa<k>EbandNN.dat`) are read too, and turned back into flux (count rate / E_central / G); the tag (PUI1, PUI2, PUI3, SW, CombPUIs) is the component.

**The key options** (in `configs/plot_settings.toml` or on the command line; the command line wins):
- **Centre:** `view = "nose"` (upwind, 255.4°/+5.1°; UPWIND, V1 and V2 marked), `"tail"` (downwind, 75.4°/−5.1°; DOWNWIND marked) or `"ecliptic"` (0°, 0°; all four marked). `center = [lon, lat]` (or `--center LON LAT`) sets any other centre. The centre goes to the middle of the map, latitude included.
- **Scaling:** `scale = S` (`--scale S`) multiplies the flux; the title then reads "Energy = 4.29 keV (Flux x S)".
- **Normalized maps:** `normalize = true` (`--normalize`) divides each map by its own maximum, taken over the cells with complete coverage (cells with missing flux are left out). The colour range is 0 to 1 (`vmin`/`vmax` change it; `[band_limits]` are not used; with `log` 0.01 to 1), the title reads "(normalized)", the colour bar "ENA flux / map maximum", and the maximum is printed at the lower left. `scale` has no effect. File names get `_norm` (e.g. `flux_band13_2010-2011_total_norm_nose.png`). In movies each frame is normalized to its own maximum, so a movie shows the change of the pattern without the change of the overall flux.
- **Colour range:** `vmin`/`vmax` for every map; otherwise `[band_limits]` per band, otherwise each map's own min and max. Values outside the range take the end colours.
  - **Default band limits** (minima 0): band 9 = 265, 10 = 140, 11 = 60, 12 = 30, 13 = 12, 14 = 0.90, 15 = 0.08, 16 = 0.025, 17 = 0.009. Bands 0–8 have no default.
  - **A range from the data** (no limits set) is rounded outwards so all five colour-bar labels are exact, e.g. 62.8–186.3 becomes 62–190.
- **Colour-bar numbers:** `legend_style` (`--legend-style`):
  - `integer`: whole numbers, e.g. 190;
  - `float`: to the nearest hundredth, e.g. 0.85;
  - `exponent`: the mantissa to the nearest hundredth, e.g. 9.00e-03;
  - `auto` (default): integer for maxima ≥ 10, float for ≥ 0.1, exponent below.

  A warning is printed when the chosen style is too coarse for the range (e.g. integer on 0–0.9).
- **Scale type:** linear by default, `log = true` or `--log` for a log scale.
- **Movies** (`movie = true` or `--movie`; time-dependent maps `flux_bandNN_<window>.h5` only): one movie per band and component, the windows in time order. Each frame is the map as drawn without `--movie`, with its window in large type at the right of the title line (`time_label = false` or `--no-time-label`: without). Every frame of a band has the same colour range: `[band_limits]` or `vmin`/`vmax` where given, otherwise one range over all the band's frames. The frames are drawn twice, first to find the area all of them need, so that every frame has the same size and the map does not move. `frame_seconds` (`--frame-seconds`, default 1) sets the time per frame, `movie_dpi` (`--movie-dpi`, default 120) the resolution. `movie_format = "gif"` (default; written with Pillow, which matplotlib installs) or `"mp4"` (`--movie-format mp4`; needs the `ffmpeg` program; H.264 at 25 frames per second, each window repeated for its time). GIF frames all use one fixed palette (192 colours evenly along the colour map and 64 greys) rather than 256 colours chosen per frame, so the colour bar is continuous and identical in every frame and a value has the same colour throughout. Output: `flux_bandNN_<component>_<view>_movie.gif` in `out` or next to the first file.

Output names: `<input stem>_<component>[_ratio][_idl]_<nose|tail|ecliptic|custom>.<format>`.

**Simulation frame.** z is along the solar rotation axis; x is 5° above the interstellar flow, in the solar equatorial plane; y completes the right-handed system.

**The rotation to ecliptic J2000** is fixed by two directions known in both frames:
- **The solar rotation axis:** inclination 7.25°, ascending node 75.76° (J2000). It is kept exactly.
- **The interstellar flow:** (26.3, 0.3, −2.3) km/s in the simulation frame; ecliptic longitude 75.4°, latitude −5.2°.
- **Agreement:** the angle between the two directions is 95.00° in one frame and 95.11° in the other, so they agree to 0.1°. The nose lands at (255.4°, +5.2°).
- **Options:** `--flow-sim`, `--flow-ecliptic` and `--solar-axis` change them.

**Compared with a plain longitude shift** of the simulation-frame map (the placement of the original IDL plots):
- the flanks move by up to ±7.2° in latitude;
- longitudes move by about 4°;
- the solar poles sit 7.3° from the ecliptic poles.

`--idl` reproduces that placement (with `--projection hammer` for the look of the original IDL plots), for comparison with existing figures.

**Cells.** The grid values are cell centred (θ = i·dθ ± dθ/2, φ = j·dφ ± dφ/2). Each cell is rotated and drawn as its own polygon, with no resampling. A 3° cap around the simulation south pole has no cell and stays blank.

**Map.**
- **Projection:** Mollweide (default) or `projection = "hammer"`.
- **Orientation:** longitude increases to the left (east to the left, as the sky seen from inside). Grid every 30°, longitude labels −180..180, plus edge meridians 4° from the map edge.
- **Missing flux:** regions with missing flux (time-dependent windows) are drawn grey with a black box per cell; the caption gives their number. `--no-boxes` hides the boxes.
- **Markers:** set by view (`markers_nose`, `markers_tail`, `markers_ecliptic`). UPWIND and DOWNWIND come from the frame; V1 and V2 from `[markers]`. `--no-markers` hides them.
- **Other settings:** `cmap` (default `jet`), `legend_style`, `common_range`, `limits` (part of the sky), `fontsize`, `figsize`, `format` (png, pdf, eps, svg), `dpi`, `out`.

`tools/plot_flux.py` (the simulation-frame quick look) is superseded by `ena_plot.py`.

## Preparing input series with Tecplot (`tools/tecplot_prep.py`)

For time-dependent runs, every MHD solution of a series has to become the four input files. `tools/tecplot_prep.py` runs the Tecplot macros of that workflow on each snapshot, either with PyTecplot (`pip install --user pytecplot`; run through `tec360-env`) or, with `--batch`, in Tecplot's own batch mode (`tec360 -b -p`; no PyTecplot). The driver needs Python 3.6 or newer and no other packages.

- **Configuration:** `configs/prep_single_ion.toml`; the macros are in `configs/macros/single_ion/`. The steps of each snapshot:
  - 1: fix the neutral solution: the program `ena_fixneutrals` (see below), which writes the fixed neutrals as a binary `.plt`;
  - 1b: Tecplot loads that file and saves it as `.szplt` (`1b_save_szplt.mcr`);
  - 2a-2d: interpolate the neutrals (2a, 2b) and the MHD solution (2c, 2d) onto the spherical grid; the four parts of the recorded macro are timed separately;
  - 3: calculate and write the neutral, plasma and boundaries files, and the streamline field: the spherical grid with the 19 streamline variables (`3_calculate_files_trace.mcr`);
  - 4: trace the streamlines through that field: the program `ena_trace` (see below), from the seeds of the recorded macro with each pole once (`configs/streamline_seeds_2AU.txt`, 1742 seeds; the macro has 30 copies of each pole).

  About 5 minutes per snapshot on a 24-core workstation with Tecplot 2025 R2, most of it the two interpolations. `configs/prep_single_ion_tecplot.toml` does steps 1 and 4 in Tecplot instead (`1_fix_neutrals.mcr`, and Tecplot's streamtraces in `3_calculate_files.mcr`): about 28 minutes per snapshot, and streamlines that are not reproducible.
- **Inputs** per snapshot: a numbered pattern, a sorted glob (e.g. `3d__var_2_t*.szplt`), a list file, or one fixed file (the spherical grid). A step can first load an input or an intermediate file into a new layout, as is done by hand, or run a program (`command = [...]`, with `{VAR}` for the file of a macro variable).
- **Files are passed as macro variables.** Every input, intermediate and output file reaches the macros as a variable (`|MHDFILE|`, `|NEUTRALSRC|`, `|NEUTRALFIXED|`, `|NEUTRALFILE|`, `|PLASMAFILE|`, `|BOUNDARYFILE|`, `|STREAMFILE|`), as an absolute path built from the starting directory, so the macros work anywhere. Before running, the driver scans the macros and stops if an absolute path is left in one.
- **Checks:** outputs are written under temporary names, checked, and only then renamed:
  - grid files: one row per cell, the expected number of columns, no NaN/Inf, the Tecplot header removed and, with `theta_index_column`/`theta_index_max`, the rows with theta index i > 29 (theta > 180 deg, where the macro grid repeats the sky) dropped;
  - the streamline file: the variable count and at least one zone.
- **Restarting:** finished snapshots are skipped, so a stopped run can be started again. `prep_index.txt` records the input files of each snapshot.
- **Without Tecplot:**
  - `--dry-run` shows what would run and checks the macros;
  - `--check-only` checks existing files.
- **Batch mode (`--batch`, `--tecplot CMD`):** the macro steps between two program steps are written into one macro in the work folder (the `$!VARSET` lines, then per step the commands that load its input and the step's macro) and run with `tec360 -b -p` (or CMD) in the work folder: one Tecplot start per snapshot, between `ena_fixneutrals` and `ena_trace`. Tecplot's messages go to `batch.log` there; a run fails on an error exit status or an error line in `batch.log`, which the driver then shows. The files are the same as with PyTecplot.
- **Folders (`--work-dir DIR`, `--stage-dir DIR`):** `--work-dir` overrides `[work] dir` (the intermediate files). With `--stage-dir` (or `[work] stage_dir`), each snapshot's inputs are copied to DIR first and Tecplot reads the copies; a fixed input (the grid) is copied once, the other copies are deleted after the snapshot. For a compute node's local disk: Tecplot reads `.szplt` files piece by piece, which is very slow from a network file system. Environment variables in both paths are expanded.
- **On NASA HECC, Aitken Rome nodes:** `tools/run_aitken_prep.pbs` (submitted from a Pleiades front end): one Rome node, `--batch`, inputs and intermediate files on the node's disk (`$TMPDIR`), the finished files to the configured output folder, Tecplot and the C++ steps on cores 0-23 (`TECPLOT_CORES`), one snapshot after the other (one Tecplot license), live log `ena_prep_<job>.log`. Measured 2026-10-26: 523 s per snapshot (input copies 23 s, `ena_fixneutrals` 14 s, Tecplot 423 s, `ena_trace` 17 s, checks about 45 s). Timing tests (`tools/run_aitken_tecplot_test.pbs`): reading the inputs from `/nobackup`, the MHD interpolation took 2711 s; from the node's disk 201 s on the whole node and 112 s on 24 cores (Tecplot starts over 1000 threads on the whole node). Build `ena_fixneutrals` and `ena_trace` on a front end: `module load gcc`, then `cmake -S enacpp -B enacpp/build_pfe -DENA_WITH_MPI=OFF -DENA_WITH_HDF5=OFF` and `cmake --build enacpp/build_pfe -j 8 --target ena_fixneutrals ena_trace`.
- **On Athena:** `tools/run_athena_prep.pbs` (batch mode), optionally with several Tecplot workers side by side, each with its own work folder. A node is charged whole, and each worker uses a Tecplot license. Status 2026-10-26: Tecplot 2025 R2 aborts at startup on the compute nodes (FIPS self-test failure; NAS: needs OpenSSL compatibility packages, not in the short term), and Tecplot 2023 R2 did not finish the MHD interpolation within 20 minutes; use the Aitken Rome nodes. `tools/run_athena_tecplot_test.pbs` times the slow steps of a Tecplot version on a compute node.
- **Tecplot versions:** Tecplot 2025 R2 interpolates quickly from the neutrals only if they are a `.szplt` file (from a binary `.plt` or ASCII `.dat` it takes hours, as does 2023 R2 from `.plt`; 2020 R1 is fast with all three). Step 1 therefore saves the fixed neutrals as `.szplt`, and step 2a reads them so.
- **Step 1 in C++:** `ena_fixneutrals` replaces Tecplot's step 1, which took about 25 minutes per snapshot, mostly to read the 18 GB neutral text file; it takes about 30 s. It reads the file in parallel through a memory map, evaluates the equations of `1_fix_neutrals.mcr` as Tecplot does (values in single precision, each equation in double precision in the written order), and writes the 33 variables Tecplot keeps. Validated 2026-10-22: the neutral input file came out byte-identical to the one of the all-Tecplot route.
- **Reproducibility (tested 2026-10-22):** Tecplot's streamtraces are not reproducible. Two runs of the manual procedure on the same snapshot gave 922 of 1582 different streamlines and 85 vs 51 streamlines stuck at the south pole; even a second trace in the same Tecplot session, on the same field, changed the streamlines. The ENA maps of such runs differ by up to 5-19% in single directions (more than 1% in 6-40% of the sky, depending on the band). Tecplot's interpolation also differs between runs, in the 7th digit; how much that still changes traced streamlines where the flow diverges is measured by tracing the fields of two runs. Hence step 4: `ena_trace` gives the same streamlines for the same field every time, on any machine and for any number of threads (two full runs: identical streamlines, bit-identical maps). `tools/compare_inputs.py` compares input files from two runs (streamline files streamline by streamline, with the seeds that never leave their seed).
- **`ena_trace` vs Tecplot's streamtraces (2010-02-06 snapshot, investigated 2026-10-23):** with identical neutral, plasma and boundaries files, the maps from the two streamline files differ in 800-1400 of the 1800 directions by more than 1% (up to 25%). Tested and found negligible: the seeds Tecplot drops, its stuck streamlines, the point spacing, and the variables written along the streamlines (equal to 1e-7 at the seeds; the paths agree to 0.01 deg out to 50 AU). The two causes:
  - **The north pole:** in the macro grid the cells around the north pole collapse to a line, and Tecplot's streamtraces stop there. Every streamline heading over the pole ends at theta = 6 deg (98-144 AU, after 107-172 points), and the seeds at theta 0 and 6 deg give none, so the polar heliosheath has no streamline points and is filled from its neighbours (1270 polar-cap cells have points only with `ena_trace`). `ena_trace` follows the flow over the pole (e.g. seeds at phi 198 pass within 0.3 deg of the pole at 110-150 AU and come down at phi 18). This gives the north-pole differences of up to 25%; the `ena_trace` maps are the faithful ones there.
  - **v_shock (`physics.vshock = "first_below"`, the default):** the TS is the entry into a region 3/6 cell, and v_shock the flow speed at the first point from there with M^2 < 0.64. That point is usually the crossing point itself (median path 0 AU), inside the MHD shock, so v_shock depends on where the points fall in it; per streamline it changes by tens of km/s with the sampling, the distributions do not (median 204-221 km/s for Tecplot and for `ena_trace` at three step sizes). Two `ena_trace` samplings of the same field (cell fraction 0.5 vs 0.89) differ by more than 1% in 420-1740 directions (max 9-16%); with `vshock = "downstream"` in 600-1300 (max 6-9%). The spread of two Tecplot runs (107-719 directions) has the same origin. The maps are therefore uncertain at this level whatever traces the streamlines; keep the tracer settings fixed within a series.
  - Diagnostic tools used (all read the files directly): `tools/seed_coverage.py` (which seeds a streamline file covers; writes the matched or de-duplicated seed list), `tools/compare_streamlines.py` (two streamline files seed by seed: stuck streamlines, point spacing, path angles, the variables at the seed and at 20/50/100 AU), `tools/compare_dumps.py` (the per-cell fields of two `ena_run --dump` runs, optionally along one line of sight), `tools/cell_visitors.py` (which streamlines put points into given cells, binned as the ENA code does), `tools/seed_paths.py` (where the streamlines of one meridian of seeds go, e.g. over a pole).

## Programs

**`ena_trace FIELD.dat STREAMLINES.dat --seeds SEEDS.txt [--cell-fraction 1] [--max-steps 750] [--substeps 4] [--threads N]`**
- The streamlines of the input preparation, traced reproducibly. FIELD.dat is the spherical grid with the streamline variables, written by Tecplot in step 3 (an ordered zone in point format, theta 0-354 deg so that it includes both poles). STREAMLINES.dat is written as Tecplot writes extracted streamtraces, so the ENA code reads it unchanged.
- Method: the field is interpolated trilinearly in the grid's (r, theta, phi) index space (phi periodic, theta 0-180 deg, so the poles are ordinary grid rows). Each streamline follows the flow direction U/|U| from its seed with the classical 4th-order Runge-Kutta method; one step is cell_fraction x min(dr, r dtheta), as Tecplot's CellFraction (at most 2 AU; Tecplot's points are about as far apart beyond 20 AU), split into `substeps` Runge-Kutta steps, with a point written after every step. A streamline ends after max_steps points (Tecplot's MaxSteps), on leaving the grid, or where |U| < 1e-6 km/s. The variables are interpolated to the points; X, Y, Z are the point itself and theta0, phi0, r0 are computed from it with the step-3 macro's formulas.
- Poles: theta and phi are computed from the position, so streamlines cross the poles and continue on the far side; Tecplot's streamtraces stop at the collapsed cells around the north pole of the macro grid (see "`ena_trace` vs Tecplot's streamtraces" above).
- Tested on analytic flows: radial outflow and uniform flow (paths exact to 1e-8 AU), flow across both poles (no deviation, none stuck), rotation across the phi = 0/360 seam (2e-4 AU after 750 steps at 300 AU). The same field always gives byte-identical output, for any thread count. A full-size field (60 x 60 x 736) takes about 11 s on 4 cores.

**`ena_fieldtopo run.toml [--snapshot N | --snapshots all|A-B|A,B,...] [--out topology.dat] [--z-top 500] [--z-bottom -500] [--ism-regions 1] [--r-sun 50] [--max-length 10000] [--cell-fraction 1] [--substeps 2] [--only-regions 3,6] [--all-radii] [--rewrite-boundaries] [--suffix _topo] [--overwrite] [--edge-ism all|nose|none] [--threads N]`**
- Magnetic field-line topology on the ENA grid, to locate the heliopause where the region codes cannot (e.g. the low-latitude heliotail, where reconnected field lines mark the mixing region). Reads the region codes and B of snapshot N (the boundaries file, or the HDF5 snapshot) of the settings file's series.
- The field line through every grid point is followed in both directions along B/|B| (4th-order Runge-Kutta, B interpolated as in `ena_trace`). A direction ends at the Sun (r <= r_sun), in the ISM (z >= z_top, z <= z_bottom, a cell of one of the `--ism-regions`, or the outer grid edge; `--ism-regions none` turns the region test off), or neither within `max_length` of path (long) or in a field null. `--edge-ism nose` counts leaving through the outer grid edge as reaching the ISM only where the line leaves at x < 0 (the nose side): down the tail a line can leave the grid while still in the heliotail, so there it ends as "grid edge (not ISM)" (end code 7) and counts like a long line; `--edge-ism none`: no edge exit counts (default `all`).
- Topology of each point: 1 closed (no end in the ISM: solar field), 2 open (Sun and ISM: reconnected, the mixing region), 3 interstellar (both ends in the ISM), 4 ISM-connected (one end in the ISM, the other long or null), 0 undetermined, -1 not traced (`--only-regions`). A summary by topology is printed.
- Output: a Tecplot ASCII file on the radii of the integration window (`--all-radii`: all), with a header, one ordered zone and the last theta row repeated at theta = 180 deg: X, Y, Z (AU, simulation frame), region, topology, the end and path length of each direction, the path length to the ISM, B_x, B_y, B_z, |B| and |B_phi|/|B| (the azimuthal share about the z axis; near 1 in the Parker spiral). Results do not depend on the thread count.
- Tested on analytic fields: rings about the z axis (closed), a radial field (open), a vertical field (interstellar), and rings inside a vertical field (closed inside, interstellar outside, ISM-connected where vertical lines run into the rings).
- Run time: lines that stay in the heliosphere run to `max_length`; on the full grid that is minutes on one Athena or Rome node.
- **Heliopause from the topology (`--rewrite-boundaries`):** traces the heliosheath cells (regions 3 and 6) at every radius and writes a copy of the boundaries file in which each of them whose field line is not closed (topology 2, 3 or 4) is region 2. The heliopause then follows the closed solar field, e.g. in the low-latitude heliotail, where lines of sight now stop at the open field. The termination shock (the 4 -> 3 entry) and every other cell are unchanged, and so is every other character of the file (a region code written as `3.000000E+00` becomes `2.000000E+00`).
  - The copy goes next to the original, with `--suffix` (default `_topo`) before the extension: `time_series/boundaries_040.txt` -> `time_series/boundaries_040_topo.txt`. Point the run at the copies, e.g. `boundaries = "time_series/boundaries_{index:03}_topo.txt"` in `[time_dependent.files]`.
  - `--snapshots all` (or `A-B`, `A,B,C`) runs over the series in one go; a copy that exists is skipped (`--overwrite` writes it again), so a stopped job can be resubmitted. The Tecplot file is written only with `--out` (with several snapshots as `OUT_NNN.dat`).
  - Use `--ism-regions none`, so that the topology does not depend on the region codes (the ISM is then |z| >= 500 AU and the grid edge), and `--edge-ism nose`. On the 2010-02-06 snapshot (grid edge 1472 AU), 1467 of the 11507 cells changed inside 800 AU reached the ISM only through the grid edge, 1411 of them down the tail with the other end long: with `--edge-ism nose` these stay region 3. ASCII boundaries files of the standard model source only (for HDF5 runs, convert the rewritten files with `ena_convert`).
  - Job script: `tools/run_athena_fieldtopo.pbs`. Tested on the analytic fields: only the region column of the region-3/6 rows that are not closed changes.

**`ena_fieldtopo_mhd MHD.plt OUT.plt [--seed-box XMIN XMAX YMIN YMAX ZMIN ZMAX] [--seed-stride N] [--points SEEDS.plt] [--no-mhd-vars] [--no-propagate] [--boundary-passes 4] [--z-top 500] [--z-bottom -500] [--hplim-ism ZMIN] [--r-sun 50] [--max-length 10000] [--cell-fraction 0.5] [--substeps 2] [--threads N]`**
- The topology of `ena_fieldtopo`, on the original BATS-R-US grid: reads the MHD solution (`.plt`, one FEBRICK zone: the cell centres and the bricks joining them; the order of the nodes does not matter) and writes `OUT.plt` on the same grid (the same nodes and bricks) with all variables of the MHD file (`--no-mhd-vars`: only X, Y, Z), then the topology, the end and path length of each direction and the path length to the ISM. `--points` also writes the seeds alone as scattered points, with the same variables (a small file).
- Each point is located in the brick that holds it (the bricks are indexed by size class and position; each thread first tries its previous brick). B is interpolated trilinearly in a brick with 8 distinct corners and by inverse-distance weights in a collapsed brick (where refinement levels meet). A point in no brick is outside the domain.
- No region codes on this grid: the ISM is z >= z_top, z <= z_bottom, the domain edge or, with `--hplim-ism ZMIN`, Hplim < 0 where |z| >= ZMIN (the heliopause variable above that height). The step is `cell_fraction` of the size of the local brick.
- Which nodes: all nodes in `--seed-box` (default: the whole domain). A field line has one topology along its length, so each traced line also labels the nodes it passes; tracing goes in rounds (every 4096-th node, then every 512-th, 64-th and 8-th of those still unlabelled, then the rest), so only part of the nodes need a line of their own. A node passed by several lines of a round takes the lowest seed's line: results do not depend on the thread count. Then boundary passes: every node labelled by a passing line that shares a brick with a node of another class is traced on its own, then those around the nodes that changed, until none change (`--boundary-passes N`, default 4). On analytic fields the result then agrees with tracing every node on its own at every node. `--no-propagate` traces each seed on its own and labels only it (every N-th node of the box with `--seed-stride N`). The output's `label source` says how each node was labelled (1 traced itself in a round, 2 by a passing line, 3 traced itself in a boundary pass, 0 not labelled).
- Memory: about 40 GB for a 224-million-cell solution with the MHD variables in the output (about 30 GB with `--no-mhd-vars`). Tested on BATS-R-US-like files with analytic fields (rings, radial, vertical, rings inside a vertical field: the results of `ena_fieldtopo`), with the nodes in block order and shuffled (identical results), and on a two-level grid with collapsed bricks at the interface (field lines cross it).
- **ENA-grid mode, heliopause from the native field:** `ena_fieldtopo_mhd --ena-config run.toml --mhd-list LIST --snapshots all --rewrite-boundaries [--suffix _topomhd] [--subcell 1] [--edge-ism nose] [--out topo.dat] [--overwrite]`. Classifies the heliosheath cells (regions 3 and 6, every radius) of the run's ENA grid by tracing through each snapshot's native MHD field, and writes the boundaries copies as `ena_fieldtopo --rewrite-boundaries` does (cells not closed -> region 2; copies named `boundaries_NNN_topomhd.txt`). This removes the ENA grid's resolution from the topology: on the 2010-02-06 snapshot the ENA-grid version opened the low-latitude tail about 100 AU too early, with lines reaching the ISM only after 4000-5000 AU of winding (leakage through the smoothed current sheet).
  - `--mhd-list LIST`: one MHD `.plt` per line in snapshot order (line 1 = snapshot 0; relative paths relative to the list's folder); `--mhd FILE` for one snapshot.
  - `--subcell n`: n x n seeds across each cell in theta and phi (default 1: the cell centre, the point every value of the boundaries file refers to); the cell becomes region 2 if more than half its seeds are not closed. The Tecplot output (`--out`) gives the closed fraction of each cell and the region after the rewrite (`region (topology)`).
  - The run configuration must read ASCII files. Memory about 20-30 GB (coordinates, B and bricks only); time per snapshot on Athena (2010-02-06, 120800 cells): 178 s with the centre seeds (133 s of it indexing the MHD file), 411 s with `--subcell 3`. Job script: `tools/run_athena_fieldtopo_mhd.pbs`.
  - `--allow-truncated` (both modes): an MHD file that ends early (an incomplete copy; the connectivity comes last, so the last bricks are missing) is used with the bricks present, and a line that reaches the gap ends there as at the domain edge. Without it, such a file fails with a message giving how many bricks are missing. Tested: a copy 734116 bytes short (22942 of 2048383 bricks) gave the same regions as the complete file.
  - A snapshot that fails (e.g. an unreadable MHD file) is reported as FAILED and the others go on; the failed ones are listed at the end (exit status 1). Reads are done in pieces of 64 MB, and a read error of the file system is retried up to 5 times.
  - `--edge-ism all|nose|none` (also in the full-grid mode): as in `ena_fieldtopo`.
  - Tested on an analytic field written both as a native file and as a boundaries file: the centre seeds agree with `ena_fieldtopo` on 99.95% of the cells (120 of 441000 differ, all at the edge of the closed region); with `--subcell 3`, 60 cells change against the centre seeds, of 3300 cells whose seeds are mixed.

**`tools/plt_header.py FILE.plt`**: prints the layout of a Tecplot binary file (variables, zones, element type, sizes, data formats, sharing) from its header, and checks the file size; seconds even for files of many GB.

**`ena_fixneutrals NEUTRAL.dat FIXED.plt [--threads N]`**
- Step 1 of the single-ion input preparation (`configs/macros/single_ion/1_fix_neutrals.mcr`) without Tecplot: reads the neutral solution (Tecplot ASCII, one FEPOINT zone of bricks, 53 variables), computes the total and per-population neutral moments and temperatures, and writes the 33 variables of Tecplot's fixed-neutral file as a Tecplot binary file (version 112). Run by `tools/tecplot_prep.py` with `configs/prep_single_ion.toml` (step 1).
- Results are the same for any thread count.

**`ena_run run.toml [--snapshot N] [--out DIR] [--dump DIR] [--bands 9,13] [--set table.key=value] [--no-legacy] [--tecplot] [--threads N]`**
- Full run for one snapshot: inputs, streamline processing, flux for every energy of each band, outputs.
- With a `[time_dependent]` section: time-dependent maps from all snapshots (see above). `--snapshot`, `--dump` and `--stream-fields` are for steady-state runs only.
- `--bands` overrides `[bands] run` in the configuration.
- `--set table.key=value` overrides any configuration value (repeatable), e.g. `--set physics.fitted_esw=false`.
- `--dump DIR` writes the checkpoints C0-C5 (about 15 MB per energy for C4).
- `--tecplot` writes the Tecplot file of the per-cell values at each band's centre energy, `OUT/tec/denstemptest_EbandNN.dat` (see Outputs).
- Prints the settings in effect and checks of the inputs:
  - r_TS,V2, and how many streamlines have a TS;
  - region-5 (undetermined) cells, lines of sight that leave and re-enter the heliosheath, and lines of sight whose heliosheath ends at a region-5 cell;
  - how many region 3/6 cells the gap fill left empty.
- `--stream-fields FILE` (testing): skips the streamline processing and takes the cell fields from a C3 file.
- **Heliosheath cells without streamline data:** the gap fill fills each empty region 3/6 cell from its neighbours, but cannot fill a part of the heliosheath that no streamline reaches. The run lists such cells after the `Gap fill` line. With `[physics] unfilled_cells = "error"` (default) a line of sight that reaches one stops the run; with `"skip"` they give no ENAs, as the supersonic pockets, and the integration continues behind them.

**`ena_load run.toml [--snapshot N] [--dump DIR]`**
- Reads one snapshot and prints value ranges, a region histogram and the streamline count.
- With `--dump`, writes C0_config, C1a_limits and C1b_inputs.

**`ena_convert run.toml OUTDIR [--first N] [--last M] [--skip-existing] [--compress L]`**
- A one-time step before the runs: converts the snapshots to `OUTDIR/snapshot_NNN.h5` (fields in internal units, all ion fluids and neutral populations of every model, plus the streamlines).
- `--first/--last` convert part of the snapshots (several jobs or a job array); `--skip-existing` resumes an interrupted conversion (files are written under a temporary name and renamed when complete); under `mpiexec -n R` the snapshots are split over the R ranks.
- `--compress L`: gzip level 0-9, default 0 (uncompressed). Uncompressed files load fastest; on the full-size synthetic grid a snapshot with its streamlines loads in 1.8 s against 4.9 s from text, while compressed files take about as long as text.
- HDF5 file locking is switched off in all programs (`HDF5_USE_FILE_LOCKING=FALSE`, unless you set it yourself): on NFS a lock makes the client re-read the whole file over the network at every open (measured: 4.3 s instead of 0.3 s per snapshot).
- Prints the configuration lines for `inputs.format = "hdf5"`. Older single-ion files (format 1) are still read; a file made for another model is refused.

## Input checks

**Grid files (plasma, neutral, boundaries):**
- Every row is placed by its i, j, k columns, in any order.
- theta/phi/r must lie within half a cell of the grid point (phi is not checked at the pole, where files write 0).
- Duplicate cells and missing cells are errors. Set `inputs.allow_incomplete = true` only for tests on file slices.
- Region codes are truncated to integers. Moscow codes are remapped: 3, 4 -> 2; 1 -> 4; 2 -> 3; anything else -> 5 (undetermined).

**Streamline file:**
- Read zone by zone (one Tecplot ZONE per streamline).

## Units inside the code
- Densities: cm^-3.
- Temperatures: K. Standard-model temperatures are halved on reading; Moscow-model temperatures are not.
- Velocities: cm/s. Each velocity column's unit is declared in the configuration (`km/s` by default).
- Magnetic field: as read from the file.

## Layout
```
include/ena/   headers: config, grid, snapshot, ascii, streamlines, regions, physics,
               stream_processing, flux, checkpoints, h5io
src/           implementations
apps/          ena_run, ena_load, ena_convert, ena_fixneutrals, ena_trace
tools/         ena_plot.py (sky maps), same_maps.py and compare_runs.py (comparing runs), tecplot_prep.py (input
               preparation), comparison tools for inputs and streamlines, run_athena*.pbs and run_aitken*.pbs (NASA jobs)
tests/         test_m2 (synthetic files in the real formats), test_physics (vs U_physics.h5)
configs/       run_reference.toml (reference single-ion inputs), slice_test.toml (file slices)
```
