# ENA code

All-sky maps of heliospheric **energetic neutral atom (ENA) flux** from 3D MHD simulations of the heliosphere, as seen from 1 AU in the energy bands of IBEX-Lo, IBEX-Hi and INCA/Ultra. For every direction on a 6° × 6° sky grid, the code sums the ENAs made by charge exchange in the heliosheath along the line of sight, and their losses on the way in.

- **Three models:** single-ion with four ion populations from the termination shock (Kornbleuth et al. 2023), an extended single-ion model for low-energies (Galli et al. 2023), and a multi-ion model (Kornbleuth et al. 2026).
- **Steady-state or time-dependent maps:** one MHD snapshot, or a series of snapshots giving maps per time window, with the ENA travel times taken into account.
- **Fast and reproducible:** C++ with OpenMP and MPI; the maps are bit-identical for any number of cores and processes.
- **Plotting:** sky maps in ecliptic J2000 coordinates (Python).
- **Input preparation:** tools that turn a series of MHD solutions into the input files, with Tecplot 360.

## Requirements

- C++17 compiler with OpenMP (e.g. GCC 9 or newer) and CMake 3.16 or newer
- HDF5 (C library)
- MPI, optional (time-dependent runs on several processes)
- Python 3 with numpy, h5py and matplotlib, for the maps and comparison tools
- Tecplot 360, only for preparing input files from MHD solutions (through PyTecplot, or in Tecplot's own batch mode)

## Quick start

Get the code and build it; all commands are run from the top folder, `ENACode`:
```
git clone https://github.com/SHIELD-DRIVE/ENACode.git
cd ENACode
cmake -S enacpp -B enacpp/build
cmake --build enacpp/build -j
enacpp/build/test_m2 enacpp/build/test_tmp
```
The test should end with `69 checks, 0 failures`.

**Sample run.** The input files of a single-ion MHD solution (Kornbleuth et al. 2023; about 1 GB) are on Zenodo, [doi:10.5281/zenodo.23675226](https://doi.org/10.5281/zenodo.23675226). Download them into `secondary_files/input/`, replacing the placeholder files there:
```
curl -L -o secondary_files/input/plasma.dat "https://zenodo.org/records/23675227/files/plasma.dat?download=1"
curl -L -o secondary_files/input/neutrals.dat "https://zenodo.org/records/23675227/files/neutrals.dat?download=1"
curl -L -o secondary_files/input/boundaries.txt "https://zenodo.org/records/23675227/files/boundaries.txt?download=1"
curl -L -o secondary_files/input/streamlines.dat "https://zenodo.org/records/23675227/files/streamlines.dat?download=1"
```
Then make the IBEX-Hi band 13 map (a few seconds) and plot it:
```
enacpp/build/ena_run enacpp/configs/all_settings.toml --bands 13
python3 enacpp/tools/ena_plot.py secondary_files/output/flux_band13.h5
```
The map is written to `secondary_files/output/flux_band13_total_nose.png`.

**Your own runs:** copy the master settings file (`cp enacpp/configs/all_settings.toml enacpp/configs/my_run.toml`), set your input files and options in the copy, and run `enacpp/build/ena_run enacpp/configs/my_run.toml`.

## Documentation

- [User guide](enacpp/docs/ENA_USER_GUIDE.md): building on a workstation or on NASA's Athena, settings, steady-state and time-dependent runs, input preparation (also on NASA's Aitken nodes), sky maps, common problems.
- [Full reference](enacpp/docs/REFERENCE.md): every program, option and input check.
- [`configs/all_settings.toml`](enacpp/configs/all_settings.toml): every setting with its default and explanation.
- [Changes from the IDL code](enacpp/CHANGES_FROM_IDL.md): the code began as a port of an IDL code; this lists every difference.

## Contact

Marc Kornbleuth, Boston University: [kmarc@bu.edu](mailto:kmarc@bu.edu)

## References

- Kornbleuth, M., et al., ApJL, 944, L47, https://doi.org/10.3847/2041-8213/acb9e0
- Galli, A. et al., ApJL, 954, L24, https://doi.org/10.3847/2041-8213/aced9b
- Kornbleuth, M., et al., ApJ, 1004, 18, https://doi.org/10.3847/1538-4357/ae67fc
- Sample input files: Kornbleuth, M., Sample Single-ion BU MHD Input Files for SHIELD ENA Model, Zenodo, https://doi.org/10.5281/zenodo.23675226
