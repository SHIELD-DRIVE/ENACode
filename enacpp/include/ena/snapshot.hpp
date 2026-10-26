#pragma once
// Data for one time snapshot of the MHD solution, on the common grid.
// Units inside the code: densities cm^-3, temperatures K, velocities cm/s,
// magnetic field as read (Gauss).

#include <cstdint>
#include <string>
#include <vector>

#include "ena/config.hpp"
#include "ena/grid.hpp"

namespace ena {

struct PlasmaFluid {         // one ion fluid (single-ion: one entry; multi-ion: SW, pu3, pu2)
  Field3D<double> n, T, vr, vtheta, vphi, vt;
};

struct NeutralPopulation {   // one neutral population (single-ion: total H only;
                             // multi-ion: total, then populations 1-4)
  Field3D<double> n, T;
  Field3D<double> vr, vt, vx, vy, vz;  // multi-ion only, as read (km/s)
};

struct FileCoordinates {     // theta/phi/r as written in the plasma file (diagnostics, C1b)
  Field3D<double> theta, phi, r;
};

struct Snapshot {
  double time = 0.0;
  std::vector<PlasmaFluid> plasma;
  std::vector<NeutralPopulation> neutrals;
  Field3D<std::int16_t> region;        // after any model remapping
  Field3D<double> Bx, By, Bz;
  Field3D<double> pbeta;               // multi-ion plasma file, last column (diagnostic only)
  FileCoordinates coords;              // filled only when requested
};

struct ReadOptions {
  bool keep_file_coordinates = false;
  bool verbose = true;
};

// Reads the plasma, neutral and boundaries files of snapshot `index`.
Snapshot read_snapshot(const RunConfig& cfg, std::size_t index, const ReadOptions& opt = {});

// Individual readers (exposed for testing).
void read_plasma_ascii(const RunConfig& cfg, const std::string& path, Snapshot& s, const ReadOptions& opt);
void read_neutral_ascii(const RunConfig& cfg, const std::string& path, Snapshot& s, const ReadOptions& opt);
// Multi-ion files: plasma with the SW, pu3 and pu2 fluids (temperatures as
// read, not halved), neutrals with total H and populations 1-4.
void read_plasma_mi_ascii(const RunConfig& cfg, const std::string& path, Snapshot& s, const ReadOptions& opt);
void read_neutral_mi_ascii(const RunConfig& cfg, const std::string& path, Snapshot& s, const ReadOptions& opt);
void read_boundaries_ascii(const RunConfig& cfg, const std::string& path, Snapshot& s, const ReadOptions& opt);

}  // namespace ena
