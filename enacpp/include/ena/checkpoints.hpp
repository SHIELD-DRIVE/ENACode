#pragma once
// Checkpoint dumps C0-C5 of the intermediate arrays (HDF5). Their layout, names
// and units (velocities in km/s here) are those of the reference dumps of the
// original IDL code, which the version 1.x archives compare with validate.py.

#include <string>
#include <utility>
#include <vector>

#include "ena/config.hpp"
#include "ena/regions.hpp"
#include "ena/snapshot.hpp"

namespace ena {

bool checkpoints_available();  // false in builds without HDF5
void write_c0_config(const RunConfig& cfg, const std::string& dir);
void write_c1_inputs(const RunConfig& cfg, const Snapshot& s, const RadialLimits& lim, const std::string& dir);

}  // namespace ena

#include "ena/flux.hpp"
#include "ena/stream_processing.hpp"

namespace ena {

void write_c2_stream_points(const StreamPoints& s, const std::string& dir);
void write_c3_stream(const StreamFields& f, const std::string& dir, const std::string& name);
void write_c4_energy(const EnergyResult& r, int band, const std::string& dir);
void write_c5_band(const BandResult& b, const std::string& dir);

// Main output: band-averaged flux maps, one file per band (<dir>/flux_bandNN.h5).
// `info` holds run diagnostics stored with the maps (e.g. r_TS,V2).
struct RunInfo {
  double r_ts_reference = 0;   // AU (0 = not computed in this run)
  std::string r_ts_note;
  long streamlines = 0, streamlines_with_ts = 0;
};
// Extra contents of a band file (time-dependent maps: window, coverage).
struct ExtraOutput {
  std::vector<std::pair<std::string, const Field2D<double>*>> maps;
  std::vector<std::pair<std::string, double>> scalars;
  std::vector<std::pair<std::string, std::string>> strings;
};
// `file_name` replaces flux_bandNN (without .h5) when given.
void write_band_flux(const RunConfig& cfg, const BandResult& b, const std::string& dir, const RunInfo& info,
                     const std::string& file_name = "", const ExtraOutput* extra = nullptr);
// Legacy count-rate text files (one value per cell, one file per component). `suffix`
// is added to the file names before .dat (time-dependent maps: the window).
void write_legacy_countrate(const RunConfig& cfg, const BandResult& b, const std::string& dir, const std::string& suffix = "");

}  // namespace ena
