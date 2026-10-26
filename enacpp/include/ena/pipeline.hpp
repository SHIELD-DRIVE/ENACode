#pragma once
// Steps shared by the steady-state run (ena_run) and the time-dependent method:
// the region checks of a snapshot and the streamline processing.

#include <cstddef>
#include <functional>
#include <string>

#include "ena/checkpoints.hpp"
#include "ena/config.hpp"
#include "ena/snapshot.hpp"
#include "ena/stream_processing.hpp"

namespace ena {

// Region checks of one snapshot: undetermined cells, lines of sight that
// re-enter the heliosheath, supersonic pockets. With `verbose` the full report
// of ena_run is printed (pockets listed in `pockets_file`, if any) and "" is
// returned; otherwise nothing is printed and a short summary is returned ("" if
// nothing is found).
std::string report_regions(const RunConfig& cfg, const Snapshot& snap, const std::string& pockets_file, bool verbose);

// Reads the streamlines of snapshot `index` and builds the cell fields (partition
// fits, extinction, binning, gap fill). `info` receives r_TS,V2 and the
// streamline counts. `on_points` (optional) sees the streamline points (C2 dump),
// `on_pass` the fields after each gap-fill pass (C3 dumps).
StreamFields build_stream_fields(const RunConfig& cfg, const Snapshot& snap, std::size_t index, RunInfo& info,
                                 bool verbose, const std::function<void(const StreamPoints&)>& on_points = nullptr,
                                 const StreamPassCallback& on_pass = nullptr);

}  // namespace ena
