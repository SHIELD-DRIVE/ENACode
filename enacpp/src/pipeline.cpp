#include "ena/pipeline.hpp"

#include <cmath>
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <stdexcept>
#include <utility>
#include <vector>

#include "ena/regions.hpp"
#include "ena/streamlines.hpp"
#ifdef ENA_HAVE_HDF5
#include "ena/h5io.hpp"
#endif

namespace ena {

namespace {
using Clock = std::chrono::steady_clock;
double since(Clock::time_point t0) { return std::chrono::duration<double>(Clock::now() - t0).count(); }
}  // namespace

std::string report_regions(const RunConfig& cfg, const Snapshot& snap, const std::string& pockets_file, bool verbose) {
  const long und = count_undetermined(cfg, snap.region);
  const ReentryStats re = count_reentries(cfg, snap.region);
  const auto pockets = heliosheath_pockets(cfg, snap.region);
  if (!verbose) {  // short summary, only if something is found
    if (!(und || re.lines || re.exits_at_undetermined || !pockets.empty())) return "";
    char b[256];
    std::snprintf(b, sizeof b, "regions: %ld undetermined cells, %ld lines of sight re-enter the heliosheath, %ld end at an "
                  "undetermined cell, %zu supersonic pockets", und, re.lines, re.exits_at_undetermined, pockets.size());
    return b;
  }
  if (und) std::printf("  WARNING: %ld undetermined (region 5) cells inside the integration window\n", und);
  if (re.lines) {
    std::printf("  %ld lines of sight leave the heliosheath and re-enter it (%ld region 3/6 cells beyond the first exit, "
                "from %.0f AU); integration stops at the first exit\n", re.lines, re.cells, re.r_first);
    std::printf("%s", reentry_details(cfg, snap.region).c_str());
  } else {
    std::printf("  No line of sight re-enters the heliosheath inside the integration window\n");
  }
  if (re.exits_at_undetermined)
    std::printf("  WARNING: %ld lines of sight leave the heliosheath into an undetermined (region 5) cell, first at "
                "theta %g, phi %g, r %g AU; integration stops there\n", re.exits_at_undetermined, re.u_theta,
                re.u_phi, re.u_r);
  // Supersonic (region 4) pockets inside the heliosheath are skipped by the
  // integration; list them so the region file can be checked.
  if (!pockets.empty()) {
    long cells = 0, nlos = 0;
    int lt = -1, lp = -1;
    for (const auto& q : pockets) {
      cells += q.k1 - q.k0 + 1;
      if (q.t != lt || q.p != lp) { ++nlos; lt = q.t; lp = q.p; }
    }
    std::printf("  WARNING: %zu supersonic (region 4) pockets inside the heliosheath on %ld lines of sight (%ld cells) "
                "are skipped: no ENAs from them, integration continues behind them. Check the regions if unexpected.\n",
                pockets.size(), nlos, cells);
    if (FILE* f = std::fopen(pockets_file.c_str(), "w")) {
      std::fprintf(f, "# Supersonic (region 4) pockets inside the heliosheath, skipped by the integration\n");
      std::fprintf(f, "# theta_deg phi_deg r_start_AU r_end_AU cells theta_index phi_index r_index_start r_index_end\n");
      for (const auto& q : pockets)
        std::fprintf(f, "%6.1f %6.1f %8.1f %8.1f %4d %3d %3d %4d %4d\n", q.t * cfg.grid.dtheta, q.p * cfg.grid.dphi,
                     cfg.grid.r(q.k0), cfg.grid.r(q.k1), q.k1 - q.k0 + 1, q.t, q.p, q.k0, q.k1);
      std::fclose(f);
    }
    const std::size_t show = std::min<std::size_t>(pockets.size(), 10);
    for (std::size_t n = 0; n < show; ++n) {
      const auto& q = pockets[n];
      std::printf("    pocket at theta %5.1f, phi %5.1f: r %g-%g AU (%d cells)\n", q.t * cfg.grid.dtheta,
                  q.p * cfg.grid.dphi, cfg.grid.r(q.k0), cfg.grid.r(q.k1), q.k1 - q.k0 + 1);
    }
    std::printf("    %s (full list: %s)\n", pockets.size() > show ? "..." : "", pockets_file.c_str());
  }
  return "";
}

StreamFields build_stream_fields(const RunConfig& cfg, const Snapshot& snap, std::size_t index, RunInfo& info,
                                 bool verbose, const std::function<void(const StreamPoints&)>& on_points,
                                 const StreamPassCallback& on_pass) {
  const SnapshotFiles& files = cfg.inputs.snapshots.at(index);
  auto t0 = Clock::now();
  StreamlineSet sl;
  if (cfg.inputs.format == "hdf5") {
#ifdef ENA_HAVE_HDF5
    sl = read_streamlines_h5(cfg.inputs.path(files.file));
#else
    throw std::runtime_error("HDF5 inputs need a build with HDF5");
#endif
  } else {
    if (files.streamlines.empty()) throw std::runtime_error("no streamline file given for this snapshot");
    sl = read_streamlines_tecplot(cfg.inputs.path(files.streamlines), int(cfg.inputs.streamlines.columns.size()));
  }
  if (verbose) std::printf("Streamlines: %zu streamlines, %zu points, read in %.1f s\n", sl.count(), sl.points(), since(t0));

  // ---- streamline processing ----
  t0 = Clock::now();
  std::string ref_note;
  const double r_ref = reference_ts_distance(cfg, snap, &ref_note);
  info.r_ts_reference = r_ref;
  info.r_ts_note = ref_note;
  if (verbose) {
    std::printf("Reference TS distance r_TS,V2: %s\n", ref_note.c_str());
    if (cfg.physics.ts_reference == "constant" && cfg.physics.ts_detection == "regions")
      std::printf("  WARNING: the TS is found from the regions but r_TS,V2 is a constant; set physics.v2_direction so "
                  "r_TS,V2 is measured with the same method (otherwise the partition fits are biased)\n");
  }
  StreamPoints pts = compute_stream_points(cfg, sl, r_ref, &snap.region);
  info.streamlines = long(pts.streamlines);
  info.streamlines_with_ts = long(pts.streamlines_with_ts);
  if (verbose) {
    std::printf("TS detection (%s): found on %zu of %zu streamlines (%zu crossings)\n", cfg.physics.ts_detection.c_str(),
                pts.streamlines_with_ts, pts.streamlines, pts.ts_crossings);
    if (cfg.extended_si())
      std::printf("Injection energy (extended single-ion): E_l = 1/2 m_p v^2, v = flow speed at the last point inside 12 AU%s\n",
                  pts.esi_no_inner_point
                      ? (" (WARNING: " + std::to_string(pts.esi_no_inner_point) +
                         " TS crossings on streamlines without a point inside 12 AU use the first point)").c_str()
                      : "");
    if (cfg.physics.ts_detection == "regions")
      std::printf("Flow speeds around the region-based TS (diagnostic):\n%s", crossing_diagnostics(cfg, pts).c_str());
  }
  if (on_points) on_points(pts);
  StreamFields stream = bin_and_fill(cfg, pts, snap.region, on_pass);
  if (verbose) {
    std::printf("Gap fill: %d passes; region 3/6 cells in the window: %ld, still empty: %ld (extinction), %ld (fractions)\n",
                stream.passes, stream.cells, stream.unfilled_cool, stream.unfilled_fractions);
    if (stream.unfilled_cool || stream.unfilled_fractions) {
      std::printf("  WARNING: region 3/6 cells left empty (no streamline data in their connected part of the heliosheath)\n");
      // where they are, and what the flux calculation does with them ([physics] unfilled_cells)
      const Grid& g = cfg.grid;
      const int rlower = g.k_of(cfg.r_lower), rupper = g.k_of(cfg.r_upper);
      int shown = 0;
      for (int t = 0; t < g.ntheta; ++t)
        for (int pp = 0; pp < g.nphi; ++pp)
          for (int k = rlower; k <= rupper; ++k) {
            const int b = snap.region(t, pp, k);
            if (b != 3 && b != 6) continue;
            const double v = stream.ntr(t, pp, k);
            if (v != 0.0 && std::isfinite(v)) continue;
            if (shown++ < 20)
              std::printf("    empty cell: theta %5.1f, phi %5.1f, r %6.1f AU (indices %d, %d, %d), region %d\n", g.theta(t),
                          g.phi(pp), g.r(k), t, pp, k, b);
          }
      if (shown > 20) std::printf("    ... %d more\n", shown - 20);
      std::printf("    [physics] unfilled_cells = \"%s\": %s\n", cfg.physics.unfilled_cells.c_str(),
                  cfg.physics.unfilled_cells == "skip"
                      ? "these cells give no ENAs; integration continues behind them"
                      : "a line of sight that reaches one of these cells stops the run with an error");
    }
    if (stream.has_u_inj)
      std::printf("Injection speed (multi-ion, carried from the upstream side of the region-based TS): region 3/6 cells "
                  "without a value: %ld%s\n", stream.unfilled_u_inj,
                  stream.unfilled_u_inj ? " (no accelerated PUIs there)" : "");
    std::printf("Streamline processing (partition fits, extinction, binning, gap fill): %.1f s\n", since(t0));
  }
  return stream;
}

}  // namespace ena
