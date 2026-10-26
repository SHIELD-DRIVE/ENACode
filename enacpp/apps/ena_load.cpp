// ena_load: reads one snapshot of a run, prints a summary and (with HDF5)
// writes the C0/C1a/C1b checkpoints (the inputs as read).
//
//   ena_load run.toml [--snapshot N] [--dump DIR] [--no-streamlines]

#include <chrono>
#include <cstdio>
#include <cstring>
#include <limits>
#include <map>
#include <string>

#include "ena/runtime.hpp"
#include "ena/checkpoints.hpp"
#include "ena/config.hpp"
#include "ena/regions.hpp"
#include "ena/snapshot.hpp"
#include "ena/streamlines.hpp"
#ifdef ENA_HAVE_HDF5
#include "ena/h5io.hpp"
#endif

using namespace ena;

static void range(const char* name, const Field3D<double>& f) {
  double lo = std::numeric_limits<double>::infinity(), hi = -lo;
  long zeros = 0;
  for (double x : f.raw()) { lo = std::min(lo, x); hi = std::max(hi, x); zeros += (x == 0.0); }
  std::printf("  %-12s min %12.5e  max %12.5e  zeros %ld\n", name, lo, hi, zeros);
}

int main(int argc, char** argv) {
  prepare_runtime();  // before any HDF5 call (file locking off, see runtime.hpp)
  if (argc < 2) {
    std::fprintf(stderr, "usage: %s run.toml [--snapshot N] [--dump DIR] [--no-streamlines]\n", argv[0]);
    return 2;
  }
  try {
    std::string cfg_path = argv[1], dump_dir;
    std::size_t snap = 0;
    bool streamlines = true;
    for (int a = 2; a < argc; ++a) {
      if (!std::strcmp(argv[a], "--snapshot") && a + 1 < argc) snap = std::stoul(argv[++a]);
      else if (!std::strcmp(argv[a], "--dump") && a + 1 < argc) dump_dir = argv[++a];
      else if (!std::strcmp(argv[a], "--no-streamlines")) streamlines = false;
      else throw std::runtime_error(std::string("unknown argument ") + argv[a]);
    }
    RunConfig cfg = load_config(cfg_path);
    const Grid& g = cfg.grid;
    std::printf("Configuration %s: model %s, grid %d x %d x %d (dtheta %g, dphi %g, dr %g, r_inner %g)\n",
                cfg_path.c_str(), cfg.model_source.c_str(), g.ntheta, g.nphi, g.nr,
                g.dtheta, g.dphi, g.dr, g.r_inner);

    auto t0 = std::chrono::steady_clock::now();
    ReadOptions opt;
    opt.keep_file_coordinates = !dump_dir.empty();
    Snapshot s = read_snapshot(cfg, snap, opt);
    double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::printf("Snapshot %zu (time %g) read in %.1f s\n", snap, s.time, secs);
    range("n_p", s.plasma[0].n);
    range("T_p", s.plasma[0].T);
    range("vr [cm/s]", s.plasma[0].vr);
    range("vt [cm/s]", s.plasma[0].vt);
    range("n_H", s.neutrals[0].n);
    std::map<int, long> hist;
    for (auto r : s.region.raw()) hist[r]++;
    std::printf("  regions:");
    for (auto& [code, n] : hist) std::printf("  %d: %ld", code, n);
    std::printf("\n");
    long und = count_undetermined(cfg, s.region);
    if (und) std::printf("  WARNING: %ld undetermined (region 5) cells inside the integration window\n", und);

    const SnapshotFiles& files = cfg.inputs.snapshots[snap];
    if (streamlines && cfg.inputs.format == "ascii" && !files.streamlines.empty()) {
      t0 = std::chrono::steady_clock::now();
      StreamlineSet sl = read_streamlines_tecplot(cfg.inputs.path(files.streamlines), int(cfg.inputs.streamlines.columns.size()));
      secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
      std::printf("Streamlines: %zu streamlines, %zu points, %zu variables, read in %.1f s\n", sl.count(), sl.points(),
                  sl.variables.size(), secs);
    }

    RadialLimits lim = compute_radial_limits(cfg, s.region);
    if (!dump_dir.empty()) {
      if (!checkpoints_available()) throw std::runtime_error("--dump needs a build with HDF5");
      write_c0_config(cfg, dump_dir);
      write_c1_inputs(cfg, s, lim, dump_dir);
      std::printf("Checkpoints C0_config, C1a_limits, C1b_inputs written to %s\n", dump_dir.c_str());
    }
  } catch (const std::exception& e) {
    std::fprintf(stderr, "error: %s\n", e.what());
    return 1;
  }
  return 0;
}
