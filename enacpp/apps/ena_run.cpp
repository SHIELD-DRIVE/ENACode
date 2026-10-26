// ena_run: steady-state ENA maps for one snapshot (single-ion, extended single-ion
// or multi-ion model), or, with a
// [time_dependent] section in the configuration, time-dependent maps from a
// sequence of snapshots (include/ena/time_dependent.hpp).
//
//   ena_run run.toml [--snapshot N] [--out DIR] [--dump DIR] [--no-legacy] [--bands 9,13] [--threads N]
//
// --out   output directory (default: output_dir from the config)
// --dump  also write the checkpoints C0-C5: intermediate arrays, e.g. the per-cell
//         streamline fields (C3) and the per-energy fluxes (C4)
// --bands overrides [bands] run from the config
// --threads N  OpenMP threads per rank (default: OMP_NUM_THREADS, or all cores);
//         the maps are identical for any thread count
// MPI (build with -DENA_WITH_MPI=ON): mpiexec -n R ena_run run.toml --threads T
//         spreads a time-dependent run over R ranks of T threads each; the maps
//         are identical for any R and T. Steady-state runs use rank 0 only.
// --set table.key=value  overrides a configuration value (repeatable), e.g.
//         --set physics.fitted_esw=false --set physics.ts_detection=mach
// --stream-fields FILE  (testing) skip the streamline processing and take the
//         cell fields from a C3 checkpoint file, e.g. C3_stream_final.h5 of --dump

#include <chrono>
#include <cstdio>
#include <cstring>
#include <functional>
#include <filesystem>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#if defined(__GLIBC__)
#include <malloc.h>
#endif

#include "ena/runtime.hpp"
#include "ena/checkpoints.hpp"
#include "ena/comm.hpp"
#include "ena/config.hpp"
#include "ena/flux.hpp"
#include "ena/parallel.hpp"
#include "ena/pipeline.hpp"
#include "ena/regions.hpp"
#include "ena/snapshot.hpp"
#include "ena/stream_processing.hpp"
#include "ena/streamlines.hpp"
#include "ena/time_dependent.hpp"
#ifdef ENA_HAVE_HDF5
#include "ena/h5io.hpp"
#endif

using namespace ena;

namespace {

using Clock = std::chrono::steady_clock;
double since(Clock::time_point t0) { return std::chrono::duration<double>(Clock::now() - t0).count(); }

char pass_name_buf[64];
std::string pass_name(int pass) {
  if (pass < 0) return "C3_stream_final";
  std::snprintf(pass_name_buf, sizeof pass_name_buf, "C3_stream_after_pass%02d", pass);
  return pass_name_buf;
}

std::vector<int> parse_bands(const std::string& s) {
  std::vector<int> b;
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, ',')) b.push_back(std::stoi(item));
  return b;
}

}  // namespace

int run_main(int argc, char** argv) {
  const bool root = Comm::world().root();  // MPI rank 0 (always true without MPI)
#if defined(__GLIBC__)
  // The flux of every energy allocates several grid-size arrays (about 10 MB
  // each). glibc serves such blocks with mmap and returns them to the system on
  // free, so every energy page-faults its arrays in again (most of the system
  // time of a run). Keep freed memory in the heap for reuse instead.
  mallopt(M_MMAP_THRESHOLD, 1 << 30);
  mallopt(M_TRIM_THRESHOLD, 1 << 30);
#endif
  if (argc < 2) {
    std::fprintf(stderr, "usage: %s run.toml [--snapshot N] [--out DIR] [--dump DIR] [--no-legacy] [--bands 9,13] [--threads N]\n",
                 argv[0]);
    return 2;
  }
  try {
    std::string cfg_path = argv[1], dump_dir, out_dir, bands_arg, stream_file;
    std::size_t snap_index = 0;
    std::vector<std::string> overrides;
    bool legacy = true, snap_given = false;
    for (int a = 2; a < argc; ++a) {
      if (!std::strcmp(argv[a], "--snapshot") && a + 1 < argc) { snap_index = std::stoul(argv[++a]); snap_given = true; }
      else if (!std::strcmp(argv[a], "--dump") && a + 1 < argc) dump_dir = argv[++a];
      else if (!std::strcmp(argv[a], "--out") && a + 1 < argc) out_dir = argv[++a];
      else if (!std::strcmp(argv[a], "--bands") && a + 1 < argc) bands_arg = argv[++a];
      else if (!std::strcmp(argv[a], "--stream-fields") && a + 1 < argc) stream_file = argv[++a];
      else if (!std::strcmp(argv[a], "--set") && a + 1 < argc) overrides.push_back(argv[++a]);
      else if (!std::strcmp(argv[a], "--no-legacy")) legacy = false;
      else if (!std::strcmp(argv[a], "--threads") && a + 1 < argc) {
        const int n = std::stoi(argv[++a]);
        if (n < 1) throw std::runtime_error("--threads needs a positive number");
        set_threads(n);
      }
      else throw std::runtime_error(std::string("unknown argument ") + argv[a]);
    }
    const auto t_start = Clock::now();
    RunConfig cfg = load_config(cfg_path, overrides);
    if (!bands_arg.empty()) cfg.bands = parse_bands(bands_arg);
    for (int b : cfg.bands)
      for (int bad : cfg.sensors.unusable)
        if (b == bad) throw std::runtime_error("band " + std::to_string(b) + " is listed as unusable in [sensors]");
    if (out_dir.empty()) out_dir = cfg.output_dir;
    std::filesystem::create_directories(out_dir);
    const bool dump = !dump_dir.empty();
    if (dump) {
      if (!checkpoints_available()) throw std::runtime_error("--dump needs a build with HDF5");
      std::filesystem::create_directories(dump_dir);
    }
    const int nranks = Comm::world().size();
    if (root) {
      std::printf("Configuration %s: model %s, populations %s\n", cfg_path.c_str(), cfg.model_source.c_str(),
                  cfg.populations.c_str());
      std::printf("  threads: %d%s", max_threads(), openmp_enabled() ? " (OpenMP)" : " (built without OpenMP)");
      if (Comm::mpi_enabled()) std::printf(" per rank, %d MPI rank%s", nranks, nranks == 1 ? "" : "s");
      std::printf("\n  %s\n", settings_summary(cfg).c_str());
      for (const std::string& n : cfg.notes) std::printf("  note: %s\n", n.c_str());
    }

    // ---- time-dependent maps ([time_dependent]) ----
    if (cfg.time_dependent.enabled) {
      if (dump || !stream_file.empty() || snap_given)
        throw std::runtime_error("--dump, --stream-fields and --snapshot are for steady-state runs; remove the "
                                 "[time_dependent] section (or run without these options)");
      TimeDependentOptions td_opt;
      td_opt.out_dir = out_dir;
      td_opt.legacy = legacy;
      run_time_dependent(cfg, td_opt);
      if (root) std::printf("Total run time %.1f s; output in %s\n", since(t_start), out_dir.c_str());
      return 0;
    }

    // Steady-state maps run on one rank: with several MPI ranks, rank 0 does the work.
    if (nranks > 1) {
      if (!root) return 0;
      std::printf("  (steady-state run: only MPI rank 0 works; the other ranks are idle)\n");
    }

    // ---- inputs ----
    auto t0 = Clock::now();
    ReadOptions opt;
    opt.keep_file_coordinates = dump;  // the file coordinates are only needed for the C1 dump
    Snapshot snap = read_snapshot(cfg, snap_index, opt);
    std::printf("Snapshot %zu (time %g) read in %.1f s\n", snap_index, snap.time, since(t0));
    report_regions(cfg, snap, out_dir + "/heliosheath_pockets.txt", true);  // prints the full report
    RadialLimits lim = compute_radial_limits(cfg, snap.region);
    if (dump) {
      write_c0_config(cfg, dump_dir);
      write_c1_inputs(cfg, snap, lim, dump_dir);
    }

    StreamFields stream;
    RunInfo info;
    if (!stream_file.empty()) {
#ifdef ENA_HAVE_HDF5
      const Grid& g = cfg.grid;
      for (auto [n, fld] : std::vector<std::pair<const char*, Field3D<double>*>>{
               {"cool", &stream.cool}, {"vx", &stream.vx}, {"vy", &stream.vy}, {"vz", &stream.vz},
               {"ntr", &stream.ntr}, {"nref", &stream.nref}, {"nen", &stream.nen}, {"etr", &stream.etr},
               {"eref", &stream.eref}, {"een", &stream.een}, {"em", &stream.em}, {"kappa_ref", &stream.kappa_ref},
               {"vshock", &stream.vshock}, {"esw", &stream.esw}, {"ptr_ct", &stream.ptr_ct}, {"ptrn_ct", &stream.ptrn_ct}})
        *fld = read_grid_field_h5(stream_file, n, g);
      std::printf("Stream fields taken from %s (streamline processing skipped)\n", stream_file.c_str());
#else
      throw std::runtime_error("--stream-fields needs a build with HDF5");
#endif
    } else {
      std::function<void(const StreamPoints&)> on_points;
      StreamPassCallback on_pass;
      if (dump) {
        on_points = [&](const StreamPoints& pts) { write_c2_stream_points(pts, dump_dir); };
        on_pass = [&](int pass, const StreamFields& f) { write_c3_stream(f, dump_dir, pass_name(pass)); };
      }
      stream = build_stream_fields(cfg, snap, snap_index, info, true, on_points, on_pass);
    }

    // ---- flux and band integration ----
    FluxContext ctx{cfg, snap, stream};
    for (int band : cfg.bands) {
      t0 = Clock::now();
      int n_energy = 0;
      BandResult b = integrate_band(
          ctx, band, [&](double) { return dump; },
          [&](const EnergyResult& r) {
            ++n_energy;
            std::printf("  band %02d  E = %.5f keV\n", band, r.energy);
            std::fflush(stdout);
            if (dump) write_c4_energy(r, band, dump_dir);
          });
      if (dump) write_c5_band(b, dump_dir);
      if (checkpoints_available()) write_band_flux(cfg, b, out_dir, info);
      if (legacy) write_legacy_countrate(cfg, b, out_dir);
      std::printf("Band %02d (%.3f-%.3f keV, %d energies) done in %.1f s\n", band, b.emin, b.emax, n_energy, since(t0));
    }
    std::printf("Total run time %.1f s; output in %s\n", since(t_start), out_dir.c_str());
  } catch (const std::exception& e) {
    // errors of the time-dependent method are agreed by all ranks; rank 0 reports them
    if (root) std::fprintf(stderr, "error: %s\n", e.what());
    return 1;
  }
  return 0;
}

int main(int argc, char** argv) {
  prepare_runtime();  // before any HDF5 call (file locking off, see runtime.hpp)
  comm_init(&argc, &argv);
  const int rc = run_main(argc, argv);
  comm_finalize();
  return rc;
}
