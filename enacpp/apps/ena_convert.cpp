// ena_convert: converts the ASCII inputs of the snapshots to one HDF5 file per
// snapshot (fields in internal units, all ion fluids and neutral populations,
// plus the streamlines). This is a one-time step before the runs: the runs then
// read the binary files (inputs.format = "hdf5"), which is much faster than
// parsing the text files, above all when many ranks read snapshots at once.
//
//   ena_convert run.toml OUTDIR [--first N] [--last M] [--skip-existing] [--compress L]
//
// --first/--last  convert only snapshots N..M (0-based, inclusive), so the
//                 conversion can be split over several jobs or job-array tasks
// --skip-existing do not convert snapshots whose output file already exists
//                 (e.g. to resume an interrupted conversion)
// --compress L    gzip level of the datasets, 0-9 (default 0: uncompressed).
//                 Uncompressed files read fastest: on the full-size synthetic
//                 grid a snapshot and its streamlines load in 1.8 s against
//                 4.9 s from the text files, while level 1 or 4 takes about as
//                 long as the text files (decompression runs on one thread).
//                 Use a level > 0 only to save disk space.
// Output files: OUTDIR/snapshot_NNN.h5 (NNN = snapshot index, from 0). The
// program prints the configuration entries to use with the converted files.
// With MPI (mpiexec -n R ena_convert ...), the selected snapshots are split into
// R contiguous blocks, one per rank, and converted at the same time.

#include <cmath>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <string>

#include "ena/runtime.hpp"
#include "ena/comm.hpp"
#include "ena/config.hpp"
#include "ena/h5io.hpp"
#include "ena/snapshot.hpp"
#include "ena/streamlines.hpp"

using namespace ena;

int run_main(int argc, char** argv) {
  const Comm& comm = Comm::world();
  if (argc < 3) {
    std::fprintf(stderr, "usage: %s run.toml OUTDIR [--first N] [--last M] [--skip-existing] [--compress L]\n", argv[0]);
    return 2;
  }
  try {
    long first = 0, last = -1;
    bool skip = false;
    int level = 0;
    for (int a = 3; a < argc; ++a) {
      if (!std::strcmp(argv[a], "--first") && a + 1 < argc) first = std::stol(argv[++a]);
      else if (!std::strcmp(argv[a], "--last") && a + 1 < argc) last = std::stol(argv[++a]);
      else if (!std::strcmp(argv[a], "--skip-existing")) skip = true;
      else if (!std::strcmp(argv[a], "--compress") && a + 1 < argc) {
        level = std::stoi(argv[++a]);
        if (level < 0 || level > 9) throw std::runtime_error("--compress needs a level from 0 to 9");
      }
      else throw std::runtime_error(std::string("unknown argument ") + argv[a]);
    }
    RunConfig cfg = load_config(argv[1]);
    if (cfg.inputs.format != "ascii") throw std::runtime_error("the configuration must describe ASCII inputs");
    const long n = long(cfg.inputs.snapshots.size());
    if (last < 0 || last >= n) last = n - 1;
    if (first < 0 || first > last) throw std::runtime_error("--first/--last select no snapshot");
    std::string out = argv[2];
    std::filesystem::create_directories(out);
    long converted = 0, skipped = 0;
    // this rank's block of the selected snapshots (all of them without MPI)
    const long nsel = last - first + 1;
    const long my_first = first + block_begin(nsel, comm.size(), comm.rank());
    const long my_last = first + block_begin(nsel, comm.size(), comm.rank() + 1) - 1;
    std::string err;
    try {
    for (long i = my_first; i <= my_last; ++i) {
      const SnapshotFiles& f = cfg.inputs.snapshots[i];
      char name[64];
      std::snprintf(name, sizeof name, "snapshot_%03ld.h5", i);
      const std::string path = out + "/" + name;
      if (skip && std::filesystem::exists(path)) { ++skipped; continue; }
      ReadOptions ro;
      ro.verbose = true;
      Snapshot s = read_snapshot(cfg, std::size_t(i), ro);
      StreamlineSet sl;
      if (!f.streamlines.empty())
        sl = read_streamlines_tecplot(cfg.inputs.path(f.streamlines), int(cfg.inputs.streamlines.columns.size()));
      // write to a temporary name first, so an interrupted conversion never leaves a partial file
      const std::string tmp = path + ".part";
      write_snapshot_h5(tmp, cfg, s, f.streamlines.empty() ? nullptr : &sl, level);
      std::filesystem::rename(tmp, path);
      ++converted;
      std::printf("snapshot %ld -> %s%s\n", i, path.c_str(),
                  comm.size() > 1 ? (" (rank " + std::to_string(comm.rank()) + ")").c_str() : "");
      std::fflush(stdout);
    }
    } catch (const std::exception& e) {
      err = e.what();
    }
    comm.agree(err);  // every rank stops if one failed
    converted = long(std::lround(comm.sum_double(double(converted))));  // totals over the ranks
    skipped = long(std::lround(comm.sum_double(double(skipped))));
    if (!comm.root()) return 0;
    std::printf("\n%ld snapshots converted, %ld skipped (already there), in %s\n", converted, skipped, out.c_str());
    std::printf("\n# To run from the converted files, set in the configuration:\n#   [inputs]\n#   format = \"hdf5\"\n"
                "#   dir = \"%s\"   (relative paths: relative to the configuration file)\n", out.c_str());
    if (cfg.time_dependent.enabled)
      std::printf("# and replace the file names in [time_dependent.files] (or the [[snapshot]] entries) by\n"
                  "#   [time_dependent.files]\n#   file = \"snapshot_{index:03}.h5\"\n#   first_index = 0\n");
    else
      for (long i = first; i <= last; ++i)
        std::printf("[[snapshot]]\ntime = %.10g\nfile = \"snapshot_%03ld.h5\"\n\n", cfg.inputs.snapshots[i].time, i);
  } catch (const std::exception& e) {
    if (comm.root()) std::fprintf(stderr, "error: %s\n", e.what());
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
