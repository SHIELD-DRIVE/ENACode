// Milestone 2 tests: configuration, readers, region handling, streamlines,
// HDF5 round trip. Uses small synthetic files written in the same formats as
// the real inputs.

#include <cmath>
#include <cstdio>
#include <fstream>
#include <functional>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <vector>
#include <algorithm>

#include "ena/comm.hpp"
#include "ena/config.hpp"
#include "ena/regions.hpp"
#include "ena/snapshot.hpp"
#include "ena/streamlines.hpp"
#include "ena/time_dependent.hpp"
#ifdef ENA_HAVE_HDF5
#include "ena/h5io.hpp"
#endif

using namespace ena;

static int failures = 0, checks = 0;
#define CHECK(cond)                                                                 \
  do {                                                                              \
    ++checks;                                                                       \
    if (!(cond)) { ++failures; std::fprintf(stderr, "FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond); } \
  } while (0)

static bool throws_with(const std::function<void()>& f, const std::string& text) {
  try { f(); } catch (const std::exception& e) {
    if (std::string(e.what()).find(text) != std::string::npos) return true;
    std::fprintf(stderr, "  (unexpected message: %s)\n", e.what());
    return false;
  }
  return false;
}

static std::string dir;
static const int NT = 4, NP = 6, NR = 5;  // 45 deg x 60 deg x 2 AU grid
static const double DT = 45, DP = 60, DR = 2, R0 = 2;

static double plasma_n(int i, int j, int k) { return 0.1 + i + 0.01 * j + 0.001 * k; }
static int region_code(int i, int j, int k) { (void)j; return k < 2 ? 4 : (k < 4 ? 3 : (i % 2 ? 2 : 1)); }

static std::string fmt(double x) { char b[32]; std::snprintf(b, sizeof b, " %.9E", x); return b; }

// theta -> r -> phi order, like the real files
static void write_grid_files(bool shuffle_plasma, bool bad_boundary_order, double region_bump = 0.0) {
  std::ofstream pl(dir + "/plasma.dat"), nh(dir + "/neutral.dat"), bd(dir + "/boundaries.txt");
  std::vector<std::string> plasma_rows;
  for (int i = 0; i < NT; ++i)
    for (int k = 0; k < NR; ++k)
      for (int j = 0; j < NP; ++j) {
        double th = i * DT + (i == 0 ? 0.0006 : 0.02), ph = (i == 0 ? 0.0 : j * DP - 0.001), r = R0 + k * DR;
        std::string row = fmt(plasma_n(i, j, k)) + fmt(2.0e5 + k) + fmt(400.0 + j) + fmt(1.0) + fmt(-2.0) + fmt(3.0 + i) +
                          fmt(i) + fmt(j) + fmt(k) + fmt(th) + fmt(ph) + fmt(r);
        plasma_rows.push_back(row);
        nh << fmt(0.1 + 0.001 * k) << fmt(8000.0) << fmt(i) << fmt(j) << fmt(k) << fmt(th) << fmt(ph) << fmt(r) << "\n";
      }
  if (shuffle_plasma) std::swap(plasma_rows[3], plasma_rows[57]);
  for (auto& r : plasma_rows) pl << r << "\n";
  for (int i = 0; i < NT; ++i)
    for (int k = 0; k < NR; ++k)
      for (int j = 0; j < NP; ++j) {
        int jj = j, kk = k;
        if (bad_boundary_order && i == 1 && k == 0) { jj = (j == 0 ? 1 : j == 1 ? 0 : j); }
        double th = i * DT + 0.01, ph = (i == 0 ? 0.0 : jj * DP), r = R0 + kk * DR;
        double code = region_code(i, jj, kk) - (region_bump > 0 && i == 2 && kk == 2 ? region_bump : 0.0);
        bd << fmt(code) << fmt(i) << fmt(jj) << fmt(kk) << fmt(th) << fmt(ph) << fmt(r) << fmt(0.1) << fmt(-0.2) << fmt(0.3 * i) << "\n";
      }
}

static std::string base_toml(const std::string& extra = "") {
  return std::string("# test configuration\n[run]\noutput_dir = \"out\"\n") +
         "[grid]\ndtheta = 45.0\ndphi = 60.0\ndr = 2.0\nr_inner = 2.0\nnr = 5\n" +
         "[integration]\nr_lower = 2.0\nr_upper = 8.0\nradial_range = 2\n" +
         "[inputs]\ndir = \"" + dir + "\"\nformat = \"ascii\"\n" + extra +
         "[[snapshot]]\ntime = 0.0\nplasma = \"plasma.dat\"\nneutral = \"neutral.dat\"\nboundaries = \"boundaries.txt\"\n"
         "streamlines = \"stream.dat\"\n";
}

static RunConfig cfg_from(const std::string& toml) { return config_from_toml(parse_toml(toml, "test"), "test"); }

static void write_streamlines() {
  std::ofstream s(dir + "/stream.dat");
  s << "TITLE     = \"\"\nVARIABLES = \"X AU\"\n";
  const char* names[] = {"Y AU", "Z AU", "U_x km/s", "U_y km/s", "U_z km/s", "nH", "UxH", "UyH", "UzH", "np", "tp",
                         "vr", "vtheta", "vphi", "vt", "theta0", "phi0", "r0"};
  for (auto n : names) s << "\"" << n << "\"\n";
  int pts[3] = {4, 3, 5};
  for (int z = 0; z < 3; ++z) {
    s << "ZONE T=\"Streamtrace\"\n STRANDID=0, SOLUTIONTIME=0\n I=" << pts[z] << ", J=1, K=1, ZONETYPE=Ordered\n"
      << " DATAPACKING=POINT\n DT=(DOUBLE DOUBLE DOUBLE SINGLE)\n";
    for (int p = 0; p < pts[z]; ++p) {
      for (int c = 0; c < 19; ++c) s << fmt(z * 100 + p + 0.01 * c);
      s << "\n";
    }
  }
}

int main(int argc, char** argv) {
  dir = argc > 1 ? argv[1] : "test_tmp";
  mkdir(dir.c_str(), 0755);

  // --- TOML ---------------------------------------------------------------
  {
    auto t = parse_toml("a = 1.5 # c\nb = \"x\"\n[s.sub]\nlist = [\"p\", \"q\",\n  \"r\"]\nflag = true\n[[snap]]\nt = 1\n[[snap]]\nt = 2\n");
    CHECK(t.get_number("a", 0) == 1.5);
    CHECK(t.get_string("b", "") == "x");
    CHECK(t.table("s") && t.table("s")->table("sub") && t.table("s")->table("sub")->get_strings("list").size() == 3);
    CHECK(t.table("s")->table("sub")->get_bool("flag", false));
    CHECK(t.table_arrays.at("snap").size() == 2 && t.table_arrays.at("snap")[1].get_number("t", 0) == 2);
    CHECK(throws_with([] { parse_toml("a = \n"); }, "cannot parse"));
  }

  // --- configuration --------------------------------------------------------
  write_grid_files(false, false);
  write_streamlines();
  RunConfig cfg = cfg_from(base_toml());
  CHECK(cfg.grid.ntheta == NT && cfg.grid.nphi == NP && cfg.grid.nr == NR);
  // version 2 defaults: fitted E_sw, TS from the regions, r_TS,V2 computed only with a V2 direction
  CHECK(cfg.physics.fitted_esw && cfg.physics.ts_detection == "regions" && cfg.physics.ts_reference == "constant");
  CHECK(cfg.physics.shock_shift == "none");
  {
    RunConfig c2 = cfg_from(base_toml() + "[physics]\nv2_direction = [90.0, 30.0]\n");
    CHECK(c2.physics.ts_reference == "computed");
    CHECK(throws_with([] { cfg_from(base_toml() + "[quirks]\nQ11 = \"fixed\"\n"); }, "version 1.0 archive"));
    std::string t = base_toml();
    t.replace(t.find("[run]\n"), 6, "[run]\nmode = \"faithful\"\n");
    CHECK(throws_with([&] { cfg_from(t); }, "no faithful"));
  }
  CHECK(throws_with([] { cfg_from("[grid]\ndtheta = 7.0\nnr = 3\n[[snapshot]]\nplasma=\"a\"\nneutral=\"b\"\nboundaries=\"c\"\n"); }, "multiple of"));

  // --- readers --------------------------------------------------------------
  {
    ReadOptions opt;
    opt.keep_file_coordinates = true;
    opt.verbose = false;
    Snapshot s = read_snapshot(cfg, 0, opt);
    bool ok_n = true, ok_T = true;
    for (int i = 0; i < NT; ++i)
      for (int j = 0; j < NP; ++j)
        for (int k = 0; k < NR; ++k) {
          ok_n &= std::fabs(s.plasma[0].n(i, j, k) - plasma_n(i, j, k)) < 1e-12;
          ok_T &= std::fabs(s.plasma[0].T(i, j, k) - 0.5 * (2.0e5 + k)) < 1e-9;  // halved (standard model)
        }
    CHECK(ok_n);
    CHECK(ok_T);
    CHECK(std::fabs(s.plasma[0].vr(1, 2, 3) - (402.0 * 1e5)) < 1e-3);  // km/s -> cm/s
    CHECK(s.region(1, 0, 4) == 2 && s.region(2, 0, 4) == 1 && s.region(0, 3, 0) == 4 && s.region(3, 5, 2) == 3);
    CHECK(s.Bz(3, 1, 1) == 0.9);                      // double precision
    CHECK(std::fabs(s.coords.theta(0, 0, 0) - 0.0006) < 1e-12);
    CHECK(std::fabs(s.neutrals[0].n(2, 3, 4) - 0.104) < 1e-12);

    // plasma rows out of order are placed by their i,j,k columns
    write_grid_files(true, false);
    Snapshot s2 = read_snapshot(cfg, 0, opt);
    CHECK(s2.plasma[0].n.raw() == s.plasma[0].n.raw());

    // boundaries rows in another order are placed by their i,j,k columns too
    write_grid_files(false, true);
    Snapshot sb = read_snapshot(cfg, 0, opt);
    CHECK(sb.region.raw() == s.region.raw() && sb.Bz.raw() == s.Bz.raw());

    // region codes are truncated (2.9999 -> 2)
    write_grid_files(false, false, 0.0001);
    Snapshot s3 = read_snapshot(cfg, 0, opt);
    CHECK(s3.region(2, 0, 2) == 2);
  }

  // coordinate mismatch, duplicates, missing cells
  {
    write_grid_files(false, false);
    {
      std::ifstream in(dir + "/neutral.dat");
      std::stringstream ss; ss << in.rdbuf();
      std::string t = ss.str();
      auto pos = t.find('\n', t.find('\n') + 1);  // after line 2
      std::string line3 = t.substr(pos + 1, t.find('\n', pos + 1) - pos - 1);
      std::ofstream out(dir + "/neutral_dup.dat");
      out << t.substr(0, pos + 1) << line3 << "\n" << t.substr(pos + 1);
    }
    RunConfig c = cfg_from(base_toml());
    c.inputs.snapshots[0].neutral = "neutral_dup.dat";
    CHECK(throws_with([&] { read_snapshot(c, 0); }, "appears twice"));

    std::ifstream in(dir + "/plasma.dat");
    std::string all((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    std::ofstream(dir + "/plasma_short.dat") << all.substr(0, all.rfind('\n', all.size() - 2) + 1);
    c = cfg_from(base_toml());
    c.inputs.snapshots[0].plasma = "plasma_short.dat";
    CHECK(throws_with([&] { read_snapshot(c, 0); }, "have no row"));
    c.inputs.allow_incomplete = true;
    ReadOptions quiet; quiet.verbose = false;
    CHECK(!throws_with([&] { read_snapshot(c, 0, quiet); }, "have no row"));

    // wrong r for k
    std::string bad = all;
    auto p = bad.find(" 2.000000000E+00\n");
    bad.replace(p, 16, " 9.000000000E+00");
    std::ofstream(dir + "/plasma_badr.dat") << bad;
    c = cfg_from(base_toml());
    c.inputs.snapshots[0].plasma = "plasma_badr.dat";
    CHECK(throws_with([&] { read_snapshot(c, 0); }, "does not match grid point"));
  }

  // --- Moscow remap and undetermined codes ------------------------------------
  {
    write_grid_files(false, false);
    std::ifstream in(dir + "/boundaries.txt");
    std::string all((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    all.replace(1, 15, "0.000000000E+00");  // first cell: Moscow code 0
    std::ofstream(dir + "/boundaries_moscow.txt") << all;
    std::string toml = base_toml();
    toml.replace(toml.find("[grid]"), 6, "[model]\nsource = \"moscow\"\n[grid]");
    RunConfig c = cfg_from(toml);
    c.inputs.snapshots[0].boundaries = "boundaries_moscow.txt";
    ReadOptions quiet; quiet.verbose = false;
    Snapshot s = read_snapshot(c, 0, quiet);
    CHECK(s.region(0, 0, 0) == 5);                 // Moscow 0 -> 5, undetermined
    CHECK(s.region(0, 1, 0) == 2);                 // Moscow 4 -> 2
    CHECK(s.region(1, 0, 4) == 3);                 // Moscow 2 -> 3
    CHECK(std::fabs(s.plasma[0].T(0, 0, 0) - 2.0e5) < 1e-9);  // Moscow: not halved
    CHECK(count_undetermined(c, s.region) == 1);
  }

  // --- radial limits (get_limits) --------------------------------------------
  {
    write_grid_files(false, false);
    ReadOptions quiet; quiet.verbose = false;
    Snapshot s = read_snapshot(cfg, 0, quiet);
    // Region 3 starts at k = 2 (r = 6 AU, the true radius): inside 60 AU, so never "found".
    RadialLimits L = compute_radial_limits(cfg, s.region);
    CHECK(L.lower(1, 1) == float(NR - 1));
    RunConfig c = cfg;
    c.grid.r_inner = 58.0;  // r = 58 + 2k: k = 2 is 62 AU > 60: found
    L = compute_radial_limits(c, s.region);
    CHECK(L.lower(1, 1) == 2.f);
    c.radial_range = 0;
    L = compute_radial_limits(c, s.region);
    CHECK(L.lower(1, 1) == 2.f && L.upper(1, 1) == 3.f);  // region 2 at k = 4
  }

  // --- streamlines -----------------------------------------------------------
  {
    StreamlineSet sl = read_streamlines_tecplot(dir + "/stream.dat", 19);
    CHECK(sl.count() == 3);
    CHECK(sl.points() == 12);
    CHECK(sl.variables.size() == 19 && sl.variables[0] == "X AU" && sl.variables[18] == "r0");
    CHECK(sl.offsets[1] == 4 && sl.offsets[2] == 7);
    CHECK(sl.gap_lines[0] == 25 && sl.gap_lines[1] == 5 && sl.gap_lines[2] == 5);
    CHECK(std::fabs(sl.at(sl.offsets[2] + 1, 3) - 201.03) < 1e-9);
    CHECK(throws_with([&] { read_streamlines_tecplot(dir + "/stream.dat", 20); }, "expected 20 numbers"));
  }

  // --- time-dependent method: configuration, windows, path survival -------------
  {
    const std::string grid = "[grid]\ndtheta = 45.0\ndphi = 60.0\ndr = 2.0\nr_inner = 2.0\nnr = 5\n";
    const std::string files = "[time_dependent.files]\nplasma = \"p_{index:03}.dat\"\nneutral = \"n.dat\"\n"
                              "boundaries = \"b_{index}.txt\"\nfirst_index = 7\n";
    RunConfig c = cfg_from(grid + "[time_dependent]\nstart = 2010.0\nend = 2012.0\ndt_months = 2\ncadence_months = 6\n" + files);
    const TimeDependentConfig& td = c.time_dependent;
    CHECK(td.enabled && td.snapshots() == 12 && std::fabs(td.dt - 2.0 / 12.0) < 1e-15);
    CHECK(c.inputs.snapshots.size() == 12 && c.inputs.snapshots[0].plasma == "p_007.dat" &&
          c.inputs.snapshots[11].plasma == "p_018.dat" && c.inputs.snapshots[3].boundaries == "b_10.txt");
    CHECK(std::fabs(c.inputs.snapshots[5].time - (2010.0 + 5.0 / 6.0)) < 1e-12);
    CHECK(td.window == td.cadence && td.output_start == 2010.0 && td.output_end == 2012.0);
    CHECK(td.coverage_extinction == 0.3 && td.coverage_missing_flux == 0.0 && !td.include_incomplete_lines);
    auto w = output_windows(td);
    CHECK(w.size() == 4 && w[1].label == "2010.5-2011" && w[3].b == 2012.0);
    // running averages: 1-yr windows every 6 months
    RunConfig c2 = cfg_from(grid + "[time_dependent]\nstart = 2010.0\nend = 2014.0\ndt = 0.5\noutput_start = 2011.0\n"
                            "cadence = 0.5\nwindow = 1.0\n" + files);
    auto w2 = output_windows(c2.time_dependent);
    CHECK(w2.size() == 5 && w2[0].label == "2011-2012" && w2[4].label == "2013-2014");
    CHECK(time_label(2020.3333333) == "2020.3333" && time_label(2021.0) == "2021");
    CHECK(throws_with([&] { cfg_from(grid + "[time_dependent]\nstart = 2010.0\nend = 2011.0\ndt = 0.3\n" + files); },
                      "whole number"));
    CHECK(throws_with([&] { cfg_from(grid + "[time_dependent]\nstart = 2010.0\nend = 2011.0\ndt = 0.5\noutput_end = 2012.0\n" + files); },
                      "inside [start, end]"));
    CHECK(throws_with([&] { cfg_from(grid + "[time_dependent]\nstart = 2010.0\nend = 2011.0\ndt = 0.5\nmonths = 2\n" + files); },
                      "unknown setting"));
    // [[snapshot]] times must be the starts of the time cells
    const std::string snap2 = "[[snapshot]]\ntime = 2010.0\nplasma=\"a\"\nneutral=\"b\"\nboundaries=\"c\"\n"
                              "[[snapshot]]\ntime = 2010.4\nplasma=\"a\"\nneutral=\"b\"\nboundaries=\"c\"\n";
    CHECK(throws_with([&] { cfg_from(grid + "[time_dependent]\nstart = 2010.0\nend = 2011.0\ndt = 0.5\n" + snap2); },
                      "time cell starts"));
    CHECK(throws_with([&] { cfg_from(grid + "[time_dependent]\nstart = 2010.0\nend = 2012.0\ndt = 0.5\n" + snap2); },
                      "4 snapshots, but 2"));

    CHECK(time_cell(-0.1, 1.0, 5) == 0 && time_cell(2.999, 1.0, 5) == 2 && time_cell(3.0, 1.0, 5) == 3 &&
          time_cell(7.0, 1.0, 5) == 4);
    // path exponent: prefix sums against a cell-by-cell sum, for many speeds and emission times
    const int nsnap = 7, k0 = 3, ncell = 40;
    std::vector<std::vector<double>> b(nsnap, std::vector<double>(ncell)), pre(nsnap, std::vector<double>(ncell - k0));
    for (int s = 0; s < nsnap; ++s) {
      double run = 0;
      for (int k = 0; k < ncell; ++k) {
        b[s][k] = (k % 5 == 2) ? 0.0 : 1.0 + 0.37 * s + 0.011 * k + 0.0001 * s * k;  // some cells without losses
        if (k >= k0) { run += b[s][k]; pre[s][k - k0] = run; }
      }
    }
    std::vector<const double*> pp;
    for (auto& v : pre) pp.push_back(v.data());
    double worst = 0;
    for (double cy : {0.001, 0.0137, 0.1, 0.25, 0.5, 1.0, 3.3, 50.0})
      for (int s = 0; s < nsnap; ++s)
        for (int k = k0; k < ncell; ++k) {
          const double te = (s + 0.5) * 0.5;
          double ref = 0;
          for (int kk = k0; kk <= k; ++kk) ref += b[time_cell(te + (k - kk) * cy, 0.5, nsnap)][kk];
          const double got = path_exponent(pp, k0, k, te, 0.5, nsnap, cy);
          worst = std::max(worst, std::fabs(got - ref) / ref);
        }
    CHECK(worst < 1e-13);
  }

  // --- MPI layout helpers (one rank here) -----------------------------------------
  {
    // block partition: contiguous, complete, sizes differ by at most one
    bool ok = true;
    for (long n : {0L, 1L, 7L, 30L, 1800L})
      for (int parts : {1, 3, 8, 13}) {
        long prev = block_begin(n, parts, 0), lo = n, hi = 0;
        ok = ok && prev == 0 && block_begin(n, parts, parts) == n;
        for (int i = 1; i <= parts; ++i) {
          const long b = block_begin(n, parts, i), size = b - prev;
          ok = ok && size >= 0;
          lo = std::min(lo, size);
          hi = std::max(hi, size);
          prev = b;
        }
        ok = ok && hi - lo <= 1;
      }
    CHECK(ok);
    const Comm& c = Comm::world();
    CHECK(c.size() == 1 && c.root());
    auto r = c.alltoallv({{1.0, 2.0}});
    CHECK(r.size() == 1 && r[0] == std::vector<double>({1.0, 2.0}));
    CHECK(throws_with([&] { c.agree("boom"); }, "boom"));
    CHECK(c.sum_double(2.5) == 2.5);
  }

#ifdef ENA_HAVE_HDF5
  // --- HDF5 snapshot round trip ------------------------------------------------
  {
    write_grid_files(false, false);
    ReadOptions quiet; quiet.verbose = false;
    Snapshot s = read_snapshot(cfg, 0, quiet);
    StreamlineSet sl = read_streamlines_tecplot(dir + "/stream.dat", 19);
    write_snapshot_h5(dir + "/snap.h5", cfg, s, &sl);
    Snapshot r = read_snapshot_h5(cfg, dir + "/snap.h5", quiet);
    CHECK(r.plasma[0].n.raw() == s.plasma[0].n.raw());
    CHECK(r.plasma[0].vt.raw() == s.plasma[0].vt.raw());
    CHECK(r.neutrals[0].T.raw() == s.neutrals[0].T.raw());
    CHECK(r.region.raw() == s.region.raw());
    CHECK(r.Bz.raw() == s.Bz.raw());
    StreamlineSet rs = read_streamlines_h5(dir + "/snap.h5");
    CHECK(rs.values == sl.values && rs.offsets == sl.offsets && rs.gap_lines == sl.gap_lines);
    RunConfig c = cfg;
    c.grid.dr = 3.0;
    CHECK(throws_with([&] { read_snapshot_h5(c, dir + "/snap.h5", quiet); }, "grid stored in the file differs"));
  }
#endif

  std::printf("%d checks, %d failures\n", checks, failures);
  return failures ? 1 : 0;
}
