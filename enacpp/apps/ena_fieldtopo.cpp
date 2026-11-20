// ena_fieldtopo: magnetic field-line topology on the ENA grid, to locate the heliopause where the
// region codes cannot (e.g. the low-latitude heliotail).
//
//   ena_fieldtopo run.toml [--snapshot N | --snapshots all|A-B|A,B,...] [--out topology.dat]
//                 [--z-top 500] [--z-bottom -500] [--ism-regions 1] [--r-sun 50] [--max-length 10000]
//                 [--cell-fraction 1] [--substeps 2] [--only-regions 3,6] [--all-radii]
//                 [--rewrite-boundaries] [--suffix _topo] [--overwrite] [--edge-ism all|nose|none]
//                 [--threads N]
//
// The field line through every grid point is followed in both directions, along B/|B|, through the
// magnetic field of the snapshot's boundaries file (Bx, By, Bz). Each direction ends where the line
//   - reaches the inner heliosphere, r <= r_sun ("Sun");
//   - reaches the interstellar medium ("ISM"): z >= z_top or z <= z_bottom, or a cell whose region
//     code is one of ism_regions, or the outer edge of the grid;
//   - or neither within max_length of path ("long"), or stops in a null of the field ("null").
// --edge-ism nose counts leaving through the outer edge as reaching the ISM only where the line leaves
// at x < 0 (the nose side); down the tail (x >= 0) a line can leave the grid while still in the
// heliotail, so there it ends as "grid edge (not ISM)" (end code 7) and counts like a long line.
// --edge-ism none: no exit through the edge counts as ISM. Default all.
// The point is then classified by its two ends:
//   1 closed:        no end in the ISM (Sun/Sun, Sun/long, long/long): solar field, inside the heliopause
//   2 open:          one end at the Sun, the other in the ISM: reconnected, the mixing region
//   3 interstellar:  both ends in the ISM
//   4 ISM-connected: one end in the ISM, the other neither at the Sun nor in the ISM (long or null)
//   0 undetermined:  both directions stopped in a field null; -1: not traced (--only-regions)
//
// Method (as ena_trace): classical 4th-order Runge-Kutta along B/|B|; one step is cell_fraction x
// min(dr, r dtheta), split into `substeps` Runge-Kutta steps; B is interpolated trilinearly in the
// grid's (r, theta, phi) index space, phi periodic; theta and phi of a point come from its position,
// so lines cross the north pole. The grid has no theta = 180 deg row: between the last row and the
// south pole the last row's values are used. Region codes are taken at the nearest grid point.
// Coordinates are the simulation frame of the MHD solution: x = r sin(theta) cos(phi),
// y = r sin(theta) sin(phi), z = r cos(theta) (AU), as X AU, Y AU, Z AU in the Tecplot macros.
//
// Output: a Tecplot ASCII file (one ordered zone, I = r, J = phi, K = theta, point format, with a
// header), on the radii of the integration window r_lower..r_upper of the settings file (or all
// radii with --all-radii), with the last theta row repeated at theta = 180 deg so that Tecplot
// closes the grid at the south pole: X, Y, Z, region, topology, the end of each direction and its
// path length, the path length to the ISM, Bx, By, Bz, |B|, and |B_phi|/|B| (the azimuthal share of
// the field about the z axis, close to 1 in the Parker spiral). The results do not depend on the
// thread count.
//
// Heliopause from the topology (--rewrite-boundaries): the heliosheath cells (regions 3 and 6) of every
// radius are traced, and a copy of the snapshot's boundaries file is written in which each of them
// whose field line is not closed (topology 2, 3 or 4: open, interstellar, ISM-connected) is region 2.
// So the heliopause becomes the boundary of the closed solar field, e.g. in the low-latitude heliotail;
// the termination shock (the 4 -> 3 entry) and all other cells are unchanged, and so is every other
// character of the file. The copy is written next to the original, with `suffix` (default "_topo")
// before the extension: time_series/boundaries_040.txt -> time_series/boundaries_040_topo.txt; point the
// run's boundaries file names at the copies (e.g. [time_dependent.files] boundaries =
// "time_series/boundaries_{index:03}_topo.txt"). --snapshots runs over several snapshots of the series
// in one go; a copy that exists already is skipped (--overwrite: written again), so an interrupted job
// can simply be restarted. The Tecplot file is then written only if --out is given (with several
// snapshots, as OUT_NNN.dat). ASCII boundaries files of the standard model source only.

#include <algorithm>
#include <array>
#include <chrono>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "ena/ascii.hpp"
#include "ena/config.hpp"
#include "ena/parallel.hpp"
#include "ena/runtime.hpp"
#include "ena/snapshot.hpp"

using namespace ena;

namespace {

constexpr double kDeg = 3.14159265358979323846 / 180.0;

enum End : int { kSun = 1, kIsmZ = 2, kIsmRegion = 3, kIsmEdge = 4, kLong = 5, kNull = 6, kEdgeNotIsm = 7 };
const char* end_name(int e) {
  switch (e) {
    case kSun: return "Sun";
    case kIsmZ: return "ISM (z limit)";
    case kIsmRegion: return "ISM (region)";
    case kIsmEdge: return "ISM (grid edge)";
    case kLong: return "long";
    case kNull: return "null";
    case kEdgeNotIsm: return "grid edge (not ISM)";
  }
  return "?";
}
bool is_ism(int e) { return e == kIsmZ || e == kIsmRegion || e == kIsmEdge; }

struct Settings {
  double z_top = 500.0, z_bottom = -500.0, r_sun = 50.0, max_length = 10000.0, cell_fraction = 1.0;
  int substeps = 2;
  std::set<int> ism_regions{1};
  // which exits through the outer grid edge count as reaching the ISM: 0 all ("all"), 1 only where the
  // line leaves at x < 0, the nose side ("nose"), 2 none ("none"); the others end as kEdgeNotIsm, which
  // counts like a long line (a heliotail line can leave the grid down the tail without being open)
  int edge_ism = 0;
};

struct Tracer {
  const Grid& g;
  const Snapshot& s;
  const Settings& o;
  double rmin, rmax;

  Tracer(const Grid& grid, const Snapshot& snap, const Settings& set)
      : g(grid), s(snap), o(set), rmin(grid.r_inner), rmax(grid.r_inner + (grid.nr - 1) * grid.dr) {}

  // fractional grid indices of x; false outside the radial range
  bool indices(const double* x, double& fr, double& ft, double& fp, double& r) const {
    r = std::sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
    if (!(r >= rmin && r <= rmax)) return false;
    const double th = std::acos(std::max(-1.0, std::min(1.0, x[2] / r))) / kDeg;
    double ph = std::atan2(x[1], x[0]) / kDeg;
    if (ph < 0) ph += 360.0;
    fr = (r - g.r_inner) / g.dr;
    ft = std::min(th / g.dtheta, double(g.ntheta - 1));  // beyond the last row: its values
    fp = ph / g.dphi;
    return true;
  }

  // B/|B| at x (trilinear in index space); false outside the grid or in a null
  bool direction(const double* x, double* d) const {
    double fr, ft, fp, r;
    if (!indices(x, fr, ft, fp, r)) return false;
    const int ir = std::min(int(fr), g.nr - 2), it = std::min(int(ft), g.ntheta - 2);
    int ip = int(fp);
    const double wr = fr - ir, wt = ft - it, wp = fp - ip;
    ip = ((ip % g.nphi) + g.nphi) % g.nphi;
    const int ip1 = (ip + 1) % g.nphi;
    const Field3D<double>* comp[3] = {&s.Bx, &s.By, &s.Bz};
    double b[3];
    for (int c = 0; c < 3; ++c) {
      const Field3D<double>& F = *comp[c];
      auto lr = [&](int t, int p) { return (1.0 - wr) * F(t, p, ir) + wr * F(t, p, ir + 1); };
      const double a0 = (1.0 - wp) * lr(it, ip) + wp * lr(it, ip1);
      const double a1 = (1.0 - wp) * lr(it + 1, ip) + wp * lr(it + 1, ip1);
      b[c] = (1.0 - wt) * a0 + wt * a1;
    }
    const double m = std::sqrt(b[0] * b[0] + b[1] * b[1] + b[2] * b[2]);
    if (!(m > 0.0) || !std::isfinite(m)) return false;
    d[0] = b[0] / m; d[1] = b[1] / m; d[2] = b[2] / m;
    return true;
  }

  // region code at the nearest grid point of x
  int region_at(const double* x) const {
    double fr, ft, fp, r;
    if (!indices(x, fr, ft, fp, r)) return -1;
    const int ir = std::min(int(std::lround(fr)), g.nr - 1), it = std::min(int(std::lround(ft)), g.ntheta - 1);
    const int ip = int(std::lround(fp)) % g.nphi;
    return s.region(it, ip, ir);
  }

  // the end of a line that leaves through the outer grid edge at x
  int edge_end(const double* x) const {
    if (o.edge_ism == 0 || (o.edge_ism == 1 && x[0] < 0.0)) return kIsmEdge;
    return kEdgeNotIsm;
  }

  // the end condition met at x, or 0
  int stop_at(const double* x) const {
    const double r = std::sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
    if (r <= o.r_sun) return kSun;
    if (x[2] >= o.z_top || x[2] <= o.z_bottom) return kIsmZ;
    if (r > rmax) return edge_end(x);
    if (!o.ism_regions.empty() && o.ism_regions.count(region_at(x))) return kIsmRegion;
    return 0;
  }

  // follow the field line from x0 in direction sign (+1 along B, -1 against): end type and path length
  std::pair<int, double> follow(const double* x0, double sign) const {
    double x[3] = {x0[0], x0[1], x0[2]};
    if (int e = stop_at(x)) return {e, 0.0};
    double length = 0.0;
    while (length < o.max_length) {
      const double r = std::sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
      const double h = sign * o.cell_fraction * std::min(g.dr, r * g.dtheta * kDeg) / o.substeps;
      double y[3] = {x[0], x[1], x[2]};
      for (int k = 0; k < o.substeps; ++k) {
        double k1[3], k2[3], k3[3], k4[3], t[3];
        bool ok = direction(y, k1);
        if (ok) { for (int i = 0; i < 3; ++i) t[i] = y[i] + 0.5 * h * k1[i]; ok = direction(t, k2); }
        if (ok) { for (int i = 0; i < 3; ++i) t[i] = y[i] + 0.5 * h * k2[i]; ok = direction(t, k3); }
        if (ok) { for (int i = 0; i < 3; ++i) t[i] = y[i] + h * k3[i]; ok = direction(t, k4); }
        if (!ok) {
          // outside the grid: outward is the ISM, inward the Sun; otherwise a null of the field
          const double ry = std::sqrt(y[0] * y[0] + y[1] * y[1] + y[2] * y[2]);
          if (ry >= rmax - 2.0 * g.dr) return {edge_end(y), length};
          if (ry <= rmin + 2.0 * g.dr) return {kSun, length};
          return {kNull, length};
        }
        for (int i = 0; i < 3; ++i) y[i] += (h / 6.0) * (((k1[i] + 2.0 * k2[i]) + 2.0 * k3[i]) + k4[i]);
        length += std::fabs(h);
      }
      x[0] = y[0]; x[1] = y[1]; x[2] = y[2];
      if (int e = stop_at(x)) return {e, length};
    }
    return {kLong, length};
  }
};

int classify(int a, int b) {
  const bool ia = is_ism(a), ib = is_ism(b), sa = a == kSun, sb = b == kSun;
  if ((ia && sb) || (ib && sa)) return 2;
  if (ia && ib) return 3;
  if (ia || ib) return 4;
  if (a == kNull && b == kNull) return 0;
  return 1;
}

std::set<int> parse_ints(const std::string& s) {
  std::set<int> out;
  std::stringstream ss(s);
  std::string tok;
  while (std::getline(ss, tok, ','))
    if (!tok.empty()) out.insert(std::stoi(tok));
  return out;
}

// "file.txt" + "_topo" -> "file_topo.txt" (the suffix goes before the extension of the file name)
std::string with_suffix(const std::string& path, const std::string& suffix) {
  const std::size_t slash = path.find_last_of('/');
  const std::size_t dot = path.find_last_of('.');
  if (dot == std::string::npos || (slash != std::string::npos && dot < slash)) return path + suffix;
  return path.substr(0, dot) + suffix + path.substr(dot);
}

// the region token rewritten from code `from` to code `to`, in the token's own number format where
// possible ("3.000000E+00" -> "2.000000E+00"), otherwise as a plain integer
std::string region_token(std::string_view tok, int from, int to) {
  std::string t(tok);
  const std::size_t p = (!t.empty() && (t[0] == '+' || t[0] == '-')) ? 1 : 0;
  if (p < t.size() && t[p] == char('0' + from) && from >= 0 && from <= 9 && to >= 0 && to <= 9) {
    t[p] = char('0' + to);
    char* end = nullptr;
    const double v = std::strtod(t.c_str(), &end);
    if (end && *end == '\0' && int(std::trunc(v)) == to) return t;
  }
  return std::to_string(to);
}

// Copies the boundaries file `in` to `out`, replacing the region code of every row whose cell has
// new_code >= 0; every other character is kept. Returns the number of rows changed.
long rewrite_boundaries(const RunConfig& cfg, const std::string& in, const std::string& out,
                        const std::vector<int>& new_code) {
  const ColumnMap& m = cfg.inputs.boundaries;
  const Grid& g = cfg.grid;
  const int ncols = int(m.columns.size());
  const int creg = m.column("region"), ci = m.column("i"), cj = m.column("j"), ck = m.column("k");
  if (creg < 0 || ci < 0 || cj < 0 || ck < 0)
    throw std::runtime_error("the boundaries column map needs region, i, j and k columns");
  const std::string tmp = out + ".part";
  FILE* fo = std::fopen(tmp.c_str(), "wb");
  if (!fo) throw std::runtime_error("cannot write " + tmp);
  long changed = 0;
  try {
    for_each_parsed_line(in, ncols, [](std::string_view) { return true; },
                         [&](std::string_view text, const double* v, int n, long line) {
      std::string row(text);
      if (n != 0) {  // blank lines are copied as they are
        if (n != ncols)
          throw std::runtime_error(in + ":" + std::to_string(line) + ": expected " + std::to_string(ncols) +
                                   " numbers, found " + (n < 0 ? std::string("a non-numeric value") : std::to_string(n)));
        const int i = int(std::lround(v[ci])), j = int(std::lround(v[cj])), k = int(std::lround(v[ck]));
        if (i < 0 || i >= g.ntheta || j < 0 || j >= g.nphi || k < 0 || k >= g.nr)
          throw std::runtime_error(in + ":" + std::to_string(line) + ": cell index outside the grid");
        const int to = new_code[g.index(i, j, k)];
        const int from = int(std::trunc(v[creg]));
        if (to >= 0 && to != from) {
          // find the region column's token: whitespace-separated, as the reader takes it
          std::size_t pos = 0, a = 0, b = 0;
          for (int c = 0; c <= creg; ++c) {
            while (pos < row.size() && std::isspace((unsigned char)row[pos])) ++pos;
            a = pos;
            while (pos < row.size() && !std::isspace((unsigned char)row[pos])) ++pos;
            b = pos;
          }
          row = row.substr(0, a) + region_token(std::string_view(row).substr(a, b - a), from, to) + row.substr(b);
          ++changed;
        }
      }
      row += '\n';
      if (std::fwrite(row.data(), 1, row.size(), fo) != row.size())
        throw std::runtime_error("error writing " + tmp + " (disk full?)");
    });
  } catch (...) {
    std::fclose(fo);
    std::remove(tmp.c_str());
    throw;
  }
  const bool ok = std::ferror(fo) == 0;
  if (std::fclose(fo) != 0 || !ok) throw std::runtime_error("error writing " + tmp + " (disk full?)");
  if (std::rename(tmp.c_str(), out.c_str()) != 0) throw std::runtime_error("cannot rename " + tmp);
  return changed;
}

// --snapshots: "all", "A-B" or "A,B,C" (indices into the configuration's snapshot list, from 0)
std::vector<std::size_t> parse_snapshots(const std::string& v, std::size_t count) {
  std::vector<std::size_t> out;
  if (v == "all") {
    for (std::size_t s = 0; s < count; ++s) out.push_back(s);
    return out;
  }
  std::stringstream ss(v);
  std::string tok;
  while (std::getline(ss, tok, ',')) {
    if (tok.empty()) continue;
    const std::size_t dash = tok.find('-');
    const std::size_t a = std::stoul(tok.substr(0, dash));
    const std::size_t b = dash == std::string::npos ? a : std::stoul(tok.substr(dash + 1));
    for (std::size_t s = a; s <= b; ++s) out.push_back(s);
  }
  for (std::size_t s : out)
    if (s >= count) throw std::runtime_error("--snapshots: index " + std::to_string(s) + " out of range (0-" +
                                             std::to_string(count - 1) + ")");
  return out;
}

}  // namespace

int main(int argc, char** argv) {
  prepare_runtime();
  std::setvbuf(stdout, nullptr, _IOLBF, 0);  // the log of a job shows each line as it is printed
  try {
    if (argc < 2) {
      std::fprintf(stderr,
                   "usage: %s run.toml [--snapshot N | --snapshots all|A-B|A,B,...] [--out topology.dat] "
                   "[--z-top 500] [--z-bottom -500] [--ism-regions 1] [--r-sun 50] [--max-length 10000] "
                   "[--cell-fraction 1] [--substeps 2] [--only-regions 3,6] [--all-radii] "
                   "[--rewrite-boundaries] [--suffix _topo] [--overwrite] [--edge-ism all|nose|none] [--threads N]\n",
                   argv[0]);
      return 2;
    }
    const std::string cfg_path = argv[1];
    std::string out, snapshots_arg, suffix = "_topo";
    std::size_t snap_index = 0;
    bool all_radii = false, rewrite = false, overwrite = false, only_given = false;
    std::set<int> only;
    Settings o;
    for (int a = 2; a < argc; ++a) {
      const std::string k = argv[a];
      auto next = [&]() -> std::string {
        if (a + 1 >= argc) throw std::runtime_error(k + " needs a value");
        return argv[++a];
      };
      if (k == "--snapshot") snap_index = std::stoul(next());
      else if (k == "--snapshots") snapshots_arg = next();
      else if (k == "--out") out = next();
      else if (k == "--z-top") o.z_top = std::stod(next());
      else if (k == "--z-bottom") o.z_bottom = std::stod(next());
      else if (k == "--ism-regions") { const std::string v = next(); o.ism_regions = v == "none" ? std::set<int>{} : parse_ints(v); }
      else if (k == "--r-sun") o.r_sun = std::stod(next());
      else if (k == "--max-length") o.max_length = std::stod(next());
      else if (k == "--cell-fraction") o.cell_fraction = std::stod(next());
      else if (k == "--substeps") o.substeps = std::stoi(next());
      else if (k == "--only-regions") { only = parse_ints(next()); only_given = true; }
      else if (k == "--all-radii") all_radii = true;
      else if (k == "--rewrite-boundaries") rewrite = true;
      else if (k == "--suffix") suffix = next();
      else if (k == "--edge-ism") {
        const std::string v = next();
        if (v == "all") o.edge_ism = 0;
        else if (v == "nose") o.edge_ism = 1;
        else if (v == "none") o.edge_ism = 2;
        else throw std::runtime_error("--edge-ism: all, nose or none");
      }
      else if (k == "--overwrite") overwrite = true;
      else if (k == "--threads") set_threads(std::stoi(next()));
      else throw std::runtime_error("unknown argument " + k);
    }
    if (!(o.cell_fraction > 0) || o.substeps < 1 || !(o.max_length > 0) || !(o.z_top > o.z_bottom))
      throw std::runtime_error("need cell_fraction > 0, substeps >= 1, max_length > 0 and z_top > z_bottom");
    if (rewrite && suffix.empty()) throw std::runtime_error("--suffix must not be empty (the originals would be overwritten)");

    const auto t_start = std::chrono::steady_clock::now();
    RunConfig cfg = load_config(cfg_path);
    const Grid& g = cfg.grid;
    if (rewrite) {
      if (cfg.inputs.format != "ascii")
        throw std::runtime_error("--rewrite-boundaries works on ASCII boundaries files (inputs.format = \"ascii\"); "
                                 "convert the rewritten files to HDF5 afterwards if needed");
      if (cfg.moscow()) throw std::runtime_error("--rewrite-boundaries: the Moscow model source is not supported");
      all_radii = true;                    // every heliosheath cell, at every radius
      if (!only_given) only = {3, 6};      // only these cells can change
    }
    std::vector<std::size_t> indices;
    if (!snapshots_arg.empty()) indices = parse_snapshots(snapshots_arg, cfg.inputs.snapshots.size());
    else {
      if (snap_index >= cfg.inputs.snapshots.size()) throw std::runtime_error("--snapshot: index out of range");
      indices.push_back(snap_index);
    }
    if (!rewrite && out.empty()) out = "topology.dat";
    std::printf("ena_fieldtopo: %s, grid %d x %d x %d, %zu snapshot(s)\n", cfg_path.c_str(), g.ntheta, g.nphi, g.nr,
                indices.size());
    std::printf("  ISM: z >= %g or z <= %g AU, the outer edge (%g AU%s)", o.z_top, o.z_bottom, g.r(g.nr - 1),
                o.edge_ism == 0 ? "" : o.edge_ism == 1 ? ", where the line leaves at x < 0 only" : ": not counted");
    if (!o.ism_regions.empty()) {
      std::printf(", regions");
      for (int r : o.ism_regions) std::printf(" %d", r);
    }
    std::printf("; Sun: r <= %g AU; at most %g AU of path per direction; step %g cell, %d substeps; %d threads\n",
                o.r_sun, o.max_length, o.cell_fraction, o.substeps, max_threads());
    if (rewrite)
      std::printf("  rewriting the boundaries files: regions 3 and 6 -> 2 where the field line is not closed; "
                  "copies named with \"%s\" before the extension\n", suffix.c_str());

    const int k0 = all_radii ? 0 : std::max(0, g.k_of(cfg.r_lower));
    const int k1 = all_radii ? g.nr : std::min(g.k_of(cfg.r_upper), g.nr);
    const int nk = k1 - k0;
    const std::size_t n = std::size_t(g.ntheta) * g.nphi * nk;

    for (std::size_t snap_index : indices) {
      const auto t0 = std::chrono::steady_clock::now();
      const std::string bpath = cfg.inputs.snapshots[snap_index].boundaries.empty()
                                    ? std::string()
                                    : cfg.inputs.path(cfg.inputs.snapshots[snap_index].boundaries);
      const std::string rpath = rewrite ? with_suffix(bpath, suffix) : std::string();
      std::string tec = out;
      if (!tec.empty() && indices.size() > 1) {
        char tag[16];
        std::snprintf(tag, sizeof tag, "_%03zu", snap_index);
        tec = with_suffix(tec, tag);
      }
      if (rewrite && !overwrite && tec.empty()) {
        if (FILE* f = std::fopen(rpath.c_str(), "rb")) {
          std::fclose(f);
          std::printf("snapshot %zu: %s exists, skipped (--overwrite writes it again)\n", snap_index, rpath.c_str());
          continue;
        }
      }
      Snapshot snap;
      ReadOptions ro;
      ro.verbose = false;
      if (cfg.inputs.format == "hdf5") {
        snap = read_snapshot(cfg, snap_index, ro);
      } else {  // the boundaries file holds all that is needed: region codes and B
        read_boundaries_ascii(cfg, bpath, snap, ro);
      }
      if (snap.Bx.empty()) throw std::runtime_error("the snapshot has no magnetic field (Bx, By, Bz)");
      std::printf("snapshot %zu (%s)\n", snap_index, cfg.inputs.format == "hdf5"
                                                         ? cfg.inputs.snapshots[snap_index].file.c_str()
                                                         : cfg.inputs.snapshots[snap_index].boundaries.c_str());

      // ---- trace every grid point of the output radii, both directions ----
      std::vector<int> topo(n, -1), end_f(n, 0), end_b(n, 0);
      std::vector<double> len_f(n, 0.0), len_b(n, 0.0);
      Tracer tr(g, snap, o);
#pragma omp parallel for schedule(dynamic, 64)
      for (long c = 0; c < long(n); ++c) {
        const int i = int(c / (std::size_t(g.nphi) * nk)), j = int((c / nk) % g.nphi), k = k0 + int(c % nk);
        if (!only.empty() && !only.count(snap.region(i, j, k))) continue;
        const double th = g.theta(i) * kDeg, ph = g.phi(j) * kDeg, r = g.r(k);
        const double x[3] = {r * std::sin(th) * std::cos(ph), r * std::sin(th) * std::sin(ph), r * std::cos(th)};
        const auto f = tr.follow(x, +1.0), b = tr.follow(x, -1.0);
        end_f[c] = f.first; len_f[c] = f.second;
        end_b[c] = b.first; len_b[c] = b.second;
        topo[c] = classify(f.first, b.first);
      }

      // ---- summary ----
      long count[6] = {0, 0, 0, 0, 0, 0}, count3[6] = {0, 0, 0, 0, 0, 0};
      for (std::size_t c = 0; c < n; ++c) {
        const int t = topo[c] + 1;  // -1..4 -> 0..5
        ++count[t];
        const int i = int(c / (std::size_t(g.nphi) * nk)), j = int((c / nk) % g.nphi), k = k0 + int(c % nk);
        const int reg = snap.region(i, j, k);
        if (reg == 3 || reg == 6) ++count3[t];
      }
      const char* names[6] = {"not traced", "undetermined", "closed", "open", "interstellar", "ISM-connected"};
      std::printf("  grid points by topology (all / region 3 and 6):\n");
      for (int t = 0; t < 6; ++t)
        if (count[t]) std::printf("    %-14s %9ld / %9ld\n", names[t], count[t], count3[t]);

      // ---- boundaries file with the heliopause from the topology ----
      if (rewrite) {
        std::vector<int> new_code(g.size(), -1);
        long to2[7] = {0, 0, 0, 0, 0, 0, 0};
        for (std::size_t c = 0; c < n; ++c) {
          const int i = int(c / (std::size_t(g.nphi) * nk)), j = int((c / nk) % g.nphi), k = k0 + int(c % nk);
          const int reg = snap.region(i, j, k);
          if ((reg == 3 || reg == 6) && topo[c] >= 2) {  // open, interstellar or ISM-connected
            new_code[g.index(i, j, k)] = 2;
            ++to2[reg];
          }
        }
        const long changed = rewrite_boundaries(cfg, bpath, rpath, new_code);
        std::printf("  wrote %s: %ld cells set to region 2 (%ld of region 3, %ld of region 6)\n", rpath.c_str(),
                    changed, to2[3], to2[6]);
      }

      // ---- Tecplot file ----
      if (!tec.empty()) {
        const std::string tmp = tec + ".part";
        FILE* fo = std::fopen(tmp.c_str(), "w");
        if (!fo) throw std::runtime_error("cannot write " + tmp);
        std::fprintf(fo, "# ena_fieldtopo: magnetic topology of snapshot %zu of %s\n", snap_index, cfg_path.c_str());
        std::fprintf(fo, "# topology: 1 closed, 2 open (Sun and ISM), 3 interstellar, 4 ISM-connected (other end "
                         "undetermined), 0 undetermined, -1 not traced; ends: 1 Sun, 2 ISM (z limit), 3 ISM (region), "
                         "4 ISM (grid edge), 5 long, 6 null, 7 grid edge not counted as ISM (--edge-ism)\n");
        std::fprintf(fo, "# ISM: z >= %g or z <= %g AU, the grid edge%s", o.z_top, o.z_bottom,
                     o.edge_ism == 0 ? "" : o.edge_ism == 1 ? " (at x < 0 only)" : " (not counted)");
        for (int r : o.ism_regions) std::fprintf(fo, ", region %d", r);
        std::fprintf(fo, "; Sun: r <= %g AU; max length %g AU\n", o.r_sun, o.max_length);
        std::fprintf(fo, "TITLE = \"field-line topology, snapshot %zu\"\n", snap_index);
        std::fprintf(fo, "VARIABLES = \"X AU\" \"Y AU\" \"Z AU\" \"region\" \"topology\" \"end along B\" "
                         "\"end against B\" \"length along B AU\" \"length against B AU\" \"length to ISM AU\" "
                         "\"B_x\" \"B_y\" \"B_z\" \"|B|\" \"|B_phi|/|B|\"\n");
        std::fprintf(fo, "ZONE T=\"topology\", I=%d, J=%d, K=%d, DATAPACKING=POINT\n", nk, g.nphi, g.ntheta + 1);
        auto row = [&](int i, int j, int k, bool pole) {
          const std::size_t c = (std::size_t(i) * g.nphi + j) * nk + (k - k0);
          const double th = (pole ? 180.0 : g.theta(i)) * kDeg, ph = g.phi(j) * kDeg, r = g.r(k);
          const double bx = snap.Bx(i, j, k), by = snap.By(i, j, k), bz = snap.Bz(i, j, k);
          const double bm = std::sqrt(bx * bx + by * by + bz * bz);
          const double bphi = std::fabs(-bx * std::sin(ph) + by * std::cos(ph));
          double to_ism = -1.0;
          if (is_ism(end_f[c])) to_ism = len_f[c];
          if (is_ism(end_b[c]) && (to_ism < 0 || len_b[c] < to_ism)) to_ism = len_b[c];
          std::fprintf(fo, "%.6e %.6e %.6e %d %d %d %d %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e\n",
                       r * std::sin(th) * std::cos(ph), r * std::sin(th) * std::sin(ph), r * std::cos(th),
                       int(snap.region(i, j, k)), topo[c], end_f[c], end_b[c], len_f[c], len_b[c], to_ism, bx, by, bz,
                       bm, bm > 0 ? bphi / bm : 0.0);
        };
        for (int i = 0; i < g.ntheta; ++i)
          for (int j = 0; j < g.nphi; ++j)
            for (int k = k0; k < k1; ++k) row(i, j, k, false);
        for (int j = 0; j < g.nphi; ++j)
          for (int k = k0; k < k1; ++k) row(g.ntheta - 1, j, k, true);
        const bool ok = std::ferror(fo) == 0;
        if (std::fclose(fo) != 0 || !ok) throw std::runtime_error("error writing " + tmp + " (disk full?)");
        if (std::rename(tmp.c_str(), tec.c_str()) != 0) throw std::runtime_error("cannot rename " + tmp);
        std::printf("  wrote %s (%zu points, radii %g-%g AU)\n", tec.c_str(), n + std::size_t(g.nphi) * nk, g.r(k0),
                    g.r(k1 - 1));
      }
      std::printf("  %.1f s\n", std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
    }
    std::printf("done in %.1f s\n",
                std::chrono::duration<double>(std::chrono::steady_clock::now() - t_start).count());
    (void)end_name;
  } catch (const std::exception& e) {
    std::fprintf(stderr, "error: %s\n", e.what());
    return 1;
  }
  return 0;
}
