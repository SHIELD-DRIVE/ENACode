// ena_fieldtopo: magnetic field-line topology on the ENA grid, to locate the heliopause where the
// region codes cannot (e.g. the low-latitude heliotail).
//
//   ena_fieldtopo run.toml [--snapshot N] [--out topology.dat] [--z-top 500] [--z-bottom -500]
//                 [--ism-regions 1] [--r-sun 50] [--max-length 10000] [--cell-fraction 1]
//                 [--substeps 2] [--only-regions 3,6] [--all-radii] [--threads N]
//
// The field line through every grid point is followed in both directions, along B/|B|, through the
// magnetic field of the snapshot's boundaries file (Bx, By, Bz). Each direction ends where the line
//   - reaches the inner heliosphere, r <= r_sun ("Sun");
//   - reaches the interstellar medium ("ISM"): z >= z_top or z <= z_bottom, or a cell whose region
//     code is one of ism_regions, or the outer edge of the grid;
//   - or neither within max_length of path ("long"), or stops in a null of the field ("null").
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

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "ena/config.hpp"
#include "ena/parallel.hpp"
#include "ena/runtime.hpp"
#include "ena/snapshot.hpp"

using namespace ena;

namespace {

constexpr double kDeg = 3.14159265358979323846 / 180.0;

enum End : int { kSun = 1, kIsmZ = 2, kIsmRegion = 3, kIsmEdge = 4, kLong = 5, kNull = 6 };
const char* end_name(int e) {
  switch (e) {
    case kSun: return "Sun";
    case kIsmZ: return "ISM (z limit)";
    case kIsmRegion: return "ISM (region)";
    case kIsmEdge: return "ISM (grid edge)";
    case kLong: return "long";
    case kNull: return "null";
  }
  return "?";
}
bool is_ism(int e) { return e == kIsmZ || e == kIsmRegion || e == kIsmEdge; }

struct Settings {
  double z_top = 500.0, z_bottom = -500.0, r_sun = 50.0, max_length = 10000.0, cell_fraction = 1.0;
  int substeps = 2;
  std::set<int> ism_regions{1};
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

  // the end condition met at x, or 0
  int stop_at(const double* x) const {
    const double r = std::sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
    if (r <= o.r_sun) return kSun;
    if (x[2] >= o.z_top || x[2] <= o.z_bottom) return kIsmZ;
    if (r > rmax) return kIsmEdge;
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
          if (ry >= rmax - 2.0 * g.dr) return {kIsmEdge, length};
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

}  // namespace

int main(int argc, char** argv) {
  prepare_runtime();
  try {
    if (argc < 2) {
      std::fprintf(stderr,
                   "usage: %s run.toml [--snapshot N] [--out topology.dat] [--z-top 500] [--z-bottom -500] "
                   "[--ism-regions 1] [--r-sun 50] [--max-length 10000] [--cell-fraction 1] [--substeps 2] "
                   "[--only-regions 3,6] [--all-radii] [--threads N]\n",
                   argv[0]);
      return 2;
    }
    const std::string cfg_path = argv[1];
    std::string out = "topology.dat";
    std::size_t snap_index = 0;
    bool all_radii = false;
    std::set<int> only;
    Settings o;
    for (int a = 2; a < argc; ++a) {
      const std::string k = argv[a];
      auto next = [&]() -> std::string {
        if (a + 1 >= argc) throw std::runtime_error(k + " needs a value");
        return argv[++a];
      };
      if (k == "--snapshot") snap_index = std::stoul(next());
      else if (k == "--out") out = next();
      else if (k == "--z-top") o.z_top = std::stod(next());
      else if (k == "--z-bottom") o.z_bottom = std::stod(next());
      else if (k == "--ism-regions") { const std::string v = next(); o.ism_regions = v == "none" ? std::set<int>{} : parse_ints(v); }
      else if (k == "--r-sun") o.r_sun = std::stod(next());
      else if (k == "--max-length") o.max_length = std::stod(next());
      else if (k == "--cell-fraction") o.cell_fraction = std::stod(next());
      else if (k == "--substeps") o.substeps = std::stoi(next());
      else if (k == "--only-regions") only = parse_ints(next());
      else if (k == "--all-radii") all_radii = true;
      else if (k == "--threads") set_threads(std::stoi(next()));
      else throw std::runtime_error("unknown argument " + k);
    }
    if (!(o.cell_fraction > 0) || o.substeps < 1 || !(o.max_length > 0) || !(o.z_top > o.z_bottom))
      throw std::runtime_error("need cell_fraction > 0, substeps >= 1, max_length > 0 and z_top > z_bottom");

    const auto t0 = std::chrono::steady_clock::now();
    RunConfig cfg = load_config(cfg_path);
    const Grid& g = cfg.grid;
    if (snap_index >= cfg.inputs.snapshots.size()) throw std::runtime_error("--snapshot: index out of range");
    Snapshot snap;
    ReadOptions ro;
    ro.verbose = false;
    if (cfg.inputs.format == "hdf5") {
      snap = read_snapshot(cfg, snap_index, ro);
    } else {  // the boundaries file holds all that is needed: region codes and B
      read_boundaries_ascii(cfg, cfg.inputs.path(cfg.inputs.snapshots[snap_index].boundaries), snap, ro);
    }
    if (snap.Bx.empty()) throw std::runtime_error("the snapshot has no magnetic field (Bx, By, Bz)");
    std::printf("ena_fieldtopo: snapshot %zu of %s, grid %d x %d x %d\n", snap_index, cfg_path.c_str(), g.ntheta,
                g.nphi, g.nr);
    std::printf("  ISM: z >= %g or z <= %g AU, the outer edge (%g AU)", o.z_top, o.z_bottom, g.r(g.nr - 1));
    if (!o.ism_regions.empty()) {
      std::printf(", regions");
      for (int r : o.ism_regions) std::printf(" %d", r);
    }
    std::printf("; Sun: r <= %g AU; at most %g AU of path per direction; step %g cell, %d substeps; %d threads\n",
                o.r_sun, o.max_length, o.cell_fraction, o.substeps, max_threads());

    // ---- trace every grid point of the output radii, both directions ----
    const int k0 = all_radii ? 0 : std::max(0, g.k_of(cfg.r_lower));
    const int k1 = all_radii ? g.nr : std::min(g.k_of(cfg.r_upper), g.nr);
    const int nk = k1 - k0;
    const std::size_t n = std::size_t(g.ntheta) * g.nphi * nk;
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

    // ---- Tecplot file ----
    const std::string tmp = out + ".part";
    FILE* fo = std::fopen(tmp.c_str(), "w");
    if (!fo) throw std::runtime_error("cannot write " + tmp);
    std::fprintf(fo, "# ena_fieldtopo: magnetic topology of snapshot %zu of %s\n", snap_index, cfg_path.c_str());
    std::fprintf(fo, "# topology: 1 closed, 2 open (Sun and ISM), 3 interstellar, 4 ISM-connected (other end "
                     "undetermined), 0 undetermined, -1 not traced; ends: 1 Sun, 2 ISM (z limit), 3 ISM (region), "
                     "4 ISM (grid edge), 5 long, 6 null\n");
    std::fprintf(fo, "# ISM: z >= %g or z <= %g AU, the grid edge", o.z_top, o.z_bottom);
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
                   int(snap.region(i, j, k)), topo[c], end_f[c], end_b[c], len_f[c], len_b[c], to_ism, bx, by, bz, bm,
                   bm > 0 ? bphi / bm : 0.0);
    };
    for (int i = 0; i < g.ntheta; ++i)
      for (int j = 0; j < g.nphi; ++j)
        for (int k = k0; k < k1; ++k) row(i, j, k, false);
    for (int j = 0; j < g.nphi; ++j)
      for (int k = k0; k < k1; ++k) row(g.ntheta - 1, j, k, true);
    const bool ok = std::ferror(fo) == 0;
    if (std::fclose(fo) != 0 || !ok) throw std::runtime_error("error writing " + tmp + " (disk full?)");
    if (std::rename(tmp.c_str(), out.c_str()) != 0) throw std::runtime_error("cannot rename " + tmp);
    std::printf("  wrote %s (%zu points, radii %g-%g AU); %.1f s in total\n", out.c_str(),
                n + std::size_t(g.nphi) * nk, g.r(k0), g.r(k1 - 1),
                std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
    (void)end_name;
  } catch (const std::exception& e) {
    std::fprintf(stderr, "error: %s\n", e.what());
    return 1;
  }
  return 0;
}
