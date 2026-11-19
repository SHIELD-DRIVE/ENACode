// ena_fieldtopo_mhd: magnetic field-line topology on the original BATS-R-US grid, from the MHD
// solution itself (a Tecplot binary .plt file), written back on that same grid. The counterpart of
// ena_fieldtopo (which works on the ENA grid), for checking it without interpolation.
//
//   ena_fieldtopo_mhd MHD.plt OUT.plt [--seed-box XMIN XMAX YMIN YMAX ZMIN ZMAX] [--seed-stride N]
//                     [--points SEEDS.plt] [--no-mhd-vars] [--no-propagate] [--boundary-passes 4] [--z-top 500] [--z-bottom -500] [--hplim-ism ZMIN]
//                     [--r-sun 50] [--max-length 10000] [--cell-fraction 0.5] [--substeps 2]
//                     [--threads N]
//
// Input: the MHD solution as BATS-R-US writes it for Tecplot: one FEBRICK zone whose nodes are the
// cell centres and whose bricks join neighbouring centres (collapsed bricks where refinement levels
// meet). The order of the nodes does not matter. Variables are found by name: X, Y, Z (AU),
// B_x, B_y, B_z and, for --hplim-ism, Hplim.
//
// Tracing, the end conditions and the classes are those of ena_fieldtopo (see there), except that
// this grid has no region codes: a direction reaches the ISM at z >= z_top or z <= z_bottom, at
// the edge of the domain, or, with --hplim-ism ZMIN, where Hplim < 0 and |z| >= ZMIN (the
// heliopause variable above ZMIN, e.g. to use the top face of the heliopause); the Sun at
// r <= r_sun. The classes: 1 closed, 2 open, 3 interstellar, 4 ISM-connected, 0 undetermined,
// -1 not traced.
//
// Interpolation: a point is located in the brick that holds it (the bricks are indexed by size
// class and position). In a brick with 8 distinct corners B is interpolated trilinearly, in a
// collapsed brick by inverse-distance weights. A point in no brick is outside the domain. The step
// is cell_fraction of the size of the brick.
//
// Which nodes are classified: all nodes inside --seed-box (default: the whole domain). A field line
// has one topology along its whole length, so each traced line also labels the nodes it passes (the
// nearest node at every step; not nodes that meet an end condition themselves). Tracing goes in
// rounds: every 4096-th node, then every 512-th, 64-th and 8-th of the nodes still unlabelled, then
// all nodes still left. The seeds of a round are fixed before it starts, and a node passed by several
// lines of a round takes the line of the lowest seed, so the result does not depend on the thread
// count. Near a boundary between classes a passing line can be up to half a cell from the node, so the
// boundary would be sharp only to about a cell; so, afterwards, every node labelled by a passing line
// that shares a brick with a node of another class is traced on its own, then the nodes around those
// that changed, until none change (--boundary-passes N, default 4; 0: none). The boundaries are then
// as exact as tracing every node. --no-propagate traces every seed on its own and labels only it
// (with --seed-stride N, every N-th node of the box); the other nodes are -1 (not labelled).
//
// Output: OUT.plt on the original grid (the same nodes and bricks): all variables of the MHD file, in
// its order (--no-mhd-vars: only X, Y, Z), then the topology, the end of each direction (1 Sun,
// 2 ISM z limit, 3 ISM Hplim, 4 ISM domain edge, 5 long, 6 null), the path length of each direction
// and to the ISM (AU; for nodes traced themselves, else -1), and the label source (1 traced from the
// node in a round, 2 labelled by a line passing through it, 3 traced from the node in a boundary pass,
// 0 not labelled). With --points, also a small file of
// the nodes traced themselves (a scattered-point zone with the same variables).

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

using Clock = std::chrono::steady_clock;
double since(Clock::time_point t) { return std::chrono::duration<double>(Clock::now() - t).count(); }

// ---------------------------------------------------------------------------------------------
// Tecplot binary (.plt, version 112): reading
// ---------------------------------------------------------------------------------------------

struct PltFile {
  std::FILE* f = nullptr;
  std::string path;
  explicit PltFile(const std::string& p) : f(std::fopen(p.c_str(), "rb")), path(p) {
    if (!f) throw std::runtime_error("cannot open " + p);
    std::setvbuf(f, nullptr, _IOFBF, 1 << 24);
  }
  ~PltFile() { if (f) std::fclose(f); }
  void read(void* p, std::size_t n) {
    if (std::fread(p, 1, n, f) != n) throw std::runtime_error(path + ": unexpected end of file");
  }
  int32_t i32() { int32_t v; read(&v, 4); return v; }
  float f32() { float v; read(&v, 4); return v; }
  double f64() { double v; read(&v, 8); return v; }
  std::string str() {
    std::string s;
    for (;;) {
      const int32_t c = i32();
      if (c == 0) return s;
      s += (c > 0 && c < 128) ? char(c) : '?';
    }
  }
  void skip(long long n) {
    if (fseeko(f, off_t(n), SEEK_CUR) != 0) throw std::runtime_error(path + ": seek failed");
  }
};

struct Mesh {
  std::vector<std::string> names;
  long long npts = 0, nelem = 0;
  std::vector<float> x, y, z, bx, by, bz, hplim;  // per node
  std::vector<int32_t> conn;                      // 8 per brick, 0-based
  std::vector<std::vector<float>> other;          // the other variables (with keep_all), in file order
  std::vector<int> role;                          // per file variable: 0-6 = x, y, z, bx, by, bz, hplim; -1 other
  // the file's variable v, as stored here (nullptr if it was not kept)
  const std::vector<float>* var(int v) const {
    const std::vector<float>* named[7] = {&x, &y, &z, &bx, &by, &bz, &hplim};
    if (role[v] >= 0) return named[role[v]]->empty() ? nullptr : named[role[v]];
    for (std::size_t k = 0, idx = 0; k < role.size(); ++k)
      if (role[k] < 0) { if (int(k) == v) return idx < other.size() ? &other[idx] : nullptr; ++idx; }
    return nullptr;
  }
};

int find_var(const std::vector<std::string>& names, std::initializer_list<const char*> prefixes) {
  for (const char* p : prefixes)
    for (std::size_t i = 0; i < names.size(); ++i) {
      const std::string& n = names[i];
      const std::size_t l = std::strlen(p);
      if (n.compare(0, l, p) == 0 && (n.size() == l || n[l] == ' ')) return int(i);
    }
  return -1;
}

Mesh read_mesh(const std::string& path, bool want_hplim, bool keep_all) {
  PltFile in(path);
  char magic[8];
  in.read(magic, 8);
  if (std::memcmp(magic, "#!TDV", 5) != 0) throw std::runtime_error(path + ": not a classic Tecplot binary file");
  if (in.i32() != 1) throw std::runtime_error(path + ": big-endian files are not supported");
  in.i32();  // file type
  in.str();  // title
  Mesh m;
  const int nv = in.i32();
  for (int v = 0; v < nv; ++v) m.names.push_back(in.str());
  int nzones = 0;
  for (;;) {
    const float marker = in.f32();
    if (marker == 357.0f) break;
    if (marker == 799.0f) { in.str(); in.i32(); in.str(); continue; }          // data-set aux
    if (marker == 899.0f) { in.i32(); in.str(); in.i32(); in.str(); continue; }  // variable aux
    if (marker != 299.0f) throw std::runtime_error(path + ": unsupported header item " + std::to_string(marker));
    if (++nzones > 1) throw std::runtime_error(path + ": more than one zone (expected the single BATS-R-US zone)");
    in.str(); in.i32(); in.i32(); in.f64(); in.i32();
    const int type = in.i32();
    if (type != 5) throw std::runtime_error(path + ": the zone is not FEBRICK (type " + std::to_string(type) + ")");
    if (in.i32() == 1)
      for (int v = 0; v < nv; ++v)
        if (in.i32() != 0) throw std::runtime_error(path + ": cell-centred variables are not supported");
    in.i32();  // raw face neighbours
    if (in.i32() != 0) { in.i32(); in.i32(); }
    m.npts = in.i32();
    m.nelem = in.i32();
    in.i32(); in.i32(); in.i32();
    while (in.i32() == 1) { in.str(); in.i32(); in.str(); }
  }
  if (nzones != 1) throw std::runtime_error(path + ": no zone");
  if (in.f32() != 299.0f) throw std::runtime_error(path + ": data section not found");
  std::vector<int> fmt(nv);
  for (int v = 0; v < nv; ++v) fmt[v] = in.i32();
  if (in.i32() == 1) for (int v = 0; v < nv; ++v) if (in.i32()) throw std::runtime_error(path + ": passive variables");
  if (in.i32() == 1) for (int v = 0; v < nv; ++v) if (in.i32() >= 0) throw std::runtime_error(path + ": shared variables");
  if (in.i32() >= 0) throw std::runtime_error(path + ": shared connectivity");
  in.skip(16LL * nv);  // min, max
  const int cx = find_var(m.names, {"X"}), cy = find_var(m.names, {"Y"}), cz = find_var(m.names, {"Z"});
  const int cbx = find_var(m.names, {"B_x", "Bx"}), cby = find_var(m.names, {"B_y", "By"}),
            cbz = find_var(m.names, {"B_z", "Bz"}), ch = find_var(m.names, {"Hplim"});
  if (cx < 0 || cy < 0 || cz < 0 || cbx < 0 || cby < 0 || cbz < 0)
    throw std::runtime_error(path + ": X, Y, Z, B_x, B_y and B_z are needed");
  if (want_hplim && ch < 0) throw std::runtime_error(path + ": --hplim-ism needs the variable Hplim");
  const std::size_t n = std::size_t(m.npts);
  m.role.assign(nv, -1);
  m.role[cx] = 0; m.role[cy] = 1; m.role[cz] = 2; m.role[cbx] = 3; m.role[cby] = 4; m.role[cbz] = 5;
  if (ch >= 0) m.role[ch] = 6;
  for (int v = 0; v < nv; ++v) {
    std::vector<float>* dst = v == cx ? &m.x : v == cy ? &m.y : v == cz ? &m.z : v == cbx ? &m.bx
                            : v == cby ? &m.by : v == cbz ? &m.bz : ((want_hplim || keep_all) && v == ch) ? &m.hplim : nullptr;
    if (!dst && keep_all) { m.other.emplace_back(); dst = &m.other.back(); }
    const std::size_t w = fmt[v] == 1 ? 4 : fmt[v] == 2 ? 8 : 0;
    if (!w) throw std::runtime_error(path + ": variable data must be float or double");
    if (!dst) { in.skip((long long)(n * w)); continue; }
    dst->resize(n);
    if (w == 4) {
      in.read(dst->data(), n * 4);
    } else {
      std::vector<double> buf(1 << 20);
      for (std::size_t i = 0; i < n; i += buf.size()) {
        const std::size_t c = std::min(buf.size(), n - i);
        in.read(buf.data(), c * 8);
        for (std::size_t k = 0; k < c; ++k) (*dst)[i + k] = float(buf[k]);
      }
    }
  }
  m.conn.resize(std::size_t(m.nelem) * 8);
  in.read(m.conn.data(), m.conn.size() * 4);
  return m;
}

// ---------------------------------------------------------------------------------------------
// field and tracing
// ---------------------------------------------------------------------------------------------

enum End : int { kSun = 1, kIsmZ = 2, kIsmHplim = 3, kIsmEdge = 4, kLong = 5, kNull = 6 };
bool is_ism(int e) { return e == kIsmZ || e == kIsmHplim || e == kIsmEdge; }

struct Settings {
  double z_top = 500.0, z_bottom = -500.0, r_sun = 50.0, max_length = 10000.0, cell_fraction = 0.5;
  double hplim_zmin = -1.0;  // < 0: off
  int substeps = 2;
};

// The bricks, indexed by size class and position: a point finds the brick that holds it without any
// assumption on the order of the nodes. Each brick goes into the table of its size class (bucket
// size = the class's size), under every bucket its box overlaps; a query tests the bricks of its
// bucket in each class, finest first. Each thread first tries the brick of its previous query.
struct Field {
  const Mesh& m;
  double gx = 1e300, gy = 1e300, gz = 1e300, ex = -1e300, ey = -1e300, ez = -1e300;  // domain box
  double smin = 0;  // the smallest brick size
  struct Level { double size; std::vector<std::pair<uint64_t, int32_t>> keys; };
  std::vector<Level> levels;
  std::vector<bool> is_collapsed;  // per brick
  long long collapsed = 0, empty = 0;

  static uint64_t key(long long i, long long j, long long k) {
    return (uint64_t(i & 0x1FFFFF) << 42) | (uint64_t(j & 0x1FFFFF) << 21) | uint64_t(k & 0x1FFFFF);
  }

  void box(std::size_t e, double* lo, double* hi) const {
    const int32_t* nd = &m.conn[e * 8];
    lo[0] = lo[1] = lo[2] = 1e300; hi[0] = hi[1] = hi[2] = -1e300;
    for (int c = 0; c < 8; ++c) {
      const double q[3] = {m.x[nd[c]], m.y[nd[c]], m.z[nd[c]]};
      for (int d = 0; d < 3; ++d) { lo[d] = std::min(lo[d], q[d]); hi[d] = std::max(hi[d], q[d]); }
    }
  }

  explicit Field(const Mesh& mesh) : m(mesh) {
    const std::size_t ne = std::size_t(m.nelem);
    std::vector<float> size(ne);
    smin = 1e300;
    for (std::size_t e = 0; e < ne; ++e) {
      double lo[3], hi[3];
      box(e, lo, hi);
      size[e] = float(std::max(hi[0] - lo[0], std::max(hi[1] - lo[1], hi[2] - lo[2])));
      if (size[e] > 0) smin = std::min(smin, double(size[e]));
      gx = std::min(gx, lo[0]); gy = std::min(gy, lo[1]); gz = std::min(gz, lo[2]);
      ex = std::max(ex, hi[0]); ey = std::max(ey, hi[1]); ez = std::max(ez, hi[2]);
    }
    if (!(smin < 1e300)) throw std::runtime_error("all bricks are collapsed to points");
    // size class of a brick: the nearest power of two times the smallest size
    auto cls = [&](float s) { return int(std::lround(std::log2(double(s) / smin))); };
    int ncls = 0;
    for (std::size_t e = 0; e < ne; ++e) if (size[e] > 0) ncls = std::max(ncls, cls(size[e]) + 1);
    levels.resize(ncls);
    is_collapsed.assign(ne, false);
    for (int c = 0; c < ncls; ++c) levels[c].size = smin * std::ldexp(1.0, c);
    for (std::size_t e = 0; e < ne; ++e) {
      if (!(size[e] > 0)) { ++empty; continue; }
      Level& L = levels[cls(size[e])];
      double lo[3], hi[3];
      box(e, lo, hi);
      const double eps = 1e-6 * L.size;
      long long i0[3], i1[3];
      const double g[3] = {gx, gy, gz};
      for (int d = 0; d < 3; ++d) {
        i0[d] = (long long)std::floor((lo[d] - g[d] + eps) / L.size);
        i1[d] = (long long)std::floor((hi[d] - g[d] - eps) / L.size);
        if (i1[d] < i0[d]) i1[d] = i0[d];
      }
      bool distinct = true;  // a collapsed brick repeats nodes
      for (int c = 0; c < 8 && distinct; ++c)
        for (int c2 = 0; c2 < c; ++c2) distinct &= m.conn[e * 8 + c] != m.conn[e * 8 + c2];
      if (!distinct) { ++collapsed; is_collapsed[e] = true; }
      for (long long i = i0[0]; i <= i1[0]; ++i)
        for (long long j = i0[1]; j <= i1[1]; ++j)
          for (long long k = i0[2]; k <= i1[2]; ++k) L.keys.push_back({key(i, j, k), int32_t(e)});
    }
    for (auto& L : levels) {
      std::sort(L.keys.begin(), L.keys.end());
      L.keys.shrink_to_fit();
    }
  }

  bool holds(std::size_t e, const double* p, double* lo, double* hi) const {
    box(e, lo, hi);
    const double tol = 1e-7 * (hi[0] - lo[0] + hi[1] - lo[1] + hi[2] - lo[2] + 1e-30);
    return p[0] >= lo[0] - tol && p[0] <= hi[0] + tol && p[1] >= lo[1] - tol && p[1] <= hi[1] + tol &&
           p[2] >= lo[2] - tol && p[2] <= hi[2] + tol && hi[0] > lo[0] - 1e-30;
  }

  // the brick holding p, or -1; its box in lo, hi. The choice does not depend on earlier queries: the
  // first holding brick in a fixed order (finest size class first, then by brick number). The previous
  // brick of the thread is reused only where it is the only possible answer: p strictly inside a brick
  // that is not collapsed (regular bricks touch only at faces).
  long long find(const double* p, double* lo, double* hi) const {
    static thread_local long long last = -1;
    if (last >= 0 && !is_collapsed[std::size_t(last)]) {
      box(std::size_t(last), lo, hi);
      bool inside = true;
      for (int d = 0; d < 3; ++d) {
        const double margin = 1e-6 * (hi[d] - lo[d]);
        inside &= p[d] > lo[d] + margin && p[d] < hi[d] - margin;
      }
      if (inside) return last;
    }
    const double g[3] = {gx, gy, gz};
    for (const auto& L : levels) {
      if (L.keys.empty()) continue;
      long long ijk[3];
      bool out = false;
      for (int d = 0; d < 3; ++d) {
        ijk[d] = (long long)std::floor((p[d] - g[d]) / L.size);
        out |= ijk[d] < 0;
      }
      if (out) continue;
      const uint64_t k = key(ijk[0], ijk[1], ijk[2]);
      auto it = std::lower_bound(L.keys.begin(), L.keys.end(), std::make_pair(k, int32_t(-1)));
      for (; it != L.keys.end() && it->first == k; ++it)
        if (holds(std::size_t(it->second), p, lo, hi)) { last = it->second; return last; }
    }
    return -1;
  }

  // B (and Hplim) at p; false outside every brick. h: the size of the brick.
  bool sample(const double* p, double* out, double* hp = nullptr, double* h = nullptr) const {
    double lo[3], hi[3];
    const long long e = find(p, lo, hi);
    if (e < 0) return false;
    const int32_t* nd = &m.conn[std::size_t(e) * 8];
    if (h) *h = std::max(hi[0] - lo[0], std::max(hi[1] - lo[1], hi[2] - lo[2]));
    int corner_of[8];
    bool used[8] = {false, false, false, false, false, false, false, false}, full = true;
    for (int c = 0; c < 8; ++c) {
      int cc = 0;
      const double q[3] = {m.x[nd[c]], m.y[nd[c]], m.z[nd[c]]};
      for (int d = 0; d < 3; ++d) if (q[d] > 0.5 * (lo[d] + hi[d])) cc |= 1 << d;
      if (used[cc]) full = false;
      used[cc] = true;
      corner_of[c] = cc;
    }
    double v[4] = {0, 0, 0, 0}, wsum = 0;
    for (int c = 0; c < 8; ++c) {
      const int32_t n = nd[c];
      double w;
      if (full) {  // trilinear in the brick's box
        w = 1.0;
        for (int d = 0; d < 3; ++d) {
          const double span = hi[d] - lo[d], t = span > 0 ? (p[d] - lo[d]) / span : 0.5;
          w *= (corner_of[c] >> d & 1) ? t : 1.0 - t;
        }
      } else {     // collapsed brick: inverse-distance weights
        const double dx = p[0] - m.x[n], dy = p[1] - m.y[n], dz = p[2] - m.z[n];
        w = 1.0 / (dx * dx + dy * dy + dz * dz + 1e-12);
      }
      v[0] += w * m.bx[n]; v[1] += w * m.by[n]; v[2] += w * m.bz[n];
      if (hp) v[3] += w * m.hplim[n];
      wsum += w;
    }
    if (!(wsum > 0)) return false;
    out[0] = v[0] / wsum; out[1] = v[1] / wsum; out[2] = v[2] / wsum;
    if (hp) *hp = v[3] / wsum;
    return true;
  }

  double cell_size(const double* p) const {
    double b[3], h = 0;
    return sample(p, b, nullptr, &h) ? h : 0.0;
  }

  // the node nearest to p among the corners of the brick holding p, or -1
  long long nearest_node(const double* p) const {
    double lo[3], hi[3];
    const long long e = find(p, lo, hi);
    if (e < 0) return -1;
    const int32_t* nd = &m.conn[std::size_t(e) * 8];
    long long best = nd[0];
    double bd = 1e300;
    for (int c = 0; c < 8; ++c) {
      const double dx = p[0] - m.x[nd[c]], dy = p[1] - m.y[nd[c]], dz = p[2] - m.z[nd[c]];
      const double d = dx * dx + dy * dy + dz * dz;
      if (d < bd || (d == bd && nd[c] < best)) { bd = d; best = nd[c]; }
    }
    return best;
  }
};

struct Tracer {
  const Field& F;
  const Settings& o;

  bool direction(const double* x, double* d) const {
    double b[3];
    if (!F.sample(x, b)) return false;
    const double s = std::sqrt(b[0] * b[0] + b[1] * b[1] + b[2] * b[2]);
    if (!(s > 0) || !std::isfinite(s)) return false;
    d[0] = b[0] / s; d[1] = b[1] / s; d[2] = b[2] / s;
    return true;
  }

  int stop_at(const double* x) const {
    const double r = std::sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
    if (r <= o.r_sun) return kSun;
    if (x[2] >= o.z_top || x[2] <= o.z_bottom) return kIsmZ;
    if (o.hplim_zmin >= 0 && std::fabs(x[2]) >= o.hplim_zmin) {
      double b[3], hp = 0;
      if (F.sample(x, b, &hp) && hp < 0) return kIsmHplim;
    }
    return 0;
  }

  // follow the line from x0 (sign +1 along B, -1 against); with `visited`, the nearest node of every
  // step is recorded (consecutive repeats once)
  std::pair<int, double> follow(const double* x0, double sign, std::vector<int32_t>* visited = nullptr) const {
    double x[3] = {x0[0], x0[1], x0[2]};
    auto mark = [&](const double* p) {
      if (!visited) return;
      const long long q = F.nearest_node(p);
      if (q < 0 || (!visited->empty() && visited->back() == q)) return;
      // a node that meets an end condition itself (e.g. beyond the z limits) is not labelled by the line:
      // traced itself, it ends at once
      const double nq[3] = {F.m.x[q], F.m.y[q], F.m.z[q]};
      if (stop_at(nq)) return;
      visited->push_back(int32_t(q));
    };
    if (int e = stop_at(x)) return {e, 0.0};
    double length = 0.0;
    while (length < o.max_length) {
      const double h0 = F.cell_size(x);
      if (!(h0 > 0)) return {kIsmEdge, length};
      const double h = sign * o.cell_fraction * h0 / o.substeps;
      double y[3] = {x[0], x[1], x[2]};
      for (int s = 0; s < o.substeps; ++s) {
        double k1[3], k2[3], k3[3], k4[3], t[3];
        bool ok = direction(y, k1);
        if (ok) { for (int i = 0; i < 3; ++i) t[i] = y[i] + 0.5 * h * k1[i]; ok = direction(t, k2); }
        if (ok) { for (int i = 0; i < 3; ++i) t[i] = y[i] + 0.5 * h * k2[i]; ok = direction(t, k3); }
        if (ok) { for (int i = 0; i < 3; ++i) t[i] = y[i] + h * k3[i]; ok = direction(t, k4); }
        if (!ok) {
          double probe[3] = {y[0] + h * k1[0], y[1] + h * k1[1], y[2] + h * k1[2]};
          if (F.cell_size(y) <= 0 || F.cell_size(probe) <= 0) {
            const double r = std::sqrt(y[0] * y[0] + y[1] * y[1] + y[2] * y[2]);
            return {r <= 2.0 * o.r_sun ? kSun : kIsmEdge, length};  // out of the domain: inner body or outer edge
          }
          return {kNull, length};
        }
        for (int i = 0; i < 3; ++i) y[i] += (h / 6.0) * (((k1[i] + 2.0 * k2[i]) + 2.0 * k3[i]) + k4[i]);
        length += std::fabs(h);
      }
      x[0] = y[0]; x[1] = y[1]; x[2] = y[2];
      mark(x);
      if (F.cell_size(x) <= 0) {
        const double r = std::sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
        return {r <= 2.0 * o.r_sun ? kSun : kIsmEdge, length};
      }
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

// ---------------------------------------------------------------------------------------------
// writing
// ---------------------------------------------------------------------------------------------

class PltWriter {
 public:
  explicit PltWriter(const std::string& path) : f_(std::fopen(path.c_str(), "wb")) {
    if (!f_) throw std::runtime_error("cannot write " + path);
    std::setvbuf(f_, nullptr, _IOFBF, 1 << 24);
  }
  ~PltWriter() { if (f_) std::fclose(f_); }
  void bytes(const void* p, std::size_t n) {
    if (n && std::fwrite(p, 1, n, f_) != n) throw std::runtime_error("write failed (disk full?)");
  }
  void i32(int32_t v) { bytes(&v, 4); }
  void f32(float v) { bytes(&v, 4); }
  void f64(double v) { bytes(&v, 8); }
  void str(const std::string& s) { for (unsigned char c : s) i32(int32_t(c)); i32(0); }
  void close() {
    if (f_ && std::fclose(f_) != 0) { f_ = nullptr; throw std::runtime_error("write failed (disk full?)"); }
    f_ = nullptr;
  }
 private:
  std::FILE* f_;
};

// one zone: FEBRICK with conn, or (conn == nullptr) an ordered zone of scattered points
void write_plt(const std::string& path, const std::string& title, const std::vector<std::string>& names,
               const std::vector<const std::vector<float>*>& cols, const std::vector<int32_t>* conn) {
  const std::size_t n = cols[0]->size();
  PltWriter w(path + ".part");
  w.bytes("#!TDV112", 8);
  w.i32(1); w.i32(0); w.str(title);
  w.i32(int32_t(names.size()));
  for (auto& s : names) w.str(s);
  w.f32(299.0f); w.str("topology"); w.i32(-1); w.i32(-1); w.f64(0.0); w.i32(-1);
  if (conn) {
    w.i32(5); w.i32(0); w.i32(0); w.i32(0);
    w.i32(int32_t(n)); w.i32(int32_t(conn->size() / 8)); w.i32(0); w.i32(0); w.i32(0);
  } else {
    w.i32(0); w.i32(0); w.i32(0); w.i32(0);
    w.i32(int32_t(n)); w.i32(1); w.i32(1);
  }
  w.i32(0);
  w.f32(357.0f);
  w.f32(299.0f);
  for (std::size_t v = 0; v < names.size(); ++v) w.i32(1);
  w.i32(0); w.i32(0); w.i32(-1);
  for (const auto* c : cols) {
    float lo = 0, hi = 0;
    if (!c->empty()) { auto mm = std::minmax_element(c->begin(), c->end()); lo = *mm.first; hi = *mm.second; }
    w.f64(lo); w.f64(hi);
  }
  for (const auto* c : cols) w.bytes(c->data(), c->size() * 4);
  if (conn) w.bytes(conn->data(), conn->size() * 4);
  w.close();
  if (std::rename((path + ".part").c_str(), path.c_str()) != 0) throw std::runtime_error("cannot rename " + path + ".part");
}

}  // namespace

int main(int argc, char** argv) {
  std::setvbuf(stdout, nullptr, _IOLBF, 0);
  try {
    std::string in, out, points;
    bool keep_vars = true, propagate = true;
    int refine = 4;  // passes over the boundaries between classes
    Settings o;
    double box[6] = {-1e300, 1e300, -1e300, 1e300, -1e300, 1e300};
    long long stride = 1;
    for (int a = 1; a < argc; ++a) {
      const std::string k = argv[a];
      auto num = [&]() -> double {
        if (a + 1 >= argc) throw std::runtime_error(k + " needs a value");
        return std::stod(argv[++a]);
      };
      if (k == "--seed-box") for (double& v : box) v = num();
      else if (k == "--seed-stride") stride = (long long)num();
      else if (k == "--points") { if (a + 1 >= argc) throw std::runtime_error("--points needs a file"); points = argv[++a]; }
      else if (k == "--no-mhd-vars") keep_vars = false;
      else if (k == "--no-propagate") propagate = false;
      else if (k == "--boundary-passes") refine = int(num());
      else if (k == "--z-top") o.z_top = num();
      else if (k == "--z-bottom") o.z_bottom = num();
      else if (k == "--hplim-ism") o.hplim_zmin = num();
      else if (k == "--r-sun") o.r_sun = num();
      else if (k == "--max-length") o.max_length = num();
      else if (k == "--cell-fraction") o.cell_fraction = num();
      else if (k == "--substeps") o.substeps = int(num());
      else if (k == "--threads") {
#ifdef _OPENMP
        omp_set_num_threads(int(num()));
#else
        num();
#endif
      } else if (k == "-h" || k == "--help") {
        std::printf("usage: ena_fieldtopo_mhd MHD.plt OUT.plt [options]; see the top of apps/ena_fieldtopo_mhd.cpp\n");
        return 0;
      } else if (in.empty()) in = k;
      else if (out.empty()) out = k;
      else throw std::runtime_error("unknown argument " + k);
    }
    if (in.empty() || out.empty())
      throw std::runtime_error("usage: ena_fieldtopo_mhd MHD.plt OUT.plt [--seed-box X0 X1 Y0 Y1 Z0 Z1] [--seed-stride N] "
                               "[--points SEEDS.plt] [--no-mhd-vars] [--no-propagate] [--boundary-passes 4] [--z-top 500] [--z-bottom -500] [--hplim-ism ZMIN] [--r-sun 50] "
                               "[--max-length 10000] [--cell-fraction 0.5] [--substeps 2] [--threads N]");
    if (stride < 1 || !(o.cell_fraction > 0) || o.substeps < 1 || !(o.z_top > o.z_bottom))
      throw std::runtime_error("need --seed-stride >= 1, cell_fraction > 0, substeps >= 1, z_top > z_bottom");
    const auto t0 = Clock::now();

    Mesh m = read_mesh(in, o.hplim_zmin >= 0, keep_vars);
    std::printf("%s: %lld nodes, %lld bricks; read in %.0f s\n", in.c_str(), m.npts, m.nelem, since(t0));
    auto t1 = Clock::now();
    Field F(m);
    std::printf("  bricks indexed in %.0f s: %zu size classes (%.4g to %.4g AU), %lld collapsed, %lld empty; domain "
                "x %.0f..%.0f, y %.0f..%.0f, z %.0f..%.0f AU\n",
                since(t1), F.levels.size(), F.levels.front().size, F.levels.back().size, F.collapsed, F.empty, F.gx,
                F.ex, F.gy, F.ey, F.gz, F.ez);
    auto in_box = [&](long long q) {
      return m.x[q] >= box[0] && m.x[q] <= box[1] && m.y[q] >= box[2] && m.y[q] <= box[3] && m.z[q] >= box[4] &&
             m.z[q] <= box[5];
    };
    long long inbox = 0;
    for (long long q = 0; q < m.npts; ++q) inbox += in_box(q);
    std::printf("  %lld nodes in the seed box; %s; ISM: z >= %g or z <= %g AU, the domain edge%s; Sun: r <= %g AU; "
                "at most %g AU per direction, step %g cell, %d substeps\n",
                inbox, propagate ? "each traced line labels the nodes it passes" : "each seed labels only itself",
                o.z_top, o.z_bottom,
                o.hplim_zmin >= 0 ? (", Hplim < 0 at |z| >= " + std::to_string(o.hplim_zmin) + " AU").c_str() : "",
                o.r_sun, o.max_length, o.cell_fraction, o.substeps);

    // per node: topology (-1 not labelled), the ends of its line, how it was labelled (1 traced from the
    // node, 2 by a line passing through it); path lengths only for nodes traced themselves
    const std::size_t n = std::size_t(m.npts);
    std::vector<int8_t> label(n, -1), endf(n, 0), endb(n, 0), source(n, 0);
    std::vector<float> LF(n, -1.0f), LB(n, -1.0f), LI(n, -1.0f);
    std::vector<uint32_t> claim(propagate ? n : 0, UINT32_MAX);
    Tracer tr{F, o};
    // rounds: each traces the not yet labelled nodes of the box whose index is a multiple of its stride;
    // the seeds of a round are fixed before it starts, and a node passed by several lines takes the line of
    // the lowest seed, so the result does not depend on the thread count
    std::vector<long long> strides;
    if (propagate) strides = {4096, 512, 64, 8, 1};
    else strides = {stride};
    t1 = Clock::now();
    for (long long st : strides) {
      std::vector<long long> seeds;
      if (propagate) {
        for (long long q = 0; q < m.npts; q += st)
          if (label[q] < 0 && in_box(q)) seeds.push_back(q);
      } else {
        long long kk = 0;
        for (long long q = 0; q < m.npts; ++q)
          if (in_box(q) && kk++ % st == 0) seeds.push_back(q);
      }
      const std::size_t ns = seeds.size();
      if (!ns) continue;
      std::vector<int8_t> topo(ns), ef(ns), eb(ns);
      std::vector<float> lf(ns), lb(ns);
      const auto tr0 = Clock::now();
      long long done = 0;
#pragma omp parallel
      {
        std::vector<int32_t> visited;
#pragma omp for schedule(dynamic, 16)
        for (long long sidx = 0; sidx < (long long)ns; ++sidx) {
          const long long q = seeds[sidx];
          const double x[3] = {m.x[q], m.y[q], m.z[q]};
          visited.clear();
          const auto f = tr.follow(x, +1.0, propagate ? &visited : nullptr);
          const auto b = tr.follow(x, -1.0, propagate ? &visited : nullptr);
          ef[sidx] = int8_t(f.first); lf[sidx] = float(f.second); eb[sidx] = int8_t(b.first); lb[sidx] = float(b.second);
          topo[sidx] = int8_t(classify(f.first, b.first));
          for (int32_t v : visited) {  // claim the passed nodes: the lowest seed wins
            uint32_t cur = __atomic_load_n(&claim[std::size_t(v)], __ATOMIC_RELAXED);
            while (uint32_t(sidx) < cur &&
                   !__atomic_compare_exchange_n(&claim[std::size_t(v)], &cur, uint32_t(sidx), true, __ATOMIC_RELAXED,
                                                __ATOMIC_RELAXED)) {
            }
          }
          long long d;
#pragma omp atomic capture
          d = ++done;
          if (d % 100000 == 0) std::printf("    %lld of %zu seeds of this round traced (%.0f s)\n", d, ns, since(tr0));
        }
      }
      // the seeds themselves, then the nodes their lines passed
      for (std::size_t sidx = 0; sidx < ns; ++sidx) {
        const std::size_t q = std::size_t(seeds[sidx]);
        label[q] = topo[sidx]; endf[q] = ef[sidx]; endb[q] = eb[sidx]; source[q] = 1;
        LF[q] = lf[sidx]; LB[q] = lb[sidx];
        double li = -1;
        if (is_ism(ef[sidx])) li = lf[sidx];
        if (is_ism(eb[sidx]) && (li < 0 || lb[sidx] < li)) li = lb[sidx];
        LI[q] = float(li);
      }
      long long passed = 0;
      if (propagate) {
#pragma omp parallel for reduction(+ : passed) schedule(static)
        for (long long q = 0; q < (long long)n; ++q) {
          const uint32_t c = claim[std::size_t(q)];
          if (c == UINT32_MAX) continue;
          claim[std::size_t(q)] = UINT32_MAX;
          if (label[q] >= 0) continue;
          label[q] = topo[c]; endf[q] = ef[c]; endb[q] = eb[c]; source[q] = 2;
          ++passed;
        }
      }
      long long labelled = 0;
      for (std::size_t q = 0; q < n; ++q) labelled += label[q] >= 0;
      std::printf("  round (every %lld-th node): %zu seeds traced in %.0f s, %lld more nodes labelled by their lines; "
                  "%lld of %lld nodes labelled (%.1f%%)\n",
                  st, ns, since(tr0), passed, labelled, (long long)n, 100.0 * labelled / n);
    }
    // boundary passes: a node labelled by a passing line that shares a brick with a node of another class is
    // traced on its own; then the same for the bricks around the nodes that changed, until none change (at
    // most `refine` passes). Away from the boundaries every line agrees, so the labels there stay.
    if (propagate && refine > 0) {
      const std::size_t ne = std::size_t(m.nelem);
      std::vector<uint8_t> changed(n, 1), cand(n, 0);  // pass 1: every brick is checked
      for (int pass = 1; pass <= refine; ++pass) {
        std::fill(cand.begin(), cand.end(), 0);
#pragma omp parallel for schedule(static)
        for (long long e = 0; e < (long long)ne; ++e) {
          const int32_t* nd = &m.conn[std::size_t(e) * 8];
          bool touched = false, mixed = false;
          for (int c = 0; c < 8; ++c) {
            touched |= changed[nd[c]] != 0;
            mixed |= label[nd[c]] != label[nd[0]];
          }
          if (!touched || !mixed) continue;
          for (int c = 0; c < 8; ++c)
            if (source[nd[c]] == 2) cand[nd[c]] = 1;  // benign race: every writer stores 1
        }
        std::vector<long long> seeds;
        for (std::size_t q = 0; q < n; ++q) if (cand[q] && in_box((long long)q)) seeds.push_back((long long)q);
        if (seeds.empty()) break;
        const auto tp0 = Clock::now();
        std::vector<int8_t> topo(seeds.size()), ef(seeds.size()), eb(seeds.size());
        std::vector<float> lf(seeds.size()), lb(seeds.size());
#pragma omp parallel for schedule(dynamic, 16)
        for (long long sidx = 0; sidx < (long long)seeds.size(); ++sidx) {
          const long long q = seeds[sidx];
          const double x[3] = {m.x[q], m.y[q], m.z[q]};
          const auto f = tr.follow(x, +1.0), b = tr.follow(x, -1.0);
          ef[sidx] = int8_t(f.first); lf[sidx] = float(f.second); eb[sidx] = int8_t(b.first); lb[sidx] = float(b.second);
          topo[sidx] = int8_t(classify(f.first, b.first));
        }
        std::fill(changed.begin(), changed.end(), 0);
        long long nchanged = 0;
        for (std::size_t sidx = 0; sidx < seeds.size(); ++sidx) {
          const std::size_t q = std::size_t(seeds[sidx]);
          if (label[q] != topo[sidx]) { changed[q] = 1; ++nchanged; }
          label[q] = topo[sidx]; endf[q] = ef[sidx]; endb[q] = eb[sidx]; source[q] = 3;
          LF[q] = lf[sidx]; LB[q] = lb[sidx];
          double li = -1;
          if (is_ism(ef[sidx])) li = lf[sidx];
          if (is_ism(eb[sidx]) && (li < 0 || lb[sidx] < li)) li = lb[sidx];
          LI[q] = float(li);
        }
        std::printf("  boundary pass %d: %zu nodes next to another class traced on their own in %.0f s; %lld changed\n",
                    pass, seeds.size(), since(tp0), nchanged);
        if (nchanged == 0) break;
      }
    }
    long long count[6] = {0, 0, 0, 0, 0, 0}, own = 0;
    for (std::size_t q = 0; q < n; ++q) { ++count[label[q] + 1]; own += source[q] == 1 || source[q] == 3; }
    const char* cname[6] = {"not labelled", "undetermined", "closed", "open", "interstellar", "ISM-connected"};
    std::printf("  done in %.0f s (%lld nodes traced themselves); nodes by topology:\n", since(t1), own);
    for (int t = 0; t < 6; ++t)
      if (count[t]) std::printf("    %-14s %12lld\n", cname[t], count[t]);

    // ---- output on the original grid ----
    std::vector<float> T(n), EF(n), EB(n), SRC(n);
    for (std::size_t q = 0; q < n; ++q) { T[q] = label[q]; EF[q] = endf[q]; EB[q] = endb[q]; SRC[q] = source[q]; }
    std::vector<long long> seeds;  // the nodes traced themselves, for --points
    for (std::size_t q = 0; q < n; ++q) if (source[q] == 1 || source[q] == 3) seeds.push_back((long long)q);
    const std::size_t ns = seeds.size();
    // the variables written: the MHD file's (all, in its order; with --no-mhd-vars only X, Y, Z), then the topology
    std::vector<std::string> names;
    std::vector<const std::vector<float>*> cols;
    for (int v = 0; v < int(m.names.size()); ++v)
      if (const std::vector<float>* c = m.var(v))
        if (keep_vars || m.role[v] <= 2) { names.push_back(m.names[v]); cols.push_back(c); }
    const std::vector<std::string> tnames = {"topology", "end along B", "end against B", "length along B AU",
                                             "length against B AU", "length to ISM AU", "label source"};
    const std::vector<const std::vector<float>*> tcols = {&T, &EF, &EB, &LF, &LB, &LI, &SRC};
    names.insert(names.end(), tnames.begin(), tnames.end());
    cols.insert(cols.end(), tcols.begin(), tcols.end());
    t1 = Clock::now();
    write_plt(out, "field-line topology of " + in, names, cols, &m.conn);
    std::printf("  wrote %s (the original grid, %zu variables) in %.0f s\n", out.c_str(), names.size(), since(t1));
    if (!points.empty()) {
      std::vector<std::vector<float>> pv(cols.size(), std::vector<float>(ns));
      for (std::size_t s = 0; s < ns; ++s) {
        const std::size_t q = std::size_t(seeds[s]);
        for (std::size_t c = 0; c < cols.size(); ++c) pv[c][s] = (*cols[c])[q];
      }
      std::vector<const std::vector<float>*> pcols;
      for (auto& v : pv) pcols.push_back(&v);
      write_plt(points, "field-line topology of the seeds of " + in, names, pcols, nullptr);
      std::printf("  wrote %s (the %zu seeds as scattered points)\n", points.c_str(), ns);
    }
    std::printf("  %.0f s in total\n", since(t0));
  } catch (const std::exception& e) {
    std::fprintf(stderr, "error: %s\n", e.what());
    return 1;
  }
  return 0;
}
