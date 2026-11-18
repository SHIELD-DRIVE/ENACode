// ena_fieldtopo_mhd: magnetic field-line topology on the original BATS-R-US grid, from the MHD
// solution itself (a Tecplot binary .plt file), written back on that same grid. The counterpart of
// ena_fieldtopo (which works on the ENA grid), for checking it without interpolation.
//
//   ena_fieldtopo_mhd MHD.plt OUT.plt [--seed-box XMIN XMAX YMIN YMAX ZMIN ZMAX] [--seed-stride N]
//                     [--points SEEDS.plt] [--z-top 500] [--z-bottom -500] [--hplim-ism ZMIN]
//                     [--r-sun 50] [--max-length 10000] [--cell-fraction 0.5] [--substeps 2]
//                     [--threads N]
//
// Input: the MHD solution as BATS-R-US writes it for Tecplot: one FEBRICK zone whose nodes are the
// cell centres, block after block (8 x 8 x 8 cells per block, as the file's BLOCKS aux data says),
// and whose bricks join neighbouring centres, also across blocks. Variables are found by name:
// X, Y, Z (AU), B_x, B_y, B_z and, for --hplim-ism, Hplim. The layout is checked when read.
//
// Tracing, the end conditions and the classes are those of ena_fieldtopo (see there), except that
// this grid has no region codes: a direction reaches the ISM at z >= z_top or z <= z_bottom, at
// the edge of the domain, or, with --hplim-ism ZMIN, where Hplim < 0 and |z| >= ZMIN (the
// heliopause variable above ZMIN, e.g. to use the top face of the heliopause); the Sun at
// r <= r_sun. The classes: 1 closed, 2 open, 3 interstellar, 4 ISM-connected, 0 undetermined,
// -1 not traced.
//
// Interpolation: a point is located in its block through one lookup table per refinement level.
// Inside the block's own lattice of cell centres B is interpolated trilinearly; between a block's
// outer centres and its neighbours' (the outer half cell of a block), the file's bricks are used:
// trilinearly in a brick whose 8 corners are distinct, by inverse-distance weighting in a
// collapsed brick (where refinement levels meet), or from the nearest centre if no brick holds
// the point.
//
// Seeds: the nodes inside --seed-box (default: a slab |y| <= 5 AU around y = 0), every N-th of them
// with --seed-stride N; the other nodes are "not traced" (-1). Each seed is traced on its own, so the
// results do not depend on the thread count.
//
// Output: OUT.plt on the original grid (the same nodes and bricks): X, Y, Z, topology, the end of
// each direction (1 Sun, 2 ISM z limit, 3 ISM Hplim, 4 ISM domain edge, 5 long, 6 null), the path
// length of each direction and to the ISM (AU). With --points, also a small file of the seeds only
// (a scattered-point zone with the same variables).

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

Mesh read_mesh(const std::string& path, bool want_hplim) {
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
  for (int v = 0; v < nv; ++v) {
    std::vector<float>* dst = v == cx ? &m.x : v == cy ? &m.y : v == cz ? &m.z : v == cbx ? &m.bx
                            : v == cby ? &m.by : v == cbz ? &m.bz : (want_hplim && v == ch) ? &m.hplim : nullptr;
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
// the BATS-R-US block structure: located by one hash table per block size
// ---------------------------------------------------------------------------------------------

constexpr int kB = 8, kBB = kB * kB * kB;

struct Blocks {
  long long n = 0;
  std::vector<double> x0, y0, z0, dx;   // centre of cell (0,0,0) and cell size, per block
  int sx = 1, sy = kB, sz = kB * kB;    // node offset of one step along x, y, z inside a block
  double gx = 0, gy = 0, gz = 0;        // lower corner of the domain
  double ex = 0, ey = 0, ez = 0;        // upper corner of the domain
  struct Level { double size; std::unordered_map<uint64_t, int32_t> map; };
  std::vector<Level> levels;            // finest first

  static uint64_t key(long long i, long long j, long long k) {
    return (uint64_t(i & 0x1FFFFF) << 42) | (uint64_t(j & 0x1FFFFF) << 21) | uint64_t(k & 0x1FFFFF);
  }

  void build(const Mesh& m) {
    if (m.npts % kBB) throw std::runtime_error("the number of nodes is not a multiple of 8 x 8 x 8: not BATS-R-US blocks");
    n = m.npts / kBB;
    // which node offsets step along x, y and z: from the first block
    auto d = [&](long long a, long long b, const std::vector<float>& c) { return double(c[b]) - double(c[a]); };
    const long long o[3] = {1, kB, kB * kB};
    for (int s = 0; s < 3; ++s) {
      const double ax = std::fabs(d(0, o[s], m.x)), ay = std::fabs(d(0, o[s], m.y)), az = std::fabs(d(0, o[s], m.z));
      if (ax > ay && ax > az) sx = int(o[s]); else if (ay > az) sy = int(o[s]); else sz = int(o[s]);
    }
    if (sx == sy || sy == sz || sx == sz) throw std::runtime_error("cannot tell the axes of the blocks' node order");
    x0.resize(n); y0.resize(n); z0.resize(n); dx.resize(n);
    gx = gy = gz = 1e300; ex = ey = ez = -1e300;
    long long bad = 0;
    for (long long b = 0; b < n; ++b) {
      const long long base = b * kBB;
      double cxmin = 1e300, cymin = 1e300, czmin = 1e300, cxmax = -1e300;
      for (int q = 0; q < kBB; ++q) {
        cxmin = std::min(cxmin, double(m.x[base + q])); cxmax = std::max(cxmax, double(m.x[base + q]));
        cymin = std::min(cymin, double(m.y[base + q])); czmin = std::min(czmin, double(m.z[base + q]));
      }
      const double h = (cxmax - cxmin) / (kB - 1);
      x0[b] = cxmin; y0[b] = cymin; z0[b] = czmin; dx[b] = h;
      // the block must be a uniform 8 x 8 x 8 lattice in the node order found above
      for (int k = 0; k < kB && !bad; k += kB - 1)
        for (int j = 0; j < kB; j += kB - 1)
          for (int i = 0; i < kB; i += kB - 1) {
            const long long q = base + i * sx + j * sy + k * sz;
            if (std::fabs(m.x[q] - (cxmin + i * h)) > 1e-3 * h || std::fabs(m.y[q] - (cymin + j * h)) > 1e-3 * h ||
                std::fabs(m.z[q] - (czmin + k * h)) > 1e-3 * h)
              ++bad;
          }
      if (bad) throw std::runtime_error("block " + std::to_string(b) + " is not a uniform 8 x 8 x 8 lattice of cell centres");
      gx = std::min(gx, cxmin - h / 2); gy = std::min(gy, cymin - h / 2); gz = std::min(gz, czmin - h / 2);
      ex = std::max(ex, cxmin + 7.5 * h); ey = std::max(ey, cymin + 7.5 * h); ez = std::max(ez, czmin + 7.5 * h);
    }
    for (long long b = 0; b < n; ++b) {
      const double s = kB * dx[b];
      Level* L = nullptr;
      for (auto& l : levels) if (std::fabs(l.size - s) < 1e-6 * s) L = &l;
      if (!L) { levels.push_back({s, {}}); L = &levels.back(); }
      const long long i = std::llround((x0[b] - dx[b] / 2 - gx) / s), j = std::llround((y0[b] - dx[b] / 2 - gy) / s),
                      k = std::llround((z0[b] - dx[b] / 2 - gz) / s);
      L->map[key(i, j, k)] = int32_t(b);
    }
    std::sort(levels.begin(), levels.end(), [](const Level& a, const Level& b) { return a.size < b.size; });
  }

  // the block holding (x, y, z), or -1
  long long find(double x, double y, double z) const {
    for (const auto& L : levels) {
      const long long i = (long long)std::floor((x - gx) / L.size), j = (long long)std::floor((y - gy) / L.size),
                      k = (long long)std::floor((z - gz) / L.size);
      if (i < 0 || j < 0 || k < 0) continue;
      auto it = L.map.find(key(i, j, k));
      if (it != L.map.end()) return it->second;
    }
    return -1;
  }
};

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

struct Field {
  const Mesh& m;
  const Blocks& B;
  std::vector<int64_t> adj_start, adj_end;  // node -> bricks: adj[adj_start[i] .. adj_end[i])
  std::vector<int32_t> adj;

  Field(const Mesh& mesh, const Blocks& blocks) : m(mesh), B(blocks) {
    const std::size_t n = std::size_t(m.npts);
    adj_start.assign(n + 1, 0);
    for (int32_t v : m.conn) ++adj_start[std::size_t(v) + 1];
    for (std::size_t i = 0; i < n; ++i) adj_start[i + 1] += adj_start[i];
    adj.resize(m.conn.size());
    adj_end.assign(adj_start.begin(), adj_start.end() - 1);
    for (std::size_t e = 0; e < std::size_t(m.nelem); ++e)
      for (int c = 0; c < 8; ++c) {
        const int32_t v = m.conn[e * 8 + c];
        bool dup = false;  // a collapsed brick repeats nodes: list it once per node
        for (int c2 = 0; c2 < c; ++c2) dup |= m.conn[e * 8 + c2] == v;
        if (!dup) adj[adj_end[v]++] = int32_t(e);
      }
  }

  // B at (x, y, z) in three components; 0 = outside the domain, 1 = ok
  int sample(const double* p, double* out, double* hp = nullptr) const {
    const long long b = B.find(p[0], p[1], p[2]);
    if (b < 0) return 0;
    const double h = B.dx[b];
    const double fi = (p[0] - B.x0[b]) / h, fj = (p[1] - B.y0[b]) / h, fk = (p[2] - B.z0[b]) / h;
    const long long base = b * kBB;
    if (fi >= 0 && fi <= kB - 1 && fj >= 0 && fj <= kB - 1 && fk >= 0 && fk <= kB - 1) {
      const int i = std::min(int(fi), kB - 2), j = std::min(int(fj), kB - 2), k = std::min(int(fk), kB - 2);
      const double wi = fi - i, wj = fj - j, wk = fk - k;
      double v[4] = {0, 0, 0, 0};
      for (int c = 0; c < 8; ++c) {
        const int di = c & 1, dj = (c >> 1) & 1, dk = (c >> 2) & 1;
        const double w = (di ? wi : 1 - wi) * (dj ? wj : 1 - wj) * (dk ? wk : 1 - wk);
        const long long q = base + (i + di) * B.sx + (j + dj) * B.sy + (k + dk) * B.sz;
        v[0] += w * m.bx[q]; v[1] += w * m.by[q]; v[2] += w * m.bz[q];
        if (hp) v[3] += w * m.hplim[q];
      }
      out[0] = v[0]; out[1] = v[1]; out[2] = v[2];
      if (hp) *hp = v[3];
      return 1;
    }
    // the outer half cell of the block: the bricks of the nearest centre
    const int i = std::min(std::max(int(std::lround(fi)), 0), kB - 1), j = std::min(std::max(int(std::lround(fj)), 0), kB - 1),
              k = std::min(std::max(int(std::lround(fk)), 0), kB - 1);
    const long long q0 = base + i * B.sx + j * B.sy + k * B.sz;
    for (int64_t a = adj_start[q0]; a < adj_end[q0]; ++a) {
      const std::size_t e = std::size_t(adj[a]);
      const int32_t* nd = &m.conn[e * 8];
      double lo[3] = {1e300, 1e300, 1e300}, hi[3] = {-1e300, -1e300, -1e300};
      for (int c = 0; c < 8; ++c) {
        const double q[3] = {m.x[nd[c]], m.y[nd[c]], m.z[nd[c]]};
        for (int d = 0; d < 3; ++d) { lo[d] = std::min(lo[d], q[d]); hi[d] = std::max(hi[d], q[d]); }
      }
      const double tol = 1e-9 * (hi[0] - lo[0] + 1.0);
      if (p[0] < lo[0] - tol || p[0] > hi[0] + tol || p[1] < lo[1] - tol || p[1] > hi[1] + tol || p[2] < lo[2] - tol ||
          p[2] > hi[2] + tol)
        continue;
      // trilinear if the 8 nodes fill the 8 corners of the brick's box, else inverse-distance weights
      int corner_of[8];
      bool used[8] = {false, false, false, false, false, false, false, false}, full = true;
      for (int c = 0; c < 8; ++c) {
        const int32_t v = nd[c];
        int cc = 0;
        for (int d = 0; d < 3; ++d) {
          const double q = d == 0 ? m.x[v] : d == 1 ? m.y[v] : m.z[v];
          if (q > 0.5 * (lo[d] + hi[d])) cc |= 1 << d;
        }
        if (used[cc]) full = false;
        used[cc] = true;
        corner_of[c] = cc;
      }
      double v[4] = {0, 0, 0, 0}, wsum = 0;
      for (int c = 0; c < 8; ++c) {
        const int32_t nn = nd[c];
        double w;
        if (full) {
          const int cc = corner_of[c];
          w = 1.0;
          for (int d = 0; d < 3; ++d) {
            const double span = hi[d] - lo[d], t = span > 0 ? (p[d] - lo[d]) / span : 0.5;
            w *= (cc >> d & 1) ? t : 1 - t;
          }
        } else {
          const double ddx = p[0] - m.x[nn], ddy = p[1] - m.y[nn], ddz = p[2] - m.z[nn];
          w = 1.0 / (ddx * ddx + ddy * ddy + ddz * ddz + 1e-12);
        }
        v[0] += w * m.bx[nn]; v[1] += w * m.by[nn]; v[2] += w * m.bz[nn];
        if (hp) v[3] += w * m.hplim[nn];
        wsum += w;
      }
      out[0] = v[0] / wsum; out[1] = v[1] / wsum; out[2] = v[2] / wsum;
      if (hp) *hp = v[3] / wsum;
      return 1;
    }
    out[0] = m.bx[q0]; out[1] = m.by[q0]; out[2] = m.bz[q0];  // no brick holds it: the nearest centre
    if (hp) *hp = m.hplim[q0];
    return 1;
  }

  double cell_size(const double* p) const {
    const long long b = B.find(p[0], p[1], p[2]);
    return b < 0 ? 0.0 : B.dx[b];
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

  std::pair<int, double> follow(const double* x0, double sign) const {
    double x[3] = {x0[0], x0[1], x0[2]};
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
    Settings o;
    double box[6] = {-1e300, 1e300, -5.0, 5.0, -1e300, 1e300};
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
                               "[--points SEEDS.plt] [--z-top 500] [--z-bottom -500] [--hplim-ism ZMIN] [--r-sun 50] "
                               "[--max-length 10000] [--cell-fraction 0.5] [--substeps 2] [--threads N]");
    if (stride < 1 || !(o.cell_fraction > 0) || o.substeps < 1 || !(o.z_top > o.z_bottom))
      throw std::runtime_error("need --seed-stride >= 1, cell_fraction > 0, substeps >= 1, z_top > z_bottom");
    const auto t0 = Clock::now();

    Mesh m = read_mesh(in, o.hplim_zmin >= 0);
    std::printf("%s: %lld nodes, %lld bricks; read in %.0f s\n", in.c_str(), m.npts, m.nelem, since(t0));
    Blocks B;
    B.build(m);
    std::printf("  %lld blocks of 8 x 8 x 8 cells, %zu refinement levels (cells %.4g to %.4g AU); domain "
                "x %.0f..%.0f, y %.0f..%.0f, z %.0f..%.0f AU\n",
                B.n, B.levels.size(), B.levels.front().size / kB, B.levels.back().size / kB, B.gx, B.ex, B.gy, B.ey,
                B.gz, B.ez);
    auto t1 = Clock::now();
    Field F(m, B);
    std::printf("  node-to-brick table built in %.0f s\n", since(t1));

    std::vector<long long> seeds;
    long long inbox = 0;
    for (long long q = 0; q < m.npts; ++q)
      if (m.x[q] >= box[0] && m.x[q] <= box[1] && m.y[q] >= box[2] && m.y[q] <= box[3] && m.z[q] >= box[4] &&
          m.z[q] <= box[5])
        if (inbox++ % stride == 0) seeds.push_back(q);
    std::printf("  seeds: %zu nodes (%lld in the box, every %lld-th); ISM: z >= %g or z <= %g AU, the domain edge%s; "
                "Sun: r <= %g AU; at most %g AU per direction, step %g cell, %d substeps\n",
                seeds.size(), inbox, stride, o.z_top, o.z_bottom,
                o.hplim_zmin >= 0 ? (", Hplim < 0 at |z| >= " + std::to_string(o.hplim_zmin) + " AU").c_str() : "",
                o.r_sun, o.max_length, o.cell_fraction, o.substeps);

    const std::size_t ns = seeds.size();
    std::vector<int> topo(ns), ef(ns), eb(ns);
    std::vector<double> lf(ns), lb(ns);
    Tracer tr{F, o};
    t1 = Clock::now();
    long long done = 0;
#pragma omp parallel for schedule(dynamic, 16)
    for (long long s = 0; s < (long long)ns; ++s) {
      const long long q = seeds[s];
      const double x[3] = {m.x[q], m.y[q], m.z[q]};
      const auto f = tr.follow(x, +1.0), b = tr.follow(x, -1.0);
      ef[s] = f.first; lf[s] = f.second; eb[s] = b.first; lb[s] = b.second;
      topo[s] = classify(f.first, b.first);
      long long d;
#pragma omp atomic capture
      d = ++done;
      if (d % 100000 == 0) std::printf("    %lld of %zu seeds traced (%.0f s)\n", d, ns, since(t1));
    }
    long long count[6] = {0, 0, 0, 0, 0, 0};
    for (int t : topo) ++count[t + 1];
    const char* cname[6] = {"not traced", "undetermined", "closed", "open", "interstellar", "ISM-connected"};
    std::printf("  traced in %.0f s; seeds by topology:\n", since(t1));
    for (int t = 1; t < 6; ++t)
      if (count[t]) std::printf("    %-14s %12lld\n", cname[t], count[t]);

    // ---- output on the original grid ----
    const std::size_t n = std::size_t(m.npts);
    std::vector<float> T(n, -1.0f), EF(n, 0.0f), EB(n, 0.0f), LF(n, 0.0f), LB(n, 0.0f), LI(n, -1.0f);
    for (std::size_t s = 0; s < ns; ++s) {
      const std::size_t q = std::size_t(seeds[s]);
      T[q] = float(topo[s]); EF[q] = float(ef[s]); EB[q] = float(eb[s]); LF[q] = float(lf[s]); LB[q] = float(lb[s]);
      double li = -1;
      if (is_ism(ef[s])) li = lf[s];
      if (is_ism(eb[s]) && (li < 0 || lb[s] < li)) li = lb[s];
      LI[q] = float(li);
    }
    const std::vector<std::string> names = {"X AU", "Y AU", "Z AU", "topology", "end along B", "end against B",
                                            "length along B AU", "length against B AU", "length to ISM AU"};
    t1 = Clock::now();
    write_plt(out, "field-line topology of " + in, names, {&m.x, &m.y, &m.z, &T, &EF, &EB, &LF, &LB, &LI}, &m.conn);
    std::printf("  wrote %s (the original grid) in %.0f s\n", out.c_str(), since(t1));
    if (!points.empty()) {
      std::vector<float> px(ns), py(ns), pz(ns), pt(ns), pef(ns), peb(ns), plf(ns), plb(ns), pli(ns);
      for (std::size_t s = 0; s < ns; ++s) {
        const std::size_t q = std::size_t(seeds[s]);
        px[s] = m.x[q]; py[s] = m.y[q]; pz[s] = m.z[q];
        pt[s] = T[q]; pef[s] = EF[q]; peb[s] = EB[q]; plf[s] = LF[q]; plb[s] = LB[q]; pli[s] = LI[q];
      }
      write_plt(points, "field-line topology of the seeds of " + in, names,
                {&px, &py, &pz, &pt, &pef, &peb, &plf, &plb, &pli}, nullptr);
      std::printf("  wrote %s (the %zu seeds as scattered points)\n", points.c_str(), ns);
    }
    std::printf("  %.0f s in total\n", since(t0));
  } catch (const std::exception& e) {
    std::fprintf(stderr, "error: %s\n", e.what());
    return 1;
  }
  return 0;
}
