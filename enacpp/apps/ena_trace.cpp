// ena_trace: the streamlines of the ENA input preparation, traced in C++ instead of Tecplot.
//
//   ena_trace FIELD.dat STREAMLINES.dat --seeds SEEDS.txt [--cell-fraction 1] [--max-steps 750]
//             [--substeps 4] [--threads N]
//
// Tecplot's streamtraces are not reproducible: two traces of the same snapshot, even in the same
// Tecplot session on the same field, give different streamlines (and different streamlines stuck at
// the south pole), which changes the ENA maps by several percent. This program traces them with a
// fixed rule, so the same field always gives the same streamlines, on any machine and for any
// number of threads.
//
// It also follows the flow over the poles. In the macro grid the cells around the north pole
// collapse to a line, and Tecplot's streamtraces stop there: on the 2010-02-06 snapshot every
// streamline that heads over the north pole ends at theta = 6 deg (98-144 AU), and the seeds at
// theta 0 and 6 deg give none. The polar heliosheath then has no streamline points, and the ENA
// code fills it from its neighbours (north-pole maps up to 25% different). Here theta and phi are
// computed from the position, so the streamlines cross the pole and continue on the far side.
//
// FIELD.dat is the spherical grid with the streamline variables, written by Tecplot in step 3
// (3_calculate_files_trace.mcr): an ordered zone in point format, the variables X, Y, Z (AU), the
// flow velocity U_x, U_y, U_z (km/s) and the other variables to carry along the streamlines.
// STREAMLINES.dat is written as Tecplot writes extracted streamtraces (one zone per streamline,
// the same variables), so the ENA code reads it unchanged.
//
// The method:
//  - The field is sampled by trilinear interpolation in the grid's (r, theta, phi) index space,
//    with phi periodic. The grid must hold theta = 0 and 180 deg (the doubled macro grid, theta
//    0-354 deg, does), so the field is defined right at both poles; theta > 180 is not used.
//  - Each streamline follows the direction of the flow, U/|U|, from its seed, with the classical
//    4th-order Runge-Kutta method. One step is cell_fraction x (the smaller of the radial cell size
//    and the polar cell length r dtheta) at the start of the step: at most 2 AU on the macro grid,
//    about the spacing of Tecplot's streamtraces with CellFraction = 1 beyond 20 AU (measured on
//    the 2010-02-06 snapshot). Each step is split into `substeps` Runge-Kutta steps; a point is
//    written after every step.
//  - A streamline ends after max_steps points (Tecplot's MaxSteps), when it leaves the grid, or
//    where the flow speed is below 1e-6 km/s.
//  - Every variable is interpolated to the points, as Tecplot does, except the position and
//    theta0, phi0, r0, which are computed from the position with the step-3 macro's formulas.

#include "ena/ascii.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kDeg = kPi / 180.0;

double seconds_since(std::chrono::steady_clock::time_point t0) {
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

int num_threads() {
#ifdef _OPENMP
  return omp_get_max_threads();
#else
  return 1;
#endif
}

bool starts_number(char c) { return (c >= '0' && c <= '9') || c == '-' || c == '+' || c == '.'; }
bool is_blank(char c) { return c == ' ' || c == '\t' || c == '\r' || c == ',' || c == '\n'; }

std::string upper(std::string s) {
  for (char& c : s) c = char(std::toupper(static_cast<unsigned char>(c)));
  return s;
}

// A read-only memory map of a whole file.
struct MappedFile {
  const char* data = nullptr;
  size_t size = 0;
  explicit MappedFile(const std::string& path) {
    int fd = open(path.c_str(), O_RDONLY);
    if (fd < 0) throw std::runtime_error("cannot open " + path);
    struct stat st {};
    if (fstat(fd, &st) != 0) { close(fd); throw std::runtime_error("cannot stat " + path); }
    size = size_t(st.st_size);
    void* p = mmap(nullptr, size, PROT_READ, MAP_PRIVATE, fd, 0);
    close(fd);
    if (p == MAP_FAILED) throw std::runtime_error("cannot map " + path);
    madvise(p, size, MADV_SEQUENTIAL);
    data = static_cast<const char*>(p);
  }
  ~MappedFile() { if (data) munmap(const_cast<char*>(data), size); }
  MappedFile(const MappedFile&) = delete;
  MappedFile& operator=(const MappedFile&) = delete;
};

std::vector<std::string> quoted(std::string_view s) {
  std::vector<std::string> out;
  size_t i = 0;
  while ((i = s.find('"', i)) != std::string_view::npos) {
    size_t j = s.find('"', i + 1);
    if (j == std::string_view::npos) break;
    out.emplace_back(s.substr(i + 1, j - i - 1));
    i = j + 1;
  }
  return out;
}

// The value of KEY=value in a ZONE header (empty if absent).
std::string zone_value(const std::string& z, const std::string& key) {
  size_t i = 0;
  while ((i = z.find(key, i)) != std::string::npos) {
    bool word_start = i == 0 || !std::isalnum(static_cast<unsigned char>(z[i - 1]));
    size_t j = i + key.size();
    while (j < z.size() && z[j] == ' ') ++j;
    if (word_start && j < z.size() && z[j] == '=') {
      ++j;
      while (j < z.size() && z[j] == ' ') ++j;
      size_t k = j;
      while (k < z.size() && !is_blank(z[k])) ++k;
      return z.substr(j, k - j);
    }
    i = j;
  }
  return "";
}

// ---------------------------------------------------------------------------------------------
// The field: an ordered zone, read and rearranged into [variable][theta][phi][r]
// ---------------------------------------------------------------------------------------------

struct Field {
  std::vector<std::string> names;
  int nv = 0;
  int nr = 0, nt = 0, np = 0;            // r, theta (0..180 deg), phi points
  double r0 = 0, dr = 0, dt = 0, dp = 0; // AU, deg, deg
  std::vector<double> v;                 // [var][it][ip][ir]
  double at(int var, int it, int ip, int ir) const {
    return v[((size_t(var) * nt + it) * np + ip) * nr + ir];
  }
  double& at(int var, int it, int ip, int ir) { return v[((size_t(var) * nt + it) * np + ip) * nr + ir]; }
};

Field read_field(const std::string& path) {
  MappedFile f(path);
  // ---- header: TITLE, VARIABLES (quoted names, possibly over several lines), ZONE ----
  Field fd;
  std::string zone;
  enum { NONE, VARS, ZONE } section = NONE;
  size_t pos = 0, data = 0;
  while (pos < f.size) {
    size_t end = pos;
    while (end < f.size && f.data[end] != '\n') ++end;
    std::string_view line(f.data + pos, end - pos);
    size_t first = line.find_first_not_of(" \t\r");
    if (first != std::string_view::npos && starts_number(line[first]) && section == ZONE) {
      data = pos;
      break;
    }
    std::string head = first == std::string_view::npos ? "" : upper(std::string(line.substr(first, 9)));
    if (head.rfind("VARIABLES", 0) == 0) {
      section = VARS;
      for (auto& s : quoted(line)) fd.names.push_back(s);
    } else if (head.rfind("ZONE", 0) == 0) {
      if (section == ZONE) throw std::runtime_error(path + ": a second ZONE before any data");
      section = ZONE;
      zone += " " + std::string(line);
    } else if (head.rfind("TITLE", 0) == 0) {
      section = NONE;
    } else if (section == VARS) {
      for (auto& s : quoted(line)) fd.names.push_back(s);
    } else if (section == ZONE) {
      zone += " " + std::string(line);
    }
    pos = end + 1;
  }
  if (!data) throw std::runtime_error(path + ": no data after the ZONE header");
  fd.nv = int(fd.names.size());
  if (fd.nv < 6) throw std::runtime_error(path + ": expected at least X, Y, Z and three velocity components");
  const std::string z = upper(zone);
  const std::string sI = zone_value(z, "I"), sJ = zone_value(z, "J"), sK = zone_value(z, "K");
  const std::string pack = zone_value(z, "DATAPACKING").empty() ? zone_value(z, "F") : zone_value(z, "DATAPACKING");
  if (sI.empty() || sJ.empty() || sK.empty()) throw std::runtime_error(path + ": the zone is not ordered (I, J, K)");
  if (!pack.empty() && pack != "POINT") throw std::runtime_error(path + ": expected DATAPACKING=POINT, got " + pack);
  const long ni = std::stol(sI), nj = std::stol(sJ), nk = std::stol(sK);
  const size_t nodes = size_t(ni) * size_t(nj) * size_t(nk);

  // ---- the numbers, in order (one point after the other) ----
  std::vector<double> raw(nodes * size_t(fd.nv));
  const char* p = f.data + data;
  const char* e = f.data + f.size;
  for (size_t n = 0; n < raw.size(); ++n) {
    while (p < e && is_blank(*p)) ++p;
    if (p < e && *p == '+') ++p;
    const char* next = ena::read_double(p, e, raw[n]);  // from_chars where the compiler has it for doubles
    if (!next)
      throw std::runtime_error(path + ": expected " + std::to_string(raw.size()) + " numbers, could read " + std::to_string(n));
    p = next;
  }
  auto node = [&](long i, long j, long k) { return &raw[((size_t(k) * nj + j) * ni + i) * fd.nv]; };

  // ---- which of I, J, K is r, theta and phi: compare neighbouring nodes near theta = 90 deg ----
  auto sph = [](const double* q, double& r, double& th, double& ph) {
    r = std::sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2]);
    th = std::acos(std::max(-1.0, std::min(1.0, q[2] / r))) / kDeg;
    ph = std::atan2(q[1], q[0]) / kDeg;
    if (ph < 0) ph += 360.0;
  };
  const long dims[3] = {ni, nj, nk};
  int role[3] = {-1, -1, -1};  // 0 = r, 1 = theta, 2 = phi
  {
    // a node away from the poles: search along each index for theta near 90 deg
    long c[3] = {ni / 2, nj / 2, nk / 2};
    double r, th, ph;
    for (int d = 0; d < 3; ++d)
      for (long a = 0; a < dims[d]; ++a) {
        long t[3] = {c[0], c[1], c[2]};
        t[d] = a;
        sph(node(t[0], t[1], t[2]), r, th, ph);
        if (std::fabs(th - 90.0) < 4.0) { c[d] = a; break; }
      }
    double r0, t0, p0;
    sph(node(c[0], c[1], c[2]), r0, t0, p0);
    for (int d = 0; d < 3; ++d) {
      if (dims[d] < 2) throw std::runtime_error(path + ": every index of the grid needs at least 2 points");
      long t[3] = {c[0], c[1], c[2]};
      t[d] = c[d] + 1 < dims[d] ? c[d] + 1 : c[d] - 1;
      sph(node(t[0], t[1], t[2]), r, th, ph);
      double dph = std::fabs(ph - p0);
      if (dph > 180) dph = 360 - dph;
      const double dr = std::fabs(r - r0) / r0, dth = std::fabs(th - t0) * kDeg, dphr = dph * kDeg;
      role[d] = (dr >= dth && dr >= dphr) ? 0 : (dth >= dphr ? 1 : 2);
    }
    if (role[0] == role[1] || role[0] == role[2] || role[1] == role[2])
      throw std::runtime_error(path + ": cannot tell which grid index is r, theta and phi");
  }
  int ir_d = 0, it_d = 0, ip_d = 0;
  for (int d = 0; d < 3; ++d) (role[d] == 0 ? ir_d : role[d] == 1 ? it_d : ip_d) = d;
  auto nodeR = [&](long a_r, long a_t, long a_p) {
    long t[3];
    t[ir_d] = a_r; t[it_d] = a_t; t[ip_d] = a_p;
    return node(t[0], t[1], t[2]);
  };
  const long NR = dims[ir_d], NT = dims[it_d], NP = dims[ip_d];
  // spacings, from nodes away from the poles
  const long tmid = [&] {
    double r, th, ph;
    for (long a = 0; a < NT; ++a) { sph(nodeR(0, a, 0), r, th, ph); if (std::fabs(th - 90.0) < 4.0) return a; }
    return NT / 2;
  }();
  double ra, ta, pa, rb, tb, pb;
  sph(nodeR(0, tmid, 0), ra, ta, pa);
  sph(nodeR(1, tmid, 0), rb, tb, pb);
  fd.r0 = ra;
  fd.dr = rb - ra;
  sph(nodeR(0, tmid, 1), rb, tb, pb);
  fd.dp = pb - pa;
  if (fd.dp < 0) fd.dp += 360.0;
  // theta of each row, from the z coordinate (also valid at the poles)
  std::vector<double> th_row(NT);
  for (long a = 0; a < NT; ++a) {
    const double* q = nodeR(NR - 1, a, 0);
    double r = std::sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2]);
    th_row[a] = std::acos(std::max(-1.0, std::min(1.0, q[2] / r))) / kDeg;
  }
  fd.dt = th_row[1] - th_row[0];
  if (!(fd.dr > 0 && fd.dt > 0 && fd.dp > 0)) throw std::runtime_error(path + ": r, theta and phi must increase with their indices");
  // checks: uniform spacing, theta from 0, phi covering the full circle
  for (long a = 0; a < NR; ++a) {
    double r, th, ph;
    sph(nodeR(a, tmid, 0), r, th, ph);
    if (std::fabs(r - (fd.r0 + a * fd.dr)) > 1e-4 * (fd.r0 + a * fd.dr))
      throw std::runtime_error(path + ": the radial spacing is not uniform");
  }
  if (std::fabs(th_row[0]) > 1e-3) throw std::runtime_error(path + ": the grid must start at theta = 0 (the north pole)");
  if (std::fabs(NP * fd.dp - 360.0) > 1e-3) throw std::runtime_error(path + ": phi must cover 0-360 deg (" +
                                                                     std::to_string(NP) + " x " + std::to_string(fd.dp) + " deg)");
  // rows used: theta 0 .. 180 deg; theta = 180 must be one of them (the doubled macro grid has it)
  long t180 = -1;
  for (long a = 0; a < NT; ++a) {
    if (std::fabs(th_row[a] - a * fd.dt) > 1e-3 && th_row[a] <= 180.0 - 1e-3)
      throw std::runtime_error(path + ": the theta spacing is not uniform");
    if (std::fabs(th_row[a] - 180.0) < 1e-3 || std::fabs(a * fd.dt - 180.0) < 1e-3) { t180 = a; break; }
  }
  if (t180 < 0) throw std::runtime_error(path + ": the grid has no row at theta = 180 deg (use the doubled theta grid, 0-354 deg)");
  fd.nr = int(NR);
  fd.nt = int(t180 + 1);
  fd.np = int(NP);
  fd.v.assign(size_t(fd.nv) * fd.nt * fd.np * fd.nr, 0.0);
  for (int it = 0; it < fd.nt; ++it)
    for (int ip = 0; ip < fd.np; ++ip)
      for (int ir = 0; ir < fd.nr; ++ir) {
        const double* q = nodeR(ir, it, ip);
        for (int var = 0; var < fd.nv; ++var) fd.at(var, it, ip, ir) = q[var];
      }
  return fd;
}

// ---------------------------------------------------------------------------------------------
// Sampling and tracing
// ---------------------------------------------------------------------------------------------

struct Tracer {
  const Field& f;
  int cx = 0, cy = 1, cz = 2, cux = 3, cuy = 4, cuz = 5;
  int ctheta0 = -1, cphi0 = -1, cr0 = -1;  // computed from the position, not interpolated
  double rmin = 0, rmax = 0;
  double cell_fraction = 1.0;
  int max_steps = 750, substeps = 4;
  static constexpr double umin = 1e-6;  // km/s

  explicit Tracer(const Field& fld) : f(fld) {
    rmin = f.r0;
    rmax = f.r0 + (f.nr - 1) * f.dr;
  }

  // Trilinear interpolation in (r, theta, phi) index space of the variables `vars` at point x.
  // Returns false outside the grid.
  bool sample(const double* x, const int* vars, int nvars, double* out) const {
    const double r = std::sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
    const double tol = 1e-9 * rmax;
    if (!(r >= rmin - tol && r <= rmax + tol)) return false;
    double th = std::acos(std::max(-1.0, std::min(1.0, x[2] / r))) / kDeg;
    double ph = std::atan2(x[1], x[0]) / kDeg;
    if (ph < 0) ph += 360.0;
    double fr = (r - f.r0) / f.dr, ft = th / f.dt, fp = ph / f.dp;
    fr = std::min(std::max(fr, 0.0), double(f.nr - 1));
    ft = std::min(std::max(ft, 0.0), double(f.nt - 1));
    int ir = std::min(int(fr), f.nr - 2), it = std::min(int(ft), f.nt - 2), ip = int(fp);
    const double wr = fr - ir, wt = ft - it, wp = fp - ip;
    ip = ((ip % f.np) + f.np) % f.np;
    const int ip1 = (ip + 1) % f.np;
    for (int n = 0; n < nvars; ++n) {
      const int var = vars[n];
      auto lerp_r = [&](int t, int p) { return (1.0 - wr) * f.at(var, t, p, ir) + wr * f.at(var, t, p, ir + 1); };
      const double a = (1.0 - wp) * lerp_r(it, ip) + wp * lerp_r(it, ip1);
      const double b = (1.0 - wp) * lerp_r(it + 1, ip) + wp * lerp_r(it + 1, ip1);
      out[n] = (1.0 - wt) * a + wt * b;
    }
    return true;
  }

  // The flow direction U/|U| at x; false outside the grid or where the flow is too slow.
  bool direction(const double* x, double* d) const {
    const int uv[3] = {cux, cuy, cuz};
    double u[3];
    if (!sample(x, uv, 3, u)) return false;
    const double s = std::sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
    if (!(s > umin)) return false;
    d[0] = u[0] / s; d[1] = u[1] / s; d[2] = u[2] / s;
    return true;
  }

  // All variables at x, as written to the streamline file.
  bool point(const double* x, std::vector<double>& row) const {
    static thread_local std::vector<int> all;
    if (int(all.size()) != f.nv) { all.resize(f.nv); for (int i = 0; i < f.nv; ++i) all[i] = i; }
    if (!sample(x, all.data(), f.nv, row.data())) return false;
    row[cx] = x[0]; row[cy] = x[1]; row[cz] = x[2];
    // theta0, phi0, r0 with the formulas of the step-3 macro (theta, phi, r of the grid)
    const double r = std::sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
    if (cr0 >= 0) row[cr0] = r;
    if (ctheta0 >= 0) row[ctheta0] = std::acos(x[2] / (r + 1e-4)) * 57.2958;
    if (cphi0 >= 0) {
      double ph = std::atan2(x[1], x[0]) * 57.2958;
      row[cphi0] = ph < 0 ? 360.0 - std::fabs(ph) : ph;
    }
    return true;
  }

  // One streamline from seed s: rows of all variables. Empty if the seed is outside the grid.
  std::vector<double> trace(const double* seed, std::string* stop) const {
    std::vector<double> rows, row(f.nv);
    double x[3] = {seed[0], seed[1], seed[2]};
    if (!point(x, row)) { *stop = "seed outside the grid"; return rows; }
    rows.insert(rows.end(), row.begin(), row.end());
    *stop = "max_steps";
    for (int n = 1; n < max_steps; ++n) {
      const double r = std::sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
      const double h = cell_fraction * std::min(f.dr, r * f.dt * kDeg) / substeps;
      double y[3] = {x[0], x[1], x[2]};
      bool ok = true;
      for (int s = 0; s < substeps && ok; ++s) {
        double k1[3], k2[3], k3[3], k4[3], t[3];
        ok = direction(y, k1);
        if (ok) { for (int i = 0; i < 3; ++i) t[i] = y[i] + 0.5 * h * k1[i]; ok = direction(t, k2); }
        if (ok) { for (int i = 0; i < 3; ++i) t[i] = y[i] + 0.5 * h * k2[i]; ok = direction(t, k3); }
        if (ok) { for (int i = 0; i < 3; ++i) t[i] = y[i] + h * k3[i]; ok = direction(t, k4); }
        if (ok)
          for (int i = 0; i < 3; ++i) y[i] += (h / 6.0) * (((k1[i] + 2.0 * k2[i]) + 2.0 * k3[i]) + k4[i]);
      }
      if (!ok || !point(y, row)) { *stop = "left the grid or reached a stagnation point"; break; }
      x[0] = y[0]; x[1] = y[1]; x[2] = y[2];
      rows.insert(rows.end(), row.begin(), row.end());
    }
    return rows;
  }
};

int find_column(const std::vector<std::string>& names, std::initializer_list<const char*> prefixes, int fallback) {
  for (const char* pre : prefixes)
    for (size_t i = 0; i < names.size(); ++i) {
      std::string n = names[i];
      size_t sp = n.find(' ');
      std::string base = sp == std::string::npos ? n : n.substr(0, sp);  // "U_x km/s" -> "U_x"
      if (upper(base) == upper(pre)) return int(i);
    }
  return fallback;
}

std::vector<double> read_seeds(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("cannot read the seed file " + path);
  std::vector<double> s;
  std::string line;
  while (std::getline(in, line)) {
    size_t c = line.find('#');
    if (c != std::string::npos) line.resize(c);
    std::istringstream is(line);
    double a, b, d;
    if (is >> a >> b >> d) { s.push_back(a); s.push_back(b); s.push_back(d); }
  }
  if (s.empty()) throw std::runtime_error(path + ": no seeds (x y z per line)");
  return s;
}

}  // namespace

int main(int argc, char** argv) {
  std::setvbuf(stdout, nullptr, _IOLBF, 0);
  try {
    std::string field_path, out_path, seeds_path;
    double cell_fraction = 1.0;
    int max_steps = 750, substeps = 4;
    for (int i = 1; i < argc; ++i) {
      std::string a = argv[i];
      auto next = [&]() -> std::string {
        if (i + 1 >= argc) throw std::runtime_error(a + " needs a value");
        return argv[++i];
      };
      if (a == "--seeds") seeds_path = next();
      else if (a == "--cell-fraction") cell_fraction = std::stod(next());
      else if (a == "--max-steps") max_steps = std::stoi(next());
      else if (a == "--substeps") substeps = std::stoi(next());
      else if (a == "--threads") {
        int n = std::stoi(next());
#ifdef _OPENMP
        omp_set_num_threads(n);
#endif
        (void)n;
      } else if (a == "-h" || a == "--help") {
        std::printf("usage: ena_trace FIELD.dat STREAMLINES.dat --seeds SEEDS.txt [--cell-fraction 1] "
                    "[--max-steps 750] [--substeps 4] [--threads N]\n");
        return 0;
      } else if (field_path.empty()) field_path = a;
      else if (out_path.empty()) out_path = a;
      else throw std::runtime_error("too many arguments: " + a);
    }
    if (field_path.empty() || out_path.empty() || seeds_path.empty())
      throw std::runtime_error("usage: ena_trace FIELD.dat STREAMLINES.dat --seeds SEEDS.txt [options]");
    if (!(cell_fraction > 0) || max_steps < 2 || substeps < 1)
      throw std::runtime_error("cell_fraction > 0, max_steps >= 2 and substeps >= 1 are required");
    {
      // the output is written under a temporary name and renamed over STREAMLINES.dat: refuse a
      // path that exists but is not a regular file (e.g. /dev/null), which the rename would replace
      struct stat st {};
      if (stat(out_path.c_str(), &st) == 0 && !S_ISREG(st.st_mode))
        throw std::runtime_error(out_path + " exists and is not a regular file; give a file name for the streamlines");
    }
    const auto t0 = std::chrono::steady_clock::now();

    Field fld = read_field(field_path);
    std::printf("%s: %d variables; grid r %.6g-%.6g AU (%d x %.6g), theta 0-180 deg (%d x %.6g), phi (%d x %.6g); %.1f s\n",
                field_path.c_str(), fld.nv, fld.r0, fld.r0 + (fld.nr - 1) * fld.dr, fld.nr, fld.dr, fld.nt, fld.dt,
                fld.np, fld.dp, seconds_since(t0));
    Tracer tr(fld);
    tr.cx = find_column(fld.names, {"X"}, 0);
    tr.cy = find_column(fld.names, {"Y"}, 1);
    tr.cz = find_column(fld.names, {"Z"}, 2);
    tr.cux = find_column(fld.names, {"U_x", "Ux"}, 3);
    tr.cuy = find_column(fld.names, {"U_y", "Uy"}, 4);
    tr.cuz = find_column(fld.names, {"U_z", "Uz"}, 5);
    tr.ctheta0 = find_column(fld.names, {"theta0"}, -1);
    tr.cphi0 = find_column(fld.names, {"phi0"}, -1);
    tr.cr0 = find_column(fld.names, {"r0"}, -1);
    tr.cell_fraction = cell_fraction;
    tr.max_steps = max_steps;
    tr.substeps = substeps;
    std::printf("  velocity: %s, %s, %s; step %.3g cell, %d Runge-Kutta substeps, at most %d points; %d threads\n",
                fld.names[tr.cux].c_str(), fld.names[tr.cuy].c_str(), fld.names[tr.cuz].c_str(), cell_fraction, substeps,
                max_steps, num_threads());

    const std::vector<double> seeds = read_seeds(seeds_path);
    const int ns = int(seeds.size() / 3);
    std::vector<std::vector<double>> lines(ns);
    std::vector<std::string> why(ns);
    // each streamline on its own: the result does not depend on the number of threads
#pragma omp parallel for schedule(dynamic)
    for (int s = 0; s < ns; ++s) lines[s] = tr.trace(&seeds[3 * size_t(s)], &why[s]);

    // ---- write as Tecplot writes extracted streamtraces ----
    long written = 0, points = 0, outside = 0, single = 0, full = 0;
    {
      std::string part = out_path + ".tmp";
      std::FILE* o = std::fopen(part.c_str(), "w");
      if (!o) throw std::runtime_error("cannot write " + out_path);
      std::setvbuf(o, nullptr, _IOFBF, 1 << 22);
      std::fprintf(o, "TITLE     = \"Streamlines traced by ena_trace\"\nVARIABLES = ");
      for (int v = 0; v < fld.nv; ++v) std::fprintf(o, v ? "\"%s\"\n" : "\"%s\"\n", fld.names[v].c_str());
      for (int s = 0; s < ns; ++s) {
        const size_t n = lines[s].size() / size_t(fld.nv);
        if (n == 0) { ++outside; continue; }
        if (n < 2) { ++single; continue; }  // stopped at the first step (flow leaving the grid, or none)
        if (int(n) == max_steps) ++full;
        std::fprintf(o, "ZONE T=\"Streamtrace\"\n STRANDID=0, SOLUTIONTIME=0\n I=%zu, J=1, K=1, ZONETYPE=Ordered\n"
                        " DATAPACKING=POINT\n DT=(", n);
        for (int v = 0; v < fld.nv; ++v) std::fprintf(o, "DOUBLE ");
        std::fprintf(o, ")\n");
        char buf[32];
        for (size_t q = 0; q < n; ++q) {
          for (int v = 0; v < fld.nv; ++v) {
            int len = std::snprintf(buf, sizeof buf, " %.9E", lines[s][q * fld.nv + v]);
            std::fwrite(buf, 1, size_t(len), o);
          }
          std::fputc('\n', o);
        }
        ++written;
        points += long(n);
      }
      if (std::fclose(o) != 0) throw std::runtime_error("write failed (disk full?)");
      if (std::rename(part.c_str(), out_path.c_str()) != 0) throw std::runtime_error("cannot rename to " + out_path);
    }
    std::printf("wrote %s: %ld streamlines (%ld points) from %d seeds; %ld reached %d points, %ld seeds outside the grid, "
                "%ld stopped at the first step (flow out of the grid, or none); %.1f s in total\n",
                out_path.c_str(), written, points, ns, full, max_steps, outside, single, seconds_since(t0));
    return written > 0 ? 0 : 1;
  } catch (const std::exception& e) {
    std::fprintf(stderr, "ena_trace: %s\n", e.what());
    return 1;
  }
}
