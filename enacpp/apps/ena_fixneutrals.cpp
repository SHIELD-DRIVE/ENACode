// ena_fixneutrals: step 1 of the single-ion input preparation, in C++ instead of Tecplot.
//
//   ena_fixneutrals NEUTRAL.dat FIXED.plt [--threads N]
//
// Reads the neutral solution (a Tecplot ASCII file: one FEPOINT zone of BRICK elements, one point
// per line, 53 variables: X, Y, Z, B, E and 11 variables for each of the 4 populations), computes
// what configs/macros/single_ion/1_fix_neutrals.mcr computes, and writes the result as a Tecplot
// binary .plt file with the same 33 variables, in the same order and with the same names.
//
// The macro's equations are evaluated the way Tecplot evaluates them: every variable is stored in
// single precision (as Tecplot stores the ASCII data and the variables it creates), and each
// equation is computed in double precision from those single-precision values, in the order the
// equation is written, then stored in single precision. The result should therefore agree with
// Tecplot's to the last digit or very nearly.
//
// The text file is read through a memory map and parsed by all threads at once (OpenMP): the
// 18 GB neutral files are read in a few minutes, mostly limited by the disk.
//
// Tecplot then loads the .plt and saves it as .szplt (macro 1b_save_szplt.mcr): Tecplot 2025 R2
// interpolates quickly only from a .szplt neutral file.

#include <algorithm>
#include <charconv>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
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

// ---------------------------------------------------------------------------------------------
// Small helpers
// ---------------------------------------------------------------------------------------------

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

// Is this character the start of a number (a data line starts with one)?
bool starts_number(char c) { return (c >= '0' && c <= '9') || c == '-' || c == '+' || c == '.'; }

bool is_blank(char c) { return c == ' ' || c == '\t' || c == '\r' || c == ','; }

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

// ---------------------------------------------------------------------------------------------
// The header of the Tecplot ASCII file
// ---------------------------------------------------------------------------------------------

struct Header {
  std::string title;
  std::vector<std::string> variables;
  long long num_points = -1;    // N=
  long long num_elements = -1;  // E=
  std::string packing;          // F= / DATAPACKING=
  std::string element_type;     // ET= / ZONETYPE=
  size_t data_start = 0;        // byte offset of the first data line
};

// The quoted strings of a line, in order (the variable names, the title).
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

std::string upper(std::string s) {
  for (char& c : s) c = char(std::toupper(static_cast<unsigned char>(c)));
  return s;
}

// The value of KEY=value in a ZONE header (empty if absent). Values end at a blank or comma.
std::string zone_value(const std::string& zone_upper, const std::string& key) {
  size_t i = 0;
  while ((i = zone_upper.find(key, i)) != std::string::npos) {
    bool word_start = i == 0 || !std::isalnum(static_cast<unsigned char>(zone_upper[i - 1]));
    size_t j = i + key.size();
    while (j < zone_upper.size() && zone_upper[j] == ' ') ++j;
    if (word_start && j < zone_upper.size() && zone_upper[j] == '=') {
      ++j;
      while (j < zone_upper.size() && zone_upper[j] == ' ') ++j;
      size_t k = j;
      while (k < zone_upper.size() && !is_blank(zone_upper[k])) ++k;
      return zone_upper.substr(j, k - j);
    }
    i = j;
  }
  return "";
}

Header read_header(const MappedFile& f) {
  Header h;
  std::string zone;  // the ZONE header text (it can span several lines)
  enum { NONE, VARIABLES, ZONE } section = NONE;
  size_t pos = 0;
  while (pos < f.size) {
    size_t end = pos;
    while (end < f.size && f.data[end] != '\n') ++end;
    std::string_view line(f.data + pos, end - pos);
    size_t first = line.find_first_not_of(" \t\r");
    if (first != std::string_view::npos && starts_number(line[first]) && section == ZONE) {
      h.data_start = pos;  // the first line of numbers after the ZONE header
      break;
    }
    std::string head = first == std::string_view::npos ? "" : upper(std::string(line.substr(first, 9)));
    if (head.rfind("TITLE", 0) == 0) {
      auto q = quoted(line);
      if (!q.empty()) h.title = q[0];
      section = NONE;
    } else if (head.rfind("VARIABLES", 0) == 0) {
      section = VARIABLES;
      for (auto& v : quoted(line)) h.variables.push_back(v);
    } else if (head.rfind("ZONE", 0) == 0) {
      section = ZONE;
      zone += " " + std::string(line);
    } else if (section == VARIABLES) {
      for (auto& v : quoted(line)) h.variables.push_back(v);
    } else if (section == ZONE) {
      zone += " " + std::string(line);
    }
    pos = end + 1;
  }
  if (h.data_start == 0) throw std::runtime_error("no data found after the ZONE header");
  std::string z = upper(zone);
  auto num = [&](const std::string& key) {
    std::string v = zone_value(z, key);
    return v.empty() ? -1LL : std::stoll(v);
  };
  h.num_points = num("N");
  if (h.num_points < 0) h.num_points = num("NODES");
  h.num_elements = num("E");
  if (h.num_elements < 0) h.num_elements = num("ELEMENTS");
  h.packing = zone_value(z, "F");
  if (h.packing.empty()) h.packing = zone_value(z, "DATAPACKING");
  h.element_type = zone_value(z, "ET");
  if (h.element_type.empty()) h.element_type = zone_value(z, "ZONETYPE");
  return h;
}

// ---------------------------------------------------------------------------------------------
// Parallel parsing of the data lines
// ---------------------------------------------------------------------------------------------

// The data lines split into byte ranges that each start at the beginning of a line, with the
// index of the first (non-blank) line of each range, so that the ranges can be parsed in parallel.
struct LineRanges {
  std::vector<size_t> begin, end;   // byte range [begin, end) of each part
  std::vector<long long> first;     // index of its first non-blank line
  long long total = 0;              // non-blank lines in all parts
};

bool blank_line(const char* p, const char* eol) {
  while (p < eol && (is_blank(*p))) ++p;
  return p == eol;
}

LineRanges split_lines(const MappedFile& f, size_t start, int parts) {
  LineRanges r;
  const size_t bytes = f.size - start;
  parts = int(std::max<size_t>(1, std::min<size_t>(size_t(parts), bytes / 4096 + 1)));
  for (int k = 0; k < parts; ++k) {
    size_t b = start + bytes / size_t(parts) * size_t(k);
    if (k > 0) {  // move to the start of the next line
      const void* q = std::memchr(f.data + b, '\n', f.size - b);
      b = q ? size_t(static_cast<const char*>(q) - f.data) + 1 : f.size;
    }
    r.begin.push_back(b);
  }
  for (int k = 0; k < parts; ++k) r.end.push_back(k + 1 < parts ? r.begin[k + 1] : f.size);
  // count the non-blank lines of each part, in parallel
  std::vector<long long> count(parts, 0);
#pragma omp parallel for schedule(dynamic)
  for (int k = 0; k < parts; ++k) {
    const char* p = f.data + r.begin[k];
    const char* e = f.data + r.end[k];
    long long n = 0;
    while (p < e) {
      const char* eol = static_cast<const char*>(std::memchr(p, '\n', size_t(e - p)));
      if (!eol) eol = e;
      if (!blank_line(p, eol)) ++n;
      p = eol + 1;
    }
    count[k] = n;
  }
  for (int k = 0; k < parts; ++k) { r.first.push_back(r.total); r.total += count[k]; }
  return r;
}

// Parses one number at p (skipping blanks and commas first); returns the position after it.
const char* parse_double(const char* p, const char* end, double& v) {
  while (p < end && is_blank(*p)) ++p;
  if (p < end && *p == '+') ++p;  // from_chars does not accept a leading '+'
  auto r = std::from_chars(p, end, v);
  if (r.ec != std::errc()) {
    // Fortran-style exponents (1.0D+05) or other forms: fall back to strtod on a copy
    char buf[64];
    size_t n = 0;
    while (p + n < end && n < sizeof(buf) - 1 && !is_blank(p[n]) && p[n] != '\n') {
      buf[n] = (p[n] == 'd' || p[n] == 'D') ? 'e' : p[n];
      ++n;
    }
    buf[n] = 0;
    char* stop = nullptr;
    v = std::strtod(buf, &stop);
    if (stop == buf) throw std::runtime_error("cannot read a number: " + std::string(buf));
    return p + (stop - buf);
  }
  return r.ptr;
}

const char* parse_int(const char* p, const char* end, long long& v) {
  while (p < end && is_blank(*p)) ++p;
  if (p < end && *p == '+') ++p;
  auto r = std::from_chars(p, end, v);
  if (r.ec != std::errc()) throw std::runtime_error("cannot read an element connection");
  return r.ptr;
}

// ---------------------------------------------------------------------------------------------
// Tecplot binary (.plt) output, version 112 (Tecplot 360 Data Format Guide, "Binary Data")
// ---------------------------------------------------------------------------------------------

class PltWriter {
 public:
  explicit PltWriter(const std::string& path) : f_(std::fopen(path.c_str(), "wb")) {
    if (!f_) throw std::runtime_error("cannot write " + path);
    std::setvbuf(f_, nullptr, _IOFBF, 1 << 24);
  }
  ~PltWriter() { if (f_) std::fclose(f_); }
  void bytes(const void* p, size_t n) {
    if (n && std::fwrite(p, 1, n, f_) != n) throw std::runtime_error("write failed (disk full?)");
  }
  void i32(int32_t v) { bytes(&v, 4); }
  void f32(float v) { bytes(&v, 4); }
  void f64(double v) { bytes(&v, 8); }
  // a string: each character as an INT32, then a 0
  void str(const std::string& s) {
    for (unsigned char c : s) i32(int32_t(c));
    i32(0);
  }
  void close() {
    if (f_ && std::fclose(f_) != 0) { f_ = nullptr; throw std::runtime_error("write failed (disk full?)"); }
    f_ = nullptr;
  }

 private:
  std::FILE* f_;
};

}  // namespace

// ---------------------------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------------------------

int main(int argc, char** argv) {
  // progress lines go out at once (the driver shows them in its log as they come)
  std::setvbuf(stdout, nullptr, _IOLBF, 0);
  try {
    std::string in, out;
    for (int i = 1; i < argc; ++i) {
      std::string a = argv[i];
      if (a == "--threads" && i + 1 < argc) {
#ifdef _OPENMP
        omp_set_num_threads(std::stoi(argv[++i]));
#else
        ++i;
#endif
      } else if (a == "-h" || a == "--help") {
        std::printf("usage: ena_fixneutrals NEUTRAL.dat FIXED.plt [--threads N]\n"
                    "Step 1 of the single-ion input preparation (1_fix_neutrals.mcr) in C++.\n");
        return 0;
      } else if (in.empty()) {
        in = a;
      } else if (out.empty()) {
        out = a;
      } else {
        throw std::runtime_error("too many arguments: " + a);
      }
    }
    if (in.empty() || out.empty()) throw std::runtime_error("usage: ena_fixneutrals NEUTRAL.dat FIXED.plt [--threads N]");
    const auto t0 = std::chrono::steady_clock::now();

    // --- the header ---------------------------------------------------------------------------
    MappedFile f(in);
    Header h = read_header(f);
    const int nv = int(h.variables.size());
    std::printf("%s: %zu bytes, %d variables, %lld points, %lld elements, %s %s; %d threads\n", in.c_str(),
                f.size, nv, h.num_points, h.num_elements, h.packing.c_str(), h.element_type.c_str(), num_threads());
    if (h.num_points <= 0 || h.num_elements <= 0) throw std::runtime_error("the ZONE header lacks N= or E=");
    if (h.packing != "FEPOINT" && h.packing != "POINT")
      throw std::runtime_error("expected point data (F=FEPOINT), got " + h.packing);
    if (h.element_type != "BRICK" && h.element_type != "FEBRICK")
      throw std::runtime_error("expected BRICK elements, got " + h.element_type);
    // The macro refers to the variables by number (v10 ... v49) and to the pressures by name
    // (pxxPop1 ...): check that this file has the layout the macro was recorded on.
    if (nv != 53) throw std::runtime_error("expected the 53 variables of the neutral file, got " + std::to_string(nv));
    const char* q[11] = {"rho", "ux", "uy", "uz", "pxx", "pyy", "pzz", "pxy", "pxz", "pyz", "ppc"};
    for (int p = 0; p < 4; ++p)
      for (int k = 0; k < 11; ++k) {
        std::string want = std::string(q[k]) + "Pop" + std::to_string(p + 1);
        if (h.variables[9 + 11 * p + k] != want)
          throw std::runtime_error("variable " + std::to_string(10 + 11 * p + k) + " is " +
                                   h.variables[9 + 11 * p + k] + ", expected " + want);
      }

    // --- the points: keep the 31 columns the macro uses ----------------------------------------
    // column (0-based) -> slot in `in_vars`: X, Y, Z, then for each population rho, ux, uy, uz,
    // pxx, pyy, pzz (v10-v16, v21-v27, v32-v38, v43-v49 in the macro's numbering)
    std::vector<int> slot(nv, -1);
    slot[0] = 0; slot[1] = 1; slot[2] = 2;
    for (int p = 0; p < 4; ++p)
      for (int k = 0; k < 7; ++k) slot[9 + 11 * p + k] = 3 + 7 * p + k;
    const int nin = 31;
    const long long np = h.num_points;
    std::vector<std::vector<float>> col(nin, std::vector<float>(size_t(np)));

    // One pass over all data lines, in parallel: lines 0 .. np-1 are the points, the next ne lines
    // the elements (8 point numbers each, 1-based).
    const long long ne = h.num_elements;
    std::vector<int32_t> conn(size_t(ne) * 8);
    LineRanges lines = split_lines(f, h.data_start, num_threads() * 16);
    std::printf("  found %lld data lines: %.1f s\n", lines.total, seconds_since(t0));
    if (lines.total < np + ne)
      throw std::runtime_error("the file has " + std::to_string(lines.total) + " data lines; the ZONE header promises " +
                               std::to_string(np) + " points + " + std::to_string(ne) + " elements");
    if (lines.total > np + ne)
      std::printf("  note: %lld lines after the elements are ignored (a second zone is not read)\n", lines.total - np - ne);
    const int nparts = int(lines.first.size());
    std::string error;
#pragma omp parallel for schedule(dynamic)
    for (int k = 0; k < nparts; ++k) {
      try {
        const char* p = f.data + lines.begin[k];
        const char* end = f.data + lines.end[k];
        long long i = lines.first[k];  // index of the next non-blank line
        while (p < end && i < np + ne) {
          const char* eol = static_cast<const char*>(std::memchr(p, '\n', size_t(end - p)));
          if (!eol) eol = end;
          if (blank_line(p, eol)) { p = eol + 1; continue; }
          if (i < np) {  // a point: 53 numbers
            for (int c = 0; c < nv; ++c) {
              while (p < eol && is_blank(*p)) ++p;
              if (p >= eol) throw std::runtime_error("point line " + std::to_string(i + 1) + " has fewer than " + std::to_string(nv) + " values");
              double v;
              p = parse_double(p, eol, v);
              if (slot[c] >= 0) col[slot[c]][size_t(i)] = float(v);  // Tecplot stores single precision
            }
            if (!blank_line(p, eol)) throw std::runtime_error("point line " + std::to_string(i + 1) + " has more than " + std::to_string(nv) + " values");
          } else {  // an element: 8 point numbers
            const long long e = i - np;
            for (int c = 0; c < 8; ++c) {
              long long v;
              p = parse_int(p, eol, v);
              if (v < 1 || v > np) throw std::runtime_error("element " + std::to_string(e + 1) + " refers to point " + std::to_string(v));
              conn[size_t(e) * 8 + c] = int32_t(v - 1);  // the binary format is 0-based
            }
          }
          p = eol + 1;
          ++i;
        }
      } catch (const std::exception& ex) {
#pragma omp critical
        error = ex.what();
      }
    }
    if (!error.empty()) throw std::runtime_error(error);
    std::printf("  read the points and elements: %.1f s\n", seconds_since(t0));

    // --- the equations of 1_fix_neutrals.mcr ------------------------------------------------------
    // Population order in the file: Pop1 = OHS, Pop2 = IHS, Pop3 = SSW, Pop4 = ISM (the renames).
    // Output: the 33 variables Tecplot writes after the macro's $!DeleteVars.
    const std::vector<std::string> out_names = {
        "X AU", "Y AU", "Z AU",
        "nH_OHS", "UxH_OHS", "UyH_OHS", "UzH_OHS", "nH_IHS", "UxH_IHS", "UyH_IHS", "UzH_IHS",
        "nH_SSW", "UxH_SSW", "UyH_SSW", "UzH_SSW", "nH_ISM", "UxH_ISM", "UyH_ISM", "UzH_ISM",
        "nH", "UxH", "UyH", "UzH", "P", "TH",
        "P_OHS", "TH_OHS", "P_IHS", "TH_IHS", "P_SSW", "TH_SSW", "P_ISM", "TH_ISM"};
    const int nout = int(out_names.size());
    // the first 19 outputs are input columns (X, Y, Z and n, ux, uy, uz of each population)
    auto in_col = [&](int p, int k) -> std::vector<float>& { return col[3 + 7 * p + k]; };  // k: 0 rho,1 ux,2 uy,3 uz,4 pxx,5 pyy,6 pzz
    std::vector<std::vector<float>> made(14, std::vector<float>(size_t(np)));  // nH UxH UyH UzH P TH and P_/TH_ x4
    long long zero_nh = 0;
    // the constants as the macro writes them, evaluated left to right as Tecplot does
    const double c_mass = 1.67e-27;
    const double c_kms2 = 1000.0 * 1000.0;  // 1000**2
    const double c_cm3 = 100.0 * 100.0 * 100.0;  // 100**3
    const double c_k = 1.38e-23;
#pragma omp parallel for schedule(static) reduction(+ : zero_nh)
    for (long long i = 0; i < np; ++i) {
      const size_t s = size_t(i);
      // values as Tecplot holds them: single precision, used in double precision
      double n[4], ux[4], uy[4], uz[4], pxx[4], pyy[4], pzz[4];
      for (int p = 0; p < 4; ++p) {
        n[p] = in_col(p, 0)[s];  ux[p] = in_col(p, 1)[s];  uy[p] = in_col(p, 2)[s];  uz[p] = in_col(p, 3)[s];
        pxx[p] = in_col(p, 4)[s]; pyy[p] = in_col(p, 5)[s]; pzz[p] = in_col(p, 6)[s];
      }
      // order in the macro's sums: OHS (v10), IHS (v21), SSW (v32), ISM (v43) = Pop1..Pop4
      // {nH}=v10+v21+v32+v43
      const float nH = float(((n[0] + n[1]) + n[2]) + n[3]);
      if (nH == 0.0f) ++zero_nh;
      // {UxH}=(v10*v11+v21*v22+v32*v33+v43*v44)/{nH}   (likewise UyH, UzH)
      auto mean = [&](const double* u) {
        return float((((n[0] * u[0] + n[1] * u[1]) + n[2] * u[2]) + n[3] * u[3]) / double(nH));
      };
      const float UxH = mean(ux), UyH = mean(uy), UzH = mean(uz);
      // {Pxx}=(v14+v25+v36+v47)*1e-9+(v10*(v11-{UxH})**2+...)*1.67e-27*(1000**2)*(100**3)
      auto ptot = [&](const double* pd, const double* u, float U) {
        const double th = (((pd[0] + pd[1]) + pd[2]) + pd[3]) * 1e-9;
        const double d0 = u[0] - double(U), d1 = u[1] - double(U), d2 = u[2] - double(U), d3 = u[3] - double(U);
        const double dyn = ((n[0] * (d0 * d0) + n[1] * (d1 * d1)) + n[2] * (d2 * d2)) + n[3] * (d3 * d3);
        return float(th + ((dyn * c_mass) * c_kms2) * c_cm3);
      };
      const float Pxx = ptot(pxx, ux, UxH), Pyy = ptot(pyy, uy, UyH), Pzz = ptot(pzz, uz, UzH);
      // {P}=1/3*({Pxx}+{Pyy}+{Pzz});  {TH}={P}/({nH}*(100**3)*1.38e-23+1e-30)
      const float P = float((1.0 / 3.0) * ((double(Pxx) + double(Pyy)) + double(Pzz)));
      const float TH = float(double(P) / (((double(nH) * c_cm3) * c_k) + 1e-30));
      made[0][s] = nH; made[1][s] = UxH; made[2][s] = UyH; made[3][s] = UzH; made[4][s] = P; made[5][s] = TH;
      // {P_OHS}=1/3*({pxxPop1}+{pyyPop1}+{pzzPop1})*1e-9;  {TH_OHS}={P_OHS}/({nH_OHS}*(100**3)*1.38e-23+1e-30)
      for (int p = 0; p < 4; ++p) {
        const float Pp = float(((1.0 / 3.0) * ((pxx[p] + pyy[p]) + pzz[p])) * 1e-9);
        const float Tp = float(double(Pp) / (((n[p] * c_cm3) * c_k) + 1e-30));
        made[6 + 2 * p][s] = Pp;
        made[7 + 2 * p][s] = Tp;
      }
    }
    if (zero_nh)
      std::printf("  note: %lld points have total nH = 0, so UxH, UyH, UzH are 0/0 there (NaN)\n", zero_nh);
    std::printf("  computed the variables: %.1f s\n", seconds_since(t0));

    // the 33 output columns, in Tecplot's order
    std::vector<const std::vector<float>*> outcol;
    outcol.push_back(&col[0]); outcol.push_back(&col[1]); outcol.push_back(&col[2]);
    for (int p : {0, 1, 2, 3})  // OHS, IHS, SSW, ISM = Pop1..Pop4
      for (int k = 0; k < 4; ++k) outcol.push_back(&in_col(p, k));
    for (int m = 0; m < 14; ++m) outcol.push_back(&made[m]);
    if (int(outcol.size()) != nout) throw std::runtime_error("internal error: output columns");

    // min and max of each variable (the binary format stores them before the data)
    std::vector<double> vmin(nout), vmax(nout);
#pragma omp parallel for schedule(dynamic)
    for (int v = 0; v < nout; ++v) {
      double lo = std::numeric_limits<double>::infinity(), hi = -lo;
      for (float x : *outcol[v])
        if (std::isfinite(x)) { lo = std::min(lo, double(x)); hi = std::max(hi, double(x)); }
      if (lo > hi) lo = hi = 0.0;
      vmin[v] = lo;
      vmax[v] = hi;
    }

    // --- write the .plt (version 112) -------------------------------------------------------------
    {
      PltWriter w(out + ".part");
      // header section
      w.bytes("#!TDV112", 8);
      w.i32(1);                      // byte order
      w.i32(0);                      // file type: full
      w.str(h.title);                // title
      w.i32(nout);                   // number of variables
      for (auto& s : out_names) w.str(s);
      w.f32(299.0f);                 // zone marker
      w.str("ZONE 001");             // zone name (as Tecplot names this zone when it loads the file)
      w.i32(-1);                     // parent zone: none
      w.i32(-1);                     // strand ID: static
      w.f64(0.0);                    // solution time
      w.i32(-1);                     // not used
      w.i32(5);                      // zone type: FEBRICK
      w.i32(0);                      // variable location: all at the nodes
      w.i32(0);                      // raw local 1-to-1 face neighbours supplied: no
      w.i32(0);                      // number of miscellaneous user-defined face neighbour connections: 0
      w.i32(int32_t(np));            // number of points
      w.i32(int32_t(ne));            // number of elements
      w.i32(0); w.i32(0); w.i32(0);  // ICellDim, JCellDim, KCellDim (unused)
      w.i32(0);                      // no auxiliary data
      w.f32(357.0f);                 // end of header
      // data section
      w.f32(299.0f);                 // zone marker
      for (int v = 0; v < nout; ++v) w.i32(1);  // all variables single precision
      w.i32(0);                      // no passive variables
      w.i32(0);                      // no variable sharing
      w.i32(-1);                     // connectivity not shared
      for (int v = 0; v < nout; ++v) { w.f64(vmin[v]); w.f64(vmax[v]); }
      for (int v = 0; v < nout; ++v) w.bytes(outcol[v]->data(), outcol[v]->size() * sizeof(float));
      w.bytes(conn.data(), conn.size() * sizeof(int32_t));
      w.close();
    }
    if (std::rename((out + ".part").c_str(), out.c_str()) != 0) throw std::runtime_error("cannot rename to " + out);
    std::printf("wrote %s: %d variables, %lld points, %lld elements; %.1f s in total\n", out.c_str(), nout, np, ne,
                seconds_since(t0));
    return 0;
  } catch (const std::exception& e) {
    std::fprintf(stderr, "ena_fixneutrals: %s\n", e.what());
    return 1;
  }
}
