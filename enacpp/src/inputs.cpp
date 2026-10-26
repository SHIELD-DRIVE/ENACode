#include <cmath>
#include <cstdio>
#include <sstream>
#include <stdexcept>

#include "ena/ascii.hpp"
#include "ena/snapshot.hpp"

#ifdef ENA_HAVE_HDF5
#include "ena/h5io.hpp"
#endif

namespace ena {

namespace {

// Places rows of a grid file into cells and keeps track of problems.
class CellPlacer {
 public:
  CellPlacer(const RunConfig& cfg, const std::string& path, const ColumnMap& cols)
      : g_(cfg.grid), path_(path), allow_incomplete_(cfg.inputs.allow_incomplete) {
    ci_ = need(cols, "i"); cj_ = need(cols, "j"); ck_ = need(cols, "k");
    ct_ = cols.column("theta"); cp_ = cols.column("phi"); cr_ = cols.column("r");
    seen_.assign(g_.size(), 0);
  }

  // Returns the cell index for a row, or throws after recording the problem.
  std::size_t place(const double* v, long line) {
    int i = int(std::lround(v[ci_])), j = int(std::lround(v[cj_])), k = int(std::lround(v[ck_]));
    if (i < 0 || i >= g_.ntheta || j < 0 || j >= g_.nphi || k < 0 || k >= g_.nr)
      fail(line, "index (" + std::to_string(i) + "," + std::to_string(j) + "," + std::to_string(k) + ") is outside the grid");
    if (ct_ >= 0 && std::fabs(v[ct_] - g_.theta(i)) > 0.5 * g_.dtheta)
      fail(line, "theta " + num(v[ct_]) + " does not match grid point " + num(g_.theta(i)) + " for i=" + std::to_string(i));
    // At the pole (i = 0) phi is undefined and files write 0 for every j.
    if (cp_ >= 0 && i != 0 && std::fabs(v[cp_] - g_.phi(j)) > 0.5 * g_.dphi)
      fail(line, "phi " + num(v[cp_]) + " does not match grid point " + num(g_.phi(j)) + " for j=" + std::to_string(j));
    if (cr_ >= 0 && std::fabs(v[cr_] - g_.r(k)) > 0.5 * g_.dr)
      fail(line, "r " + num(v[cr_]) + " does not match grid point " + num(g_.r(k)) + " for k=" + std::to_string(k));
    std::size_t idx = g_.index(i, j, k);
    if (seen_[idx]) fail(line, "cell (" + std::to_string(i) + "," + std::to_string(j) + "," + std::to_string(k) + ") appears twice");
    seen_[idx] = 1;
    ++placed_;
    last_i_ = i; last_j_ = j; last_k_ = k;
    return idx;
  }

  void finish(bool verbose) const {
    std::size_t missing = g_.size() - placed_;
    if (missing == 0) return;
    std::string msg = path_ + ": " + std::to_string(missing) + " of " + std::to_string(g_.size()) + " grid cells have no row";
    if (!allow_incomplete_) throw std::runtime_error(msg + " (set inputs.allow_incomplete = true only for testing on file slices)");
    if (verbose) std::fprintf(stderr, "warning: %s; they are left as 0\n", msg.c_str());
  }

  int i() const { return last_i_; }
  int j() const { return last_j_; }
  int k() const { return last_k_; }

 private:
  static int need(const ColumnMap& m, const char* role) {
    int c = m.column(role);
    if (c < 0) throw std::runtime_error(std::string("column map lacks required column '") + role + "'");
    return c;
  }
  static std::string num(double x) { std::ostringstream o; o << x; return o.str(); }
  [[noreturn]] void fail(long line, const std::string& what) const {
    throw std::runtime_error(path_ + ":" + std::to_string(line) + ": " + what);
  }
  const Grid& g_;
  std::string path_;
  bool allow_incomplete_;
  int ci_, cj_, ck_, ct_, cp_, cr_;
  std::vector<unsigned char> seen_;
  std::size_t placed_ = 0;
  int last_i_ = -1, last_j_ = -1, last_k_ = -1;
};

int col(const ColumnMap& m, const char* role) {
  int c = m.column(role);
  if (c < 0) throw std::runtime_error(std::string("column map lacks required column '") + role + "'");
  return c;
}

}  // namespace

void read_plasma_ascii(const RunConfig& cfg, const std::string& path, Snapshot& s, const ReadOptions& opt) {
  const ColumnMap& m = cfg.inputs.plasma;
  const Grid& g = cfg.grid;
  s.plasma.assign(1, PlasmaFluid{});
  PlasmaFluid& f = s.plasma[0];
  for (auto* fld : {&f.n, &f.T, &f.vr, &f.vtheta, &f.vphi, &f.vt}) fld->resize(g);
  if (opt.keep_file_coordinates) for (auto* fld : {&s.coords.theta, &s.coords.phi, &s.coords.r}) fld->resize(g);

  const int cn = col(m, "n"), cT = col(m, "T"), cvr = col(m, "vr"), cvth = col(m, "vtheta"),
            cvph = col(m, "vphi"), cvt = col(m, "vt");
  const int ct = m.column("theta"), cp = m.column("phi"), cr = m.column("r");
  const double fr = m.velocity_factor("vr"), fth = m.velocity_factor("vtheta"),
               fph = m.velocity_factor("vphi"), ft = m.velocity_factor("vt");
  // The MHD temperature is halved, except for the Moscow model (moscow_temp).
  const double tfac = cfg.moscow() ? 1.0 : 0.5;

  CellPlacer placer(cfg, path, m);
  for_each_row(path, int(m.columns.size()), [&](const double* v, long line) {
    std::size_t c = placer.place(v, line);
    f.n.raw()[c] = v[cn];
    f.T.raw()[c] = v[cT] * tfac;
    f.vr.raw()[c] = v[cvr] * fr;
    f.vtheta.raw()[c] = v[cvth] * fth;
    f.vphi.raw()[c] = v[cvph] * fph;
    f.vt.raw()[c] = v[cvt] * ft;
    if (opt.keep_file_coordinates) {
      if (ct >= 0) s.coords.theta.raw()[c] = v[ct];
      if (cp >= 0) s.coords.phi.raw()[c] = v[cp];
      if (cr >= 0) s.coords.r.raw()[c] = v[cr];
    }
  });
  placer.finish(opt.verbose);
}

void read_plasma_mi_ascii(const RunConfig& cfg, const std::string& path, Snapshot& s, const ReadOptions& opt) {
  const ColumnMap& m = cfg.inputs.plasma;
  const Grid& g = cfg.grid;
  s.plasma.assign(3, PlasmaFluid{});
  for (auto& f : s.plasma)
    for (auto* fld : {&f.n, &f.T, &f.vr, &f.vtheta, &f.vphi, &f.vt}) fld->resize(g);
  s.pbeta.resize(g);
  if (opt.keep_file_coordinates) for (auto* fld : {&s.coords.theta, &s.coords.phi, &s.coords.r}) fld->resize(g);
  struct Cols { int n, T, vr, vth, vph, vt; double fr, fth, fph, ft; };
  Cols cs[3];
  const char* suffix[3] = {"", "_pu3", "_pu2"};  // SW, pu3, pu2
  for (int f = 0; f < 3; ++f) {
    auto nm = [&](const char* base) { return std::string(base) + suffix[f]; };
    cs[f] = {col(m, nm("n").c_str()), col(m, nm("T").c_str()), col(m, nm("vr").c_str()), col(m, nm("vtheta").c_str()),
             col(m, nm("vphi").c_str()), col(m, nm("vt").c_str()), m.velocity_factor(nm("vr")),
             m.velocity_factor(nm("vtheta")), m.velocity_factor(nm("vphi")), m.velocity_factor(nm("vt"))};
  }
  const int ct = m.column("theta"), cp = m.column("phi"), cr = m.column("r"), cb = m.column("pbeta");
  CellPlacer placer(cfg, path, m);
  for_each_row(path, int(m.columns.size()), [&](const double* v, long line) {
    std::size_t c = placer.place(v, line);
    for (int f = 0; f < 3; ++f) {
      PlasmaFluid& pf = s.plasma[f];
      pf.n.raw()[c] = v[cs[f].n];
      pf.T.raw()[c] = v[cs[f].T];  // the multi-ion code does not halve the temperatures
      pf.vr.raw()[c] = v[cs[f].vr] * cs[f].fr;
      pf.vtheta.raw()[c] = v[cs[f].vth] * cs[f].fth;
      pf.vphi.raw()[c] = v[cs[f].vph] * cs[f].fph;
      pf.vt.raw()[c] = v[cs[f].vt] * cs[f].ft;
    }
    if (cb >= 0) s.pbeta.raw()[c] = v[cb];
    if (opt.keep_file_coordinates) {
      if (ct >= 0) s.coords.theta.raw()[c] = v[ct];
      if (cp >= 0) s.coords.phi.raw()[c] = v[cp];
      if (cr >= 0) s.coords.r.raw()[c] = v[cr];
    }
  });
  placer.finish(opt.verbose);
}

void read_neutral_mi_ascii(const RunConfig& cfg, const std::string& path, Snapshot& s, const ReadOptions& opt) {
  const ColumnMap& m = cfg.inputs.neutral;
  s.neutrals.assign(5, NeutralPopulation{});
  int cn[5], cvr[5], cvt[5], cT[5], cvx[5], cvy[5], cvz[5];
  for (int q = 0; q < 5; ++q) {  // q = 0: total H; 1-4: populations (pristine ISM, OHS, IHS, supersonic SW)
    std::string suf = q == 0 ? "" : std::to_string(q);
    NeutralPopulation& p = s.neutrals[q];
    for (auto* fld : {&p.n, &p.T, &p.vr, &p.vt, &p.vx, &p.vy, &p.vz}) fld->resize(cfg.grid);
    cn[q] = col(m, ("n" + suf).c_str()); cvr[q] = col(m, ("vr" + suf).c_str()); cvt[q] = col(m, ("vt" + suf).c_str());
    cT[q] = col(m, ("T" + suf).c_str()); cvx[q] = col(m, ("vx" + suf).c_str()); cvy[q] = col(m, ("vy" + suf).c_str());
    cvz[q] = col(m, ("vz" + suf).c_str());
  }
  CellPlacer placer(cfg, path, m);
  for_each_row(path, int(m.columns.size()), [&](const double* v, long line) {
    std::size_t c = placer.place(v, line);
    for (int q = 0; q < 5; ++q) {
      NeutralPopulation& p = s.neutrals[q];
      p.n.raw()[c] = v[cn[q]]; p.vr.raw()[c] = v[cvr[q]]; p.vt.raw()[c] = v[cvt[q]]; p.T.raw()[c] = v[cT[q]];
      p.vx.raw()[c] = v[cvx[q]]; p.vy.raw()[c] = v[cvy[q]]; p.vz.raw()[c] = v[cvz[q]];
    }
  });
  placer.finish(opt.verbose);
}

void read_neutral_ascii(const RunConfig& cfg, const std::string& path, Snapshot& s, const ReadOptions& opt) {
  const ColumnMap& m = cfg.inputs.neutral;
  s.neutrals.assign(1, NeutralPopulation{});
  NeutralPopulation& p = s.neutrals[0];
  p.n.resize(cfg.grid);
  p.T.resize(cfg.grid);
  const int cn = col(m, "n"), cT = col(m, "T");
  CellPlacer placer(cfg, path, m);
  for_each_row(path, int(m.columns.size()), [&](const double* v, long line) {
    std::size_t c = placer.place(v, line);
    p.n.raw()[c] = v[cn];
    p.T.raw()[c] = v[cT];
  });
  placer.finish(opt.verbose);
}

void read_boundaries_ascii(const RunConfig& cfg, const std::string& path, Snapshot& s, const ReadOptions& opt) {
  const ColumnMap& m = cfg.inputs.boundaries;
  const Grid& g = cfg.grid;
  s.region.resize(g, 0);
  for (auto* fld : {&s.Bx, &s.By, &s.Bz}) fld->resize(g);
  const int creg = col(m, "region"), cbx = m.column("Bx"), cby = m.column("By"), cbz = m.column("Bz");
  CellPlacer placer(cfg, path, m);
  long unexpected_moscow = 0;
  for_each_row(path, int(m.columns.size()), [&](const double* v, long line) {
    // Rows are placed by their i, j, k columns, so any row order works.
    std::size_t c = placer.place(v, line);
    auto rd = [&](int c) { return v[c]; };  // double precision
    // Region codes are written as floats: truncated.
    int code = int(std::trunc(rd(creg)));
    if (cfg.moscow()) {
      int mapped = 5;  // undetermined, caught by the region-5 check
      if (code == 3 || code == 4) mapped = 2;
      else if (code == 1) mapped = 4;
      else if (code == 2) mapped = 3;
      else if (code != 0) ++unexpected_moscow;
      code = mapped;
    }
    s.region.raw()[c] = std::int16_t(code);
    if (cbx >= 0) s.Bx.raw()[c] = rd(cbx);
    if (cby >= 0) s.By.raw()[c] = rd(cby);
    if (cbz >= 0) s.Bz.raw()[c] = rd(cbz);
  });
  placer.finish(opt.verbose);
  if (cfg.multi_ion()) {
    // Multi-ion model: a region-3 cell directly beyond a region-1/2 cell (along
    // r, already relabelled) takes that region, so heliosheath stretches
    // re-entered beyond the heliopause count as ISM / outer heliosheath.
    for (int i = 0; i < g.ntheta; ++i)
      for (int j = 0; j < g.nphi; ++j)
        for (int k = 1; k < g.nr; ++k) {
          std::int16_t& b = s.region(i, j, k);
          const std::int16_t before = s.region(i, j, k - 1);
          if (b == 3 && (before == 1 || before == 2)) b = before;
        }
  }
  if (unexpected_moscow && opt.verbose)
    std::fprintf(stderr, "warning: %s: %ld Moscow region codes outside 0-4 (set to 5, undetermined)\n", path.c_str(), unexpected_moscow);
}

Snapshot read_snapshot(const RunConfig& cfg, std::size_t index, const ReadOptions& opt) {
  if (index >= cfg.inputs.snapshots.size()) throw std::runtime_error("snapshot index out of range");
  const SnapshotFiles& f = cfg.inputs.snapshots[index];
  Snapshot s;
  s.time = f.time;
  if (cfg.inputs.format == "hdf5") {
#ifdef ENA_HAVE_HDF5
    s = read_snapshot_h5(cfg, cfg.inputs.path(f.file), opt);
    s.time = f.time;
    return s;
#else
    throw std::runtime_error("inputs.format = \"hdf5\" needs a build with HDF5");
#endif
  }
  if (cfg.multi_ion()) {
    read_plasma_mi_ascii(cfg, cfg.inputs.path(f.plasma), s, opt);
    read_neutral_mi_ascii(cfg, cfg.inputs.path(f.neutral), s, opt);
  } else if (cfg.extended_si()) {  // single-ion plasma file, 41-column neutral file
    read_plasma_ascii(cfg, cfg.inputs.path(f.plasma), s, opt);
    read_neutral_mi_ascii(cfg, cfg.inputs.path(f.neutral), s, opt);
  } else {
    read_plasma_ascii(cfg, cfg.inputs.path(f.plasma), s, opt);
    read_neutral_ascii(cfg, cfg.inputs.path(f.neutral), s, opt);
  }
  read_boundaries_ascii(cfg, cfg.inputs.path(f.boundaries), s, opt);
  return s;
}

}  // namespace ena
