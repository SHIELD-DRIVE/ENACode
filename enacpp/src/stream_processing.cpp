#include "ena/stream_processing.hpp"

#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <string>
#include <tuple>

#include "ena/physics.hpp"

namespace ena {

namespace {

const double NaN = std::numeric_limits<double>::quiet_NaN();

// Nearest integer, halves away from zero; NaN gives the most negative 32-bit
// integer, so such points are never binned.
std::int32_t round_index(double x) {
  if (std::isnan(x)) return INT32_MIN;
  return std::int32_t(std::lround(x));
}

bool fin(double x) { return std::isfinite(x); }

struct Fit { double a, b, c; };
// Kornbleuth et al. (2023) Eqs. 6-13 (Table 1 values).
const Fit kFits[10] = {
    {-4.65496, 2.46160, -1.71480},    // ntr
    {-2.27104, 0.395963, -2.65926},   // nref
    {8.10281, -3.06275, -2.52573},    // nen
    {0.101026, 0.472896, -0.579818},  // etr
    {-15.8367, 3.91212, -1.71480},    // eref
    {2.73671, -0.477155, -0.223144},  // een
    {18.7905, -10.5773, 3.21888},     // em
    {-13.5408, 6.52448, 2.01490},     // kappa_ref (delta)
    {0.0, 0.0, 0.0},                  // vshock: not a fit
    {2.73671, -0.477155, -3.21888},   // esw
};

}  // namespace

double reference_ts_distance(const RunConfig& cfg, const Snapshot& snap, std::string* note) {
  const PhysicsConfig& ph = cfg.physics;
  if (ph.ts_reference == "constant") {
    if (note) {
      char buf[64];
      std::snprintf(buf, sizeof buf, "constant %g AU", ph.ts_reference_r);
      *note = buf;
    }
    return ph.ts_reference_r;
  }
  const Grid& g = cfg.grid;
  const PlasmaFluid& pl = snap.plasma.at(0);
  const bool traj = !ph.v2_trajectory_file.empty();
  if (traj && ph.v2_r.size() < 2)
    throw std::runtime_error("physics.v2_trajectory_file was not read (use load_config)");
  if (!traj && ph.v2_direction.size() != 2) throw std::runtime_error("no V2 path for ts_reference = \"computed\"");

  // Direction of the V2 path at radius r: interpolated in the trajectory file.
  // Beyond its last row (97 AU in configs/v2_trajectory_sim.txt) theta and phi
  // are extrapolated linearly in r, with a least-squares line through the rows
  // of the last 5 AU of the file (the multi-ion TS along V2 lies beyond
  // 97 AU). Before the first row there is no path (false).
  double ext_th0 = 0, ext_th1 = 0, ext_ph0 = 0, ext_ph1 = 0;  // theta = th0 + th1 * (r - r_last), same for phi
  if (traj) {
    const auto& R = ph.v2_r;
    const double rl = R.back();
    double sw = 0, sx = 0, sxx = 0, sth = 0, sxth = 0, sph = 0, sxph = 0;
    for (std::size_t i = 0; i < R.size(); ++i) {
      if (R[i] < rl - 5.0) continue;
      const double x = R[i] - rl;
      sw += 1; sx += x; sxx += x * x;
      sth += ph.v2_theta[i]; sxth += x * ph.v2_theta[i];
      sph += ph.v2_phi[i]; sxph += x * ph.v2_phi[i];
    }
    const double det = sw * sxx - sx * sx;
    if (sw >= 2 && det > 0) {
      ext_th1 = (sw * sxth - sx * sth) / det; ext_th0 = (sth - ext_th1 * sx) / sw;
      ext_ph1 = (sw * sxph - sx * sph) / det; ext_ph0 = (sph - ext_ph1 * sx) / sw;
    } else {
      ext_th0 = ph.v2_theta.back(); ext_ph0 = ph.v2_phi.back();
    }
  }
  auto direction = [&](double r, double& th, double& phi) {
    if (!traj) {
      th = ph.v2_direction[0];
      phi = ph.v2_direction[1];
      return true;
    }
    const auto& R = ph.v2_r;
    if (r < R.front()) return false;
    if (r > R.back()) {
      const double x = r - R.back();
      th = std::min(180.0, std::max(0.0, ext_th0 + ext_th1 * x));
      phi = ext_ph0 + ext_ph1 * x;
      return true;
    }
    std::size_t i = std::upper_bound(R.begin(), R.end(), r) - R.begin();
    if (i >= R.size()) i = R.size() - 1;
    const double w = (r - R[i - 1]) / (R[i] - R[i - 1]);
    th = ph.v2_theta[i - 1] + w * (ph.v2_theta[i] - ph.v2_theta[i - 1]);
    phi = ph.v2_phi[i - 1] + w * (ph.v2_phi[i] - ph.v2_phi[i - 1]);
    return true;
  };
  auto pindex = [&](int p) { return ((p % g.nphi) + g.nphi) % g.nphi; };
  // M^2 at (theta, phi) on grid radius k: bilinear over the four nearest lines of sight
  auto m2_at = [&](double th, double phi, int k) {
    const double ft = std::min(std::max(th / g.dtheta, 0.0), double(g.ntheta - 1));
    const int t0 = int(std::floor(ft)), t1 = std::min(t0 + 1, g.ntheta - 1);
    const double wt = ft - t0;
    const double fp = std::fmod(std::fmod(phi, 360.0) + 360.0, 360.0) / g.dphi;
    const int p0 = pindex(int(std::floor(fp))), p1 = pindex(p0 + 1);
    const double wp = fp - std::floor(fp);
    auto m2 = [&](int t, int p) {
      double v2 = pl.vr(t, p, k) * pl.vr(t, p, k) + pl.vt(t, p, k) * pl.vt(t, p, k);
      double cs2 = 5.0 * (2.0 * pl.T(t, p, k)) * cgs::kb / (3.0 * cgs::mp);
      return v2 / cs2;
    };
    return (1 - wt) * ((1 - wp) * m2(t0, p0) + wp * m2(t0, p1)) + wt * ((1 - wp) * m2(t1, p0) + wp * m2(t1, p1));
  };
  // region at (theta, phi) on grid radius k: nearest line of sight
  auto region_at = [&](double th, double phi, int k) {
    const int t = std::min(std::max(int(std::lround(th / g.dtheta)), 0), g.ntheta - 1);
    const int p = pindex(int(std::lround(std::fmod(std::fmod(phi, 360.0) + 360.0, 360.0) / g.dphi)));
    return int(snap.region(t, p, k));
  };
  const bool by_regions = ph.ts_detection == "regions" || ph.ts_detection == "moscow_column";
  double r_ts = std::numeric_limits<double>::quiet_NaN(), th_ts = 0, phi_ts = 0;
  for (int k = 1; k < g.nr && std::isnan(r_ts); ++k) {
    double th0, ph0, th1, ph1;
    if (!direction(g.r(k - 1), th0, ph0) || !direction(g.r(k), th1, ph1)) continue;
    if (by_regions) {  // region change 4 -> 3/6, midway between the two grid radii
      const int a = region_at(th0, ph0, k - 1), b = region_at(th1, ph1, k);
      if (a == 4 && (b == 3 || b == 6)) r_ts = g.r(k - 1) + 0.5 * g.dr;
    } else {  // first crossing of M^2 = 0.64, linear in r
      const double thr = 0.64, a = m2_at(th0, ph0, k - 1), b = m2_at(th1, ph1, k);
      if (a > thr && b <= thr) r_ts = g.r(k - 1) + (a - thr) / (a - b) * g.dr;
    }
    if (!std::isnan(r_ts)) direction(r_ts, th_ts, phi_ts);
  }
  if (std::isnan(r_ts)) {
    char buf[200];
    if (traj)
      std::snprintf(buf, sizeof buf, "ts_reference = \"computed\": no TS found along the V2 trajectory (from %g AU outward)",
                    ph.v2_r.front());
    else
      std::snprintf(buf, sizeof buf, "ts_reference = \"computed\": no TS found along theta %g, phi %g",
                    ph.v2_direction[0], ph.v2_direction[1]);
    throw std::runtime_error(buf);
  }
  if (note) {
    char buf[320];
    std::snprintf(buf, sizeof buf, "computed %.2f AU along the V2 %s (theta %.2f, phi %.2f there%s), %s method",
                  r_ts, traj ? "trajectory" : "direction", th_ts, std::fmod(phi_ts + 360.0, 360.0),
                  traj && r_ts > ph.v2_r.back() ? "; beyond the last row of the trajectory file: extrapolated from its last 5 AU" : "",
                  ph.ts_detection.c_str());
    *note = buf;
  }
  return r_ts;
}

std::string crossing_diagnostics(const RunConfig& cfg, const StreamPoints& s) {
  // For each region-based TS crossing (row i), walk down the same streamline:
  // path length (sum of straight segments, AU) to the first point with
  // M^2 < 0.64, and the flow speed at the crossing, vshock_offset_cells past
  // it, and at that first subsonic point.
  const double vsh_dist = cfg.physics.vshock_offset_cells * cfg.grid.dr;
  std::vector<double> d_sub, v_cross, v_down, v_sub;
  long never = 0;
  for (std::size_t i : s.ts_rows) {
    double path = 0, dsub = -1, vdown = -1, vsub = -1;
    for (std::size_t j = i; j < s.num && std::isfinite(s.r[j]); ++j) {
      if (j > i) {
        const double dx = s.x[j] - s.x[j - 1], dy = s.y[j] - s.y[j - 1], dz = s.z[j] - s.z[j - 1];
        path += std::sqrt(dx * dx + dy * dy + dz * dz);
      }
      const double ex = s.x[j] - s.x[i], ey = s.y[j] - s.y[i], ez = s.z[j] - s.z[i];
      if (vdown < 0 && std::sqrt(ex * ex + ey * ey + ez * ez) >= vsh_dist) vdown = s.v[j];
      const bool sub = cfg.multi_ion() ? s.v[j] < cfg.mi.vshock_vmin_kms * 1e5 : s.m2[j] < 0.64;
      if (dsub < 0 && sub) { dsub = path; vsub = s.v[j]; }
      if (dsub >= 0 && vdown >= 0) break;
    }
    v_cross.push_back(s.v[i] / 1e5);
    if (vdown >= 0) v_down.push_back(vdown / 1e5);
    if (dsub >= 0) { d_sub.push_back(dsub); v_sub.push_back(vsub / 1e5); } else ++never;
  }
  auto pct = [](std::vector<double> v) {
    if (v.empty()) return std::string("none");
    std::sort(v.begin(), v.end());
    auto q = [&](double f) { return v[std::min(v.size() - 1, std::size_t(f * (v.size() - 1) + 0.5))]; };
    char b[96];
    std::snprintf(b, sizeof b, "p10 %.1f, median %.1f, p90 %.1f", q(0.1), q(0.5), q(0.9));
    return std::string(b);
  };
  std::string out = "  TS crossings (regions): " + std::to_string(s.ts_rows.size()) + "\n";
  const std::string crit = cfg.multi_ion() ? "pu3 speed < " + std::to_string(int(cfg.mi.vshock_vmin_kms)) + " km/s" : "M^2 < 0.64";
  out += "  path from the crossing to " + crit + " [AU]: " + pct(d_sub) + " (never on " + std::to_string(never) + ")\n";
  out += "  flow speed at the crossing [km/s]: " + pct(v_cross) + "\n";
  out += "  flow speed at the first point with " + crit + " (v_shock \"first_below\") [km/s]: " + pct(v_sub) + "\n";
  out += "  flow speed " + std::to_string(int(cfg.physics.vshock_offset_cells)) +
         " cells past the crossing (v_shock \"downstream\") [km/s]: " + pct(v_down) + "\n";
  return out;
}

StreamPoints compute_stream_points(const RunConfig& cfg, const StreamlineSet& sl, double r_ref,
                                   const Field3D<std::int16_t>* region) {
  const ColumnMap& m = cfg.inputs.streamlines;
  auto colx = [&](const char* role) {
    int c = m.column(role);
    if (c < 0) throw std::runtime_error(std::string("streamline column map lacks '") + role + "'");
    return c;
  };
  const int cx = colx("x"), cy = colx("y"), cz = colx("z"), cux = colx("ux"), cuy = colx("uy"), cuz = colx("uz"),
            cnh = colx("nH"), cnp = colx("np"), ctp = colx("tp"), cvr = colx("vr"), cvt = colx("vt"),
            cth = colx("theta"), cr = colx("r");
  const int creg = cfg.moscow() ? colx("region") : -1;
  const double vf = 1e5;  // km/s -> cm/s
  const bool moscow = cfg.moscow();
  // Multi-ion streamlines also carry the pu3 fluid: its speed is the flow speed
  // used for the TS and v_shock, and the sound speed uses the SW + pu3
  // temperature.
  const bool mi = cfg.multi_ion();
  const int cnp3 = mi ? colx("n_pu3") : -1, ctp3 = mi ? colx("T_pu3") : -1, cvr3 = mi ? colx("vr_pu3") : -1,
            cvt3 = mi ? colx("vt_pu3") : -1;
  const Grid& g = cfg.grid;

  // ---- rows, with NaN rows for the Tecplot header lines between zones ----
  StreamPoints s;
  std::size_t num = sl.points();
  for (std::size_t z = 1; z < sl.count(); ++z) num += std::size_t(sl.gap_lines[z]);
  s.num = num;
  for (auto* v : {&s.x, &s.y, &s.z, &s.ux, &s.uy, &s.uz, &s.nh, &s.np, &s.tp, &s.vr, &s.vt, &s.theta, &s.phi, &s.r, &s.region,
                  &s.v, &s.cs, &s.m2})
    v->assign(num, NaN);
  if (mi) for (auto* v : {&s.npu3, &s.tpu3, &s.vrpu3, &s.vtpu3, &s.tssw}) v->assign(num, NaN);
  std::size_t row = 0;
  std::vector<std::int32_t> line_of(num, -1);  // streamline of each row (-1: header row)
  s.streamlines = sl.count();
  for (std::size_t z = 0; z < sl.count(); ++z) {
    if (z > 0) row += std::size_t(sl.gap_lines[z]);
    for (std::size_t q = sl.offsets[z]; q < sl.offsets[z + 1]; ++q, ++row) {
      line_of[row] = std::int32_t(z);
      s.x[row] = sl.at(q, cx);
      s.y[row] = sl.at(q, cy);
      s.z[row] = sl.at(q, cz);
      s.ux[row] = sl.at(q, cux) * vf;
      s.uy[row] = sl.at(q, cuy) * vf;
      s.uz[row] = sl.at(q, cuz) * vf;
      s.nh[row] = sl.at(q, cnh);
      s.np[row] = sl.at(q, cnp);
      // standard single-ion model: halved, as the grid files; Moscow and multi-ion: as read
      s.tp[row] = (moscow || mi) ? sl.at(q, ctp) : sl.at(q, ctp) / 2.0;
      if (mi) {
        s.npu3[row] = sl.at(q, cnp3);
        s.tpu3[row] = sl.at(q, ctp3);
        s.vrpu3[row] = sl.at(q, cvr3) * vf;
        s.vtpu3[row] = sl.at(q, cvt3) * vf;
      }
      s.vr[row] = sl.at(q, cvr) * vf;
      s.vt[row] = sl.at(q, cvt) * vf;
      s.theta[row] = sl.at(q, cth);
      s.r[row] = sl.at(q, cr);
      if (creg >= 0) s.region[row] = sl.at(q, creg);
    }
  }

  // ---- point quantities ----
  // (independent per point: parallel)
  const double deg_per_rad = 180.0 / cgs::pi;
#pragma omp parallel for schedule(static)
  for (std::size_t i = 0; i < num; ++i) {
    double ph = std::atan2(s.y[i], s.x[i]) * deg_per_rad;
    if (ph < 0) ph = 360.0 - std::fabs(ph);
    s.phi[i] = ph;
    if (mi) {
      s.v[i] = std::sqrt(s.vrpu3[i] * s.vrpu3[i] + s.vtpu3[i] * s.vtpu3[i]);
      s.tssw[i] = (s.np[i] * s.tp[i] + s.npu3[i] * s.tpu3[i]) / (s.np[i] + s.npu3[i]);
      s.cs[i] = std::sqrt((5.0 * s.tssw[i] * cgs::kb) / (3.0 * cgs::mp));
    } else {
      s.v[i] = std::sqrt(s.vr[i] * s.vr[i] + s.vt[i] * s.vt[i]);
      double t2 = (2.0 * s.tp[i]) * 5.0;  // 5 * (T_p + T_e), with T_e = T_p
      s.cs[i] = std::sqrt((t2 * cgs::kb) / (3.0 * cgs::mp));
    }
    s.m2[i] = (s.v[i] * s.v[i]) / (s.cs[i] * s.cs[i]);
  }

  // ---- path length and optical depth between consecutive points ----
  // Spherical line element with r in cm in every term and the phi step from
  // the phi difference.
  const std::size_t n1 = num ? num - 1 : 0;
  for (auto* v : {&s.ds, &s.tau, &s.ext, &s.ntr0, &s.nref0, &s.nen0, &s.etr0, &s.eref0, &s.een0, &s.em0, &s.kappa_ref0,
                  &s.vshock0, &s.esw0})
    v->assign(n1, 0.0);
  const double deg = cgs::pi / 180.0;
#pragma omp parallel for schedule(static)
  for (std::size_t i = 0; i < n1; ++i) {
    double dr0 = (s.r[i + 1] - s.r[i]) * cgs::au_cm;
    double dth = (s.theta[i + 1] - s.theta[i]) * deg;
    double d = s.phi[i + 1] - s.phi[i];
    if (d > 180) d -= 360;
    if (d < -180) d += 360;
    double dph = d * deg;
    double r1 = s.r[i + 1];
    double term2 = (r1 * r1) * (dth * dth);
    double sn = std::sin(s.theta[i + 1] * deg);
    double term3 = ((r1 * r1) * (sn * sn)) * (dph * dph);
    term2 *= cgs::au_cm * cgs::au_cm;
    term3 *= cgs::au_cm * cgs::au_cm;
    s.ds[i] = std::sqrt((dr0 * dr0 + term2) + term3);
    s.tau[i] = -(s.nh[i] / s.v[i + 1]) * s.ds[i];
  }

  // ---- grid indices of each point: nearest grid point, radius from r_inner ----
  s.r_in.resize(num);
  s.t_in.resize(num);
  s.p_in.resize(num);
#pragma omp parallel for schedule(static)
  for (std::size_t i = 0; i < num; ++i) {
    s.r_in[i] = round_index((s.r[i] - g.r_inner) / g.dr);
    s.t_in[i] = round_index(s.theta[i] / g.dtheta);
    std::int32_t pi = round_index(s.phi[i] / g.dphi);
    s.p_in[i] = (pi == INT32_MIN) ? pi : ((pi % g.nphi) + g.nphi) % g.nphi;
  }

  // ---- termination shock: partition fits, carried downstream ----
  std::vector<double>* fr[10] = {&s.ntr0, &s.nref0, &s.nen0, &s.etr0, &s.eref0, &s.een0, &s.em0, &s.kappa_ref0, &s.vshock0, &s.esw0};
  // Multi-ion fits at the region crossing, r' = ln(x) with x capped at ts_cap:
  //  "compressed" (default): x = ((r - r_TS,V2)*ts_compress + r_TS,V2)/r_TS,V2,
  //   anchored at r_TS,V2 and compressed by ts_compress, for the hybrid
  //   model whose fits cover a narrower range of TS distances;
  //  "ratio": x = r / r_TS,V2, strictly the MI solution's TS distances.
  // Polar blend: beyond 45 deg from the equator the fits are blended linearly
  // toward fixed polar values. It is applied at the TS point only at
  // theta <= 25 or >= 155 deg (between 25 and 155 the fit stays unblended).
  auto mi_set_fits = [&](std::size_t i, double rnorm) {
    const MultiIonConfig& mc = cfg.mi;
    const double th = s.theta[i];
    const double x = mc.ts_scaling == "ratio" ? s.r[i] / rnorm : ((s.r[i] - rnorm) * mc.ts_compress + rnorm) / rnorm;
    double rts = x >= 1.0 ? std::log(x) : 0.0;
    if (x > mc.ts_cap) rts = std::log(mc.ts_cap);
    const double r2 = rts * rts;
    const double polar[10] = {0.165, 0.085, 0.15, 0.72, 0.07, 0.7, 2.0, 18.0, 0.0, 0.04};
    for (int f = 0; f < 10; ++f) {
      if (f == 8) continue;
      const double base = std::exp(kFits[f].a * r2 + kFits[f].b * rts + kFits[f].c);
      double v = base;
      if (th >= 155.0) v = base + (polar[f] - base) * (th - 135.0) / 45.0;
      else if (th <= 25.0) v = base + (polar[f] - base) * (45.0 - th) / 45.0;
      (*fr[f])[i] = v;
    }
  };
  auto set_fits = [&](std::size_t i, double rnorm) {
    if (mi) {  // multi-ion scaling and polar blend (below)
      mi_set_fits(i, rnorm);
      return;
    }
    double rts0 = s.r[i] / rnorm;
    double rts = rts0 >= 1.0 ? std::log(rts0) : 0.0;
    if (rts0 > 1.62) rts = std::log(1.62);
    double r2 = rts * rts;
    for (int f = 0; f < 10; ++f) {
      if (f == 8) continue;
      double arg = ((kFits[f].a * r2) + (kFits[f].b * rts)) + kFits[f].c;
      (*fr[f])[i] = std::exp(arg);
    }
  };
  auto copy_prev = [&](std::size_t i) {
    for (int f = 0; f < 10; ++f) if (f != 8) (*fr[f])[i] = (*fr[f])[i - 1];
    if (s.ntr0[i] > 1 || s.nref0[i] > 1 || s.etr0[i] > 1 || s.eref0[i] > 1)
      throw std::runtime_error("streamline processing: a population fraction exceeds 1 at streamline row " + std::to_string(i));
  };
  const double m064 = 0.64, m081 = 0.81;
  std::vector<double> reg = s.region;
  const std::string& method = cfg.physics.ts_detection;
  s.r_ts_reference = r_ref;
  std::vector<char> line_has_ts(sl.count(), 0);
  auto mark_ts = [&](std::size_t i) {
    ++s.ts_crossings;
    if (line_of[i] >= 0) line_has_ts[line_of[i]] = 1;
  };
  // Region of the cell a point is binned to (0 outside the grid).
  std::vector<std::int16_t> preg;
  if (method == "regions") {
    if (!region) throw std::runtime_error("ts_detection = \"regions\" needs the region array");
    preg.assign(num, 0);
    for (std::size_t i = 0; i < num; ++i) {
      if (!fin(s.r[i]) || s.t_in[i] == INT32_MIN || s.p_in[i] == INT32_MIN) continue;
      long kk = std::lround((s.r[i] - g.r_inner) / g.dr);
      int t = s.t_in[i], pp = ((s.p_in[i] % g.nphi) + g.nphi) % g.nphi;
      if (t < 0 || t >= g.ntheta || kk < 0 || kk >= g.nr) continue;
      preg[i] = (*region)(t, pp, int(kk));
    }
  }
  // Shock-frame speed for the regions method (physics.vshock):
  //  downstream_speed: the flow speed at the first point at least
  //   vshock_offset_cells grid cells (straight-line distance) past the crossing
  //   point i, on the same streamline (its last point if it ends first);
  //  first_below_speed: the flow speed at the first point from i on with
  //   M^2 < 0.64 (the Mach criterion, started from the region crossing), or the
  //   downstream speed if the streamline never gets there.
  const double vsh_dist = cfg.physics.vshock_offset_cells * g.dr;
  auto downstream_speed = [&](std::size_t i) {
    std::size_t j = i;
    while (j + 1 < num && fin(s.r[j + 1])) {
      const double dx = s.x[j] - s.x[i], dy = s.y[j] - s.y[i], dz = s.z[j] - s.z[i];
      if (std::sqrt(dx * dx + dy * dy + dz * dz) >= vsh_dist) break;
      ++j;
    }
    return s.v[j];
  };
  // first point below the subsonic criterion: M^2 < 0.64 (single-ion), or
  // pu3 speed < vshock_vmin_kms (multi-ion; the pu3 flow lags the region
  // boundary by a median 4.9 AU in the MI solution, while the SW fluid slows there)
  const double mi_vmin = cfg.mi.vshock_vmin_kms * 1e5;
  auto below = [&](std::size_t j) { return mi ? s.v[j] < mi_vmin : s.m2[j] < m064; };
  auto first_below_speed = [&](std::size_t i) {
    for (std::size_t j = i; j < num && fin(s.r[j]); ++j)
      if (below(j)) return s.v[j];
    return downstream_speed(i);
  };
  // Extended single-ion model: the DSA injection energy is E_l = 1/2 m_p v10^2,
  // with v10 the flow speed at the last point of the same streamline inside
  // 12 AU, replacing the E_l fit at the TS point; the other fits are
  // unchanged. v10 is reset at the start of every streamline; a streamline
  // with no point inside 12 AU uses its first point, and is counted in
  // StreamPoints::esi_no_inner_point.
  const bool esi = cfg.extended_si();
  double v10 = NaN;
  bool v10_inner = false;
  auto esi_injection = [&](std::size_t i) {
    if (!esi) return;
    if (!v10_inner) ++s.esi_no_inner_point;
    s.een0[i] = 0.5 * cgs::mp * v10 * v10 / cgs::kev_erg;
  };
  for (std::size_t i = 1; i + 1 < num; ++i) {
    const bool okr = fin(s.r[i]) && fin(s.r[i - 1]);
    if (esi) {
      if (fin(s.r[i]) && !fin(s.r[i - 1])) { v10 = s.v[i]; v10_inner = s.r[i] < 12.0; }  // first point of a streamline
      else if (i == 1 && fin(s.r[0])) { v10 = s.v[0]; v10_inner = s.r[0] < 12.0; }
      if (fin(s.r[i]) && s.r[i] < 12.0) { v10 = s.v[i]; v10_inner = true; }
    }
    if (method == "regions") {
      // The TS is where the streamline enters region 3/6 from region 4.
      // The fractions and v_shock (see physics.vshock) are then carried along
      // the rest of the streamline,
      // whatever cells it passes (only region 3/6 cells are binned), until a
      // new 4 -> 3/6 crossing.
      const bool down = preg[i] == 3 || preg[i] == 6;
      if (down && preg[i - 1] == 4 && okr) {
        set_fits(i, r_ref);
        esi_injection(i);
        mark_ts(i);
        s.ts_rows.push_back(i);
        s.vshock0[i] = cfg.physics.vshock == "first_below" ? first_below_speed(i)
                       : cfg.physics.vshock == "crossing"  ? s.v[i]  // version 2.0, for comparison only
                                                          : downstream_speed(i);
      } else if (okr) {
        copy_prev(i);
        s.vshock0[i] = s.vshock0[i - 1];
      }
    } else if (method == "moscow_column") {
      if (reg[i - 1] == 2 && reg[i] < 2 && reg[i] > 1) reg[i] = 2;
      if (reg[i] >= 1 && reg[i] < 2 && reg[i + 1] == 2 && s.r[i] < 200) { set_fits(i, r_ref); mark_ts(i); }
      else if (reg[i] == 2 && okr) copy_prev(i);
      if (reg[i] == 2 && okr) {
        if (s.vshock0[i - 1] == 0.0) s.vshock0[i] = s.v[i];
        else if (s.vshock0[i - 1] > 0.0) s.vshock0[i] = s.vshock0[i - 1];
      }
    } else {  // "mach": first point with M^2 < 0.64 after a supersonic point, at r < 130 AU
      if (s.m2[i - 1] > m064 && s.m2[i] < m064 && s.r[i] < 130) { set_fits(i, r_ref); esi_injection(i); mark_ts(i); }
      else if (s.m2[i] < m064 && okr) copy_prev(i);
      if (s.m2[i] < m064 && okr) {
        if (s.vshock0[i - 1] == 0.0) s.vshock0[i] = s.v[i];
        else if (s.vshock0[i - 1] > 0.0) s.vshock0[i] = s.vshock0[i - 1];
      }
    }
  }
  // ---- multi-ion, injection = "streamline": upstream SW speed at the
  //      region crossing, carried along the rest of the streamline ----
  if (mi && cfg.mi.injection == "streamline") {
    s.u_inj0.assign(num, 0.0);
    for (std::size_t i = 1; i < num; ++i) {
      if (!(fin(s.r[i]) && fin(s.r[i - 1]))) continue;  // first point of a streamline: nothing carried in
      if ((preg[i] == 3 || preg[i] == 6) && preg[i - 1] == 4)
        s.u_inj0[i] = std::sqrt(s.vr[i - 1] * s.vr[i - 1] + s.vt[i - 1] * s.vt[i - 1]);  // SW fluid
      else
        s.u_inj0[i] = s.u_inj0[i - 1];
    }
  }
  // ---- extinction: optical depth summed from the termination shock ----
  for (std::size_t z = 0; z < line_has_ts.size(); ++z) s.streamlines_with_ts += line_has_ts[z];
  for (std::size_t i = 1; i + 1 < num; ++i) {
    if (method == "regions") {
      if (preg[i] == 4) s.tau[i] = 0.0;  // upstream of the TS
    } else if (method == "moscow_column") {
      if (reg[i] >= 1 && reg[i] < 2 && s.r[i] < 200) s.tau[i] = 0.0;
    } else if (s.m2[i] > m081) {
      s.tau[i] = 0.0;
    }
    s.ext[i] = (fin(s.r[i]) && fin(s.r[i - 1])) ? s.tau[i] + s.ext[i - 1] : 0.0;
  }
  return s;
}

// ===========================================================================
// Binning and gap filling
// ===========================================================================
StreamFields bin_stream_points(const RunConfig& cfg, const StreamPoints& s, const Field3D<std::int16_t>& region) {
  const Grid& g = cfg.grid;
  const int nt = g.ntheta, np = g.nphi, nr = g.nr;
  StreamFields o;
  for (auto* f : {&o.cool, &o.vx, &o.vy, &o.vz, &o.ntr, &o.nref, &o.nen, &o.etr, &o.eref, &o.een, &o.em, &o.kappa_ref,
                  &o.vshock, &o.esw, &o.ptr_ct, &o.ptrn_ct})
    f->resize(g, 0.0);
  const int rlower = g.k_of(cfg.r_lower), rupper = g.k_of(cfg.r_upper);
  if (rlower < 0 || rupper >= nr || rlower > rupper) throw std::runtime_error("integration window outside the grid");
  auto in36 = [&](int t, int pp, int k) { int b = region(t, pp, k); return b == 3 || b == 6; };

  Field3D<double>* fr[10] = {&o.ntr, &o.nref, &o.etr, &o.eref, &o.nen, &o.een, &o.em, &o.kappa_ref, &o.vshock, &o.esw};
  const std::vector<double>* fr0[10] = {&s.ntr0, &s.nref0, &s.etr0, &s.eref0, &s.nen0, &s.een0, &s.em0, &s.kappa_ref0,
                                        &s.vshock0, &s.esw0};

  // ---- average the points in each region 3/6 cell of the window ----
  std::vector<std::vector<std::size_t>> bucket(g.size());
  for (std::size_t q = 0; q < s.num; ++q) {
    int t = s.t_in[q], pp = s.p_in[q], k = s.r_in[q];
    if (t >= 0 && t < nt && pp >= 0 && pp < np && k >= rlower && k <= rupper) bucket[g.index(t, pp, k)].push_back(q);
  }
  // ext and the fractions have num-1 entries: the last point uses entry num-2.
  const std::size_t last = s.num >= 2 ? s.num - 2 : 0;
  auto clip = [&](std::size_t q) { return q > last ? last : q; };
  for (int t = 0; t < nt; ++t)
    for (int pp = 0; pp < np; ++pp)
      for (int k = rlower; k <= rupper; ++k) {
        if (!in36(t, pp, k)) continue;
        const auto& b = bucket[g.index(t, pp, k)];
        if (b.empty()) continue;
        double sum = 0, sx = 0, sy = 0, sz = 0; int n = 0;
        for (auto q : b) {
          double e = s.ext[clip(q)];
          if (!fin(e)) continue;
          sum += e;
          sx += s.ux[q]; sy += s.uy[q]; sz += s.uz[q];
          ++n;
        }
        if (n > 0) {
          o.ptr_ct(t, pp, k) = n;
          o.cool(t, pp, k) = sum / n;
          o.vx(t, pp, k) = sx / n;
          o.vy(t, pp, k) = sy / n;
          o.vz(t, pp, k) = sz / n;
        }
        double fs[10] = {0}; int nn = 0;
        for (auto q : b) {
          double v = s.ntr0[clip(q)];
          if (!fin(v) || v == 0.0) continue;
          for (int f = 0; f < 10; ++f) fs[f] += (*fr0[f])[clip(q)];
          ++nn;
        }
        if (nn > 0) {
          o.ptrn_ct(t, pp, k) = nn;
          for (int f = 0; f < 10; ++f) (*fr[f])(t, pp, k) = fs[f] / nn;
        }
      }
  return o;
}

namespace {

// Region 3/6 cells inside the integration window that are still empty.
void count_unfilled(const RunConfig& cfg, StreamFields& o, const Field3D<std::int16_t>& region) {
  const Grid& g = cfg.grid;
  const int rlower = g.k_of(cfg.r_lower), rupper = g.k_of(cfg.r_upper);
  o.cells = o.unfilled_cool = o.unfilled_fractions = 0;
  for (int t = 0; t < g.ntheta; ++t)
    for (int pp = 0; pp < g.nphi; ++pp)
      for (int k = rlower; k <= rupper; ++k) {
        int b = region(t, pp, k);
        if (b != 3 && b != 6) continue;
        ++o.cells;
        if (!(o.cool(t, pp, k) != 0.0 && fin(o.cool(t, pp, k)))) ++o.unfilled_cool;
        if (!(o.ntr(t, pp, k) != 0.0 && fin(o.ntr(t, pp, k)))) ++o.unfilled_fractions;
      }
}

// Every pass reads the previous pass's values (Jacobi), so the
// result does not depend on the visiting order. All quantities wrap in phi;
// theta index 0 is the pole, a single point whose neighbours are the whole of
// row 1; the last theta row is one-sided. Empty cells take the mean of their
// non-empty neighbours: the 8 angular neighbours at the same radius for the
// extinction (cool), or, if none of those has a value, the neighbours at
// r -/+ 1 (so cells in rings without streamline data, e.g. near the
// heliopause, are filled too); and the 26 neighbours including r -/+ 1 for the
// population fractions (all ten divided by the count).
// Passes continue until no cell changes or gap_fill_max_passes is reached.
void fill_jacobi(const RunConfig& cfg, StreamFields& o, const Field3D<std::int16_t>& region,
                const StreamPassCallback& on_pass) {
  const Grid& g = cfg.grid;
  const int nt = g.ntheta, np = g.nphi;
  const int rlower = g.k_of(cfg.r_lower), rupper = g.k_of(cfg.r_upper);
  auto in36 = [&](int t, int pp, int k) { int b = region(t, pp, k); return b == 3 || b == 6; };
  auto has = [](double v) { return v != 0.0 && fin(v); };
  Field3D<double>* fr[10] = {&o.ntr, &o.nref, &o.etr, &o.eref, &o.nen, &o.een, &o.em, &o.kappa_ref, &o.vshock, &o.esw};
  // cells outside region 3/6 hold no value
  for (int t = 0; t < nt; ++t)
    for (int pp = 0; pp < np; ++pp)
      for (int k = rlower; k <= rupper; ++k)
        if (!in36(t, pp, k)) o.cool(t, pp, k) = 0.0;

  // angular neighbours (theta, phi) of a line of sight
  auto neighbours = [&](int t, int pp, std::vector<std::pair<int, int>>& out) {
    out.clear();
    if (t == 0) {
      for (int p2 = 0; p2 < np; ++p2) if (nt > 1) out.push_back({1, p2});
      return;
    }
    for (int dt = -1; dt <= 1; ++dt) {
      int tt = t + dt;
      if (tt < 0 || tt >= nt) continue;
      if (tt == 0) { out.push_back({0, pp}); continue; }  // the pole, once
      for (int dp = -1; dp <= 1; ++dp) {
        if (dt == 0 && dp == 0) continue;
        out.push_back({tt, ((pp + dp) % np + np) % np});
      }
    }
  };

  const int max_passes = cfg.physics.gap_fill_max_passes;
  int pass = 1;
  for (; pass <= max_passes; ++pass) {
    if (on_pass) on_pass(pass - 1, o);
    const Field3D<double> cool0 = o.cool;
    std::vector<Field3D<double>> fr0;
    fr0.reserve(10);
    for (auto* f : fr) fr0.push_back(*f);
    long changed = 0;
    // Jacobi pass: every cell reads only the copies of the previous pass and
    // writes only itself, so lines of sight are filled in parallel with the
    // same result for any thread count.
#pragma omp parallel for schedule(dynamic, 8) reduction(+ : changed)
    for (long los = 0; los < long(nt) * np; ++los) {
        const int t = int(los / np), pp = int(los % np);
        std::vector<std::pair<int, int>> nb;
        neighbours(t, pp, nb);
        for (int k = rlower; k <= rupper; ++k) {
          if (!in36(t, pp, k)) continue;
          if (!has(cool0(t, pp, k))) {
            double sum = 0; int n = 0;
            for (auto [tt, p2] : nb) if (has(cool0(tt, p2, k))) { sum += cool0(tt, p2, k); ++n; }
            if (!n) {  // no data at this radius nearby (e.g. near the heliopause): use r -/+ 1
              for (int kk : {k - 1, k + 1}) {
                if (kk < rlower || kk > rupper) continue;
                if (has(cool0(t, pp, kk))) { sum += cool0(t, pp, kk); ++n; }
                for (auto [tt, p2] : nb) if (has(cool0(tt, p2, kk))) { sum += cool0(tt, p2, kk); ++n; }
              }
            }
            if (n) { o.cool(t, pp, k) = sum / n; ++changed; }
          }
          if (!has(fr0[0](t, pp, k))) {
            double fs[10] = {0}; int n = 0;
            for (int kk = std::max(k - 1, rlower); kk <= std::min(k + 1, rupper); ++kk) {
              auto take = [&](int tt, int p2) {
                if (!has(fr0[0](tt, p2, kk))) return;
                for (int f = 0; f < 10; ++f) fs[f] += fr0[f](tt, p2, kk);
                ++n;
              };
              for (auto [tt, p2] : nb) take(tt, p2);
              if (kk != k) take(t, pp);
            }
            if (n) {
              for (int f = 0; f < 10; ++f) (*fr[f])(t, pp, k) = fs[f] / n;
              ++changed;
            }
          }
        }
      }
    if (!changed) break;
  }
  o.passes = std::min(pass, max_passes);
  for (int k = 0; k < g.nr; ++k) {  // diagnostic velocities at the pole: those of phi index 1
    o.vx(0, 0, k) = o.vx(0, 1, k);
    o.vy(0, 0, k) = o.vy(0, 1, k);
    o.vz(0, 0, k) = o.vz(0, 1, k);
  }
}

}  // namespace

void fill_stream_fields(const RunConfig& cfg, StreamFields& o, const Field3D<std::int16_t>& region,
                        const StreamPassCallback& on_pass) {
  fill_jacobi(cfg, o, region, on_pass);
  count_unfilled(cfg, o, region);
  if (on_pass) on_pass(-1, o);
}

namespace {
// Multi-ion, injection = "streamline": average u_inj0 over the points in
// each region 3/6 cell that have a value, then fill empty cells from their
// neighbours (Jacobi passes: the 8 angular neighbours at the same radius, else
// also those at r -/+ 1; the pole is one point whose neighbours are the whole
// of theta row 1), as the other cell fields.
void bin_fill_injection(const RunConfig& cfg, const StreamPoints& s, const Field3D<std::int16_t>& region, StreamFields& o) {
  const Grid& g = cfg.grid;
  const int nt = g.ntheta, np = g.nphi;
  const int rlower = g.k_of(cfg.r_lower), rupper = g.k_of(cfg.r_upper);
  auto in36 = [&](int t, int pp, int k) { int b = region(t, pp, k); return b == 3 || b == 6; };
  o.u_inj.resize(g, 0.0);
  o.has_u_inj = true;
  std::vector<double> sum(g.size(), 0.0);
  std::vector<int> cnt(g.size(), 0);
  for (std::size_t q = 0; q < s.num && q < s.u_inj0.size(); ++q) {
    const double u = s.u_inj0[q];
    if (!(u > 0)) continue;
    const int t = s.t_in[q], pp = s.p_in[q], k = s.r_in[q];
    if (t < 0 || t >= nt || pp < 0 || pp >= np || k < rlower || k > rupper || !in36(t, pp, k)) continue;
    sum[g.index(t, pp, k)] += u;
    ++cnt[g.index(t, pp, k)];
  }
  for (std::size_t c = 0; c < g.size(); ++c)
    if (cnt[c]) o.u_inj.raw()[c] = sum[c] / cnt[c];
  auto neighbours = [&](int t, int pp, std::vector<std::pair<int, int>>& out) {
    out.clear();
    if (t == 0) { for (int p2 = 0; p2 < np; ++p2) if (nt > 1) out.push_back({1, p2}); return; }
    for (int dt = -1; dt <= 1; ++dt) {
      const int tt = t + dt;
      if (tt < 0 || tt >= nt) continue;
      if (tt == 0) { out.push_back({0, pp}); continue; }
      for (int dp = -1; dp <= 1; ++dp) {
        if (dt == 0 && dp == 0) continue;
        out.push_back({tt, ((pp + dp) % np + np) % np});
      }
    }
  };
  std::vector<std::pair<int, int>> nb;
  for (int pass = 0; pass < cfg.physics.gap_fill_max_passes; ++pass) {
    const Field3D<double> prev = o.u_inj;
    long changed = 0;
    for (int t = 0; t < nt; ++t)
      for (int pp = 0; pp < np; ++pp) {
        neighbours(t, pp, nb);
        for (int k = rlower; k <= rupper; ++k) {
          if (!in36(t, pp, k) || prev(t, pp, k) > 0) continue;
          double sm = 0; int n = 0;
          for (auto [tt, p2] : nb) if (prev(tt, p2, k) > 0) { sm += prev(tt, p2, k); ++n; }
          if (!n)
            for (int kk : {k - 1, k + 1}) {
              if (kk < rlower || kk > rupper) continue;
              if (prev(t, pp, kk) > 0) { sm += prev(t, pp, kk); ++n; }
              for (auto [tt, p2] : nb) if (prev(tt, p2, kk) > 0) { sm += prev(tt, p2, kk); ++n; }
            }
          if (n) { o.u_inj(t, pp, k) = sm / n; ++changed; }
        }
      }
    if (!changed) break;
  }
  o.unfilled_u_inj = 0;
  for (int t = 0; t < nt; ++t)
    for (int pp = 0; pp < np; ++pp)
      for (int k = rlower; k <= rupper; ++k)
        if (in36(t, pp, k) && !(o.u_inj(t, pp, k) > 0)) ++o.unfilled_u_inj;
}
}  // namespace

StreamFields bin_and_fill(const RunConfig& cfg, const StreamPoints& s, const Field3D<std::int16_t>& region,
                          const StreamPassCallback& on_pass) {
  StreamFields o = bin_stream_points(cfg, s, region);
  if (cfg.multi_ion() && cfg.mi.injection == "streamline") bin_fill_injection(cfg, s, region, o);
  fill_stream_fields(cfg, o, region, on_pass);
  return o;
}

}  // namespace ena
