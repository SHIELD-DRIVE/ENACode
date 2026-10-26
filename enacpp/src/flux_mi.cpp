// Flux of the multi-ion model (Kornbleuth+2026).
//
// Three ion fluids: SW; pu3 (PUIs picked up in the supersonic wind), split at
// the termination shock into transmitted, reflected and DSA-accelerated parts
// with the partition fits carried along the streamlines; pu2 (PUIs created in
// the heliosheath). In the reconnection region ("INCA belt") the pu2 and SW
// slots hold two kappa populations instead. Unlike the single-ion model, the
// densities come straight from the MHD fluids (no streamline extinction).
//
// Notes:
//  - double precision; integration stops where the line of sight first leaves
//    the heliosheath, or at an undetermined (region 5) cell;
//  - q = 5 for the accelerated PUIs everywhere;
//  - E_l from the upstream SW speed carried along the streamline (default);
//  - fits and v_shock from the regions-based TS (stream_processing);
//  - the kappa of the pu2 slot is set in every cell;
//  - E_SW from the fits when physics.fitted_esw (as the single-ion model).

#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>

#include "ena/flux.hpp"
#include "ena/parallel.hpp"
#include "ena/regions.hpp"

namespace ena {

namespace {

using namespace cgs;
constexpr int kMiComponents = 6;  // transmitted, reflected, accelerated, pu2, solar_wind, total

[[noreturn]] void mi_error(const std::string& what, int t, int pp, int k, double energy) {
  std::ostringstream o;
  o << what << " at theta index " << t << ", phi index " << pp << ", r index " << k << ", energy " << energy << " keV";
  throw PhysicsError(o.str());
}

}  // namespace

EnergyResult make_flux_mi(const FluxContext& fc, double E, int band, bool diag) {
  (void)band;
  const RunConfig& cfg = fc.cfg;
  const Grid& g = cfg.grid;
  const PhysicsConfig& phys = cfg.physics;
  const MultiIonConfig& mi = cfg.mi;
  const Snapshot& snap = fc.snap;
  const StreamFields& st = fc.stream;
  if (snap.plasma.size() != 3 || snap.neutrals.empty()) throw std::runtime_error("multi-ion make_flux needs 3 ion fluids");
  if (phys.shock_shift != "none") throw std::runtime_error("physics.shock_shift is not implemented for the multi-ion model");
  const PlasmaFluid &sw = snap.plasma[0], &pu3 = snap.plasma[1], &pu2 = snap.plasma[2];
  const NeutralPopulation& nh = snap.neutrals[0];
  const int nr = g.nr;
  auto wrapk = [&](int k) { return k < 0 ? k + nr : k; };
  PhysicsSettings ps;
  ps.survival = phys.survival;
  ps.dr_au = g.dr;

  EnergyResult R;
  R.energy = E;
  R.vena = energy_velocity(E);
  R.sig_cx = xsection(E);
  R.flux.assign(kMiComponents, Field2D<double>(g, 0.0));
  R.has_diagnostics = diag;
  R.multi_ion = true;

  // Streamline extinction: only for the cooling-length diagnostic (lc); the
  // multi-ion densities are taken from the MHD fluids directly.
  Field3D<double> ext(g);
#pragma omp parallel for schedule(static)
  for (std::size_t i = 0; i < g.size(); ++i)
    ext.raw()[i] = phys.extinction ? std::exp((st.cool.raw()[i] * R.vena) * R.sig_cx) : 1.0;
  BetField bet = ionization_rate(fc, E);

  const std::size_t ncell = g.size();
  if (diag) {
    // flux_test1..3: transmitted, reflected, accelerated per cell; flux_test_sw: pu2 + SW per cell
    for (auto* f : {&R.v_plas, &R.flux_test1, &R.flux_test2, &R.flux_test3, &R.flux_test_sw}) f->resize(g, 0.0);
    R.shock_strength.resize(g, 0.0);
    R.pops_n.assign(ncell * 5, 0.0);
    R.pops_t.assign(ncell * 5, 0.0);
    R.pops_kappa.assign(ncell * 5, 0.0);
    R.zank_par.assign(ncell * 4, 0.0);
    R.cell_par.assign(ncell * 4, 0.0);
    R.belt_mask.resize(g, 0);
  }
  R.lc.resize(g, 0.f);

  const int rlower = g.k_of(cfg.r_lower);
  const double e_target = 1.0 / std::exp(1.0);
  const double kap = phys.kappa, kap_pui = phys.kappa_pui, kap_ref = phys.kappa_ref, kap_oh = phys.kappa_oh;
  const double c1em10 = 1e-10;
  const bool inj_streamline = mi.injection == "streamline";
  if (inj_streamline && !st.has_u_inj)
    throw std::runtime_error("multi_ion.injection = \"streamline\" needs the streamline processing (not --stream-fields)");
  if (mi.injection_k < 0 || mi.injection_k >= nr) throw std::runtime_error("multi_ion.injection_k outside the grid");
  auto up = [&](int t, int pp, int k) {  // SW bulk speed (cm/s)
    k = wrapk(k);
    return std::sqrt(sw.vr(t, pp, k) * sw.vr(t, pp, k) + sw.vt(t, pp, k) * sw.vt(t, pp, k));
  };

  // Lines of sight are independent work items: each writes only its own
  // pixels and cells, so the result does not depend on the thread count.
  ParallelErrors perr;
  const long nlos = long(g.ntheta) * g.nphi;
#pragma omp parallel for schedule(dynamic, 4)
  for (long los = 0; los < nlos; ++los) {
    const int t = int(los / g.nphi), pp = int(los % g.nphi);
    perr.capture(std::size_t(los), [&] {
        {  // cooling length: index where ext is closest to 1/e
          double best = std::numeric_limits<double>::infinity();
          int bi = 0;
          for (int k = 0; k < nr; ++k) {
            const double d = std::fabs(e_target - ext(t, pp, k));
            if (d < best) { best = d; bi = k; }
          }
          R.lc(t, pp) = float(bi);
        }
        if (fc.cells) fc.cells->clear_los(g, t, pp);
        double acc[5] = {0, 0, 0, 0, 0};
        // Integration stops where the line of sight first leaves the heliosheath.
        const int k_end = heliosheath_exit(cfg, snap.region, t, pp) - 1;
        for (int k = rlower; k <= k_end; ++k) {
          const int b = snap.region(t, pp, k);
          const int prev = snap.region(t, pp, wrapk(k - 1));
          if (b != 3 && b != 6) continue;           // only heliosheath cells produce ENAs (radial_range = 0)
          if (prev == 1 || prev == 2) continue;     // cannot occur before the first exit; kept for safety

          // ---- plasma: three fluids ----
          const double nsw0 = sw.n(t, pp, k), tsw0 = sw.T(t, pp, k), ur = sw.vr(t, pp, k), ut = sw.vt(t, pp, k);
          const double npu3 = pu3.n(t, pp, k), tpu3 = pu3.T(t, pp, k), urpu3 = pu3.vr(t, pp, k), utpu3 = pu3.vt(t, pp, k);
          const double npu2 = pu2.n(t, pp, k), tpu2 = pu2.T(t, pp, k), urpu2 = pu2.vr(t, pp, k), utpu2 = pu2.vt(t, pp, k);
          const double np = nsw0 + npu3 + npu2;                                   // total ion density
          double tp = (tsw0 * nsw0 + tpu3 * npu3 + tpu2 * npu2) / np;             // density-weighted temperature
          const double nn = nh.n(t, pp, k);
          const double urh = nh.vr(t, pp, k) * 1e5;                               // neutral H radial speed, km/s -> cm/s
          // Alfven energy of the pu2 fluid relative to its thermal energy (reconnection-region test)
          const double bx = snap.Bx(t, pp, k), by = snap.By(t, pp, k), bz = snap.Bz(t, pp, k);
          const double btot = std::sqrt(bx * bx + by * by + bz * bz);              // nT
          const double va = (btot * 1e-5) / std::sqrt(4.0 * pi * npu2 * mp);     // cm/s (B in G)
          const double mva2_kev = mp * va * va / kev_erg;
          const double ratio = mva2_kev / (tpu2 * kb / kev_erg);

          // ---- pu3 partition at the termination shock (fits carried along the streamline) ----
          // A cell the gap fill left empty has no partition (all fractions 0). "skip" drops the whole
          // cell, its other populations included.
          if (!(st.ntr(t, pp, k) != 0.0 && std::isfinite(st.ntr(t, pp, k)))) {
            if (phys.unfilled_cells == "skip") continue;
            mi_error("heliosheath cell without streamline data (the gap fill left it empty, see the Gap fill line); set [physics] unfilled_cells = \"skip\" to give such cells no ENAs", t, pp, k, E);
          }
          const double ntr1 = st.ntr(t, pp, k), nref1 = st.nref(t, pp, k), nen1 = st.nen(t, pp, k),
                       etr1 = st.etr(t, pp, k), eref1 = st.eref(t, pp, k), een1 = st.een(t, pp, k),
                       em1 = st.em(t, pp, k), esw1 = st.esw(t, pp, k);
          const double kappaEn = st.kappa_ref(t, pp, k);  // spectral index delta of the DSA spectrum
          const double vshock = st.vshock(t, pp, k);
          const double npu3_frac = ntr1 + nref1;
          const double ntr_frac = ntr1 / npu3_frac, nref_frac = nref1 * (1.0 - nen1) / npu3_frac,
                       nen_frac = nref1 * nen1 / npu3_frac;
          const double nsw_frac = 1.0 - (ntr1 + nref1);
          const double esw_used = phys.fitted_esw ? esw1 : 0.04;
          const double ttr_frac = etr1 / ntr1, tref_frac = eref1 / (nref1 * (1.0 - nen1)), tsw_frac = esw_used / nsw_frac;
          // DSA injection energy E_l = 0.5 m_p u^2; E_m scaled by the same factor (fit ratio kept)
          const double u = inj_streamline ? st.u_inj(t, pp, k) : up(t, pp, mi.injection_k);
          const double een = 0.5 * mp * u * u / kev_erg;
          const double em = een1 > 0 ? em1 / een1 * een : 0.0;

          double np1 = ntr_frac * npu3, np2 = nref_frac * npu3, np3 = nen_frac * npu3, np4 = npu2, nsw = nsw0;
          double tp1 = tp * ttr_frac, tp2 = tp * tref_frac, tp4 = tpu2, tsw = tp * tsw_frac;
          for (double x : {tp1, np1, tp2, np2, np3, np4, tp4, tsw, nsw})
            if (!(x > 0)) mi_error("non-positive or NaN population density/temperature", t, pp, k, E);
          double kappaIHS = kap_oh;  // pu2 slot
          double kappa0 = b == 3 ? kap : kap_oh;  // SW slot

          // ---- reconnection region ("INCA belt") ----
          // Grid coordinates: at the pole (theta = 0) x = 0, so the polar row is
          // excluded.
          bool belt = false;
          if (mi.belt && b == 3) {
            const double r = g.r(k), th = g.theta(t), ph = g.phi(pp);
            const double x0 = r * std::cos(ph * pi / 180.0) * std::sin(th * pi / 180.0);
            const bool r_ok = mi.belt_r_max <= 0 || r <= mi.belt_r_max;
            const bool phi_ok = mi.belt_phi_min >= mi.belt_phi_max || ph < mi.belt_phi_min || ph > mi.belt_phi_max;
            if (x0 > 0.0 && r_ok && phi_ok && ratio >= mi.belt_ratio_min && mva2_kev < mi.belt_mva2_max) {
              belt = true;
              nsw = mi.belt_accel_frac * np;               // accelerated ions (kappa belt_kappa_sw)
              np4 = mi.belt_pu2_frac * np;                 // non-accelerated ions (kappa belt_kappa_pu2)
              tp = tpu2;
              tp4 = mi.belt_t_pu2_factor * tp;
              tsw = mi.belt_t_sw_factor * tp;
              kappaIHS = mi.belt_kappa_pu2;
              kappa0 = mi.belt_kappa_sw;
            }
          }

          // ---- ENA flux from each population ----
          const double betk = bet(t, pp, k);
          auto cdf = [&](double T, double vr_, double vt_, double n, double kappa) {
            FluxInputs in{E, T, vr_, vt_, nn, n, kappa, betk};
            return calc_diff_flux_mi(in, urh, ps);
          };
          auto zank = [&](double vr_, double vt_) {
            FluxInputs in{E, 0.0, vr_, vt_, nn, np3, 0.0, betk};
            return calc_diff_flux_zank(in, kappaEn, een, em, vshock, ps);  // q = 5
          };
          // In the inner heliosheath each population moves with its own fluid
          // (pu3 parts with pu3, pu2 with pu2, SW with SW); in region 6 (open
          // heliosheath) all use the SW velocity.
          const bool ihs = b == 3;
          const double j0 = cdf(tp1, ihs ? urpu3 : ur, ihs ? utpu3 : ut, np1, kap_pui);
          const double j1 = cdf(tp2, ihs ? urpu3 : ur, ihs ? utpu3 : ut, np2, kap_ref);
          const double j2 = zank(ihs ? urpu3 : ur, ihs ? utpu3 : ut);
          const double j3 = cdf(tp4, ihs ? urpu2 : ur, ihs ? utpu2 : ut, np4, kappaIHS);
          const double j4 = cdf(tsw, ur, ut, nsw, kappa0);
          acc[0] += j0; acc[1] += j1; acc[2] += j2; acc[3] += j3; acc[4] += j4;
          if (fc.cells) {
            const double jj[5] = {j0, j1, j2, j3, j4};
            for (int m = 0; m < 5; ++m) fc.cells->set(g, m, t, pp, k, jj[m]);
          }

          if (diag) {
            const std::size_t c = g.index(t, pp, k);
            const double pn[5] = {np1, np2, np3, np4, nsw}, pt[5] = {tp1, tp2, c1em10, tp4, tsw},
                         pk[5] = {kap_pui, kap_ref, kappaEn, kappaIHS, kappa0};
            for (int q = 0; q < 5; ++q) { R.pops_n[c * 5 + q] = pn[q]; R.pops_t[c * 5 + q] = pt[q]; R.pops_kappa[c * 5 + q] = pk[q]; }
            const double zp[4] = {een, em, vshock, 5.0}, cp[4] = {np, tp, ratio, mva2_kev};
            for (int q = 0; q < 4; ++q) { R.zank_par[c * 4 + q] = zp[q]; R.cell_par[c * 4 + q] = cp[q]; }
            if (belt) R.belt_mask(t, pp, k) = 1;
            R.v_plas(t, pp, k) = plasma_frame(E, ur, ut);
            R.flux_test1(t, pp, k) = j0;
            R.flux_test2(t, pp, k) = j1;
            R.flux_test3(t, pp, k) = j2;
            R.flux_test_sw(t, pp, k) = j3 + j4;
            if (b == 3 && prev == 4) R.shock_strength(t, pp) = up(t, pp, k - 2) / up(t, pp, k);
          }
        }
        for (int m = 0; m < 5; ++m) R.flux[m](t, pp) = acc[m];
        R.flux[5](t, pp) = acc[0] + acc[1] + acc[2] + acc[3] + acc[4];
    });
  }
  perr.rethrow();
  if (diag) {
    R.ext = std::move(ext);
    R.bet = std::move(bet);
  }
  return R;
}

}  // namespace ena
