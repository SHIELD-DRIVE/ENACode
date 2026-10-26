// Flux of the extended single-ion model.
//
// The single-ion termination-shock partition (Kornbleuth+2023: transmitted,
// reflected and DSA-accelerated PUIs and the solar wind, fits carried along the
// streamlines, densities reduced by the streamline extinction), with:
//  - four neutral populations (41-column neutral file); ENA production uses
//    the total H density;
//  - the charge-exchange cross section at the thermally averaged speed of the
//    parent ion relative to the neutral bulk flow (vrel, Heerikhuisen+2006),
//    for every population except the accelerated one;
//  - component 3, PUIs created in the heliosheath: density n_p - n1 - n2 - n3
//    (the ions removed by the extinction plus the SW share),
//    temperature 1/2 m_p |u_p - u_H|^2 / k from the ion-neutral relative flow,
//    kappa = kappa_oh;
//  - components 5 and 6 (Ne4, Neu): neutral populations 1 and 2 themselves, a
//    Maxwellian in their bulk frame at the inner boundary, density x 0.639307
//    (photoionization inside it);
//  - the DSA injection energy E_l = 1/2 m_p v^2 with v the flow speed inside
//    12 AU on the streamline (stream_processing).
//
// Notes:
//  - region-6 cells get the region-3 populations; their kappas and an
//    accelerated fraction of the ions are set in [region6] (defaults: the
//    region-3 treatment with the SW at kappa_oh);
//  - survival losses are counted in the heliosheath, regions 3 and 6
//    (ionization_rate), not in region 2;
//  - neutral population 1/2 velocities converted km/s -> cm/s (x 1e5) like
//    the total H (checked on the input file);
//  - no reconnection-region ("INCA belt") test populations;
//  - as in the single-ion model: double precision, integration stops at the
//    first exit from the heliosheath (region-4 pockets skipped with a
//    warning), fits from the regions-based TS, v_shock "first_below", fitted
//    E_sw, complete gap fill, sub-band centres.

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
constexpr int kEsiComponents = 9;  // transmitted, reflected, accelerated, heliosheath_pui, solar_wind, ne4, neu,
                                    // region6_accelerated, total
constexpr int kEsiSlots = 6;        // populations in the C4 dump: the five above (without ne4/neu) + region-6 accelerated

[[noreturn]] void esi_error(const std::string& what, int t, int pp, int k, double energy) {
  std::ostringstream o;
  o << what << " at theta index " << t << ", phi index " << pp << ", r index " << k << ", energy " << energy << " keV";
  throw PhysicsError(o.str());
}

}  // namespace

EnergyResult make_flux_esi(const FluxContext& fc, double E, int band, bool diag) {
  (void)band;
  const RunConfig& cfg = fc.cfg;
  const Grid& g = cfg.grid;
  const PhysicsConfig& phys = cfg.physics;
  const Snapshot& snap = fc.snap;
  const StreamFields& st = fc.stream;
  if (snap.plasma.size() != 1 || snap.neutrals.size() != 5)
    throw std::runtime_error("extended single-ion make_flux needs one ion fluid and five neutral populations");
  if (phys.shock_shift != "none") throw std::runtime_error("physics.shock_shift is not implemented for the extended single-ion model");
  const PlasmaFluid& pl = snap.plasma[0];
  const NeutralPopulation &nh = snap.neutrals[0], &nh1 = snap.neutrals[1], &nh2 = snap.neutrals[2];
  const int nr = g.nr;
  auto wrapk = [&](int k) { return k < 0 ? k + nr : k; };
  PhysicsSettings ps;
  ps.survival = phys.survival;
  ps.dr_au = g.dr;

  EnergyResult R;
  R.energy = E;
  R.vena = energy_velocity(E);
  R.sig_cx = xsection(E);
  R.flux.assign(kEsiComponents, Field2D<double>(g, 0.0));
  R.has_diagnostics = diag;
  R.extended_si = true;

  // Energy-dependent extinction along streamlines (Kornbleuth+2021 Eq. 4).
  Field3D<double> ext(g);
#pragma omp parallel for schedule(static)
  for (std::size_t i = 0; i < g.size(); ++i)
    ext.raw()[i] = phys.extinction ? std::exp((st.cool.raw()[i] * R.vena) * R.sig_cx) : 1.0;
  BetField bet = ionization_rate(fc, E);

  const std::size_t ncell = g.size();
  if (diag) {
    // flux_test1..3: transmitted, reflected, accelerated per cell; flux_test_sw: heliosheath PUIs + SW per cell
    for (auto* f : {&R.v_plas, &R.flux_test1, &R.flux_test2, &R.flux_test3, &R.flux_test_sw}) f->resize(g, 0.0);
    R.shock_strength.resize(g, 0.0);
    R.pops_n.assign(ncell * kEsiSlots, 0.0);
    R.pops_t.assign(ncell * kEsiSlots, 0.0);
    R.pops_kappa.assign(ncell * kEsiSlots, 0.0);
    R.zank_par.assign(ncell * 3, 0.0);
    R.cell_par.assign(ncell * 4, 0.0);
  }
  R.lc.resize(g, 0.f);

  const int rlower = g.k_of(cfg.r_lower);
  const double e_target = 1.0 / std::exp(1.0);
  const double kap = phys.kappa, kap_pui = phys.kappa_pui, kap_ref = phys.kappa_ref;
  const double kap_hs = phys.kappa_heliosheath_pui;  // heliosheath PUIs in region 3 (default: kappa_oh)
  const Region6Config& reg6 = cfg.region6;
  const double acc_f = reg6.accelerated_fraction, acc_e = reg6.accelerated_energy_share;
  const bool acc_pooled = reg6.accelerated_from == "pooled";
  const double c1em10 = 1e-10;
  const double photo = 0.639307;  // photoionization of the neutrals inside the inner boundary
  auto up = [&](int t, int pp, int k) {  // flow speed (cm/s)
    k = wrapk(k);
    return std::sqrt(pl.vr(t, pp, k) * pl.vr(t, pp, k) + pl.vt(t, pp, k) * pl.vt(t, pp, k));
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
        double acc[8] = {0, 0, 0, 0, 0, 0, 0, 0};

        // ---- neutral populations 1 and 2 at the inner boundary (k = 0), velocities km/s -> cm/s ----
        acc[5] = calc_diff_flux_neutrals(E, nh1.T(t, pp, 0), nh1.vr(t, pp, 0) * 1e5, nh1.vt(t, pp, 0) * 1e5,
                                         nh1.n(t, pp, 0) * photo);
        acc[6] = calc_diff_flux_neutrals(E, nh2.T(t, pp, 0), nh2.vr(t, pp, 0) * 1e5, nh2.vt(t, pp, 0) * 1e5,
                                         nh2.n(t, pp, 0) * photo);
        if (fc.cells) {  // the neutral components come from the inner boundary (k = 0)
          fc.cells->set(g, 5, t, pp, 0, acc[5]);
          fc.cells->set(g, 6, t, pp, 0, acc[6]);
        }

        // Integration stops where the line of sight first leaves the heliosheath.
        const int k_end = heliosheath_exit(cfg, snap.region, t, pp) - 1;
        for (int k = rlower; k <= k_end; ++k) {
          const int b = snap.region(t, pp, k);
          const int prev = snap.region(t, pp, wrapk(k - 1));
          if (b != 3 && b != 6) continue;        // only heliosheath cells produce ENAs (radial_range = 0)
          if (prev == 1 || prev == 2) continue;  // cannot occur before the first exit; kept for safety

          double& e = ext(t, pp, k);
          if (std::isnan(e) || e == 0.0) e = c1em10;
          const double np = pl.n(t, pp, k), ur = pl.vr(t, pp, k), ut = pl.vt(t, pp, k);
          double tp = pl.T(t, pp, k);
          const double nn = nh.n(t, pp, k);
          const double urh = nh.vr(t, pp, k) * 1e5, uth = nh.vt(t, pp, k) * 1e5;  // total H bulk velocity, km/s -> cm/s

          // ---- termination-shock partition (Kornbleuth+2023), carried along the streamline ----
          // A cell the gap fill left empty has no partition (all fractions 0).
          if (!(st.ntr(t, pp, k) != 0.0 && std::isfinite(st.ntr(t, pp, k)))) {
            if (phys.unfilled_cells == "skip") continue;
            esi_error("heliosheath cell without streamline data (the gap fill left it empty, see the Gap fill line); set [physics] unfilled_cells = \"skip\" to give such cells no ENAs", t, pp, k, E);
          }
          const double ntr = st.ntr(t, pp, k), nref = st.nref(t, pp, k), nen = st.nen(t, pp, k),
                       etr = st.etr(t, pp, k), eref = st.eref(t, pp, k), esw = st.esw(t, pp, k);
          const double een = st.een(t, pp, k);        // E_l = 1/2 m_p v(< 12 AU)^2
          const double em = st.em(t, pp, k);
          const double kappaEn = st.kappa_ref(t, pp, k);  // spectral index delta of the DSA spectrum
          const double vshock = st.vshock(t, pp, k);
          const double esw_used = phys.fitted_esw ? esw : 0.04;
          const double nsw_frac = 1.0 - ntr - nref;
          const double ttr_frac = etr / ntr;
          const double tref_frac = eref / (nref * (1.0 - nen));
          const double tsw_frac = esw_used / nsw_frac;
          const double ten_frac = (1.0 - etr - eref - esw_used) / (nen * nref);
          if (phys.cold_electrons) tp = 2.0 * tp / (1.0 + tsw_frac);
          double np1 = ntr * np * e;                          // transmitted PUIs
          double np2 = nref * (1.0 - nen) * np * e;           // reflected PUIs
          double nsw = nsw_frac * np * e;                     // solar wind
          double np3 = nref * nen * np * e;                   // DSA-accelerated PUIs
          if (np3 <= 0) np3 = c1em10;
          const double rest = np - np1 - np2 - np3;
          double np4 = rest > 0 ? rest : c1em10;              // heliosheath PUIs (includes the SW share)
          double tp1 = tp * ttr_frac, tp2 = tp * tref_frac, tsw = tp * tsw_frac;
          const double tp3 = tp * ten_frac;
          double tp4 = 0.5 * ((ur - urh) * (ur - urh) + (ut - uth) * (ut - uth)) * mp / kb;
          for (double x : {tp1, np1, tp2, np2, tp3, np3, np4, tp4, tsw, nsw})
            if (!(x > 0)) esi_error("non-positive or NaN population density/temperature", t, pp, k, E);

          // kappas: region 3 or region 6 ([region6], defaults = region 3; SW kappa_oh)
          const bool r6 = b == 6;
          const double k_tr = r6 ? reg6.kappa_transmitted : kap_pui, k_ref = r6 ? reg6.kappa_reflected : kap_ref,
                       k_hs = r6 ? reg6.kappa_heliosheath_pui : kap_hs, k_sw = r6 ? reg6.kappa_sw : kap;
          // Region 6: a fraction f of the transmitted, reflected, heliosheath PUIs and SW ions is
          // accelerated, taking a fraction e of their pressure; density and pressure are conserved
          // (Region6Config in config.hpp). The DSA PUIs are not changed.
          // accelerated parts: (n, T) per parent population (transmitted, reflected, heliosheath PUIs, SW),
          // or one pooled population (slot 0)
          double nacc_i[4] = {0, 0, 0, 0}, tacc_i[4] = {0, 0, 0, 0};
          double nacc = 0.0, tacc = 0.0;  // totals for the dump: density, pressure-weighted temperature
          if (r6 && acc_f > 0.0) {
            const double n0[4] = {np1, np2, np4, nsw}, t0[4] = {tp1, tp2, tp4, tsw};
            if (acc_pooled) {
              double ntot = 0, ptot = 0;
              for (int q = 0; q < 4; ++q) { ntot += n0[q]; ptot += n0[q] * t0[q]; }
              nacc_i[0] = acc_f * ntot;
              tacc_i[0] = acc_e * ptot / nacc_i[0];
            } else {
              for (int q = 0; q < 4; ++q) { nacc_i[q] = acc_f * n0[q]; tacc_i[q] = acc_e / acc_f * t0[q]; }
            }
            double pacc = 0.0;
            for (int q = 0; q < 4; ++q) { nacc += nacc_i[q]; pacc += nacc_i[q] * tacc_i[q]; }
            tacc = pacc / nacc;
            const double tfac = (1.0 - acc_e) / (1.0 - acc_f);
            np1 *= 1.0 - acc_f; np2 *= 1.0 - acc_f; np4 *= 1.0 - acc_f; nsw *= 1.0 - acc_f;
            tp1 *= tfac; tp2 *= tfac; tp4 *= tfac; tsw *= tfac;
          }

          // ---- ENA flux from each population ----
          const double betk = bet(t, pp, k);
          auto cdf = [&](double T, double n, double kappa) {
            FluxInputs in{E, T, ur, ut, nn, n, kappa, betk};
            return calc_diff_flux_mi(in, urh, ps);  // cross section at vrel(T, -v_ENA, u_H,r)
          };
          FluxInputs i3{E, 0.0, ur, ut, nn, np3, 0.0, betk};
          const double j0 = cdf(tp1, np1, k_tr);
          const double j1 = cdf(tp2, np2, k_ref);
          const double j2 = calc_diff_flux_zank(i3, kappaEn, een, em, vshock, ps);  // q = 5
          const double j3 = cdf(tp4, np4, k_hs);
          const double j4 = cdf(tsw, nsw, r6 ? k_sw : kap);
          double j7 = 0.0;
          for (int q = 0; q < 4; ++q)
            if (nacc_i[q] > 0.0) j7 += cdf(tacc_i[q], nacc_i[q], reg6.accelerated_kappa);
          acc[0] += j0; acc[1] += j1; acc[2] += j2; acc[3] += j3; acc[4] += j4; acc[7] += j7;
          if (fc.cells) {
            const double jj[8] = {j0, j1, j2, j3, j4, 0.0, 0.0, j7};
            for (int m : {0, 1, 2, 3, 4, 7}) fc.cells->set(g, m, t, pp, k, jj[m]);
          }

          if (diag) {
            const std::size_t c = g.index(t, pp, k);
            const double pn[kEsiSlots] = {np1, np2, np3, np4, nsw, nacc}, pt[kEsiSlots] = {tp1, tp2, tp3, tp4, tsw, tacc},
                         pk[kEsiSlots] = {k_tr, k_ref, kappaEn, k_hs, r6 ? k_sw : kap, nacc > 0 ? reg6.accelerated_kappa : 0.0};
            for (int q = 0; q < kEsiSlots; ++q) {
              R.pops_n[c * kEsiSlots + q] = pn[q]; R.pops_t[c * kEsiSlots + q] = pt[q]; R.pops_kappa[c * kEsiSlots + q] = pk[q];
            }
            const double zp[3] = {een, em, vshock}, cp[4] = {np, tp, urh, uth};
            for (int q = 0; q < 3; ++q) R.zank_par[c * 3 + q] = zp[q];
            for (int q = 0; q < 4; ++q) R.cell_par[c * 4 + q] = cp[q];
            R.v_plas(t, pp, k) = plasma_frame(E, ur, ut);
            R.flux_test1(t, pp, k) = j0;
            R.flux_test2(t, pp, k) = j1;
            R.flux_test3(t, pp, k) = j2;
            R.flux_test_sw(t, pp, k) = j3 + j4 + j7;
            if (b == 3 && prev == 4) R.shock_strength(t, pp) = up(t, pp, k - 2) / up(t, pp, k);
          }
        }
        double tot = 0.0;
        for (int m = 0; m < 8; ++m) { R.flux[m](t, pp) = acc[m]; tot += acc[m]; }
        R.flux[8](t, pp) = tot;  // CombPUIs: all components, including the neutral ones
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
