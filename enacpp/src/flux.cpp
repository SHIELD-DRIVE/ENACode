#include "ena/flux.hpp"
#include "ena/parallel.hpp"

#include "ena/regions.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace ena {

const ComponentSet& components(const RunConfig& cfg) {
  static const ComponentSet si{{"transmitted", "reflected", "accelerated", "solar_wind", "total"},
                               {"PUI1", "PUI2", "PUI3", "SW", "CombPUIs"}};
  static const ComponentSet mi{{"transmitted", "reflected", "accelerated", "pu2", "solar_wind", "total"},
                               {"PUI1", "PUI2", "PUI3", "PUI4", "SW", "CombPUIs"}};
  // extended single-ion: heliosheath_pui = PUIs created in the heliosheath (n_p - n1 - n2 - n3, T from the
  // ion-neutral relative flow); ne4/neu = neutral populations 1 and 2 at the inner boundary (file tags Ne4, Neu)
  // region6_accelerated = the accelerated fraction of the region-6 ions ([region6]; 0 by default)
  static const ComponentSet esi{{"transmitted", "reflected", "accelerated", "heliosheath_pui", "solar_wind", "ne4", "neu",
                                 "region6_accelerated", "total"},
                                {"PUI1", "PUI2", "PUI3", "PUI4", "SW", "Ne4", "Neu", "Acc6", "CombPUIs"}};
  return cfg.multi_ion() ? mi : cfg.extended_si() ? esi : si;
}

namespace {

using namespace cgs;

[[noreturn]] void cell_error(const std::string& what, int t, int pp, int k, double energy) {
  std::ostringstream o;
  o << what << " at theta index " << t << ", phi index " << pp << ", r index " << k << ", energy " << energy << " keV";
  throw PhysicsError(o.str());
}

}  // namespace

// ---------------------------------------------------------------------------
int survival_first_index(const RunConfig& cfg) {
  return std::max(0, int(std::lround((cfg.physics.survival_observer_r - cfg.grid.r_inner) / cfg.grid.dr)));
}

double loss_increment(const PlasmaFluid& pl, double E, double dr, int t, int pp, int k) {
  double vp = plasma_frame2(E, pl.vr(t, pp, k), pl.vt(t, pp, k), pl.T(t, pp, k));
  double erel = 0.5 * mp * (vp * vp);  // erg
  if (std::isnan(erel)) cell_error("ionization_rate: E_rel is NaN", t, pp, k, E);
  double sig = xsection(erel / kev_erg);  // xsection takes keV
  return ((sig * vp) * pl.n(t, pp, k)) * dr;
}

BetField ionization_rate(const FluxContext& fc, double E) {
  const Grid& g = fc.cfg.grid;
  BetField bet(g, 0);
  if (!fc.cfg.physics.survival) return bet;
  const PlasmaFluid& pl = fc.snap.plasma.at(0);
  const int k0 = survival_first_index(fc.cfg);
  const double dr = g.dr;
  ParallelErrors perr;
  const long nlos = long(g.ntheta) * g.nphi;
#pragma omp parallel for schedule(static)
  for (long los = 0; los < nlos; ++los) {
    const int t = int(los / g.nphi), pp = int(los % g.nphi);
    perr.capture(std::size_t(los), [&] {
      for (int k = k0; k < g.nr; ++k) {
        // Losses are counted in the heliosheath, regions 3 and 6, and the
        // accumulated exponent is carried along the line of sight through
        // every other cell.
        const int code = fc.snap.region(t, pp, k);
        if (code != 3 && code != 6) {
          if (k > 0) bet(t, pp, k) = bet(t, pp, k - 1);
          continue;
        }
        const double b = loss_increment(pl, E, dr, t, pp, k);
        double prev = k > 0 ? double(bet(t, pp, k - 1)) : 0.0;
        bet(t, pp, k) = BetField::value_type(b + prev);
      }
    });
  }
  perr.rethrow();
  return bet;
}

// ---------------------------------------------------------------------------
EnergyResult make_flux(const FluxContext& fc, double E, int band, bool diag) {
  if (fc.cfg.multi_ion()) return make_flux_mi(fc, E, band, diag);
  if (fc.cfg.extended_si()) return make_flux_esi(fc, E, band, diag);
  (void)band;
  const RunConfig& cfg = fc.cfg;
  const Grid& g = cfg.grid;
  const PhysicsConfig& phys = cfg.physics;
  const Snapshot& snap = fc.snap;
  const PlasmaFluid& pl = snap.plasma.at(0);
  const NeutralPopulation& nh = snap.neutrals.at(0);
  const StreamFields& st = fc.stream;
  if (cfg.radial_range != 2)
    throw std::runtime_error("radial_range " + std::to_string(cfg.radial_range) + " is not implemented in this model version");
  if (cfg.populations != "four_population_ts")
    throw std::runtime_error("populations.model '" + cfg.populations + "' is not implemented yet");
  PhysicsSettings ps;
  ps.survival = phys.survival;
  ps.dr_au = g.dr;
  auto wrapk = [&](int k) { return k < 0 ? k + g.nr : k; };
  auto reg = [&](int t, int pp, int k) { return snap.region(t, pp, wrapk(k)); };

  EnergyResult R;
  R.energy = E;
  R.vena = energy_velocity(E);
  R.sig_cx = xsection(E);
  R.flux.assign(5, Field2D<double>(g, 0.0));
  R.has_diagnostics = diag;

  // Energy-dependent extinction along streamlines (Kornbleuth+2021 Eq. 4).
  Field3D<double> ext(g);
#pragma omp parallel for schedule(static)
  for (std::size_t i = 0; i < g.size(); ++i)
    ext.raw()[i] = phys.extinction ? std::exp((st.cool.raw()[i] * R.vena) * R.sig_cx) : 1.0;

  BetField bet = ionization_rate(fc, E);

  if (diag) {
    for (auto* f : {&R.v_plas, &R.f_test1, &R.dens_test1, &R.temp_test1, &R.dens_test2, &R.temp_test2, &R.dens_test3,
                    &R.temp_test3, &R.dens_test_sw, &R.temp_test_sw, &R.temp_test_p, &R.flux_test1, &R.flux_test2,
                    &R.flux_test3, &R.flux_test_sw})
      f->resize(g, 0.0);
    R.lc_test.resize(g, 0.f);
    R.shock_strength.resize(g, 0.0);
  }
  R.lc.resize(g, 0.f);

  const int rlower = g.k_of(cfg.r_lower);
  const double e_target = 1.0 / std::exp(1.0);
  const double kap = phys.kappa, kap_pui = phys.kappa_pui, kap_ref = phys.kappa_ref, kap_ism = phys.kappa_ism,
               kap_oh = phys.kappa_oh;
  const double esw_const = 0.04, c1em10 = 1e-10;
  const bool shift_scalar = phys.shock_shift == "scalar", shift_vector = phys.shock_shift == "vector";

  // Lines of sight are independent work items: each writes only its own
  // pixels and cells, so the result does not depend on the thread count.
  ParallelErrors perr;
  const long nlos = long(g.ntheta) * g.nphi;
#pragma omp parallel for schedule(dynamic, 4)
  for (long los = 0; los < nlos; ++los) {
    const int t = int(los / g.nphi), pp = int(los % g.nphi);
    perr.capture(std::size_t(los), [&] {
        // cooling length: index where ext is closest to 1/e
        {
          double best = std::numeric_limits<double>::infinity();
          int bi = 0;
          for (int k = 0; k < g.nr; ++k) {
            double d = std::fabs(e_target - ext(t, pp, k));
            if (d < best) { best = d; bi = k; }
          }
          R.lc(t, pp) = float(bi);
        }
        if (fc.cells) fc.cells->clear_los(g, t, pp);
        double kappa0 = 1.0, kappaEn = kap;
        std::array<double, 4> acc{0, 0, 0, 0};
        // Integration stops where the line of sight first leaves the heliosheath.
        const int k_end = heliosheath_exit(cfg, snap.region, t, pp) - 1;
        for (int k = rlower; k <= k_end; ++k) {
          const int b = snap.region(t, pp, k);
          if (!(b == 3 || b == 2 || b == 1 || b == 6 || b == 4)) continue;
          if (b == 3) kappa0 = kap;
          if (b == 6) kappa0 = kap_oh;
          if (b < 3) kappa0 = kap_ism;
          if (b == 4 || b == 5) kappa0 = 1.0;
          if (b != 3) kappaEn = kap;
          double& e = ext(t, pp, k);
          if (std::isnan(e) || e == 0.0) e = c1em10;
          const double np = pl.n(t, pp, k), nn = nh.n(t, pp, k), ur = pl.vr(t, pp, k), ut = pl.vt(t, pp, k);
          double tp = pl.T(t, pp, k);
          const int prev = reg(t, pp, k - 1);

          const bool populate = (b == 3 || b == 6) && prev != 2 && prev != 1;
          if (!populate) continue;
          // Termination-shock partition (Kornbleuth+2023), carried along streamlines.
          // A cell the gap fill left empty has no termination-shock partition (all fractions 0).
          if (!(st.ntr(t, pp, k) != 0.0 && std::isfinite(st.ntr(t, pp, k)))) {
            if (phys.unfilled_cells == "skip") continue;
            cell_error("heliosheath cell without streamline data (the gap fill left it empty, see the Gap fill line); set [physics] unfilled_cells = \"skip\" to give such cells no ENAs", t, pp, k, E);
          }
          const double ntr = st.ntr(t, pp, k), nref = st.nref(t, pp, k), nen = st.nen(t, pp, k),
                       etr = st.etr(t, pp, k), eref = st.eref(t, pp, k), een = st.een(t, pp, k), em = st.em(t, pp, k),
                       vshock = st.vshock(t, pp, k), esw = st.esw(t, pp, k);
          kappaEn = st.kappa_ref(t, pp, k);
          const double esw_used = phys.fitted_esw ? esw : esw_const;
          const double nsw_frac = 1.0 - ntr - nref;
          const double ttr_frac = etr / ntr;
          const double tref_frac = eref / (nref * (1.0 - nen));
          const double tsw_frac = esw_used / nsw_frac;
          const double ten_frac = (((1.0 - etr) - eref) - esw_used) / (nen * nref);
          if (phys.cold_electrons) tp = (2.0 * tp) / (1.0 + tsw_frac);
          const double ex = e;
          const double np1 = (ntr * np) * ex;
          const double np2 = ((nref * (1.0 - nen)) * np) * ex;
          const double nsw = (nsw_frac * np) * ex;
          double np3 = ((nref * nen) * np) * ex;
          if (np3 <= 0) np3 = c1em10;
          const double tp1 = tp * ttr_frac, tp2 = tp * tref_frac, tp3 = tp * ten_frac, tsw = tp * tsw_frac;
          for (double x : {tp1, np1, tp2, np2, tp3, np3, tsw, nsw})
            if (!(x > 0)) cell_error("non-positive or NaN population density/temperature", t, pp, k, E);

          FluxInputs base{E, 0, ur, ut, nn, 0, 0, double(bet(t, pp, k))};
          FluxInputs i1 = base; i1.T = tp1; i1.np = np1; i1.kappa = kap_pui;
          FluxInputs i2 = base; i2.T = tp2; i2.np = np2; i2.kappa = kap_ref;
          FluxInputs i3 = base; i3.np = np3;
          FluxInputs isw = base; isw.T = tsw; isw.np = nsw; isw.kappa = kappa0;
          if (shift_scalar) {  // option: shock frame, as for the accelerated PUIs
            i1.vshift = i2.vshift = isw.vshift = vshock;
          } else if (shift_vector) {  // option: frame moving at v_shock along the flow
            const double u = std::sqrt(ur * ur + ut * ut);
            if (u > 0) {
              const double fac = 1.0 - vshock / u;
              for (FluxInputs* in : {&i1, &i2, &isw}) { in->ur = ur * fac; in->ut = ut * fac; }
            }
          }
          const double j1 = calc_diff_flux(i1, ps);
          const double j2 = calc_diff_flux(i2, ps);
          const double j3 = calc_diff_flux_zank(i3, kappaEn, een, em, vshock, ps);
          const double jsw = calc_diff_flux(isw, ps);
          // region 3 and region 6 (open heliosheath) both contribute
          acc[0] = j1 + acc[0];
          acc[1] = j2 + acc[1];
          acc[2] = j3 + acc[2];
          acc[3] = jsw + acc[3];
          if (fc.cells) {
            fc.cells->set(g, 0, t, pp, k, j1); fc.cells->set(g, 1, t, pp, k, j2);
            fc.cells->set(g, 2, t, pp, k, j3); fc.cells->set(g, 3, t, pp, k, jsw);
          }
          if (diag) {
            R.dens_test1(t, pp, k) = np1; R.temp_test1(t, pp, k) = tp1;
            R.dens_test2(t, pp, k) = np2; R.temp_test2(t, pp, k) = tp2;
            R.dens_test3(t, pp, k) = np3; R.temp_test3(t, pp, k) = tp3;
            R.dens_test_sw(t, pp, k) = nsw; R.temp_test_sw(t, pp, k) = tsw;
            R.temp_test_p(t, pp, k) = tp;
            const double vpl = plasma_frame(E, ur, ut);
            R.v_plas(t, pp, k) = vpl;
            R.f_test1(t, pp, k) = maxwellian_fn(tp1, np1, vpl);
            R.lc_test(t, pp, k) = (k == int(R.lc(t, pp))) ? 1.f : 0.f;
            if (b == 3 && prev == 4) {
              auto up = [&](int kk) { kk = wrapk(kk); return std::sqrt(pl.vr(t, pp, kk) * pl.vr(t, pp, kk) + pl.vt(t, pp, kk) * pl.vt(t, pp, kk)); };
              R.shock_strength(t, pp) = up(k - 2) / up(k);
            }
            R.flux_test1(t, pp, k) = j1; R.flux_test2(t, pp, k) = j2;
            R.flux_test3(t, pp, k) = j3; R.flux_test_sw(t, pp, k) = jsw;
          }
        }
        for (int m = 0; m < 4; ++m) R.flux[m](t, pp) = acc[m];
        R.flux[4](t, pp) = ((acc[0] + acc[1]) + acc[2]) + acc[3];
    });
  }
  perr.rethrow();
  if (diag) {
    R.ext = std::move(ext);
    R.bet = std::move(bet);
  }
  return R;
}

// ---------------------------------------------------------------------------
std::vector<double> band_energies(const RunConfig& cfg, int band) {
  const SensorTable& s = cfg.sensors;
  if (s.central_only(band)) return {s.ecen.at(band)};
  std::vector<double> e;
  const double emin = s.emin.at(band), emax = s.emax.at(band);
  // sub-band centres
  for (int i = 0; i < cfg.eres; ++i) e.push_back(emin + (i + 0.5) * (emax - emin) / cfg.eres);
  return e;
}

double band_energy_weight(const RunConfig& cfg, int band, int i, double E) {
  const SensorTable& s = cfg.sensors;
  const double eg = E * s.g.at(band);
  if (s.central_only(band)) return eg;
  auto wit = s.weights.find(band);
  if (wit == s.weights.end()) return eg / double(cfg.eres);
  double tot = 0;
  for (double w : wit->second) tot += w;
  return eg * wit->second.at(i) / tot;
}

BandAccumulator::BandAccumulator(const RunConfig& cfg, int band) : cfg_(cfg) {
  const SensorTable& s = cfg.sensors;
  BandResult& B = B_;
  B.band = band;
  B.eres = cfg.eres;
  B.emin = s.emin.at(band);
  B.emax = s.emax.at(band);
  B.ecen = s.ecen.at(band);
  B.g = s.g.at(band);
  B.denergy = (B.emax - B.emin) / double(cfg.eres);
  B.countrate.assign(components(cfg).size(), Field2D<double>(cfg.grid, 0.0));
  central_ = s.central_only(band);
  auto wit = s.weights.find(band);
  if (wit != s.weights.end()) {
    double tot = 0;
    for (double w : wit->second) tot = tot + w;
    for (double w : wit->second) wt_.push_back(w / tot);
    if (int(wt_.size()) != cfg.eres)
      throw std::runtime_error("band " + std::to_string(band) + ": number of weights differs from bands.eres");
  }
}

void BandAccumulator::add(int i, double E, const MapSet& flux) {
  BandResult& B = B_;
  const double span = B.emax - B.emin;
  for (std::size_t m = 0; m < B.countrate.size(); ++m) {
    auto& cr = B.countrate[m].raw();
    const auto& fx = flux[m].raw();
    for (std::size_t q = 0; q < cr.size(); ++q) {
      double v = (fx[q] * E) * B.g;
      // central-only bands: one energy, countrate = flux E G; else the sub-band
      // weight (instrument response where given, else dE/span)
      if (!central_) v = wt_.empty() ? (v * B.denergy) / span : v * wt_.at(i);
      cr[q] += v;
    }
  }
}

BandResult BandAccumulator::finish() {
  BandResult& B = B_;
  const int ncomp = int(B.countrate.size());
  B.flux.assign(ncomp, Field2D<double>(cfg_.grid, 0.0));
  for (int m = 0; m < ncomp; ++m) {
    for (std::size_t i = 0; i < B.flux[m].raw().size(); ++i) B.flux[m].raw()[i] = (B.countrate[m].raw()[i] / B.ecen) / B.g;
  }
  return std::move(B_);
}

BandResult integrate_band(const FluxContext& fc, int band, const std::function<bool(double)>& want_diag,
                          const std::function<void(const EnergyResult&)>& on_energy) {
  BandAccumulator acc(fc.cfg, band);
  const std::vector<double> energies = band_energies(fc.cfg, band);
  for (std::size_t i = 0; i < energies.size(); ++i) {
    const double E = energies[i];
    EnergyResult r = make_flux(fc, E, band, want_diag && want_diag(E));
    if (on_energy) on_energy(r);
    acc.add(int(i), E, r.flux);
  }
  return acc.finish();
}

}  // namespace ena
