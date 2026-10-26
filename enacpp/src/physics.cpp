#include "ena/physics.hpp"

#include <cmath>
#include <sstream>

namespace ena {

using namespace cgs;

namespace {
std::string str(double x) { std::ostringstream o; o.precision(10); o << x; return o.str(); }
}  // namespace

double energy_velocity(double e) { return std::sqrt(2.0 * e * kev_erg / mp); }

double xsection(double ea) {
  if (!(ea >= 0)) throw PhysicsError("xsection: invalid energy " + str(ea));
  const double a1 = 4.15, a2 = 0.531, a3 = 67.3;
  double s1 = a1 - a2 * std::log(ea);
  double s2 = std::pow(1.0 - std::exp(-a3 / ea), 4.5);
  double sigma = (s1 * s1) * s2;
  sigma *= 1e-16;
  if (std::isnan(sigma)) throw PhysicsError("xsection: NaN for energy " + str(ea));
  return sigma;
}

double plasma_frame(double e, double ur, double ut) {
  double v = energy_velocity(e);
  double a = v + ur;
  return std::sqrt(a * a + ut * ut);
}

double plasma_frame2(double e, double ur, double ut, double T) {
  double vp = plasma_frame(e, ur, ut);
  double vth = std::sqrt(2.0 * kb * T / mp);
  double w = vp / vth;
  double v1 = std::exp(-(w * w)) / std::sqrt(pi);
  double v2 = (w + 1.0 / (2.0 * w)) * std::erf(w);
  return vth * (v1 + v2);
}

double maxwellian_fn(double T, double n, double v) {
  double twopi = 2.0 * pi;
  double f1 = n * std::pow(mp / ((twopi * kb) * T), 1.5);
  double f2 = std::exp(((-1.0 * mp) * (v * v)) / ((2.0 * kb) * T));
  double fm = f1 * f2;
  if (std::isnan(fm)) throw PhysicsError("maxwellian_fn: NaN (T=" + str(T) + ", n=" + str(n) + ", v=" + str(v) + ")");
  return fm;
}

double kappa_fn(double T, double n, double v, double k) {
  double a = 2.0 * (k - 1.5);
  a = a / k;
  double w = std::pow(((a * kb) * T) / mp, 0.5);
  double pi15 = std::pow(pi, 1.5);
  double f1 = n / (w * w * w) / pi15;
  double f2 = (std::tgamma(k + 1.0) / std::tgamma(k - 0.5)) / std::pow(k, 1.5);
  double f3 = std::pow(1.0 + (v * v) / k / (w * w), -1.0 - k);
  double fk = (f1 * f2) * f3;
  if (std::isnan(fk)) throw PhysicsError("kappa_fn: NaN (kappa=" + str(k) + ", T=" + str(T) + ")");
  return fk;
}

double zank_fn(double v, double n, double del, double el, double em, double q) {
  const double vl = std::sqrt(2.0 * el * kev_erg / mp);
  // No accelerated PUIs below the lower energy limit E_l of the DSA spectrum
  // (een, from the TS fits), in the shock frame; this also covers negative
  // shock-frame speeds.
  if (v <= vl) return 0.0;
  const double vm = std::sqrt(2.0 * em * kev_erg / mp);
  const double fourpi = 4.0 * pi;
  const double c = ((n * (del - 3.0)) * (q - 3.0)) / ((fourpi * std::pow(vl, 3.0)) * (q - del));
  double fz;
  if (v < vm) fz = c * (std::pow(v / vl, -del) - std::pow(v / vl, -q));
  else if (v >= vm) fz = (c * (std::pow(vm / vl, q - del) - 1.0)) * std::pow(v / vl, -q);
  else throw PhysicsError("zank_fn: speed is NaN");
  if (vl == 0.0) fz = 0.0;
  if (fz < 0) fz = 0.0;
  return fz;
}

double survival_prob(double e, double bet) {
  double va = energy_velocity(e);
  double S = std::exp(((-1.0 * bet) * au_cm) / va);
  if (std::isnan(S)) throw PhysicsError("SurvivalProb: NaN");
  return S;
}

double flux_constants(double e, double nn, double bet, const PhysicsSettings& s) {
  double L = s.dr_au * au_cm;
  double sigma = xsection(e);
  double fc = (sigma * L) * nn;
  if (s.survival) fc *= survival_prob(e, bet);
  if (std::isnan(fc)) throw PhysicsError("fluxconstants: NaN");
  return fc;
}

namespace {
double finish_flux(double e, double f, double fc, const char* who) {
  double j = (((2.0 * e * kev_erg) * f) * fc) / (mp * mp);
  j *= kev_erg;  // per erg -> per keV
  if (std::isnan(j)) throw PhysicsError(std::string(who) + ": NaN flux");
  return j;
}
}  // namespace

double calc_diff_flux(const FluxInputs& in, const PhysicsSettings& s) {
  if (std::isnan(in.energy) || std::isnan(in.nn) || std::isnan(in.bet)) throw PhysicsError("calc_diff_flux: NaN input");
  double v = plasma_frame(in.energy, in.ur, in.ut);
  if (in.vshift != 0.0) v -= in.vshift;
  double f = (in.kappa <= 1.5) ? maxwellian_fn(in.T, in.np, v) : kappa_fn(in.T, in.np, v, in.kappa);
  double fc = flux_constants(in.energy, in.nn, in.bet, s);
  return finish_flux(in.energy, f, fc, "calc_diff_flux");
}

double calc_diff_flux_zank(const FluxInputs& in, double del, double el, double em, double vshock,
                           const PhysicsSettings& s) {
  double v = plasma_frame(in.energy, in.ur, in.ut) - vshock;
  double f = zank_fn(v, in.np, del, el, em);
  double fc = flux_constants(in.energy, in.nn, in.bet, s);
  return finish_flux(in.energy, f, fc, "calc_diff_flux_zank");
}

// ---- Multi-ion model ----

double vrel(double T, double vp, double vh) {
  const double vth = std::sqrt(2.0 * kb * T / mp);
  const double w = std::fabs(vp - vh) / vth;
  const double v1 = std::exp(-(w * w)) / std::sqrt(pi);
  const double v2 = (w + 1.0 / (2.0 * w)) * std::erf(w);
  return vth * (v1 + v2);
}

double flux_constants_vdiff(double e, double nn, double bet, double vdiff, const PhysicsSettings& s) {
  const double L = s.dr_au * au_cm;
  const double sigma = vdiff != 0.0 ? xsection(0.5 * mp * vdiff * vdiff / kev_erg) : xsection(e);
  double fc = (sigma * L) * nn;
  if (s.survival) fc *= survival_prob(e, bet);  // survival uses the ENA energy
  if (std::isnan(fc)) throw PhysicsError("fluxconstants: NaN");
  return fc;
}

double calc_diff_flux_mi(const FluxInputs& in, double urh, const PhysicsSettings& s) {
  if (std::isnan(in.energy) || std::isnan(in.nn) || std::isnan(in.bet)) throw PhysicsError("calc_diff_flux: NaN input");
  if (!(in.T > 0) || !(in.np > 0)) throw PhysicsError("calc_diff_flux: T or n <= 0");
  const double v = plasma_frame(in.energy, in.ur, in.ut);
  const double f = (in.kappa <= 1.5) ? maxwellian_fn(in.T, in.np, v) : kappa_fn(in.T, in.np, v, in.kappa);
  const double vdiff = std::isnan(urh) ? 0.0 : vrel(in.T, -energy_velocity(in.energy), urh);
  const double fc = flux_constants_vdiff(in.energy, in.nn, in.bet, vdiff, s);
  return finish_flux(in.energy, f, fc, "calc_diff_flux");
}

// ---- Extended single-ion model ----

double calc_diff_flux_neutrals(double e, double Tn, double ur, double ut, double nn) {
  const double v = plasma_frame(e, ur, ut);
  const double f = maxwellian_fn(Tn, nn, v);
  double j = ((2.0 * e * kev_erg) * f) / (mp * mp);
  j *= kev_erg;  // per erg -> per keV
  if (std::isnan(j)) throw PhysicsError("calc_diff_flux_neutrals: NaN flux");
  return j;
}

}  // namespace ena
