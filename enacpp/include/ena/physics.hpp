#pragma once
// Physics functions of the ENA model, in double precision.

#include <cmath>
#include <stdexcept>
#include <string>

namespace ena {

namespace cgs {
constexpr double mp = 1.66e-24;            // g (value used by the IDL code)
constexpr double kb = 1.3806503e-16;       // erg/K
constexpr double kev_erg = 1.60217646e-9;  // erg per keV
constexpr double au_cm = 1.49598e13;       // cm
constexpr double pi = 3.14159265358979323846;
}  // namespace cgs

// Parameters of one population's flux evaluation.
struct FluxInputs {
  double energy;        // keV
  double T;             // K
  double ur, ut;        // plasma bulk velocity, radial and tangential (cm/s)
  double nn;            // neutral H density (cm^-3)
  double np;            // ion density (cm^-3)
  double kappa;         // kappa <= 1.5 -> Maxwellian
  double bet;           // survival exponent from ionization_rate
  double vshift = 0.0;  // shock_shift = "scalar" option: speed subtracted from |v_ENA + u| (cm/s)
};

struct PhysicsSettings {
  bool survival = true;  // survival probability on
  double dr_au = 2.0;    // line-of-sight step L_los (AU)
};

// energy_velocity: sqrt(2E/m_p), E in keV.
double energy_velocity(double e_kev);

// xsection: Lindsay & Stebbings (2005) charge-exchange cross-section (cm^2), E in keV.
double xsection(double e_kev);

// plasma_frame: speed of the parent ion in the plasma frame for an inward ENA.
double plasma_frame(double e_kev, double ur, double ut);

// plasma_frame2: thermally averaged relative speed (Kornbleuth+2021 Eq. 8).
double plasma_frame2(double e_kev, double ur, double ut, double T);

double maxwellian_fn(double T, double n, double v);
double kappa_fn(double T, double n, double v, double kappa);
// Zank-Wang distribution (accelerated PUIs); q = 5 in the single-ion model.
// Zero below the lower energy limit E_l (shock-frame speed v <= v_l), which
// also covers negative speeds.
double zank_fn(double v, double n, double delta, double el_kev, double em_kev, double q = 5.0);

// SurvivalProb: exp(-bet * AU / v_ENA).
double survival_prob(double e_kev, double bet);

// fluxconstants: sigma(E) * L * n_H * S.
double flux_constants(double e_kev, double nn, double bet, const PhysicsSettings& s);

// calc_diff_flux: ENA flux produced in one cell, cm^-2 s^-1 sr^-1 keV^-1.
double calc_diff_flux(const FluxInputs& in, const PhysicsSettings& s);

// calc_diff_flux_zank: accelerated population, evaluated in the shock frame.
double calc_diff_flux_zank(const FluxInputs& in, double delta, double el_kev, double em_kev, double vshock,
                           const PhysicsSettings& s);

// ---- Multi-ion model ----
// vrel: thermally averaged speed of the parent ion relative to the neutrals
// (Heerikhuisen+2006); vplasma is the ENA velocity (-v_ENA), vr_h the neutral
// radial bulk speed (cm/s).
double vrel(double T, double vplasma, double vr_h);
// fluxconstants with the cross section at the relative proton-neutral energy
// 0.5 m vdiff^2 (vdiff != 0), else at the ENA energy.
double flux_constants_vdiff(double e_kev, double nn, double bet, double vdiff, const PhysicsSettings& s);
// calc_diff_flux of the multi-ion code: vdiff = vrel(T, -v_ENA, urh) (urh in
// cm/s; NaN urh: vdiff = 0).
double calc_diff_flux_mi(const FluxInputs& in, double urh, const PhysicsSettings& s);

// ---- Extended single-ion model ----
// calc_diff_flux_neutrals: flux of a neutral population itself (components
// Ne4/Neu): a Maxwellian of temperature Tn and density nn in the frame of the
// neutral bulk velocity (ur, ut in cm/s), j = 2 E f / m_p^2 per keV. No cross
// section, path length or survival factor.
double calc_diff_flux_neutrals(double e_kev, double Tn, double ur, double ut, double nn);

// Thrown on an invalid physical state (e.g. a NaN or a negative energy).
struct PhysicsError : std::runtime_error {
  using std::runtime_error::runtime_error;
};

}  // namespace ena
