#pragma once
// ENA flux along radial lines of sight (make_flux, ionization_rate) and the
// band integration (countrate), single-ion four-population model
// (Kornbleuth+2023), radial_range = 2 (termination shock outward).

#include <array>
#include <cstdint>
#include <functional>
#include <string>
#include <vector>

#include "ena/config.hpp"
#include "ena/grid.hpp"
#include "ena/physics.hpp"
#include "ena/snapshot.hpp"
#include "ena/stream_processing.hpp"

namespace ena {

// Output components. Single-ion: transmitted, reflected, accelerated,
// solar_wind, total (file tags PUI1, PUI2, PUI3, SW, CombPUIs). Multi-ion:
// transmitted, reflected, accelerated, pu2, solar_wind, total (PUI1, PUI2,
// PUI3, PUI4, SW, CombPUIs). Extended single-ion: transmitted, reflected, accelerated,
// heliosheath_pui, solar_wind, ne4, neu, total (PUI1, PUI2, PUI3, PUI4, SW,
// Ne4, Neu, CombPUIs), plus region6_accelerated (Acc6), the accelerated
// fraction of the region-6 ions ([region6], 0 by default), before the total.
// The last component is always the total.
struct ComponentSet {
  std::vector<std::string> names;  // dataset names in flux_bandNN.h5
  std::vector<std::string> tags;   // legacy count-rate file tags
  int size() const { return int(names.size()); }
};
const ComponentSet& components(const RunConfig& cfg);
using MapSet = std::vector<Field2D<double>>;

// Survival exponent along each line of sight.
using BetField = Field3D<double>;

// Flux map for one energy, plus the per-cell arrays dumped at checkpoint C4.
struct EnergyResult {
  double energy = 0;                                 // keV
  double vena = 0, sig_cx = 0;
  MapSet flux;  // per component, cm^-2 s^-1 sr^-1 keV^-1
  // diagnostics (filled when requested)
  bool has_diagnostics = false;
  BetField bet;
  Field3D<double> ext, v_plas, f_test1;
  Field3D<double> dens_test1, temp_test1, dens_test2, temp_test2, dens_test3, temp_test3, dens_test_sw, temp_test_sw,
      temp_test_p;
  Field3D<double> flux_test1, flux_test2, flux_test3, flux_test_sw;
  Field3D<float> lc_test;
  Field2D<float> lc;                                 // cooling-length index per line of sight
  Field2D<double> shock_strength;
  // multi-ion diagnostics, per cell [theta][phi][r][slot]
  bool multi_ion = false;
  std::vector<double> pops_n, pops_t, pops_kappa;  // 5 slots: transmitted, reflected, DSA, pu2, SW (as passed to the flux calls)
  std::vector<double> zank_par, cell_par;          // 4 slots: E_l, E_m, v_shock, q; n, T, m vA^2/kT_pu2, m vA^2 (keV)
  Field3D<std::int16_t> belt_mask;                 // reconnection-region cells
  // extended single-ion: pops_* as above (slot 3 = heliosheath PUIs), zank_par 3 slots (E_l, E_m,
  // v_shock), cell_par 4 slots (n, T after the cold-electron correction, u_H,r, u_H,t)
  bool extended_si = false;
};

// One band: the count-rate map (legacy file layout) and the band-averaged flux.
struct BandResult {
  int band = 0;
  MapSet countrate;  // flux * E * G summed with the band weights
  MapSet flux;       // countrate / (Ecentral * G)
  double emin = 0, emax = 0, ecen = 0, g = 0, denergy = 0;
  int eres = 10;
};

// Per-cell ENA production, for the time-dependent method: every make_flux
// variant stores the flux each cell adds to each component (all components
// except the total), [component][theta][phi][r]. Cells that add nothing hold 0.
// The values include the survival factor only if physics.survival is on; the
// time-dependent method runs make_flux with survival off and applies the
// survival along the ENA's path itself.
struct CellProduction {
  int ncomp = 0;            // components(cfg).size() - 1
  std::vector<double> j;    // ncomp * grid size
  void resize(const Grid& g, int n) { ncomp = n; j.assign(std::size_t(n) * g.size(), 0.0); }
  double* comp(int m) { return j.data() + std::size_t(m) * (j.size() / std::size_t(ncomp)); }
  const double* comp(int m) const { return j.data() + std::size_t(m) * (j.size() / std::size_t(ncomp)); }
  // zero one line of sight (called by make_flux before it fills that line of sight)
  void clear_los(const Grid& g, int t, int p) {
    for (int m = 0; m < ncomp; ++m) {
      double* c = comp(m) + g.index(t, p, 0);
      for (int k = 0; k < g.nr; ++k) c[k] = 0.0;
    }
  }
  void set(const Grid& g, int m, int t, int p, int k, double v) { comp(m)[g.index(t, p, k)] = v; }
};

struct FluxContext {
  const RunConfig& cfg;
  const Snapshot& snap;
  const StreamFields& stream;
  CellProduction* cells = nullptr;  // time-dependent method: per-cell production (see above)
};

// ionization_rate: survival exponent beta along each line of sight, summed
// over region 3 and 6 cells from physics.survival_observer_r outward.
BetField ionization_rate(const FluxContext& ctx, double energy);
// The loss term one heliosheath cell adds to beta: sigma(E_rel) v_rel n dr, with
// the relative proton-neutral speed (Kornbleuth+2021 Eq. 8); pl is the ion
// fluid used for the losses (plasma[0]).
double loss_increment(const PlasmaFluid& pl, double energy, double dr, int t, int p, int k);
// First radial index of the survival losses (physics.survival_observer_r).
int survival_first_index(const RunConfig& cfg);

// make_flux for one energy. `diagnostics` keeps the per-cell arrays.
EnergyResult make_flux(const FluxContext& ctx, double energy, int band, bool diagnostics);
// Multi-ion make_flux (src/flux_mi.cpp); make_flux calls it for that model.
EnergyResult make_flux_mi(const FluxContext& ctx, double energy, int band, bool diagnostics);
// Extended single-ion make_flux; called by make_flux.
EnergyResult make_flux_esi(const FluxContext& ctx, double energy, int band, bool diagnostics);

// countrate: samples the centre of each of the bands.eres sub-bands (or the
// central energy for central-only bands). `on_energy` sees each energy's result
// (e.g. to write checkpoint C4) before it is accumulated.
BandResult integrate_band(const FluxContext& ctx, int band, const std::function<bool(double energy)>& want_diagnostics,
                          const std::function<void(const EnergyResult&)>& on_energy);

// The energies a band is evaluated at, in order (integrate_band).
std::vector<double> band_energies(const RunConfig& cfg, int band);
// Weight of energy i (band_energies order) in the band count rate, per unit flux:
// E G for central-only bands, else E G times the sub-band weight (instrument
// response where given, else 1/eres). The band count rate is sum_i weight_i flux_i.
double band_energy_weight(const RunConfig& cfg, int band, int i, double energy);
// Band sum of integrate_band: add the flux maps of energy i (in band_energies
// order), then finish() gives the band result. Used by integrate_band and by
// the time-dependent method, so both combine energies identically.
class BandAccumulator {
 public:
  BandAccumulator(const RunConfig& cfg, int band);
  void add(int i, double energy, const MapSet& flux);
  BandResult finish();

 private:
  const RunConfig& cfg_;
  BandResult B_;
  std::vector<double> wt_;
  bool central_ = false;
};

}  // namespace ena
