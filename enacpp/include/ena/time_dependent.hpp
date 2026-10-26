#pragma once
// Time-dependent ENA maps (milestone 8): ENA maps per energy band and per
// output time window from a sequence of snapshots ([time_dependent] in the
// configuration, TimeDependentConfig in config.hpp). Guided by the
// time-dependent flux of Zirnstein et al. (2015, Sec. 2.2.4, Eq. 2), not a copy.
//
// Method
//  - Time cells: snapshot s holds, unchanged, for [T0 + s dt, T0 + (s+1) dt).
//  - Production: per cell, energy and component, from the steady-state model of
//    that snapshot (make_flux with the survival factor off; CellProduction).
//  - Travel: ENAs move radially inward at v = (2E/m_p)^1/2 (gravity and
//    radiation pressure neglected). From radius r_k they reach the observer at
//    r_obs after tau_k = (r_k - r_obs)/v.
//  - Survival along the path: an ENA made in cell k at the centre of time cell s
//    passes cell k' < k at t = T0 + (s + 1/2) dt + (r_k - r_k')/v. The loss
//    term of cell k' (loss_increment: regions 3 and 6, from
//    physics.survival_observer_r outward, as the steady state) is taken from the
//    snapshot of that time; after the last snapshot, from the last one. The
//    exponent is summed along the path with per-snapshot prefix sums, so a path
//    costs one difference per snapshot it crosses. Identical snapshots give the
//    steady-state exponent.
//  - Arrival: ENAs made in cell k during time cell s arrive during
//    [T0 + s dt + tau_k, T0 + (s+1) dt + tau_k), at the steady flux rate of that
//    cell (production x survival). Each output window [a, b) receives the part
//    that overlaps it, divided by b - a: the window-averaged flux.
//  - Coverage: see TimeDependentConfig.
//
// Order of work: snapshots from the last to the first, so the loss terms a path
// needs (the snapshot of emission and later ones) are always available when a
// snapshot's production is integrated. Every sum has a fixed order (snapshots
// last to first, cells along the line of sight, energies in band order), so the
// maps do not depend on the thread count.

#include <string>
#include <vector>

#include "ena/config.hpp"

namespace ena {

struct OutputWindow {
  double a = 0, b = 0;  // years
  std::string label;    // e.g. "2020-2021", "2020.5-2021"
};
std::vector<OutputWindow> output_windows(const TimeDependentConfig& td);
// Shortest decimal form of a time (up to 4 decimals): 2020, 2020.5, 2020.3333.
std::string time_label(double t);

// Snapshot whose time cell contains time t_rel (years after T0), clamped to
// [0, nsnap - 1].
int time_cell(double t_rel, double dt, int nsnap);

// Survival exponent (beta, as ionization_rate) of an ENA made in cell k at
// t_emit (years after T0). prefix[s][i] is the running sum of the loss terms of
// snapshot s from cell k0 to cell k0 + i along this line of sight; cell_years
// is the time to cross one cell, dr/v. Cells k0 .. k are summed (k >= k0).
double path_exponent(const std::vector<const double*>& prefix, int k0, int k, double t_emit, double dt, int nsnap,
                     double cell_years);

struct TimeDependentOptions {
  std::string out_dir;
  bool legacy = true;  // also write the legacy count-rate .dat files
};
// Runs the method for cfg.bands and writes one file per band and window:
// <out_dir>/flux_bandNN_<window>.h5 (flux maps and coverage maps), the legacy
// .dat files with _<window> before .dat, and <out_dir>/time_dependent_summary.txt.
void run_time_dependent(const RunConfig& cfg, const TimeDependentOptions& opt);

}  // namespace ena
