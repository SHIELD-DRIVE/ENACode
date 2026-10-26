#pragma once

#include <string>
#include <vector>
// Region codes and the radial limits of each line of sight.

#include <cstdint>

#include "ena/config.hpp"
#include "ena/grid.hpp"

namespace ena {

enum Region : std::int16_t {
  kISM = 1,           // beyond the bow shock
  kOuterHS = 2,       // heliopause to bow shock
  kInnerHS = 3,       // termination shock to heliopause
  kSupersonicSW = 4,
  kUndetermined = 5,  // deliberately unset: reported as a warning
  kOpenHS = 6,
};

struct RadialLimits {
  Field2D<float> lower, upper;  // radial indices (float)
};

// get_limits: radial index limits per line of sight for cfg.radial_range.
// Diagnostic output only; the flux integration uses the region codes directly.
RadialLimits compute_radial_limits(const RunConfig& cfg, const Field3D<std::int16_t>& region);

// First radial index, inside the integration window, at which a line of sight
// leaves the heliosheath (region 3/6) after first entering it; rupper + 1 if
// it never leaves (or never enters). Integration stops there.
// Leaving means entering region 1/2 (beyond the heliopause) or any region
// other than 3, 4 and 6 (e.g. 5, undetermined). A supersonic (region 4) pocket
// between heliosheath cells, where a wavy termination shock is crossed again,
// is not an exit: it produces no ENAs and integration continues behind it
// (see heliosheath_pockets). Region 4 cells after the
// heliosheath with no heliosheath cell behind them are harmless (no ENAs).
int heliosheath_exit(const RunConfig& cfg, const Field3D<std::int16_t>& region, int t, int p);

// A supersonic (region 4) stretch inside the heliosheath of one line of sight,
// before its exit: skipped by the integration and reported, so the region
// file can be checked.
struct HeliosheathPocket {
  int t, p;        // line of sight (indices)
  int k0, k1;      // first and last radial index of the pocket
};
std::vector<HeliosheathPocket> heliosheath_pockets(const RunConfig& cfg, const Field3D<std::int16_t>& region);

// Lines of sight that leave the heliosheath and re-enter it (e.g. 3 -> 2 -> 3)
// inside the integration window, and the region 3/6 cells beyond the first exit.
struct ReentryStats {
  long lines = 0, cells = 0;
  double r_first = 0;  // smallest radius (AU) of a re-entered cell
  // lines of sight whose heliosheath ends at an undetermined (region 5) cell,
  // and the first of them (theta, phi in degrees, r in AU)
  long exits_at_undetermined = 0;
  double u_theta = 0, u_phi = 0, u_r = 0;
};
ReentryStats count_reentries(const RunConfig& cfg, const Field3D<std::int16_t>& region);
// Diagnostic for the lines of sight that leave the heliosheath and re-enter it:
// how many exit into each region, and the lines that lose the most heliosheath
// cells (region sequence around the first exit). Printed by ena_run.
std::string reentry_details(const RunConfig& cfg, const Field3D<std::int16_t>& region, int max_lines = 8);

// Number of region-5 (undetermined) cells inside the integration window.
long count_undetermined(const RunConfig& cfg, const Field3D<std::int16_t>& region);

}  // namespace ena
