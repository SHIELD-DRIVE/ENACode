#include "ena/regions.hpp"

#include <algorithm>
#include <cstdio>
#include <string>
#include <vector>

#include <algorithm>
#include <cmath>
#include <limits>

namespace ena {

RadialLimits compute_radial_limits(const RunConfig& cfg, const Field3D<std::int16_t>& region) {
  const Grid& g = cfg.grid;
  RadialLimits L;
  L.lower.resize(g, 0.f);
  L.upper.resize(g, 0.f);
  auto radius = [&](int k) { return g.r(k); };  // true radius
  for (int i = 0; i < g.ntheta; ++i)
    for (int j = 0; j < g.nphi; ++j) {
      bool begin = false, end = false;
      const std::int16_t* reg = region.los(i, j);
      for (int k = 0; k < g.nr; ++k) {
        if (cfg.radial_range == 0) {  // TS to HP
          if (!begin && reg[k] == 3 && radius(k) > 60.0) { L.lower(i, j) = float(k); begin = true; }
          if (!end && reg[k] == 2) { L.upper(i, j) = float(k - 1); end = true; }
          if (end && reg[k] == 3) end = false;
        } else if (cfg.radial_range == 1) {  // HP to outer boundary
          if (!begin && reg[k] == 2 && radius(k) > 80.0) { L.lower(i, j) = float(k); begin = true; }
          if (begin && reg[k] == 3) begin = false;
          if (!end && reg[k] == 1) { L.upper(i, j) = float(k - 1); end = true; }
        } else if (cfg.radial_range == 2) {  // TS to outer boundary
          if (!begin && reg[k] == 3 && radius(k) > 60.0) { L.lower(i, j) = float(k); begin = true; }
        }
      }
      if (!begin) L.lower(i, j) = float(g.nr - 1);
      if (!end) L.upper(i, j) = float(g.nr - 1);
    }
  return L;
}

int heliosheath_exit(const RunConfig& cfg, const Field3D<std::int16_t>& region, int t, int p) {
  const Grid& g = cfg.grid;
  const int k0 = std::max(0, g.k_of(cfg.r_lower)), k1 = std::min(g.nr - 1, g.k_of(cfg.r_upper));
  bool inside = false;
  for (int k = k0; k <= k1; ++k) {
    const int b = region(t, p, k);
    const bool hs = b == kInnerHS || b == kOpenHS;
    if (hs) inside = true;
    else if (inside && b != kSupersonicSW) return k;  // region-4 pockets are skipped, not exits
  }
  return k1 + 1;
}

std::vector<HeliosheathPocket> heliosheath_pockets(const RunConfig& cfg, const Field3D<std::int16_t>& region) {
  const Grid& g = cfg.grid;
  const int k0 = std::max(0, g.k_of(cfg.r_lower));
  std::vector<HeliosheathPocket> out;
  for (int t = 0; t < g.ntheta; ++t)
    for (int p = 0; p < g.nphi; ++p) {
      const int kx = heliosheath_exit(cfg, region, t, p);  // pockets lie before the exit
      auto hs = [&](int k) { return region(t, p, k) == kInnerHS || region(t, p, k) == kOpenHS; };
      bool inside = false;
      for (int k = k0; k < kx; ++k) {
        if (hs(k)) { inside = true; continue; }
        if (!inside || region(t, p, k) != kSupersonicSW) continue;
        int e = k;
        while (e + 1 < kx && region(t, p, e + 1) == kSupersonicSW) ++e;
        if (e + 1 < kx && hs(e + 1)) out.push_back({t, p, k, e});  // region 4 with heliosheath behind it
        k = e;
      }
    }
  return out;
}

ReentryStats count_reentries(const RunConfig& cfg, const Field3D<std::int16_t>& region) {
  const Grid& g = cfg.grid;
  const int k1 = std::min(g.nr - 1, g.k_of(cfg.r_upper));
  ReentryStats st;
  st.r_first = std::numeric_limits<double>::infinity();
  for (int t = 0; t < g.ntheta; ++t)
    for (int p = 0; p < g.nphi; ++p) {
      long n = 0;
      const int kx = heliosheath_exit(cfg, region, t, p);
      if (kx <= k1 && region(t, p, kx) == kUndetermined) {
        if (!st.exits_at_undetermined) { st.u_theta = t * g.dtheta; st.u_phi = p * g.dphi; st.u_r = g.r(kx); }
        ++st.exits_at_undetermined;
      }
      for (int k = kx; k <= k1; ++k)
        if (region(t, p, k) == kInnerHS || region(t, p, k) == kOpenHS) {
          ++n;
          st.r_first = std::min(st.r_first, g.r(k));
        }
      if (n) { ++st.lines; st.cells += n; }
    }
  if (!st.lines) st.r_first = 0;
  return st;
}

long count_undetermined(const RunConfig& cfg, const Field3D<std::int16_t>& region) {
  const Grid& g = cfg.grid;
  int k0 = std::max(0, g.k_of(cfg.r_lower)), k1 = std::min(g.nr - 1, g.k_of(cfg.r_upper));
  long n = 0;
  for (int i = 0; i < g.ntheta; ++i)
    for (int j = 0; j < g.nphi; ++j)
      for (int k = k0; k <= k1; ++k)
        if (region(i, j, k) == kUndetermined) ++n;
  return n;
}

std::string reentry_details(const RunConfig& cfg, const Field3D<std::int16_t>& region, int max_lines) {
  const Grid& g = cfg.grid;
  const int k0 = std::max(0, g.k_of(cfg.r_lower)), k1 = std::min(g.nr - 1, g.k_of(cfg.r_upper));
  struct Line { int t, p, kx, before, after; };
  std::vector<Line> lines;
  long by_region[8] = {0, 0, 0, 0, 0, 0, 0, 0};
  for (int t = 0; t < g.ntheta; ++t)
    for (int p = 0; p < g.nphi; ++p) {
      const int kx = heliosheath_exit(cfg, region, t, p);
      if (kx > k1) continue;
      int before = 0, after = 0;
      for (int k = k0; k < kx; ++k) before += region(t, p, k) == kInnerHS || region(t, p, k) == kOpenHS;
      for (int k = kx; k <= k1; ++k) after += region(t, p, k) == kInnerHS || region(t, p, k) == kOpenHS;
      if (!after) continue;
      const int rx = region(t, p, kx);
      ++by_region[(rx >= 0 && rx < 8) ? rx : 7];
      lines.push_back({t, p, kx, before, after});
    }
  std::string out;
  char buf[256];
  if (lines.empty()) return out;
  std::snprintf(buf, sizeof buf, "  first exit into region: 1: %ld, 2: %ld, 4: %ld, 5: %ld, other: %ld\n", by_region[1],
                by_region[2], by_region[4], by_region[5], by_region[0] + by_region[3] + by_region[6] + by_region[7]);
  out += buf;
  std::sort(lines.begin(), lines.end(), [](const Line& a, const Line& b) { return a.after > b.after; });
  out += "  lines losing the most heliosheath cells (theta, phi [deg]; exit radius; cells before/after the exit; regions from 4 cells before the exit):\n";
  for (int n = 0; n < int(lines.size()) && n < max_lines; ++n) {
    const Line& L = lines[n];
    std::string seq;  // run-length region sequence around the exit, e.g. 3x4 4x2 3x60
    int k = std::max(k0, L.kx - 4);
    while (k <= k1 && seq.size() < 90) {
      int r = region(L.t, L.p, k), m = 0;
      while (k <= k1 && region(L.t, L.p, k) == r) { ++m; ++k; }
      seq += std::to_string(r) + "x" + std::to_string(m) + " ";
    }
    std::snprintf(buf, sizeof buf, "    theta %5.1f, phi %5.1f: exit at %g AU into region %d; %d / %d cells; %s\n",
                  L.t * g.dtheta, L.p * g.dphi, g.r(L.kx), int(region(L.t, L.p, L.kx)), L.before, L.after, seq.c_str());
    out += buf;
  }
  return out;
}

}  // namespace ena
