// Time-dependent ENA maps (milestone 8). The method is described in
// include/ena/time_dependent.hpp, the settings in config.hpp
// (TimeDependentConfig).

#include "ena/time_dependent.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <vector>

#include "ena/checkpoints.hpp"
#include "ena/comm.hpp"
#include "ena/flux.hpp"
#include "ena/parallel.hpp"
#include "ena/physics.hpp"
#include "ena/pipeline.hpp"
#include "ena/regions.hpp"
#include "ena/snapshot.hpp"

namespace ena {

namespace {

using Clock = std::chrono::steady_clock;
double since(Clock::time_point t0) { return std::chrono::duration<double>(Clock::now() - t0).count(); }

constexpr double kYear = 3.15576e7;  // s (Julian year)

std::string band_label(int band) {
  char b[16];
  std::snprintf(b, sizeof b, "%02d", band);
  return b;
}

// One energy of one band (the unit of work of the production step).
struct EnergyTask {
  int band = 0, index = 0;  // band, position in band_energies
  double energy = 0;        // keV
  double v_cm = 0;          // ENA speed (cm/s)
  double v_au_yr = 0;       // ENA speed (AU/yr)
};

}  // namespace

// ---------------------------------------------------------------------------
std::string time_label(double t) {
  char b[64];
  std::snprintf(b, sizeof b, "%.4f", t);
  std::string s = b;
  while (!s.empty() && s.back() == '0') s.pop_back();
  if (!s.empty() && s.back() == '.') s.pop_back();
  if (s == "-0") s = "0";
  return s;
}

std::vector<OutputWindow> output_windows(const TimeDependentConfig& td) {
  std::vector<OutputWindow> w;
  const double eps = 1e-9 * std::max(1.0, std::fabs(td.output_end - td.output_start));
  for (long n = 0;; ++n) {
    const double a = td.output_start + double(n) * td.cadence, b = a + td.window;
    if (b > td.output_end + eps) break;
    w.push_back({a, b, time_label(a) + "-" + time_label(b)});
  }
  return w;
}

int time_cell(double t_rel, double dt, int nsnap) {
  const double x = std::floor(t_rel / dt);
  if (!(x > 0)) return 0;  // also NaN
  if (x >= double(nsnap - 1)) return nsnap - 1;
  return int(x);
}

double path_exponent(const std::vector<const double*>& prefix, int k0, int k, double t_emit, double dt, int nsnap,
                     double cell_years) {
  // Snapshot used for cell k' = k - d: the one covering t_emit + d * cell_years.
  // It does not decrease with d, so the path splits into stretches of constant
  // snapshot; each adds prefix[s][end] - prefix[s][start - 1].
  auto snap_at = [&](int kk) { return time_cell(t_emit + double(k - kk) * cell_years, dt, nsnap); };
  double beta = 0.0;
  int kk = k;
  while (kk >= k0) {
    const int s = snap_at(kk);
    int kl = k0;  // lowest cell of this stretch
    if (s < nsnap - 1) {
      // estimate from d < ((s + 1) dt - t_emit)/cell_years, then correct for rounding
      const double x = (double(s + 1) * dt - t_emit) / cell_years;
      const double dmax = std::ceil(x) - 1.0;
      kl = dmax >= double(k - k0) ? k0 : std::max(k0, k - int(std::max(dmax, 0.0)));
      if (kl > kk) kl = kk;
      while (kl > k0 && snap_at(kl - 1) == s) --kl;
      while (snap_at(kl) != s) ++kl;
    }
    const double* p = prefix[s];
    beta += p[kk - k0] - (kl > k0 ? p[kl - 1 - k0] : 0.0);
    kk = kl - 1;
  }
  return beta;
}

// ---------------------------------------------------------------------------
void run_time_dependent(const RunConfig& cfg, const TimeDependentOptions& opt) {
  const TimeDependentConfig& td = cfg.time_dependent;
  if (!td.enabled) throw std::runtime_error("run_time_dependent needs a [time_dependent] section");
  const Grid& g = cfg.grid;
  const int nsnap = td.snapshots();
  const std::vector<OutputWindow> windows = output_windows(td);
  const int nwin = int(windows.size());
  const ComponentSet& cs = components(cfg);
  const int ncomp = cs.size();   // including the total (last)
  const int nprod = ncomp - 1;   // components with a per-cell production
  const long nlos = long(g.ntheta) * g.nphi;
  const bool survival = cfg.physics.survival;
  // losses: cells k0 .. k1 (k1: end of the integration window, the last cell that can emit)
  const int k0 = survival_first_index(cfg);
  const int k1 = std::min(g.nr - 1, g.k_of(cfg.r_upper));
  const int nloss = survival ? std::max(0, k1 - k0 + 1) : 0;
  const double r_obs = td.observer_r_au;

  // Production is evaluated without the survival factor; survival is applied along the path.
  RunConfig cfg_prod = cfg;
  cfg_prod.physics.survival = false;

  // ---- energies, in band order ----
  std::vector<EnergyTask> tasks;
  for (int band : cfg.bands) {
    const std::vector<double> es = band_energies(cfg, band);
    for (std::size_t i = 0; i < es.size(); ++i) {
      EnergyTask e;
      e.band = band;
      e.index = int(i);
      e.energy = es[i];
      e.v_cm = energy_velocity(es[i]);
      e.v_au_yr = e.v_cm * kYear / cgs::au_cm;
      tasks.push_back(e);
    }
  }
  const int ntask = int(tasks.size());

  // ---- layout: snapshots and lines of sight over the MPI ranks ----
  // Production phase: rank r reads the snapshots [s_begin(r), s_begin(r+1)) (contiguous
  // blocks), computes their loss sums and per-cell production, and sends each line of
  // sight's share to the rank that owns it. Integration phase: rank r integrates the lines
  // of sight [l_begin(r), l_begin(r+1)) over all snapshots. Each line of sight is summed in
  // the same order on one rank whatever the layout, so the maps are identical for any
  // number of ranks and threads. Without MPI there is one rank and the exchange is a copy.
  const Comm& comm = Comm::world();
  const int nranks = comm.size(), me = comm.rank();
  const bool root = comm.root();
  auto s_begin = [&](int r) { return int(block_begin(nsnap, nranks, r)); };
  auto l_begin = [&](int r) { return long(block_begin(nlos, nranks, r)); };
  const long l0 = l_begin(me), nlos_own = l_begin(me + 1) - l0;
  int rounds = 0;  // production rounds: the largest snapshot block
  for (int r = 0; r < nranks; ++r) rounds = std::max(rounds, s_begin(r + 1) - s_begin(r));
  // snapshot handled by rank r in round j (from its last snapshot down), or -1
  auto snap_of = [&](int r, int j) { const int hi = s_begin(r + 1); return hi - 1 - j >= s_begin(r) ? hi - 1 - j : -1; };

  // ---- every input file must exist before the (long) run starts ----
  {
    std::vector<std::string> missing;
    for (int s = 0; s < nsnap; ++s) {
      const SnapshotFiles& f = cfg.inputs.snapshots[s];
      const std::vector<std::string> names =
          cfg.inputs.format == "hdf5" ? std::vector<std::string>{f.file}
                                      : std::vector<std::string>{f.plasma, f.neutral, f.boundaries, f.streamlines};
      for (const std::string& n : names)
        if (n.empty() || !std::filesystem::exists(cfg.inputs.path(n)))
          missing.push_back("snapshot " + std::to_string(s) + ": " + (n.empty() ? "(no streamline file given)" : cfg.inputs.path(n)));
    }
    std::string msg;
    if (!missing.empty()) {
      msg = std::to_string(missing.size()) + " input files are missing:";
      for (std::size_t i = 0; i < std::min<std::size_t>(missing.size(), 10); ++i) msg += "\n  " + missing[i];
      if (missing.size() > 10) msg += "\n  ...";
    }
    comm.agree(msg);
  }

  // ---- energy groups: what a rank stores for all snapshots must fit in memory_gb ----
  // Per energy and rank: the loss sums of its lines of sight (N x lines x cells) and their
  // per-cell production (N x emitting cells x (components + 1)). The number of emitting
  // cells is estimated from the first snapshot (heliosheath cells in the integration
  // window, + 25%), read once by rank 0 for this.
  double hs_cells = 0.0;
  {
    std::string err;
    if (root) {
      try {
        ReadOptions ro;
        ro.verbose = false;
        Snapshot s0 = read_snapshot(cfg, 0, ro);
        const int rl = g.k_of(cfg.r_lower);
        for (int t = 0; t < g.ntheta; ++t)
          for (int pp = 0; pp < g.nphi; ++pp) {
            const int kx = heliosheath_exit(cfg, s0.region, t, pp);
            for (int k = std::max(rl, 0); k < kx; ++k)
              if (s0.region(t, pp, k) == 3 || s0.region(t, pp, k) == 6) hs_cells += 1.0;
            hs_cells += 1.0;  // a component at the inner boundary (extended single-ion neutrals)
          }
      } catch (const std::exception& e) {
        err = e.what();
      }
    }
    comm.agree(err);
    hs_cells = comm.max_double(hs_cells);  // from rank 0
  }
  const long nlos_max = l_begin(1) - l_begin(0);  // the largest line-of-sight block (the first)
  const double cells_own = 1.25 * hs_cells * double(nlos_max) / double(nlos);
  const double per_energy_bytes =
      8.0 * double(nsnap) * (double(nlos_max) * double(nloss) + cells_own * double(nprod + 1)) + 8.0 * cells_own;
  int group = ntask;
  if (per_energy_bytes > 0) {
    group = int(std::floor(td.memory_gb * 1e9 / per_energy_bytes));
    if (group < 1)
      throw std::runtime_error("time_dependent.memory_gb = " + std::to_string(td.memory_gb) + " GB is too small: one energy needs " +
                               std::to_string(per_energy_bytes / 1e9) + " GB per rank (use more ranks or more memory)");
    group = std::min(group, ntask);
  }
  const int ngroups = (ntask + group - 1) / group;

  if (root) {
    std::printf("Time-dependent maps: %d snapshots from %s to %s (dt = %.6g yr = %.4g months); observer at %g AU\n", nsnap,
                time_label(td.start).c_str(), time_label(td.end).c_str(), td.dt, td.dt * 12.0, r_obs);
    std::printf("  %d output windows of %.6g yr every %.6g yr, from %s to %s:", nwin, td.window, td.cadence,
                time_label(windows.front().a).c_str(), time_label(windows.back().b).c_str());
    for (int w = 0; w < std::min(nwin, 8); ++w) std::printf(" %s", windows[w].label.c_str());
    std::printf("%s\n", nwin > 8 ? " ..." : "");
    std::printf("  %d energies in %zu bands; survival %s\n", ntask, cfg.bands.size(),
                survival ? "along the ENA path" : "off");
    std::printf("  layout: %d MPI rank%s x %d threads; %d snapshots and %ld lines of sight per rank (at most); "
                "stored per rank: %.2f GB per energy, %d energ%s per pass, %d pass%s over the snapshots\n",
                nranks, nranks == 1 ? "" : "s", max_threads(), rounds, nlos_max, per_energy_bytes / 1e9, group,
                group == 1 ? "y" : "ies", ngroups, ngroups == 1 ? "" : "es");
    std::printf("  coverage: incomplete where the extinction of the missing heliosheath cells is >= %g", td.coverage_extinction);
    if (td.coverage_missing_flux > 0)
      std::printf(" or the missing share of the window's band flux is >= %g", td.coverage_missing_flux);
    std::printf("; their flux is %s\n", td.include_incomplete_lines ? "kept and flagged" : "set to 0 and flagged");
    std::fflush(stdout);
  }

  // ---- accumulators of this rank's lines of sight, per energy ----
  // acc_own[e]: [window][component][own line of sight]; cov_own[e]: [quantity][window][own line]
  // with quantity 0: max extinction of the missing cells, 1: missing flux, 2: total flux
  // (window-averaged, first snapshot)
  std::vector<std::vector<double>> acc_own(ntask, std::vector<double>(std::size_t(nwin) * ncomp * nlos_own, 0.0));
  std::vector<std::vector<double>> cov_own(ntask, std::vector<double>(std::size_t(3) * nwin * nlos_own, 0.0));
  auto acc_at = [&](int e, int w, int m, long li) -> double& {
    return acc_own[e][(std::size_t(w) * ncomp + m) * nlos_own + li];
  };
  auto cov_at = [&](int e, int q, int w, long li) -> double& {
    return cov_own[e][(std::size_t(q) * nwin + w) * nlos_own + li];
  };

  CellProduction cells;
  cells.resize(g, nprod);
  const auto t_all = Clock::now();
  double t_prod = 0, t_exch = 0, t_int = 0;

  for (int gi = 0; gi < ngroups; ++gi) {
    const int e_first = gi * group, e_last = std::min(ntask, e_first + group);  // [e_first, e_last)
    const int ng = e_last - e_first;
    if (ngroups > 1 && root) std::printf("Energy group %d of %d: energies %d-%d\n", gi + 1, ngroups, e_first + 1, e_last);

    // Data of this rank's lines of sight, for every snapshot and energy of the group:
    //  P[e][s]: running loss sums, [own line][k - k0];
    //  off[e][s]: start of each own line's cells in rec[e][s] (nlos_own + 1 entries);
    //  rec[e][s]: per emitting cell (k ascending): k, the production of each component,
    //  and, for the first snapshot, the streamline extinction (coverage test).
    std::vector<std::vector<std::vector<double>>> P(ng, std::vector<std::vector<double>>(nsnap));
    std::vector<std::vector<std::vector<long>>> off(ng, std::vector<std::vector<long>>(nsnap));
    std::vector<std::vector<std::vector<double>>> rec(ng, std::vector<std::vector<double>>(nsnap));
    auto rec_len = [&](int s) { return 1 + nprod + (s == 0 ? 1 : 0); };

    // ---- production phase: rounds over the snapshots of each rank ----
    for (int j = 0; j < rounds; ++j) {
      const int s = snap_of(me, j);
      std::vector<std::vector<double>> send(nranks);
      std::string err, line;
      const auto t_p = Clock::now();
      if (s >= 0) {
        try {
          const auto t_snap = Clock::now();
          ReadOptions ro;
          ro.verbose = false;
          Snapshot snap = read_snapshot(cfg, std::size_t(s), ro);
          const double t_read = since(t_snap);
          // region checks (as ena_run), summarized on the snapshot's line
          const std::string regions = gi == 0 ? report_regions(cfg, snap, "", false) : "";

          // loss sums of this snapshot, every energy of the group, all lines of sight
          std::vector<std::vector<double>> Ps(ng);
          if (survival) {
            const PlasmaFluid& pl = snap.plasma.at(0);
            for (int e = 0; e < ng; ++e) {
              std::vector<double>& Pe = Ps[e];
              Pe.assign(std::size_t(nlos) * nloss, 0.0);
              const double E = tasks[e_first + e].energy;
              ParallelErrors perr;
#pragma omp parallel for schedule(static)
              for (long los = 0; los < nlos; ++los) {
                const int t = int(los / g.nphi), pp = int(los % g.nphi);
                perr.capture(std::size_t(los), [&] {
                  double* p = Pe.data() + std::size_t(los) * nloss;
                  double prev = 0.0;
                  for (int k = k0; k <= k1; ++k) {
                    const int code = snap.region(t, pp, k);
                    // same sum as ionization_rate: regions 3 and 6, carried through other cells
                    if (code == 3 || code == 6) prev = loss_increment(pl, E, g.dr, t, pp, k) + prev;
                    p[k - k0] = prev;
                  }
                });
              }
              perr.rethrow();
            }
          }

          // streamline fields and production; pack each destination's lines of sight
          RunInfo info;
          const auto t_stream = Clock::now();
          StreamFields stream = build_stream_fields(cfg_prod, snap, std::size_t(s), info, false);
          const double t_st = since(t_stream);
          const auto t_flux = Clock::now();
          FluxContext ctx{cfg_prod, snap, stream, &cells};
          for (int e = 0; e < ng; ++e) {
            const EnergyTask& task = tasks[e_first + e];
            make_flux(ctx, task.energy, task.band, false);  // fills `cells`; its maps are not used
            const double sig_cx = xsection(task.energy);
            for (int d = 0; d < nranks; ++d) {
              std::vector<double>& b = send[d];
              const long a0 = l_begin(d), a1 = l_begin(d + 1);
              if (survival)
                b.insert(b.end(), Ps[e].begin() + std::size_t(a0) * nloss, Ps[e].begin() + std::size_t(a1) * nloss);
              for (long los = a0; los < a1; ++los) {
                const int t = int(los / g.nphi), pp = int(los % g.nphi);
                const std::size_t count_at = b.size();
                b.push_back(0.0);
                long count = 0;
                for (int k = 0; k < g.nr; ++k) {
                  bool any = false;
                  for (int m = 0; m < nprod; ++m) any = any || cells.comp(m)[g.index(t, pp, k)] != 0.0;
                  if (!any) continue;
                  if (survival && k > k1)
                    throw std::runtime_error("time-dependent method: emitting cell beyond the integration window");
                  b.push_back(double(k));
                  for (int m = 0; m < nprod; ++m) b.push_back(cells.comp(m)[g.index(t, pp, k)]);
                  if (s == 0) {  // streamline extinction, for the coverage test (as make_flux)
                    double ext = cfg.physics.extinction ? std::exp((stream.cool(t, pp, k) * task.v_cm) * sig_cx) : 1.0;
                    if (std::isnan(ext) || ext == 0.0) ext = 1e-10;
                    b.push_back(ext);
                  }
                  ++count;
                }
                b[count_at] = double(count);
              }
            }
          }
          char buf[512];
          const double t_cell0 = td.dt * double(s);
          std::snprintf(buf, sizeof buf, "  snapshot %d (%s-%s): read %.1f s, streamlines %.1f s, flux %.1f s%s%s%s", s,
                        time_label(td.start + t_cell0).c_str(), time_label(td.start + t_cell0 + td.dt).c_str(), t_read,
                        t_st, since(t_flux), nranks > 1 ? (" (rank " + std::to_string(me) + ")").c_str() : "",
                        regions.empty() ? "" : "; ", regions.c_str());
          line = buf;
        } catch (const std::exception& ex) {
          err = "snapshot " + std::to_string(s) + ": " + ex.what();
        }
      }
      t_prod += since(t_p);
      comm.agree(err);  // all ranks stop together if any snapshot failed

      // ---- exchange: every line of sight's data to the rank that integrates it ----
      const auto t_x = Clock::now();
      double send_mb = 0;
      for (const auto& v : send) send_mb += 8e-6 * double(v.size());
      std::vector<std::vector<double>> recv = comm.alltoallv(std::move(send));
      const double t_a2a = since(t_x);
      for (int src = 0; src < nranks; ++src) {
        const int ss = snap_of(src, j);
        if (ss < 0) continue;
        std::vector<double>& b = recv[src];
        std::size_t pos = 0;
        const int rl = rec_len(ss);
        for (int e = 0; e < ng; ++e) {
          if (survival) {
            P[e][ss].assign(b.begin() + pos, b.begin() + pos + std::size_t(nlos_own) * nloss);
            pos += std::size_t(nlos_own) * nloss;
          }
          std::vector<long>& o = off[e][ss];
          std::vector<double>& r = rec[e][ss];
          o.assign(nlos_own + 1, 0);
          r.clear();
          for (long li = 0; li < nlos_own; ++li) {
            const long count = long(b[pos++]);
            o[li] = long(r.size());
            r.insert(r.end(), b.begin() + pos, b.begin() + pos + std::size_t(count) * rl);
            pos += std::size_t(count) * rl;
          }
          o[nlos_own] = long(r.size());
        }
        if (pos != b.size()) throw std::logic_error("time-dependent exchange: unexpected message length");
        recv[src].clear();
        recv[src].shrink_to_fit();
      }
      t_exch += since(t_x);
      if (std::getenv("ENA_DEBUG_EXCHANGE"))
        std::fprintf(stderr, "rank %d round %d: sends %.0f MB, alltoallv %.2f s, total exchange %.2f s\n", me, j, send_mb, t_a2a, since(t_x));
      // log lines of this round, rank 0 prints them (highest snapshot first)
      const std::vector<std::string> lines = comm.gather_text(line);
      if (root) {
        for (int r = nranks - 1; r >= 0; --r)
          if (!lines[r].empty()) std::printf("%s\n", lines[r].c_str());
        std::fflush(stdout);
      }
    }

    // ---- integration phase: this rank's lines of sight over all snapshots ----
    const auto t_i = Clock::now();
    std::string err;
    for (int e = 0; e < ng && err.empty(); ++e) {
      const EnergyTask& task = tasks[e_first + e];
      const double E = task.energy;
      const double cell_years = g.dr / task.v_au_yr;
      ParallelErrors perr;
#pragma omp parallel for schedule(dynamic, 4)
      for (long li = 0; li < nlos_own; ++li) {
        perr.capture(std::size_t(li), [&] {
          std::vector<const double*> pl_los;
          if (survival)
            for (int q = 0; q < nsnap; ++q) pl_los.push_back(P[e][q].data() + std::size_t(li) * nloss);
          std::vector<double> part(std::size_t(nwin) * ncomp);
          for (int s = nsnap - 1; s >= 0; --s) {  // snapshot order of version 2.7 (identical sums)
            const double t_cell0 = td.dt * double(s);  // start of this time cell, years after T0
            const double t_emit = t_cell0 + 0.5 * td.dt;
            const bool first = s == 0;
            const int rl = rec_len(s);
            std::fill(part.begin(), part.end(), 0.0);
            // first snapshot: per emitting cell, its radius, extinction and flux (all components,
            // including the heliosheath-created PUIs, at this energy)
            std::vector<double> ck_r, ck_tau, ck_ext, ck_flux;
            const double* rr = rec[e][s].data();
            for (long c = off[e][s][li]; c < off[e][s][li + 1]; c += rl) {
              const int k = int(rr[c]);
              const double* j = rr + c + 1;
              double jtot = 0.0;
              for (int m = 0; m < nprod; ++m) jtot += j[m];
              double S = 1.0;
              if (survival && k >= k0) {
                const double beta = path_exponent(pl_los, k0, k, t_emit, td.dt, nsnap, cell_years);
                S = survival_prob(E, beta);
              }
              // arrival interval of this cell's ENAs (absolute years)
              const double tau = (g.r(k) - r_obs) / task.v_au_yr;
              const double lo = td.start + t_cell0 + tau, hi = lo + td.dt;
              int w_lo = int(std::floor((lo - td.output_start - td.window) / td.cadence));
              int w_hi = int(std::floor((hi - td.output_start) / td.cadence));
              w_lo = std::max(w_lo, 0);
              w_hi = std::min(w_hi, nwin - 1);
              for (int w = w_lo; w <= w_hi; ++w) {
                const double a = windows[w].a, b = windows[w].b;
                const double ov = std::min(hi, b) - std::max(lo, a);
                if (!(ov > 0)) continue;
                const double f = ov / (b - a);
                double* pw = part.data() + std::size_t(w) * ncomp;
                double tot = 0.0;
                for (int m = 0; m < nprod; ++m) {
                  const double v = j[m] * S;
                  pw[m] += v * f;
                  tot += v;
                }
                pw[ncomp - 1] += tot * f;
              }
              if (first) {
                ck_r.push_back(g.r(k));
                ck_tau.push_back(tau);
                ck_ext.push_back(j[nprod]);
                ck_flux.push_back(jtot * S);
              }
            }
            for (int w = 0; w < nwin; ++w)
              for (int m = 0; m < ncomp; ++m) acc_at(e_first + e, w, m, li) += part[std::size_t(w) * ncomp + m];
            if (first) {
              // Coverage. ENAs from cell c reaching the observer at time t were made at t - tau_c, so
              // they are missing (made before T0) while t < T0 + tau_c.
              //  - extinction test: cells missing at the window start a (r_c > r_cov), max extinction;
              //  - missing flux: each cell's flux times the fraction of the window it is missing,
              //    (T0 + tau_c - a)/(b - a) clipped to [0, 1], summed; the total flux alongside.
              double ftot = 0.0;
              for (double x : ck_flux) ftot += x;
              for (int w = 0; w < nwin; ++w) {
                const double a = windows[w].a, b = windows[w].b;
                const double r_cov = r_obs + task.v_au_yr * (a - td.start);
                double emax = 0.0, fmiss = 0.0;
                for (std::size_t c = 0; c < ck_r.size(); ++c) {
                  if (ck_r[c] > r_cov) emax = std::max(emax, ck_ext[c]);
                  const double frac = std::min(1.0, std::max(0.0, (td.start + ck_tau[c] - a) / (b - a)));
                  fmiss += ck_flux[c] * frac;
                }
                cov_at(e_first + e, 0, w, li) = emax;
                cov_at(e_first + e, 1, w, li) = fmiss;
                cov_at(e_first + e, 2, w, li) = ftot;
              }
            }
          }
        });
      }
      try {
        perr.rethrow();
      } catch (const std::exception& ex) {
        err = ex.what();
      }
    }
    comm.agree(err);
    t_int += since(t_i);
  }

  // ---- collect every rank's lines of sight on rank 0 ----
  std::vector<std::vector<MapSet>> acc;
  std::vector<std::vector<Field2D<double>>> cov_ext, cov_miss, cov_tot;
  if (root) {
    acc.assign(ntask, std::vector<MapSet>(nwin, MapSet(ncomp, Field2D<double>(g, 0.0))));
    cov_ext.assign(ntask, std::vector<Field2D<double>>(nwin, Field2D<double>(g, 0.0)));
    cov_miss = cov_ext;
    cov_tot = cov_ext;
  }
  for (int e = 0; e < ntask; ++e) {
    std::vector<double> mine = acc_own[e];
    mine.insert(mine.end(), cov_own[e].begin(), cov_own[e].end());
    acc_own[e].clear();
    acc_own[e].shrink_to_fit();
    const std::vector<std::vector<double>> all = comm.gather(mine);
    if (!root) continue;
    for (int r = 0; r < nranks; ++r) {
      const long a0 = l_begin(r), n = l_begin(r + 1) - a0;
      const std::vector<double>& b = all[r];
      for (int w = 0; w < nwin; ++w) {
        for (int m = 0; m < ncomp; ++m)
          for (long li = 0; li < n; ++li) acc[e][w][m].raw()[a0 + li] = b[(std::size_t(w) * ncomp + m) * n + li];
        const std::size_t c0 = std::size_t(nwin) * ncomp * n;
        for (long li = 0; li < n; ++li) {
          cov_ext[e][w].raw()[a0 + li] = b[c0 + (std::size_t(0) * nwin + w) * n + li];
          cov_miss[e][w].raw()[a0 + li] = b[c0 + (std::size_t(1) * nwin + w) * n + li];
          cov_tot[e][w].raw()[a0 + li] = b[c0 + (std::size_t(2) * nwin + w) * n + li];
        }
      }
    }
  }
  t_prod = comm.max_double(t_prod);
  t_exch = comm.max_double(t_exch);
  t_int = comm.max_double(t_int);
  if (!root) return;  // rank 0 writes the outputs
  std::printf("Phases (slowest rank): production %.1f s, exchange %.1f s, integration %.1f s\n", t_prod, t_exch, t_int);

  // ---- band sums, coverage and output ----
  std::filesystem::create_directories(opt.out_dir);
  const std::string summary_path = opt.out_dir + "/time_dependent_summary.txt";
  FILE* sum = std::fopen(summary_path.c_str(), "w");
  if (!sum) throw std::runtime_error("cannot write " + summary_path);
  std::fprintf(sum, "# Time-dependent ENA maps, ENA C++ %s\n", kVersion);
  std::fprintf(sum, "# input: %d snapshots, %s to %s, dt = %.8g yr; observer at %g AU\n", nsnap, time_label(td.start).c_str(),
               time_label(td.end).c_str(), td.dt, r_obs);
  std::fprintf(sum, "# coverage: a line of sight is incomplete where the extinction of the heliosheath cells beyond r_cov "
                    "is >= %g (extinction test)", td.coverage_extinction);
  if (td.coverage_missing_flux > 0)
    std::fprintf(sum, " or the share of the window's band flux made before the first snapshot is >= %g (missing-flux test)",
                 td.coverage_missing_flux);
  std::fprintf(sum, "\n# r_cov = observer + v (window start - input start), at the lowest energy of the band\n");
  std::fprintf(sum, "# band window_start window_end r_cov_AU incomplete_lines lines by_extinction by_missing_flux "
                    "max_extinction_missing max_missing_flux_fraction\n");
  long n_incomplete_maps = 0;
  for (int band : cfg.bands) {
    for (int w = 0; w < nwin; ++w) {
      BandAccumulator ba(cfg, band);
      Field2D<double> complete(g, 1.0), max_ext(g, 0.0), miss(g, 0.0), miss_num(g, 0.0), miss_den(g, 0.0);
      double v_min = 0.0;
      for (int e = 0; e < ntask; ++e) {
        if (tasks[e].band != band) continue;
        ba.add(tasks[e].index, tasks[e].energy, acc[e][w]);
        v_min = v_min == 0.0 ? tasks[e].v_au_yr : std::min(v_min, tasks[e].v_au_yr);
        // the missing share of the band: missing and total flux combined with the band weights of
        // the count rate, so each energy counts as much as it does in the band map
        const double wE = band_energy_weight(cfg, band, tasks[e].index, tasks[e].energy);
        for (std::size_t q = 0; q < max_ext.raw().size(); ++q) {
          max_ext.raw()[q] = std::max(max_ext.raw()[q], cov_ext[e][w].raw()[q]);
          miss_num.raw()[q] += wE * cov_miss[e][w].raw()[q];
          miss_den.raw()[q] += wE * cov_tot[e][w].raw()[q];
        }
      }
      for (std::size_t q = 0; q < miss.raw().size(); ++q)
        miss.raw()[q] = miss_den.raw()[q] > 0 ? miss_num.raw()[q] / miss_den.raw()[q] : 0.0;
      long incomplete = 0, by_ext = 0, by_flux = 0;
      double worst_ext = 0, worst_frac = 0;
      for (std::size_t q = 0; q < complete.raw().size(); ++q) {
        const bool bad_ext = max_ext.raw()[q] >= td.coverage_extinction;
        const bool bad_flux = td.coverage_missing_flux > 0 && miss.raw()[q] >= td.coverage_missing_flux;
        by_ext += bad_ext;
        by_flux += bad_flux;
        if (bad_ext || bad_flux) { complete.raw()[q] = 0.0; ++incomplete; }
        worst_ext = std::max(worst_ext, max_ext.raw()[q]);
        worst_frac = std::max(worst_frac, miss.raw()[q]);
      }
      const double r_cov = r_obs + v_min * (windows[w].a - td.start);
      BandResult B = ba.finish();
      // Regions with missing flux: set to 0 unless include_incomplete_lines (flagged either way)
      if (!td.include_incomplete_lines && incomplete)
        for (std::size_t q = 0; q < complete.raw().size(); ++q)
          if (complete.raw()[q] == 0.0)
            for (std::size_t m = 0; m < B.flux.size(); ++m) B.flux[m].raw()[q] = B.countrate[m].raw()[q] = 0.0;
      const std::string name = "flux_band" + band_label(band) + "_" + windows[w].label;
      ExtraOutput ex;
      ex.maps = {{"coverage_complete", &complete}, {"coverage_max_extinction_missing", &max_ext},
                 {"coverage_missing_flux_fraction", &miss}};
      ex.scalars = {{"window_start_yr", windows[w].a}, {"window_end_yr", windows[w].b}, {"input_start_yr", td.start},
                    {"input_end_yr", td.end}, {"input_dt_yr", td.dt}, {"observer_r_au", r_obs},
                    {"coverage_extinction", td.coverage_extinction}, {"coverage_missing_flux", td.coverage_missing_flux},
                    {"r_coverage_au", r_cov}, {"incomplete_lines", double(incomplete)},
                    {"incomplete_lines_extinction", double(by_ext)}, {"incomplete_lines_missing_flux", double(by_flux)}};
      ex.scalars.push_back({"include_incomplete_lines", td.include_incomplete_lines ? 1.0 : 0.0});
      ex.strings = {{"method", "time-dependent: production per snapshot (time cells of dt), survival along the ENA "
                               "path, flux averaged over the arrival window"},
                    {"coverage_note", "coverage_complete = 0 where ENAs made before the first snapshot are missing and "
                                      "the extinction of the missing heliosheath cells (max: "
                                      "coverage_max_extinction_missing) is >= coverage_extinction; "
                                      "coverage_missing_flux_fraction: the share of the window-averaged band flux "
                                      "(all components, including the heliosheath-created PUIs) made before the first "
                                      "snapshot, estimated from the first snapshot; also incomplete where it is >= "
                                      "coverage_missing_flux (if > 0). These lines of sight are regions with missing "
                                      "flux: all maps are 0 there unless include_incomplete_lines = 1"}};
      if (checkpoints_available()) write_band_flux(cfg, B, opt.out_dir, RunInfo{}, name, &ex);
      if (opt.legacy) {
        write_legacy_countrate(cfg, B, opt.out_dir, "_" + windows[w].label);
        // flag file in the layout of the legacy count-rate files: 1 = missing flux, 0 = complete
        char kappa[32];
        std::snprintf(kappa, sizeof kappa, "%4.2f", cfg.physics.kappa);
        const std::string fp = opt.out_dir + "/Opher_multiion_withneutralsMissingFluxFlagKappa" + kappa + "Eband" +
                               band_label(band) + "_" + windows[w].label + ".dat";
        FILE* ff = std::fopen(fp.c_str(), "w");
        if (!ff) throw std::runtime_error("cannot write " + fp);
        for (int t = 0; t < g.ntheta; ++t)
          for (int p = 0; p < g.nphi; ++p) std::fprintf(ff, "%d\n", complete(t, p) == 0.0 ? 1 : 0);
        std::fclose(ff);
      }
      std::fprintf(sum, "%3d %10s %10s %9.1f %6ld %6ld %6ld %6ld %10.4g %10.4g\n", band, time_label(windows[w].a).c_str(),
                   time_label(windows[w].b).c_str(), r_cov, incomplete, nlos, by_ext, by_flux, worst_ext, worst_frac);
      if (incomplete) {
        ++n_incomplete_maps;
        std::printf("WARNING (coverage): band %02d, window %s: %ld of %ld lines of sight are INCOMPLETE "
                    "(ENAs made before the first snapshot, %s, are missing).\n"
                    "    extinction test: %ld lines; the extinction of the heliosheath cells beyond %.0f AU reaches %.3g "
                    "(threshold %g)\n"
                    "    missing-flux test: %ld lines; up to %.3g of the window's band flux is missing (threshold %g%s)\n"
                    "    %s. See coverage_complete in %s.h5\n",
                    band, windows[w].label.c_str(), incomplete, nlos, time_label(td.start).c_str(), by_ext, r_cov, worst_ext,
                    td.coverage_extinction, by_flux, worst_frac, td.coverage_missing_flux,
                    td.coverage_missing_flux > 0 ? "" : ", test off",
                    td.include_incomplete_lines ? "Their (incomplete) flux is kept (include_incomplete_lines = true)"
                                                : "Their flux is set to 0 in every component map",
                    name.c_str());
      }
    }
  }
  std::fclose(sum);
  std::printf("Time-dependent maps: %zu bands x %d windows written in %s (%ld maps with incomplete coverage); %.1f s\n",
              cfg.bands.size(), nwin, opt.out_dir.c_str(), n_incomplete_maps, since(t_all));
  std::printf("Coverage summary: %s\n", summary_path.c_str());
}

}  // namespace ena
