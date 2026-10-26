#include "ena/checkpoints.hpp"

#include <cstdio>
#include <stdexcept>

#ifdef ENA_HAVE_HDF5
#include "ena/h5io.hpp"
#endif

namespace ena {

namespace {
[[maybe_unused]] std::string band_label(int band) { char b[16]; std::snprintf(b, sizeof b, "%02d", band); return b; }
[[maybe_unused]] std::string energy_label(double e) { char b[32]; std::snprintf(b, sizeof b, "%.5f", e); return b; }
}  // namespace

#ifdef ENA_HAVE_HDF5
bool checkpoints_available() { return true; }

void write_c0_config(const RunConfig& cfg, const std::string& dir) {
  CheckpointWriter w(dir, "C0_config");
  const Grid& g = cfg.grid;
  // Types as in the reference dumps: single-precision floats for the grid parameters.
  w.scalar_float("dtheta", float(g.dtheta));
  w.scalar_float("dphi", float(g.dphi));
  w.scalar_float("dr", float(g.dr));
  w.scalar_float("ntheta", float(g.ntheta));
  w.scalar_float("nphi", float(g.nphi));
  w.scalar_float("r_inner", float(g.r_inner));
  w.scalar_float("r_lower", float(cfg.r_lower));
  w.scalar_float("r_upper", float(cfg.r_upper));
  w.scalar_int("radial_range", cfg.radial_range);
  w.scalar_int("moscow_regs", cfg.moscow() ? 1 : 0);
  w.scalar_int("moscow_temp", cfg.moscow() ? 1 : 0);
}

void write_c1_inputs(const RunConfig& cfg, const Snapshot& s, const RadialLimits& lim, const std::string& dir) {
  {
    CheckpointWriter w(dir, "C1a_limits");
    w.field("boundaries_arr", s.region);
    w.field("r_limits_lower", lim.lower);
    w.field("r_limits_upper", lim.upper);
    w.field("bx", s.Bx);
    w.field("by", s.By);
    w.field("bz", s.Bz);
  }
  {
    CheckpointWriter w(dir, "C1b_inputs");
    const PlasmaFluid& p = s.plasma.at(0);
    auto kms = [](const Field3D<double>& f) {  // velocities in km/s, as in the reference dumps
      Field3D<double> out = f;
      for (double& x : out.raw()) x /= 1e5;
      return out;
    };
    w.field("density", p.n);
    w.field("temp", p.T);
    w.field("vr", kms(p.vr));
    w.field("vtheta", kms(p.vtheta));
    w.field("vphi", kms(p.vphi));
    w.field("vt", kms(p.vt));
    if (!s.coords.theta.empty()) {
      w.field("theta", s.coords.theta);
      w.field("phi", s.coords.phi);
      w.field("r", s.coords.r);
    }
    if (s.neutrals.size() == 5) {  // multi-ion / extended single-ion: total H and populations 1-4, [theta][phi][r][population]
      const Grid& g = cfg.grid;
      auto stack = [&](Field3D<double> NeutralPopulation::*m) {
        std::vector<double> out(g.size() * 5);
        for (std::size_t c = 0; c < g.size(); ++c)
          for (std::size_t q = 0; q < 5; ++q) out[c * 5 + q] = (s.neutrals[q].*m).raw()[c];
        return out;
      };
      w.field4("nh_density", stack(&NeutralPopulation::n), g, 5);
      w.field4("nh_vr", stack(&NeutralPopulation::vr), g, 5);
      w.field4("nh_vt", stack(&NeutralPopulation::vt), g, 5);
      w.field4("nh_temp", stack(&NeutralPopulation::T), g, 5);
    } else {
      w.field("nh_density", s.neutrals.at(0).n);
      w.field("nh_temp", s.neutrals.at(0).T);
    }
  }
}


void write_c2_stream_points(const StreamPoints& s, const std::string& dir) {
  CheckpointWriter w(dir, "C2_stream_points");
  auto dimn = [](std::size_t n) { return "n" + std::to_string(n); };
  for (auto [name, vec] : std::vector<std::pair<const char*, const std::vector<double>*>>{
           {"x", &s.x}, {"y", &s.y}, {"z", &s.z}, {"ux", &s.ux}, {"uy", &s.uy}, {"uz", &s.uz}, {"nh", &s.nh}, {"np", &s.np},
           {"tp", &s.tp}, {"vr", &s.vr}, {"vt", &s.vt}, {"theta", &s.theta}, {"phi", &s.phi}, {"r", &s.r}, {"v", &s.v},
           {"cs", &s.cs}, {"m2", &s.m2}, {"ds", &s.ds}, {"tau", &s.tau}, {"ext", &s.ext}, {"ntr0", &s.ntr0},
           {"nref0", &s.nref0}, {"nen0", &s.nen0}, {"etr0", &s.etr0}, {"eref0", &s.eref0}, {"een0", &s.een0},
           {"em0", &s.em0}, {"kappa_ref0", &s.kappa_ref0}, {"vshock0", &s.vshock0}, {"esw0", &s.esw0}})
    w.vector(name, *vec, dimn(vec->size()));
  if (!s.tssw.empty())  // multi-ion
    for (auto [name, vec] : std::vector<std::pair<const char*, const std::vector<double>*>>{
             {"npu3", &s.npu3}, {"tpu3", &s.tpu3}, {"vrpu3", &s.vrpu3}, {"vtpu3", &s.vtpu3}, {"tssw", &s.tssw},
             {"u_inj0", &s.u_inj0}})
      if (!vec->empty()) w.vector(name, *vec, dimn(vec->size()));
  for (auto [name, vec] : std::vector<std::pair<const char*, const std::vector<std::int32_t>*>>{
           {"r_in", &s.r_in}, {"t_in", &s.t_in}, {"p_in", &s.p_in}}) {
    std::vector<long> l(vec->begin(), vec->end());
    w.vector(name, l, dimn(l.size()));
  }
}

void write_c3_stream(const StreamFields& f, const std::string& dir, const std::string& name) {
  CheckpointWriter w(dir, name);
  for (auto [n, fld] : std::vector<std::pair<const char*, const Field3D<double>*>>{
           {"cool", &f.cool}, {"vx", &f.vx}, {"vy", &f.vy}, {"vz", &f.vz}, {"ntr", &f.ntr}, {"nref", &f.nref},
           {"nen", &f.nen}, {"etr", &f.etr}, {"eref", &f.eref}, {"een", &f.een}, {"em", &f.em},
           {"kappa_ref", &f.kappa_ref}, {"vshock", &f.vshock}, {"esw", &f.esw}, {"ptr_ct", &f.ptr_ct},
           {"ptrn_ct", &f.ptrn_ct}})
    w.field(n, *fld);
  if (f.has_u_inj) w.field("u_inj", f.u_inj);
}

namespace {
// [component][theta][phi], written with dim_names ['n5','theta','phi'].
void write_components(CheckpointWriter& w, const std::string& name, const MapSet& c) {
  const int nc = int(c.size()), nt = c[0].ntheta(), np = c[0].nphi();
  Field3D<double> stack;
  Grid g; g.ntheta = nc; g.nphi = nt; g.nr = np;
  stack.resize(g);
  for (int m = 0; m < nc; ++m)
    for (int t = 0; t < nt; ++t)
      for (int p = 0; p < np; ++p) stack(m, t, p) = c[m](t, p);
  w.field_dims(name, stack, {"n" + std::to_string(nc), "theta", "phi"});
}
}  // namespace

void write_c4_energy(const EnergyResult& r, int band, const std::string& dir) {
  CheckpointWriter w(dir, "C4_eband" + band_label(band) + "_E" + energy_label(r.energy));
  w.scalar_float("energy", float(r.energy));
  w.scalar_int("eband", band);
  w.scalar("vena", r.vena);
  w.scalar("sig_cx", r.sig_cx);
  write_components(w, "flux", r.flux);
  if (!r.has_diagnostics) return;
  if (r.extended_si) {  // extended single-ion: populations as passed to the flux calls, per cell
    // slots: transmitted, reflected, DSA, heliosheath PUIs, SW, region-6 accelerated ions
    Grid g; g.ntheta = r.ext.ntheta(); g.nphi = r.ext.nphi(); g.nr = r.ext.nr();
    w.field("bet", r.bet);
    w.field("ext", r.ext);
    w.field("lc", r.lc);
    w.field4("pops_n", r.pops_n, g, 6);
    w.field4("pops_t", r.pops_t, g, 6);
    w.field4("pops_kappa", r.pops_kappa, g, 6);
    w.field4("zank_par", r.zank_par, g, 3);
    w.field4("cell_par", r.cell_par, g, 4);
    for (auto [n, f] : std::vector<std::pair<const char*, const Field3D<double>*>>{
             {"v_plas", &r.v_plas}, {"flux_test1", &r.flux_test1}, {"flux_test2", &r.flux_test2},
             {"flux_test3", &r.flux_test3}, {"flux_test_sw", &r.flux_test_sw}})
      w.field(n, *f);
    w.field("shock_strength", r.shock_strength);
    return;
  }
  if (r.multi_ion) {  // multi-ion: populations as passed to the flux calls, per cell
    Grid g; g.ntheta = r.ext.ntheta(); g.nphi = r.ext.nphi(); g.nr = r.ext.nr();
    w.field("bet", r.bet);
    w.field("ext", r.ext);
    w.field("lc", r.lc);
    w.field4("pops_n", r.pops_n, g, 5);
    w.field4("pops_t", r.pops_t, g, 5);
    w.field4("pops_kappa", r.pops_kappa, g, 5);
    w.field4("zank_par", r.zank_par, g, 4);
    w.field4("cell_par", r.cell_par, g, 4);
    w.field("belt_mask", r.belt_mask);
    for (auto [n, f] : std::vector<std::pair<const char*, const Field3D<double>*>>{
             {"v_plas", &r.v_plas}, {"flux_test1", &r.flux_test1}, {"flux_test2", &r.flux_test2},
             {"flux_test3", &r.flux_test3}, {"flux_test_sw", &r.flux_test_sw}})
      w.field(n, *f);
    w.field("shock_strength", r.shock_strength);
    return;
  }
  w.field("bet", r.bet);
  w.field("ext", r.ext);
  w.field("lc", r.lc);
  w.field("lc_test", r.lc_test);
  for (auto [n, f] : std::vector<std::pair<const char*, const Field3D<double>*>>{
           {"v_plas", &r.v_plas}, {"dens_test1", &r.dens_test1}, {"temp_test1", &r.temp_test1},
           {"dens_test2", &r.dens_test2}, {"temp_test2", &r.temp_test2}, {"dens_test3", &r.dens_test3},
           {"temp_test3", &r.temp_test3}, {"dens_test_sw", &r.dens_test_sw}, {"temp_test_sw", &r.temp_test_sw},
           {"temp_test_p", &r.temp_test_p}, {"f_test1", &r.f_test1}, {"flux_test1", &r.flux_test1},
           {"flux_test2", &r.flux_test2}, {"flux_test3", &r.flux_test3}, {"flux_test_sw", &r.flux_test_sw}})
    w.field(n, *f);
  w.field("shock_strength", r.shock_strength);
}

void write_c5_band(const BandResult& b, const std::string& dir) {
  CheckpointWriter w(dir, "C5_eband" + band_label(b.band));
  write_components(w, "countrate_map", b.countrate);
  w.scalar_int("eband", b.band);
  w.scalar_float("eres", float(b.eres));
  w.scalar_float("emin", float(b.emin));
  w.scalar_float("emax", float(b.emax));
  w.scalar_float("ecen", float(b.ecen));
  w.scalar_float("g", float(b.g));
  w.scalar_float("denergy", float(b.denergy));
}

void write_band_flux(const RunConfig& cfg, const BandResult& b, const std::string& dir, const RunInfo& info,
                     const std::string& file_name, const ExtraOutput* extra) {
  CheckpointWriter w(dir, file_name.empty() ? "flux_band" + band_label(b.band) : file_name);
  const ComponentSet& cs = components(cfg);
  for (int m = 0; m < cs.size(); ++m) w.field(cs.names[m], b.flux[m]);
  w.string("model", cfg.populations);
  w.scalar_int("band", b.band);
  w.scalar("emin_kev", b.emin);
  w.scalar("emax_kev", b.emax);
  w.scalar("ecentral_kev", b.ecen);
  w.scalar("geometric_factor", b.g);
  w.string("units", "ENA flux, cm^-2 s^-1 sr^-1 keV^-1, averaged over the band");
  w.string("version", kVersion);
  w.string("settings", settings_summary(cfg));
  if (info.r_ts_reference > 0) {
    w.scalar("r_ts_v2_au", info.r_ts_reference);
    w.string("r_ts_v2_note", info.r_ts_note);
    w.scalar_int("streamlines", info.streamlines);
    w.scalar_int("streamlines_with_ts", info.streamlines_with_ts);
  }
  w.string("grid_note", "theta = i*dtheta (i = 0: north pole), phi = j*dphi");
  if (extra) {
    for (auto& [n, f] : extra->maps) w.field(n, *f);
    for (auto& [n, v] : extra->scalars) w.scalar(n, v);
    for (auto& [n, v] : extra->strings) w.string(n, v);
  }
}

#else
bool checkpoints_available() { return false; }
void write_c0_config(const RunConfig&, const std::string&) { throw std::runtime_error("built without HDF5"); }
void write_c1_inputs(const RunConfig&, const Snapshot&, const RadialLimits&, const std::string&) {
  throw std::runtime_error("built without HDF5");
}
void write_c2_stream_points(const StreamPoints&, const std::string&) { throw std::runtime_error("built without HDF5"); }
void write_c3_stream(const StreamFields&, const std::string&, const std::string&) {
  throw std::runtime_error("built without HDF5");
}
void write_c4_energy(const EnergyResult&, int, const std::string&) { throw std::runtime_error("built without HDF5"); }
void write_c5_band(const BandResult&, const std::string&) { throw std::runtime_error("built without HDF5"); }
void write_band_flux(const RunConfig&, const BandResult&, const std::string&, const RunInfo&, const std::string&,
                     const ExtraOutput*) {
  throw std::runtime_error("built without HDF5");
}
#endif

void write_legacy_countrate(const RunConfig& cfg, const BandResult& b, const std::string& dir, const std::string& suffix) {
  char kappa[32];
  std::snprintf(kappa, sizeof kappa, "%4.2f", cfg.physics.kappa);  // f4.2, as in the legacy file names
  const ComponentSet& cs = components(cfg);
  for (int m = 0; m < cs.size(); ++m) {
    std::string path = dir + "/Opher_multiion_withneutrals" + cs.tags[m] + "CountRateKappa" + kappa + "Eband" +
                       band_label(b.band) + suffix + ".dat";
    FILE* f = std::fopen(path.c_str(), "w");
    if (!f) throw std::runtime_error("cannot write " + path);
    const auto& c = b.countrate[m];
    for (int t = 0; t < c.ntheta(); ++t)
      for (int p = 0; p < c.nphi(); ++p) std::fprintf(f, "%24.16e\n", c(t, p));
    std::fclose(f);
  }
}

}  // namespace ena
