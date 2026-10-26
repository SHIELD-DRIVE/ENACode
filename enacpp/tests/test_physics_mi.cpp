// Compares the multi-ion physics functions with IDL's dump_unit_tests_mi output
// (U_physics_mi.h5): test_physics_mi PATH/U_physics_mi.h5 [rtol]
// Version 2 computes in double precision, so it agrees with IDL (which mixes
// single precision in) to about 1e-7; the default tolerance is 1e-5.
// Also takes the extended single-ion U_physics_esi.h5 (same inputs plus
// calc_diff_flux_neutrals, which is then checked too).
#include <hdf5.h>

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "ena/physics.hpp"

using namespace ena;

static std::vector<double> rd(hid_t f, const char* name) {
  hid_t d = H5Dopen2(f, name, H5P_DEFAULT);
  if (d < 0) { std::fprintf(stderr, "missing %s\n", name); std::exit(2); }
  hid_t s = H5Dget_space(d);
  hssize_t n = H5Sget_simple_extent_npoints(s);
  std::vector<double> v(n);
  H5Dread(d, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, v.data());
  H5Sclose(s); H5Dclose(d);
  return v;
}

int main(int argc, char** argv) {
  if (argc < 2) { std::fprintf(stderr, "usage: %s U_physics_mi.h5 [rtol]\n", argv[0]); return 2; }
  const double rtol = argc > 2 ? std::atof(argv[2]) : 1e-5;
  hid_t f = H5Fopen(argv[1], H5F_ACC_RDONLY, H5P_DEFAULT);
  if (f < 0) return 2;
  auto e = rd(f, "in_e"), tp = rd(f, "in_tp"), np = rd(f, "in_np"), nn = rd(f, "in_nn"), ur = rd(f, "in_ur"),
       un = rd(f, "in_un"), urh = rd(f, "in_urh"), vdf = rd(f, "in_vdiff"), bet = rd(f, "in_bet"),
       kk = rd(f, "in_kappa_k"), kap = rd(f, "in_kappa"), del = rd(f, "in_del"), el = rd(f, "in_el"),
       em = rd(f, "in_em"), qr = rd(f, "in_q_r"), vz = rd(f, "in_vz"), vsh = rd(f, "in_vshock");
  std::vector<std::pair<std::string, std::vector<double>>> ref = {
      {"vrel", rd(f, "vrel")}, {"zank_fn", rd(f, "zank_fn")}, {"fluxconstants", rd(f, "fluxconstants")},
      {"calc_diff_flux", rd(f, "calc_diff_flux")}, {"calc_diff_flux_zank", rd(f, "calc_diff_flux_zank")}};
  std::vector<double> en, tn, urn, unn;
  if (H5Lexists(f, "calc_diff_flux_neutrals", H5P_DEFAULT) > 0) {  // U_physics_esi.h5
    ref.push_back({"calc_diff_flux_neutrals", rd(f, "calc_diff_flux_neutrals")});
    en = rd(f, "in_en"); tn = rd(f, "in_tn"); urn = rd(f, "in_urn"); unn = rd(f, "in_unn");
  }
  H5Fclose(f);

  PhysicsSettings s;  // survival on, dr = 2 AU
  int bad = 0;
  for (auto& [name, r] : ref) {
    double worst = 0; std::size_t wi = 0; int nan_mismatch = 0, used = 0;
    for (std::size_t i = 0; i < e.size(); ++i) {
      FluxInputs in{e[i], tp[i], ur[i], un[i], nn[i], np[i], kap[i], bet[i]};
      double c;
      if (name == "vrel") c = vrel(tp[i], -energy_velocity(e[i]), urh[i]);
      else if (name == "zank_fn") c = zank_fn(vz[i], np[i], del[i], el[i], em[i], qr[i] > 1.0 ? qr[i] : 5.0);
      else if (name == "fluxconstants") c = flux_constants_vdiff(e[i], nn[i], bet[i], vdf[i], s);
      else if (name == "calc_diff_flux") c = calc_diff_flux_mi(in, urh[i], s);
      else if (name == "calc_diff_flux_neutrals") c = calc_diff_flux_neutrals(en[i], tn[i], urn[i], unn[i], nn[i]);
      else {
        if (qr[i] > 1.0) continue;  // q = 5 is used: compare the q = 5 rows only
        c = calc_diff_flux_zank(in, del[i], std::max(el[i], 0.3), em[i], vsh[i], s);
      }
      ++used;
      if (std::isnan(c) || std::isnan(r[i])) { nan_mismatch += std::isnan(c) != std::isnan(r[i]); continue; }
      const double scale = std::max(std::fabs(c), std::fabs(r[i]));
      const double rel = scale > 0 ? std::fabs(c - r[i]) / scale : 0.0;
      if (rel > worst) { worst = rel; wi = i; }
    }
    const bool ok = worst <= rtol && nan_mismatch == 0;
    bad += !ok;
    std::printf("%-20s %3d values, max rel diff %.3e at %3zu%s  %s\n", name.c_str(), used, worst, wi,
                nan_mismatch ? "  NaN MISMATCH" : "", ok ? "ok" : "FAIL");
  }
  return bad ? 1 : 0;
}
