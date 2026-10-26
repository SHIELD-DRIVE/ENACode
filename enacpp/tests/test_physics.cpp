// Compares the C++ physics functions with IDL's dump_unit_tests output
// (U_physics.h5): test_physics PATH/U_physics.h5 [rtol]
// The C++ code is double precision throughout, while IDL uses single-precision
// literals and !pi, so agreement is to about 1e-7.
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
  if (argc < 2) { std::fprintf(stderr, "usage: %s U_physics.h5 [rtol]\n", argv[0]); return 2; }
  double rtol = argc > 2 ? std::atof(argv[2]) : 1e-5;
  hid_t f = H5Fopen(argv[1], H5F_ACC_RDONLY, H5P_DEFAULT);
  if (f < 0) return 2;
  auto e = rd(f, "in_e"), tp = rd(f, "in_tp"), np = rd(f, "in_np"), nn = rd(f, "in_nn"), ur = rd(f, "in_ur"),
       un = rd(f, "in_un"), bet = rd(f, "in_bet"), kk = rd(f, "in_kappa_k"), kap = rd(f, "in_kappa"),
       del = rd(f, "in_del"), el = rd(f, "in_el"), em = rd(f, "in_em"), vz = rd(f, "in_vz"), vsh = rd(f, "in_vshock");
  std::vector<std::pair<std::string, std::vector<double>>> ref = {
      {"energy_velocity", rd(f, "energy_velocity")}, {"xsection", rd(f, "xsection")}, {"plasma_frame", rd(f, "plasma_frame")},
      {"plasma_frame2", rd(f, "plasma_frame2")}, {"maxwellian_fn", rd(f, "maxwellian_fn")}, {"kappa_fn", rd(f, "kappa_fn")},
      {"zank_fn", rd(f, "zank_fn")}, {"survivalprob", rd(f, "survivalprob")}, {"fluxconstants", rd(f, "fluxconstants")},
      {"calc_diff_flux", rd(f, "calc_diff_flux")}, {"calc_diff_flux_zank", rd(f, "calc_diff_flux_zank")}};
  H5Fclose(f);

  PhysicsSettings s;  // survival on, dr = 2 AU
  int bad = 0;
  for (auto& [name, r] : ref) {
    double worst = 0; std::size_t wi = 0;
    for (std::size_t i = 0; i < e.size(); ++i) {
      double vp = plasma_frame(e[i], ur[i], un[i]);
      double c;
      if (name == "energy_velocity") c = energy_velocity(e[i]);
      else if (name == "xsection") c = xsection(e[i]);
      else if (name == "plasma_frame") c = vp;
      else if (name == "plasma_frame2") c = plasma_frame2(e[i], ur[i], un[i], tp[i]);
      else if (name == "maxwellian_fn") c = maxwellian_fn(tp[i], np[i], vp);
      else if (name == "kappa_fn") c = kappa_fn(tp[i], np[i], vp, kk[i]);
      else if (name == "zank_fn") c = zank_fn(vz[i], np[i], del[i], el[i], em[i]);
      else if (name == "survivalprob") c = survival_prob(e[i], bet[i]);
      else if (name == "fluxconstants") c = flux_constants(e[i], nn[i], bet[i], s);
      else {
        FluxInputs in{e[i], tp[i], ur[i], un[i], nn[i], np[i], kap[i], bet[i]};
        c = name == "calc_diff_flux" ? calc_diff_flux(in, s) : calc_diff_flux_zank(in, del[i], el[i], em[i], vsh[i], s);
      }
      double scale = std::max(std::fabs(c), std::fabs(r[i]));
      double rel = scale > 0 ? std::fabs(c - r[i]) / scale : 0.0;
      if (rel > worst) { worst = rel; wi = i; }
    }
    bool ok = worst <= rtol;
    bad += !ok;
    std::printf("%-20s max rel diff %.3e at %3zu  %s\n", name.c_str(), worst, wi, ok ? "ok" : "FAIL");
  }
  return bad ? 1 : 0;
}
