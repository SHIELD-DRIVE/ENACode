#pragma once
// Run configuration, read from a TOML file (a subset of TOML: tables,
// arrays of tables, strings, numbers, booleans and flat arrays).

#include <cmath>
#include <map>
#include <string>
#include <vector>

#include "ena/grid.hpp"

namespace ena {

// ---------------------------------------------------------------------------
// Minimal TOML value tree
// ---------------------------------------------------------------------------
struct TomlValue {
  enum Kind { None, String, Number, Bool, Array } kind = None;
  std::string str;
  double num = 0.0;
  bool boolean = false;
  std::vector<TomlValue> arr;
};

struct TomlTable {
  std::map<std::string, TomlValue> values;
  std::map<std::string, TomlTable> tables;
  std::map<std::string, std::vector<TomlTable>> table_arrays;

  bool has(const std::string& key) const { return values.count(key) != 0; }
  const TomlTable* table(const std::string& name) const;
  double get_number(const std::string& key, double fallback) const;
  double require_number(const std::string& key) const;
  std::string get_string(const std::string& key, const std::string& fallback) const;
  bool get_bool(const std::string& key, bool fallback) const;
  std::vector<std::string> get_strings(const std::string& key) const;
};

TomlTable parse_toml(const std::string& text, const std::string& source_name = "<string>");
TomlTable parse_toml_file(const std::string& path);

// Code version, written into the output files. Version 1.0 (archived) has the
// faithful IDL mode and per-quirk switches; from 2.0 the accepted fixes are
// built in (CHANGES_FROM_IDL.md).
constexpr const char* kVersion = "2.8";

// ---------------------------------------------------------------------------
// Input description
// ---------------------------------------------------------------------------
struct ColumnMap {
  std::vector<std::string> columns;       // role of each column, by position
  std::vector<std::string> units;         // optional, same length ("" = default)
  int column(const std::string& role) const;  // -1 if absent
  double velocity_factor(const std::string& role) const;  // to cm/s
};

struct SnapshotFiles {
  double time = 0.0;       // steady state: a label; time-dependent: start of the snapshot's time cell (years)
  bool has_time = false;   // 'time' given in the configuration
  std::string plasma, neutral, boundaries, streamlines;  // ASCII inputs
  std::string file;                                      // HDF5 input (one file per snapshot)
};

struct InputSpec {
  std::string dir = ".";
  std::string format = "ascii";  // "ascii" or "hdf5"
  bool allow_incomplete = false;  // permit missing grid cells (testing on file slices)
  ColumnMap plasma, neutral, boundaries, streamlines;
  std::vector<SnapshotFiles> snapshots;
  std::string path(const std::string& file) const;
};

// Instrument energy bands (IBEX-Lo 0-7, IBEX-Hi 8-13, INCA 14-17).
struct SensorTable {
  std::vector<double> emin, emax, ecen, g;   // keV, keV, keV, geometric factor
  std::map<int, std::vector<double>> weights;  // band -> sub-band response weights
  int central_first = 14, central_last = 20;   // bands evaluated at the central energy only
  bool cg_lo = false;                          // IBEX-Lo (bands 0-7) at the CG-corrected central energy only (multi-ion)
  std::vector<int> unusable = {8};             // faulty bands, rejected
  bool central_only(int band) const { return (band >= central_first && band <= central_last) || (cg_lo && band < 8); }
  static SensorTable single_ion_default();
  static SensorTable multi_ion_default();      // IBEX-Lo CG tail centres, Ultra/INCA bands 14-20
};

// Multi-ion model (Kornbleuth+2026): three ion
// fluids (SW; pu3, PUIs picked up in the supersonic wind; pu2, PUIs created in
// the heliosheath). The version 1.1 archive reproduces the original IDL code.
// Extended single-ion model: the ions of region 6 (open heliosheath), set
// separately from region 3 ([region6] in the configuration). Every default
// reproduces the region-3 treatment (version 2.4 maps).
//  - kappas of the four Maxwellian/kappa populations (kappa <= 1.5: Maxwellian;
//    at fixed temperature a kappa changes the shape, not the pressure);
//  - accelerated ions: a fraction accelerated_fraction (f) of the transmitted
//    PUIs, reflected PUIs, heliosheath PUIs and SW (the same f for each) is
//    accelerated, taking a fraction accelerated_energy_share (e) of their
//    pressure, with kappa accelerated_kappa. The remaining ions keep the rest:
//    n_i -> (1 - f) n_i, T_i -> T_i (1 - e)/(1 - f). Density and pressure are
//    conserved; e = f keeps the temperatures, e > f makes the accelerated ions
//    hotter. accelerated_from sets how the accelerated ions are formed:
//     "each_population" (default): each population keeps its own accelerated
//      part, n = f n_i at T = (e/f) T_i, so the accelerated ions are always
//      e/f times hotter than the ions they come from;
//     "pooled": one population, n = f sum n_i at T = e sum n_i T_i / (f sum n_i),
//      i.e. e/f times the density-weighted mean temperature. In the heliosheath
//      most ions are the cool SW and heliosheath PUIs, so this can be cooler
//      than the transmitted/reflected PUIs and lower the keV flux.
//    Both are written as one component, region6_accelerated. The DSA-accelerated
//    PUIs are not changed.
struct Region6Config {
  double kappa_transmitted = -1.0;      // default: physics.kappa_pui
  double kappa_reflected = -1.0;        // default: physics.kappa_ref
  double kappa_heliosheath_pui = -1.0;  // default: physics.kappa_oh
  double kappa_sw = -1.0;               // default: physics.kappa_oh
  double accelerated_fraction = 0.0;
  double accelerated_energy_share = 0.0;
  double accelerated_kappa = 1.6;
  std::string accelerated_from = "each_population";  // or "pooled"
};

struct MultiIonConfig {
  // Reconnection region ("INCA belt"): region-3 cells with x > 0,
  // r <= belt_r_max, m vA^2 / kT_pu2 >= belt_ratio_min, m vA^2 < belt_mva2_max
  // keV, and phi < belt_phi_min or phi > belt_phi_max. There the pu2 and SW
  // slots are replaced by two kappa populations: belt_pu2_frac and
  // belt_accel_frac of the total density, at belt_t_pu2_factor and
  // belt_t_sw_factor times T_pu2. The radius and longitude limits are free
  // parameters: belt_r_max <= 0 removes the radius limit,
  // belt_phi_min >= belt_phi_max the longitude limit.
  bool belt = true;
  double belt_r_max = 550.0;
  double belt_phi_min = 30.0, belt_phi_max = 330.0;
  double belt_ratio_min = 1.0, belt_mva2_max = 0.5;
  double belt_accel_frac = 0.00087, belt_pu2_frac = 0.013;
  double belt_t_pu2_factor = 53.0, belt_t_sw_factor = 177.0;
  double belt_kappa_pu2 = 8.6, belt_kappa_sw = 8.15;
  // DSA injection energy E_l = 0.5 m_p u^2:
  //  "line_of_sight" (default): u = SW speed at radial index injection_k
  //   (4: 10 AU) on the cell's own line of sight. E_l is then about 1-2 keV,
  //   in line with the range of Kornbleuth+2023 (about 0.8-1.2 keV).
  //  "streamline" (option): u = SW speed at the last region-4 point
  //   before the streamline enters region 3/6, carried along that streamline.
  //   The wind there is slowed by PUI mass loading (E_l about 0.5-0.8 keV):
  //   accelerated flux x0.44-0.49 at 3-8 keV, a new component at 0.5-0.9 keV.
  std::string injection = "line_of_sight";
  int injection_k = 4;
  // TS partition fits (set where the streamline enters region 3/6 from region
  // 4, as for the single-ion model). r' = ln(x), x capped at ts_cap:
  //  "compressed" (default): x = ((r - r_TS,V2)*ts_compress + r_TS,V2)/r_TS,V2,
  //   tailored to the hybrid model: it compresses the MI solution's wider range of TS
  //   distances onto the range the fits were derived for.
  //  "ratio": x = r / r_TS,V2, strictly the MI solution's TS distances (for
  //   an MI solution whose tail TS is closer to the hybrid model's).
  std::string ts_scaling = "compressed";
  double ts_compress = 38.0 / 90.0;
  double ts_cap = 1.62;
  // v_shock = "first_below" (physics.vshock): the pu3 speed at the first point
  // past the region crossing below this speed.
  double vshock_vmin_kms = 215.0;
};

struct PhysicsConfig {
  bool survival = true;          // ENA survival probability
  bool extinction = true;
  bool cold_electrons = true;    // cold-electron correction of the ion temperature
  double kappa = 1.0, kappa_pui = 1.0, kappa_ref = 1.0, kappa_ism = 1.0, kappa_oh = 1.0;
  // Extended single-ion: kappa of the heliosheath-created PUIs in region 3.
  // Default (unset): kappa_oh (kappaIHS = kappaOH).
  double kappa_heliosheath_pui = -1.0;
  double survival_observer_r = 100.0;  // AU; beta is accumulated from here outward
  bool fitted_esw = true;        // fitted E_sw (false: the constant 0.04)
  // Option: shock-frame shift for the transmitted, reflected and SW
  // Maxwellians: "none" (default), "scalar" (|v_ENA + u| - v_shock, as for the
  // accelerated PUIs) or "vector" (u replaced by u - v_shock u/|u|).
  std::string shock_shift = "none";
  // How the termination shock is found along each streamline:
  // "regions" (default: entry into region 3/6 from region 4 in boundaries.txt),
  // "mach" (the first point with M^2 < 0.64 after a supersonic point, at
  // r < 130 AU) or "moscow_column" (the region column of the Moscow model).
  std::string ts_detection;      // empty = "regions"
  // Reference TS distance r_TS,V2 used to normalize the partition fits:
  // "computed" (the first TS crossing along the Voyager 2 path, found with the
  // ts_detection method; the default when a V2 path is given) or "constant"
  // (ts_reference_r; model default 78.5 AU, 82 AU Moscow).
  std::string ts_reference;
  double ts_reference_r = 0.0;   // 0 = model default
  // The V2 path, simulation frame, degrees (theta from the north pole, as the
  // grid): a fixed direction [theta, phi], or a trajectory file with columns
  // r [AU], theta, phi (path relative to the configuration file).
  std::vector<double> v2_direction;
  std::string v2_trajectory_file;
  std::vector<double> v2_r, v2_theta, v2_phi;  // trajectory, read by load_config
  int gap_fill_max_passes = 200; // gap fill: passes stop when nothing changes, or here
  // v_shock (shock-frame speed of the accelerated PUIs) with ts_detection =
  // "regions". The point where the streamline enters region 3/6 can still be
  // upstream of the velocity jump (points go to the nearest grid node, and the
  // MHD shock is spread over a few cells), so its speed is not used.
  //  "first_below" (default): the Mach criterion, started from the region
  //   crossing: the flow speed at the first point from the
  //   crossing on with M^2 < 0.64; if the streamline never gets there, the
  //   "downstream" speed.
  //  "downstream": the flow speed at the first point at least
  //   vshock_offset_cells grid cells (x dr, straight-line distance) past the
  //   crossing (version 2.1).
  //  "crossing": the flow speed at the crossing point itself (version 2.0).
  //   Kept only to compare with version 2.0; it can give the upstream speed.
  std::string vshock = "first_below";
  double vshock_offset_cells = 2.0;
  // Heliosheath (region 3/6) cells that the streamline gap fill leaves empty: their connected
  // part of the heliosheath holds no streamline data, so the termination-shock partition is
  // unknown there. "error" (default): stop at the first such cell on a line of sight.
  // "skip": the cell gives no ENAs (as the supersonic pockets), integration continues behind it.
  std::string unfilled_cells = "error";
};

// Time-dependent maps ([time_dependent], milestone 8). The input snapshots are
// time cells: snapshot s (s = 0 .. N-1) holds for [start + s dt, start + (s+1) dt),
// N = (end - start)/dt, with no interpolation in time. ENAs travel radially at
// v = (2E/m_p)^1/2 to an observer at observer_r_au; their production is taken
// from the snapshot of the emission time, their losses from the snapshots of the
// times they pass each cell. Maps are written for the output windows
// [output_start + w cadence, output_start + w cadence + window) up to output_end:
// the flux arriving at the observer, averaged over the window. All times in
// years (dt, cadence and window can also be given in months).
//
// Coverage: ENAs made before the first snapshot are missing from the first
// windows. A line of sight is flagged incomplete (a region with missing flux) in
// a window by either test:
//  - extinction test (default limit 0.3, about the cooling length): ENAs arriving at the window start a from beyond
//    r_cov = observer_r_au + v (a - start) are missing; incomplete if the
//    streamline extinction e^-tau of those heliosheath cells still reaches
//    coverage_extinction;
//  - missing-flux test: the extinction does not reduce the heliosheath-created
//    PUIs, so the missing flux itself is estimated, from the first snapshot: each
//    cell's flux (all components, at each energy, with its survival) times the
//    fraction of the window during which its ENAs were made before the first
//    snapshot, summed, over the window-averaged flux of all cells; the energies
//    are combined with the band count-rate weights. Incomplete if this share
//    reaches coverage_missing_flux (default 0: test off; the share is always
//    written).
// Flagged lines of sight: by default their flux is set to 0 in every component
// map of that band and window (include_incomplete_lines = false);
// include_incomplete_lines = true keeps their (incomplete) flux. Either way the
// output marks them (coverage_complete = 0, and a legacy flag file), so plots
// can box them as regions to be evaluated further.
struct TimeDependentConfig {
  bool enabled = false;        // a [time_dependent] section switches it on, unless it has enabled = false
  double start = 0, end = 0, dt = 0;                    // input range and step (years)
  double output_start = 0, output_end = 0, cadence = 0, window = 0;  // output windows (years)
  double observer_r_au = 1.0;
  double coverage_extinction = 0.3;
  double coverage_missing_flux = 0.0;
  bool include_incomplete_lines = false;
  // Memory for the stored survival terms (GB); if all energies do not fit, the
  // energies are processed in groups and the snapshots are read once per group.
  double memory_gb = 8.0;
  int snapshots() const { return enabled ? int(std::lround((end - start) / dt)) : 0; }
};

struct RunConfig {
  std::string source_path;
  std::string output_dir = "output";
  std::string model_source = "standard";  // "standard" | "moscow"
  // "four_population_ts" (single ion, Kornbleuth+2023), "extended_si" or
  // "multi_ion_inca" (Kornbleuth+2026)
  std::string populations = "four_population_ts";
  Grid grid;
  double r_lower = 2.0, r_upper = 800.0;
  int radial_range = 2;
  InputSpec inputs;
  PhysicsConfig physics;
  SensorTable sensors = SensorTable::single_ion_default();
  MultiIonConfig mi;
  Region6Config region6;  // extended single-ion only
  TimeDependentConfig time_dependent;
  // Notes from reading the configuration (e.g. settings skipped because they belong
  // to another model); ena_run prints them.
  std::vector<std::string> notes;
  std::vector<int> bands = {13};
  int eres = 10;

  bool moscow() const { return model_source == "moscow"; }
  bool multi_ion() const { return populations == "multi_ion_inca"; }
  // Extended single-ion model: the
  // single-ion partition with four neutral populations, the ion-neutral relative
  // speed in the cross section, heliosheath-created PUIs and the neutral
  // components.
  bool extended_si() const { return populations == "extended_si"; }
};

// `overrides` are "table.key=value" strings (e.g. "physics.fitted_esw=false")
// applied on top of the file.
RunConfig load_config(const std::string& path, const std::vector<std::string>& overrides = {});
// One line with the code version and the physics options, for logs and output
// files.
std::string settings_summary(const RunConfig& cfg);

RunConfig config_from_toml(const TomlTable& t, const std::string& source = "<toml>");

}  // namespace ena
