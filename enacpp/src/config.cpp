#include "ena/config.hpp"

#include "ena/ascii.hpp"

#include <cctype>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace ena {

// ===========================================================================
// Minimal TOML parser
// ===========================================================================
namespace {

struct Cursor {
  const std::string& s;
  std::size_t p = 0;
  int line = 1;
  std::string src;
  [[noreturn]] void fail(const std::string& msg) const {
    throw std::runtime_error(src + ":" + std::to_string(line) + ": " + msg);
  }
  bool eof() const { return p >= s.size(); }
  char peek() const { return eof() ? '\0' : s[p]; }
  void skip_ws() {  // spaces/tabs only
    while (!eof() && (s[p] == ' ' || s[p] == '\t' || s[p] == '\r')) ++p;
  }
  void skip_comment() {
    if (peek() == '#') while (!eof() && s[p] != '\n') ++p;
  }
  void skip_ws_nl_comments() {  // inside arrays
    for (;;) {
      skip_ws();
      if (peek() == '#') { skip_comment(); continue; }
      if (peek() == '\n') { ++p; ++line; continue; }
      break;
    }
  }
};

std::string parse_key(Cursor& c) {
  c.skip_ws();
  std::string key;
  if (c.peek() == '"') {
    ++c.p;
    while (!c.eof() && c.peek() != '"') key += c.s[c.p++];
    if (c.eof()) c.fail("unterminated quoted key");
    ++c.p;
  } else {
    while (!c.eof() && (std::isalnum((unsigned char)c.peek()) || c.peek() == '_' || c.peek() == '-' || c.peek() == '.'))
      key += c.s[c.p++];
  }
  if (key.empty()) c.fail("expected a key");
  return key;
}

TomlValue parse_value(Cursor& c) {
  c.skip_ws();
  TomlValue v;
  char ch = c.peek();
  if (ch == '"' || ch == '\'') {
    char q = ch;
    ++c.p;
    v.kind = TomlValue::String;
    while (!c.eof() && c.peek() != q) {
      if (c.peek() == '\n') c.fail("newline in string");
      if (q == '"' && c.peek() == '\\' && c.p + 1 < c.s.size()) {
        char e = c.s[c.p + 1];
        v.str += (e == 'n') ? '\n' : (e == 't') ? '\t' : e;
        c.p += 2;
        continue;
      }
      v.str += c.s[c.p++];
    }
    if (c.eof()) c.fail("unterminated string");
    ++c.p;
    return v;
  }
  if (ch == '[') {
    ++c.p;
    v.kind = TomlValue::Array;
    for (;;) {
      c.skip_ws_nl_comments();
      if (c.peek() == ']') { ++c.p; break; }
      v.arr.push_back(parse_value(c));
      c.skip_ws_nl_comments();
      if (c.peek() == ',') { ++c.p; continue; }
      if (c.peek() == ']') { ++c.p; break; }
      c.fail("expected ',' or ']' in array");
    }
    return v;
  }
  std::string tok;
  while (!c.eof() && !std::isspace((unsigned char)c.peek()) && c.peek() != ',' && c.peek() != ']' && c.peek() != '#')
    tok += c.s[c.p++];
  if (tok == "true" || tok == "false") {
    v.kind = TomlValue::Bool;
    v.boolean = (tok == "true");
    return v;
  }
  std::string clean;
  for (char t : tok) if (t != '_') clean += t;
  char* end = nullptr;
  v.num = std::strtod(clean.c_str(), &end);
  if (clean.empty() || *end != '\0') c.fail("cannot parse value '" + tok + "'");
  v.kind = TomlValue::Number;
  return v;
}

std::vector<std::string> split_dotted(const std::string& name) {
  std::vector<std::string> parts;
  std::string cur;
  for (char ch : name) {
    if (ch == '.') { parts.push_back(cur); cur.clear(); }
    else cur += ch;
  }
  parts.push_back(cur);
  return parts;
}

}  // namespace

TomlTable parse_toml(const std::string& text, const std::string& source_name) {
  TomlTable root;
  TomlTable* current = &root;
  Cursor c{text, 0, 1, source_name};
  while (!c.eof()) {
    c.skip_ws();
    if (c.peek() == '#') { c.skip_comment(); continue; }
    if (c.peek() == '\n') { ++c.p; ++c.line; continue; }
    if (c.eof()) break;
    if (c.peek() == '[') {
      bool array_of_tables = (c.p + 1 < text.size() && text[c.p + 1] == '[');
      c.p += array_of_tables ? 2 : 1;
      std::string name = parse_key(c);
      c.skip_ws();
      for (int n = 0; n < (array_of_tables ? 2 : 1); ++n) {
        if (c.peek() != ']') c.fail("expected ']'");
        ++c.p;
      }
      auto parts = split_dotted(name);
      TomlTable* t = &root;
      for (std::size_t i = 0; i + 1 < parts.size(); ++i) {
        auto it = t->table_arrays.find(parts[i]);
        t = (it != t->table_arrays.end() && !it->second.empty()) ? &it->second.back() : &t->tables[parts[i]];
      }
      if (array_of_tables) {
        auto& vec = t->table_arrays[parts.back()];
        vec.emplace_back();
        current = &vec.back();
      } else {
        current = &t->tables[parts.back()];
      }
    } else {
      std::string key = parse_key(c);
      c.skip_ws();
      if (c.peek() != '=') c.fail("expected '=' after key '" + key + "'");
      ++c.p;
      TomlValue v = parse_value(c);
      auto parts = split_dotted(key);
      TomlTable* t = current;
      for (std::size_t i = 0; i + 1 < parts.size(); ++i) t = &t->tables[parts[i]];
      t->values[parts.back()] = v;
    }
    c.skip_ws();
    c.skip_comment();
    if (!c.eof() && c.peek() != '\n') c.fail("unexpected text after value");
  }
  return root;
}

TomlTable parse_toml_file(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("cannot open configuration file " + path);
  std::stringstream ss;
  ss << in.rdbuf();
  return parse_toml(ss.str(), path);
}

const TomlTable* TomlTable::table(const std::string& name) const {
  auto it = tables.find(name);
  return it == tables.end() ? nullptr : &it->second;
}
double TomlTable::get_number(const std::string& key, double fallback) const {
  auto it = values.find(key);
  if (it == values.end()) return fallback;
  if (it->second.kind != TomlValue::Number) throw std::runtime_error("'" + key + "' must be a number");
  return it->second.num;
}
double TomlTable::require_number(const std::string& key) const {
  if (!has(key)) throw std::runtime_error("missing required setting '" + key + "'");
  return get_number(key, 0.0);
}
std::string TomlTable::get_string(const std::string& key, const std::string& fallback) const {
  auto it = values.find(key);
  if (it == values.end()) return fallback;
  if (it->second.kind != TomlValue::String) throw std::runtime_error("'" + key + "' must be a string");
  return it->second.str;
}
bool TomlTable::get_bool(const std::string& key, bool fallback) const {
  auto it = values.find(key);
  if (it == values.end()) return fallback;
  if (it->second.kind != TomlValue::Bool) throw std::runtime_error("'" + key + "' must be true or false");
  return it->second.boolean;
}
std::vector<std::string> TomlTable::get_strings(const std::string& key) const {
  std::vector<std::string> out;
  auto it = values.find(key);
  if (it == values.end()) return out;
  if (it->second.kind != TomlValue::Array) throw std::runtime_error("'" + key + "' must be an array");
  for (auto& v : it->second.arr) {
    if (v.kind != TomlValue::String) throw std::runtime_error("'" + key + "' must contain strings");
    out.push_back(v.str);
  }
  return out;
}

// ===========================================================================
// Columns, grid
// ===========================================================================

int ColumnMap::column(const std::string& role) const {
  for (std::size_t i = 0; i < columns.size(); ++i)
    if (columns[i] == role) return int(i);
  return -1;
}

double ColumnMap::velocity_factor(const std::string& role) const {
  int c = column(role);
  std::string u = (c >= 0 && std::size_t(c) < units.size()) ? units[c] : "";
  if (u.empty() || u == "km/s") return 1e5;
  if (u == "m/s") return 1e2;
  if (u == "cm/s") return 1.0;
  throw std::runtime_error("unknown velocity unit '" + u + "' for column '" + role + "'");
}

int Grid::k_of(double radius) const { return int(std::lround((radius - r_inner) / dr)); }

void Grid::validate() const {
  if (dtheta <= 0 || dphi <= 0 || dr <= 0) throw std::runtime_error("grid resolutions must be positive");
  if (std::fabs(180.0 / dtheta - std::round(180.0 / dtheta)) > 1e-9)
    throw std::runtime_error("180 must be a multiple of grid.dtheta");
  if (std::fabs(360.0 / dphi - std::round(360.0 / dphi)) > 1e-9)
    throw std::runtime_error("360 must be a multiple of grid.dphi");
  if (nr <= 0) throw std::runtime_error("grid.nr must be positive");
}

std::string InputSpec::path(const std::string& file) const {
  if (file.empty() || file[0] == '/' || dir.empty() || dir == ".") return file;
  return dir + "/" + file;
}

// ===========================================================================
// RunConfig
// ===========================================================================
namespace {

ColumnMap read_columns(const TomlTable* t, const std::vector<std::string>& defaults) {
  ColumnMap m;
  m.columns = defaults;
  if (t) {
    auto cols = t->get_strings("columns");
    if (!cols.empty()) m.columns = cols;
    m.units = t->get_strings("units");
    if (!m.units.empty() && m.units.size() != m.columns.size())
      throw std::runtime_error("'units' must have one entry per column");
  }
  return m;
}

}  // namespace

SensorTable SensorTable::single_ion_default() {
  SensorTable s;
  s.emin = {0.009, 0.016, 0.031, 0.061, 0.118, 0.271, 0.545, 1.142, 0.38, 0.52, 0.84, 1.36, 1.99, 3.13, 5.2, 13.5, 24., 35.};
  s.emax = {0.020, 0.038, 0.073, 0.143, 0.276, 0.631, 1.27, 2.66, 0.59, 0.95, 1.55, 2.5, 3.75, 6.0, 13.5, 24., 35., 55.};
  s.ecen = {0.014, 0.027, 0.052, 0.102, 0.197, 0.451, 0.908, 1.903, 0.45, 0.71, 1.11, 1.74, 2.73, 4.29,
            8.37854, 18.0, 28.9828, 43.8748};
  s.g = {2.6E-5, 5.3E-5, 8.1E-5, 9.1E-5, 9.0E-5, 1.0E-4, 1.9E-4, 2.7E-4, 1.3E-4, 4.1E-4, 7.5E-4, 1.3E-3, 2.4E-3, 4.5E-3,
         5.0E-3, 5.5E-3, 6.0E-3, 6.5E-3};
  s.weights[9] = {0.0013277909, 0.0016871014, 0.0019651474, 0.0022037048, 0.0022894956, 0.0022060268, 0.0020821493,
                  0.001934731, 0.0016461171, 0.0012819526};
  s.weights[10] = {0.0008472022, 0.0011271987, 0.0013325584, 0.0014113985, 0.0013964271, 0.0013630371, 0.0012635717,
                   0.0011411818, 0.0009836919, 0.000787456};
  s.weights[11] = {0.0005482951, 0.0007248594, 0.0008530854, 0.0008994173, 0.0008982827, 0.000832246, 0.0007757703,
                   0.000709315, 0.0006125453, 0.0005138144};
  s.weights[12] = {0.0003588889, 0.0004682034, 0.0005517594, 0.0005910363, 0.0005914153, 0.0005547477, 0.000513061,
                   0.0004612037, 0.0004018451, 0.0003388693};
  s.weights[13] = {0.0002324185, 0.0002991092, 0.000347635, 0.0003659699, 0.0003612539, 0.000345767, 0.0003205023,
                   0.0002890332, 0.0002574077, 0.000213177};
  return s;
}

SensorTable SensorTable::multi_ion_default() {
  SensorTable s = single_ion_default();
  // IBEX-Lo: widened bands 1-3 and CG-corrected central energies (tail, Galli+2022)
  s.emax[1] = 0.07; s.emax[2] = 0.1; s.emax[3] = 0.2;
  const double lo_c[8] = {0, 0.048, 0.082, 0.147, 0.260, 0.512, 0.975, 1.968};
  for (int i = 0; i < 8; ++i) s.ecen[i] = lo_c[i];
  s.cg_lo = true;
  // Ultra: INCA bands plus V2 energies (bands 14-20; the table's 8th slot is empty)
  const double umin[7] = {5.2, 13.5, 24., 35., 28., 60., 70.}, umax[7] = {13.5, 24., 35., 55., 43., 80., 90.},
               ucen[7] = {8.37854, 18.0, 28.9828, 43.8748, 33.4664, 70., 80.},
               ug[7] = {5.0E-3, 5.5E-3, 6.0E-3, 6.5E-3, 7.0E-3, 7.5E-3, 8.0E-3};
  s.emin.resize(14); s.emax.resize(14); s.ecen.resize(14); s.g.resize(14);
  for (int i = 0; i < 7; ++i) { s.emin.push_back(umin[i]); s.emax.push_back(umax[i]); s.ecen.push_back(ucen[i]); s.g.push_back(ug[i]); }
  return s;
}

RunConfig config_from_toml(const TomlTable& t, const std::string& source) {
  RunConfig c;
  c.source_path = source;

  if (auto* run = t.table("run")) {
    c.output_dir = run->get_string("output_dir", c.output_dir);
    const std::string mode = run->get_string("mode", "fixed");
    if (mode != "fixed")
      throw std::runtime_error("run.mode = \"" + mode + "\": this version has no faithful mode; use the version 1.0 "
                               "archive to reproduce the results of the original IDL code");
  }
  if (t.table("quirks"))
    throw std::runtime_error("[quirks]: this version has the fixes built in (see CHANGES_FROM_IDL.md); "
                             "per-quirk switches exist only in the version 1.0 archive");

  if (auto* m = t.table("model")) {
    c.model_source = m->get_string("source", c.model_source);
    c.populations = m->get_string("populations", c.populations);
  }
  if (c.model_source != "standard" && c.model_source != "moscow")
    throw std::runtime_error("model.source must be \"standard\" or \"moscow\"");
  if (c.populations != "four_population_ts" && c.populations != "multi_ion_inca" && c.populations != "extended_si")
    throw std::runtime_error("model.populations must be \"four_population_ts\", \"extended_si\" or \"multi_ion_inca\"");
  if (c.multi_ion() || c.extended_si()) {
    if (c.moscow()) throw std::runtime_error("the multi-ion and extended single-ion models need model.source = \"standard\"");
    // The extended single-ion model uses the multi-ion sensor table
    // (IBEX-Lo at the CG-corrected central energies, Ultra/INCA).
    c.sensors = SensorTable::multi_ion_default();
    c.radial_range = 0;  // TS to HP (the range of both codes)
  }

  const TomlTable* g = t.table("grid");
  if (!g) throw std::runtime_error("missing [grid] section");
  c.grid.dtheta = g->get_number("dtheta", c.grid.dtheta);
  c.grid.dphi = g->get_number("dphi", c.grid.dphi);
  c.grid.dr = g->get_number("dr", c.grid.dr);
  c.grid.r_inner = g->get_number("r_inner", c.grid.r_inner);
  c.grid.nr = int(g->require_number("nr"));
  c.grid.validate();
  c.grid.ntheta = int(std::lround(180.0 / c.grid.dtheta));
  c.grid.nphi = int(std::lround(360.0 / c.grid.dphi));

  if (auto* in = t.table("integration")) {
    c.r_lower = in->get_number("r_lower", c.r_lower);
    c.r_upper = in->get_number("r_upper", c.r_upper);
    c.radial_range = int(in->get_number("radial_range", c.radial_range));
  }
  if ((c.multi_ion() || c.extended_si()) && c.radial_range != 0)
    throw std::runtime_error("the multi-ion and extended single-ion models are implemented for radial_range = 0 (TS to HP) only");
  if (t.table("legacy"))
    std::fprintf(stderr, "warning: [legacy] is ignored (ri_opher is no longer used)\n");

  const TomlTable* in = t.table("inputs");
  c.inputs.dir = in ? in->get_string("dir", ".") : ".";
  c.inputs.format = in ? in->get_string("format", "ascii") : "ascii";
  c.inputs.allow_incomplete = in ? in->get_bool("allow_incomplete", false) : false;
  auto sub = [&](const char* n) { return in ? in->table(n) : nullptr; };
  if (c.multi_ion()) {
    // multi-ion files: 25 plasma columns (SW, pu3, pu2), 41 neutral columns (total H and populations 1-4)
    c.inputs.plasma = read_columns(sub("plasma"), {"n", "T", "n_pu3", "T_pu3", "n_pu2", "T_pu2", "vr", "vtheta", "vphi", "vt",
                                                   "vr_pu3", "vtheta_pu3", "vphi_pu3", "vt_pu3", "vr_pu2", "vtheta_pu2",
                                                   "vphi_pu2", "vt_pu2", "i", "j", "k", "theta", "phi", "r", "pbeta"});
    c.inputs.neutral = read_columns(sub("neutral"), {"n", "vr", "vt", "T", "i", "j", "k", "theta", "phi", "r",
                                                     "n1", "vr1", "vt1", "T1", "n2", "vr2", "vt2", "T2", "n3", "vr3", "vt3", "T3",
                                                     "n4", "vr4", "vt4", "T4", "vx4", "vy4", "vz4", "vx3", "vy3", "vz3",
                                                     "vx2", "vy2", "vz2", "vx1", "vy1", "vz1", "vx", "vy", "vz"});
  } else {
    c.inputs.plasma = read_columns(sub("plasma"), {"n", "T", "vr", "vtheta", "vphi", "vt", "i", "j", "k", "theta", "phi", "r"});
    if (c.extended_si())  // nhpops = 1: the 41-column neutral file of the multi-ion model
      c.inputs.neutral = read_columns(sub("neutral"), {"n", "vr", "vt", "T", "i", "j", "k", "theta", "phi", "r",
                                                       "n1", "vr1", "vt1", "T1", "n2", "vr2", "vt2", "T2", "n3", "vr3", "vt3", "T3",
                                                       "n4", "vr4", "vt4", "T4", "vx4", "vy4", "vz4", "vx3", "vy3", "vz3",
                                                       "vx2", "vy2", "vz2", "vx1", "vy1", "vz1", "vx", "vy", "vz"});
    else
      c.inputs.neutral = read_columns(sub("neutral"), {"n", "T", "i", "j", "k", "theta", "phi", "r"});
  }
  c.inputs.boundaries = read_columns(sub("boundaries"), {"region", "i", "j", "k", "theta", "phi", "r", "Bx", "By", "Bz"});
  std::vector<std::string> sl = {"x", "y", "z", "ux", "uy", "uz", "nH", "uxH", "uyH", "uzH",
                                 "np", "tp", "vr", "vtheta", "vphi", "vt", "theta", "phi", "r"};
  if (c.moscow()) sl.push_back("region");
  if (c.multi_ion())  // multi-ion streamlines: SW and pu3 quantities
    sl = {"x", "y", "z", "ux", "uy", "uz", "nH", "uxH", "uyH", "uzH", "np", "tp", "vr", "vtheta", "vphi", "vt",
          "n_pu3", "T_pu3", "vr_pu3", "vtheta_pu3", "vphi_pu3", "vt_pu3", "theta", "phi", "r"};
  c.inputs.streamlines = read_columns(sub("streamlines"), sl);

  if (auto* ph = t.table("physics")) {
    PhysicsConfig& p = c.physics;
    p.survival = ph->get_bool("survival", p.survival);
    p.extinction = ph->get_bool("extinction", p.extinction);
    p.cold_electrons = ph->get_bool("cold_electrons", p.cold_electrons);
    p.kappa = ph->get_number("kappa", p.kappa);
    p.kappa_pui = ph->get_number("kappa_pui", p.kappa_pui);
    p.kappa_ref = ph->get_number("kappa_ref", p.kappa_ref);
    p.kappa_ism = ph->get_number("kappa_ism", p.kappa_ism);
    p.kappa_oh = ph->get_number("kappa_oh", p.kappa_oh);
    p.kappa_heliosheath_pui = ph->get_number("kappa_heliosheath_pui", p.kappa_heliosheath_pui);
    p.survival_observer_r = ph->get_number("survival_observer_r", p.survival_observer_r);
    p.unfilled_cells = ph->get_string("unfilled_cells", p.unfilled_cells);
    if (p.unfilled_cells != "error" && p.unfilled_cells != "skip")
      throw std::runtime_error("[physics] unfilled_cells must be \"error\" or \"skip\", not \"" + p.unfilled_cells + "\"");
    p.fitted_esw = ph->get_bool("fitted_esw", p.fitted_esw);
    p.shock_shift = ph->get_string("shock_shift", p.shock_shift);
    p.ts_detection = ph->get_string("ts_detection", p.ts_detection);
    p.ts_reference = ph->get_string("ts_reference", "");
    p.ts_reference_r = ph->get_number("ts_reference_r", p.ts_reference_r);
    p.gap_fill_max_passes = int(ph->get_number("gap_fill_max_passes", p.gap_fill_max_passes));
    p.vshock_offset_cells = ph->get_number("vshock_offset_cells", p.vshock_offset_cells);
    p.vshock = ph->get_string("vshock", p.vshock);
    if (p.vshock != "first_below" && p.vshock != "downstream" && p.vshock != "crossing")
      throw std::runtime_error("physics.vshock must be \"first_below\", \"downstream\" or \"crossing\"");
    if (p.vshock_offset_cells < 0) throw std::runtime_error("physics.vshock_offset_cells must be >= 0");
    p.v2_trajectory_file = ph->get_string("v2_trajectory_file", "");
    auto it = ph->values.find("v2_direction");
    if (it != ph->values.end()) {
      if (it->second.kind != TomlValue::Array || it->second.arr.size() != 2)
        throw std::runtime_error("physics.v2_direction must be [theta, phi] in degrees");
      p.v2_direction.clear();
      for (auto& v : it->second.arr) p.v2_direction.push_back(v.num);
    }
  }
  {
    PhysicsConfig& p = c.physics;
    if (p.ts_detection.empty()) p.ts_detection = "regions";
    // r_TS,V2 is computed from the model when the V2 direction is known.
    const bool have_v2 = p.v2_direction.size() == 2 || !p.v2_trajectory_file.empty();
    if (p.v2_direction.size() == 2 && !p.v2_trajectory_file.empty())
      throw std::runtime_error("physics: give either v2_direction or v2_trajectory_file, not both");
    if (p.ts_reference.empty()) p.ts_reference = have_v2 ? "computed" : "constant";
    if (p.ts_reference_r == 0.0) p.ts_reference_r = c.moscow() ? 82.0 : 78.5;
    if (p.shock_shift != "none" && p.shock_shift != "scalar" && p.shock_shift != "vector")
      throw std::runtime_error("physics.shock_shift must be \"none\", \"scalar\" or \"vector\"");
    if (p.ts_detection != "mach" && p.ts_detection != "moscow_column" && p.ts_detection != "regions")
      throw std::runtime_error("physics.ts_detection must be \"mach\", \"moscow_column\" or \"regions\"");
    if (p.ts_detection == "moscow_column" && !c.moscow())
      throw std::runtime_error("physics.ts_detection = \"moscow_column\" needs model.source = \"moscow\"");
    if (p.ts_reference != "constant" && p.ts_reference != "computed")
      throw std::runtime_error("physics.ts_reference must be \"constant\" or \"computed\"");
    if (p.ts_reference == "computed" && !have_v2)
      throw std::runtime_error("physics.ts_reference = \"computed\" needs physics.v2_trajectory_file or physics.v2_direction");
    if (p.gap_fill_max_passes < 1) throw std::runtime_error("physics.gap_fill_max_passes must be >= 1");
  }
  {
    PhysicsConfig& p = c.physics;
    if (p.kappa_heliosheath_pui >= 0 && !c.extended_si()) {
      // a master configuration can keep the setting of another model: skipped, with a note
      c.notes.push_back("physics.kappa_heliosheath_pui is skipped (extended single-ion model only)");
      p.kappa_heliosheath_pui = -1.0;
    }
    if (p.kappa_heliosheath_pui < 0) p.kappa_heliosheath_pui = p.kappa_oh;
    Region6Config& r6 = c.region6;
    const TomlTable* t6 = t.table("region6");
    if (t6 && !c.extended_si()) {
      c.notes.push_back("[region6] is skipped (extended single-ion model only)");
      t6 = nullptr;
    }
    if (t6) {
      for (auto& [key, v] : t6->values) {
        (void)v;
        static const char* known[] = {"kappa_transmitted", "kappa_reflected", "kappa_heliosheath_pui", "kappa_sw",
                                      "accelerated_fraction", "accelerated_energy_share", "accelerated_kappa",
                                      "accelerated_from"};
        bool ok = false;
        for (const char* k : known) ok |= key == k;
        if (!ok) throw std::runtime_error("[region6]: unknown setting '" + key + "'");
      }
      r6.kappa_transmitted = t6->get_number("kappa_transmitted", r6.kappa_transmitted);
      r6.kappa_reflected = t6->get_number("kappa_reflected", r6.kappa_reflected);
      r6.kappa_heliosheath_pui = t6->get_number("kappa_heliosheath_pui", r6.kappa_heliosheath_pui);
      r6.kappa_sw = t6->get_number("kappa_sw", r6.kappa_sw);
      r6.accelerated_fraction = t6->get_number("accelerated_fraction", r6.accelerated_fraction);
      r6.accelerated_energy_share = t6->get_number("accelerated_energy_share", r6.accelerated_energy_share);
      r6.accelerated_kappa = t6->get_number("accelerated_kappa", r6.accelerated_kappa);
      r6.accelerated_from = t6->get_string("accelerated_from", r6.accelerated_from);
    }
    if (r6.kappa_transmitted < 0) r6.kappa_transmitted = p.kappa_pui;
    if (r6.kappa_reflected < 0) r6.kappa_reflected = p.kappa_ref;
    if (r6.kappa_heliosheath_pui < 0) r6.kappa_heliosheath_pui = p.kappa_oh;
    if (r6.kappa_sw < 0) r6.kappa_sw = p.kappa_oh;
    const double f = r6.accelerated_fraction, e = r6.accelerated_energy_share;
    if (!(f >= 0 && f < 1)) throw std::runtime_error("region6.accelerated_fraction must be in [0, 1)");
    if (!(e >= 0 && e < 1)) throw std::runtime_error("region6.accelerated_energy_share must be in [0, 1)");
    if (f > 0 && !(e > 0)) throw std::runtime_error("region6.accelerated_fraction > 0 needs accelerated_energy_share > 0");
    if (f == 0 && e > 0) throw std::runtime_error("region6.accelerated_energy_share > 0 needs accelerated_fraction > 0");
    if (r6.accelerated_from != "each_population" && r6.accelerated_from != "pooled")
      throw std::runtime_error("region6.accelerated_from must be \"each_population\" or \"pooled\"");
  }
  const TomlTable* mi = t.table("multi_ion");
  if (mi && !c.multi_ion()) {
    c.notes.push_back("[multi_ion] is skipped (multi-ion model only)");
    mi = nullptr;
  }
  if (mi) {
    MultiIonConfig& m = c.mi;
    m.belt = mi->get_bool("belt", m.belt);
    for (auto [key, val] : std::vector<std::pair<const char*, double*>>{
             {"belt_r_max", &m.belt_r_max}, {"belt_phi_min", &m.belt_phi_min}, {"belt_phi_max", &m.belt_phi_max},
             {"belt_ratio_min", &m.belt_ratio_min}, {"belt_mva2_max", &m.belt_mva2_max},
             {"belt_accel_frac", &m.belt_accel_frac}, {"belt_pu2_frac", &m.belt_pu2_frac},
             {"belt_t_pu2_factor", &m.belt_t_pu2_factor}, {"belt_t_sw_factor", &m.belt_t_sw_factor},
             {"belt_kappa_pu2", &m.belt_kappa_pu2}, {"belt_kappa_sw", &m.belt_kappa_sw},
             {"ts_compress", &m.ts_compress}, {"ts_cap", &m.ts_cap}, {"vshock_vmin_kms", &m.vshock_vmin_kms}})
      *val = mi->get_number(key, *val);
    m.injection = mi->get_string("injection", m.injection);
    m.injection_k = int(mi->get_number("injection_k", m.injection_k));
    m.ts_scaling = mi->get_string("ts_scaling", m.ts_scaling);
    if (m.injection != "streamline" && m.injection != "line_of_sight")
      throw std::runtime_error("multi_ion.injection must be \"streamline\" or \"line_of_sight\"");
    if (m.ts_scaling != "compressed" && m.ts_scaling != "ratio")
      throw std::runtime_error("multi_ion.ts_scaling must be \"compressed\" or \"ratio\"");
    if (!(m.ts_compress > 0) || !(m.ts_cap >= 1)) throw std::runtime_error("multi_ion: ts_compress must be > 0, ts_cap >= 1");
  }
  if (c.multi_ion() && c.physics.ts_detection != "regions")
    throw std::runtime_error("the multi-ion model needs physics.ts_detection = \"regions\"");
  auto numbers = [](const TomlTable& tt, const std::string& key) {
    std::vector<double> out;
    auto it = tt.values.find(key);
    if (it == tt.values.end()) return out;
    if (it->second.kind != TomlValue::Array) throw std::runtime_error("'" + key + "' must be an array of numbers");
    for (auto& v : it->second.arr) {
      if (v.kind != TomlValue::Number) throw std::runtime_error("'" + key + "' must contain numbers");
      out.push_back(v.num);
    }
    return out;
  };
  if (auto* b = t.table("bands")) {
    auto run = numbers(*b, "run");
    if (!run.empty()) { c.bands.clear(); for (double x : run) c.bands.push_back(int(x)); }
    c.eres = int(b->get_number("eres", c.eres));
  }
  if (auto* s = t.table("sensors")) {
    for (auto [key, vec] : {std::pair<const char*, std::vector<double>*>{"emin", &c.sensors.emin}, {"emax", &c.sensors.emax},
                           {"ecen", &c.sensors.ecen}, {"g", &c.sensors.g}}) {
      auto v = numbers(*s, key);
      if (!v.empty()) *vec = v;
    }
    c.sensors.central_first = int(s->get_number("central_first", c.sensors.central_first));
    c.sensors.central_last = int(s->get_number("central_last", c.sensors.central_last));
    if (auto* w = s->table("weights"))
      for (auto& [band, val] : w->values) c.sensors.weights[std::stoi(band)] = numbers(*w, band);
  }
  const std::size_t nb = c.sensors.emin.size();
  if (c.sensors.emax.size() != nb || c.sensors.ecen.size() != nb || c.sensors.g.size() != nb)
    throw std::runtime_error("sensors: emin, emax, ecen and g must have the same length");
  for (int b : c.bands) {
    if (b < 0 || std::size_t(b) >= nb) throw std::runtime_error("band " + std::to_string(b) + " is not in the sensor table");
    for (int u : c.sensors.unusable)
      if (b == u) throw std::runtime_error("band " + std::to_string(b) + " is marked unusable (faulty IBEX-Hi channel)");
    if (!(c.sensors.emax[b] > c.sensors.emin[b]) && !c.sensors.central_only(b))
      throw std::runtime_error("band " + std::to_string(b) + " has no energy range");
  }
  if (c.eres < 1) throw std::runtime_error("bands.eres must be >= 1");

  // ---- time-dependent maps ----
  TimeDependentConfig& td = c.time_dependent;
  const TomlTable* tdt = t.table("time_dependent");
  // enabled = false keeps the section in the file (e.g. a master configuration) without using it
  if (tdt && !tdt->get_bool("enabled", true)) tdt = nullptr;
  if (tdt) {
    static const char* known[] = {"start", "end", "dt", "dt_months", "output_start", "output_end", "cadence",
                                  "cadence_months", "window", "window_months", "observer_r_au", "coverage_extinction",
                                  "coverage_missing_flux", "memory_gb", "include_incomplete_lines", "enabled"};
    for (auto& [key, v] : tdt->values) {
      bool ok = false;
      for (const char* k : known) ok = ok || key == k;
      if (!ok) throw std::runtime_error("[time_dependent]: unknown setting '" + key + "'");
    }
    // a duration in years, or in months with the _months key (one of the two)
    auto years = [&](const std::string& key, double fallback) {
      const bool y = tdt->has(key), m = tdt->has(key + "_months");
      if (y && m) throw std::runtime_error("time_dependent: give " + key + " or " + key + "_months, not both");
      if (m) return tdt->require_number(key + "_months") / 12.0;
      return y ? tdt->require_number(key) : fallback;
    };
    td.enabled = true;
    td.start = tdt->require_number("start");
    td.end = tdt->require_number("end");
    td.dt = years("dt", 0.0);
    if (!(td.end > td.start)) throw std::runtime_error("time_dependent: end must be after start");
    if (!(td.dt > 0)) throw std::runtime_error("time_dependent: dt (or dt_months) must be positive");
    // N = (end - start)/dt must be a whole number; dt is then taken as exactly
    // (end - start)/N, so e.g. dt = 0.1667 (2 months) is accepted
    const double nf = (td.end - td.start) / td.dt;
    const long n = std::lround(nf);
    if (n < 1 || std::fabs(nf - double(n)) > 1e-3)
      throw std::runtime_error("time_dependent: (end - start)/dt must be a whole number of snapshots");
    td.dt = (td.end - td.start) / double(n);
    td.output_start = tdt->get_number("output_start", td.start);
    td.output_end = tdt->get_number("output_end", td.end);
    td.cadence = years("cadence", 1.0);
    td.window = years("window", td.cadence);
    td.observer_r_au = tdt->get_number("observer_r_au", td.observer_r_au);
    td.coverage_extinction = tdt->get_number("coverage_extinction", td.coverage_extinction);
    td.coverage_missing_flux = tdt->get_number("coverage_missing_flux", td.coverage_missing_flux);
    td.memory_gb = tdt->get_number("memory_gb", td.memory_gb);
    td.include_incomplete_lines = tdt->get_bool("include_incomplete_lines", td.include_incomplete_lines);
    const double eps = 1e-9 * (td.end - td.start);
    if (td.output_start < td.start - eps || td.output_end > td.end + eps)
      throw std::runtime_error("time_dependent: the output range must lie inside [start, end] (ENAs arriving after "
                               "end were partly made after the last snapshot)");
    if (!(td.cadence > 0) || !(td.window > 0)) throw std::runtime_error("time_dependent: cadence and window must be positive");
    if (td.output_start + td.window > td.output_end + eps)
      throw std::runtime_error("time_dependent: no output window fits between output_start and output_end");
    if (!(td.observer_r_au >= 0)) throw std::runtime_error("time_dependent.observer_r_au must be >= 0");
    if (!(td.coverage_extinction > 0 && td.coverage_extinction < 1))
      throw std::runtime_error("time_dependent.coverage_extinction must be in (0, 1)");
    if (!(td.coverage_missing_flux >= 0 && td.coverage_missing_flux < 1))
      throw std::runtime_error("time_dependent.coverage_missing_flux must be in [0, 1) (0: not used)");
    if (!(td.memory_gb > 0)) throw std::runtime_error("time_dependent.memory_gb must be positive");
  }

  // ---- snapshots: [[snapshot]] entries, or (time-dependent) file name patterns ----
  auto check_files = [&](const SnapshotFiles& f, const std::string& where) {
    if (c.inputs.format != "ascii" && c.inputs.format != "hdf5")
      throw std::runtime_error("inputs.format must be \"ascii\" or \"hdf5\"");
    if (c.inputs.format == "ascii" && (f.plasma.empty() || f.neutral.empty() || f.boundaries.empty()))
      throw std::runtime_error(where + " needs plasma, neutral and boundaries files");
    if (c.inputs.format == "hdf5" && f.file.empty())
      throw std::runtime_error("with inputs.format = \"hdf5\" " + where + " needs 'file'");
  };
  auto it = t.table_arrays.find("snapshot");
  const bool listed = it != t.table_arrays.end() && !it->second.empty();
  const TomlTable* pat = tdt ? tdt->table("files") : nullptr;
  // Time-dependent runs take [time_dependent.files] when given (any [[snapshot]] entries, e.g. the
  // steady-state snapshot of a master configuration, are then not used); otherwise [[snapshot]].
  if (listed && pat) c.notes.push_back("[[snapshot]] entries are not used: the time-dependent run takes [time_dependent.files]");
  if (listed && !pat) {
    for (auto& s : it->second) {
      SnapshotFiles f;
      f.has_time = s.has("time");
      f.time = s.get_number("time", 0.0);
      f.plasma = s.get_string("plasma", "");
      f.neutral = s.get_string("neutral", "");
      f.boundaries = s.get_string("boundaries", "");
      f.streamlines = s.get_string("streamlines", "");
      f.file = s.get_string("file", "");
      check_files(f, "each [[snapshot]]");
      c.inputs.snapshots.push_back(f);
    }
  } else if (pat) {
    // [time_dependent.files]: one name pattern per input; {index} is replaced by
    // the snapshot number first_index + s, {index:03} pads it with zeros to 3 digits.
    for (auto& [key, v] : pat->values)
      if (key != "plasma" && key != "neutral" && key != "boundaries" && key != "streamlines" && key != "file" &&
          key != "first_index")
        throw std::runtime_error("[time_dependent.files]: unknown setting '" + key + "'");
    const long first = long(pat->get_number("first_index", 0));
    auto expand = [&](const std::string& p, long index) {
      std::string out;
      for (std::size_t i = 0; i < p.size(); ++i) {
        if (p.compare(i, 6, "{index") == 0) {
          const std::size_t close = p.find('}', i);
          if (close == std::string::npos) throw std::runtime_error("[time_dependent.files]: unclosed {index in '" + p + "'");
          const std::string spec = p.substr(i + 6, close - i - 6);
          int width = 0;
          if (!spec.empty()) {
            if (spec.size() < 2 || spec[0] != ':' || spec.find_first_not_of("0123456789", 1) != std::string::npos)
              throw std::runtime_error("[time_dependent.files]: use {index} or {index:0N} in '" + p + "'");
            width = std::stoi(spec.substr(1));
          }
          std::string num = std::to_string(index);
          if (int(num.size()) < width) num = std::string(width - num.size(), '0') + num;
          out += num;
          i = close;
        } else {
          out += p[i];
        }
      }
      return out;
    };
    for (int s = 0; s < td.snapshots(); ++s) {
      SnapshotFiles f;
      for (auto [key, dst] : {std::pair<const char*, std::string*>{"plasma", &f.plasma}, {"neutral", &f.neutral},
                              {"boundaries", &f.boundaries}, {"streamlines", &f.streamlines}, {"file", &f.file}})
        *dst = expand(pat->get_string(key, ""), first + s);
      check_files(f, "[time_dependent.files]");
      c.inputs.snapshots.push_back(f);
    }
  } else {
    throw std::runtime_error("no [[snapshot]] entries: at least one is required");
  }
  if (td.enabled) {
    // snapshot s is the time cell starting at start + s dt; 'time', if given, must be that start
    const int n = td.snapshots();
    if (int(c.inputs.snapshots.size()) != n)
      throw std::runtime_error("time_dependent: (end - start)/dt = " + std::to_string(n) + " snapshots, but " +
                               std::to_string(c.inputs.snapshots.size()) + " are given");
    for (int s = 0; s < n; ++s) {
      SnapshotFiles& f = c.inputs.snapshots[s];
      const double t0 = td.start + s * td.dt;
      if (f.has_time && std::fabs(f.time - t0) > 1e-3 * td.dt)
        throw std::runtime_error("time_dependent: snapshot " + std::to_string(s) + " has time " + std::to_string(f.time) +
                                 ", but its time cell starts at start + " + std::to_string(s) + " dt = " +
                                 std::to_string(t0));
      f.time = t0;
      f.has_time = true;
    }
  }
  for (std::size_t i = 1; i < c.inputs.snapshots.size(); ++i)
    if (c.inputs.snapshots[i].time <= c.inputs.snapshots[i - 1].time)
      throw std::runtime_error("[[snapshot]] times must be strictly increasing");
  return c;
}

std::string settings_summary(const RunConfig& c) {
  const PhysicsConfig& p = c.physics;
  char buf[320];
  std::snprintf(buf, sizeof buf,
                "ENA C++ %s; fitted_esw %s, shock_shift %s, ts_detection %s%s, ts_reference %s, survival %s, "
                "extinction %s, cold_electrons %s",
                kVersion, p.fitted_esw ? "on" : "off", p.shock_shift.c_str(), p.ts_detection.c_str(),
                p.ts_detection == "regions" ? (" (v_shock " + p.vshock + ")").c_str() : "",
                p.ts_reference.c_str(), p.survival ? "on" : "off", p.extinction ? "on" : "off",
                p.cold_electrons ? "on" : "off");
  if (c.multi_ion()) {
    const MultiIonConfig& m = c.mi;
    char mb[320];
    std::snprintf(mb, sizeof mb,
                  "; multi-ion: TS fits %s scaling, injection %s, reconnection region %s (r <= %g AU, phi < %g or > %g deg)",
                  m.ts_scaling.c_str(), m.injection.c_str(), m.belt ? "on" : "off", m.belt_r_max, m.belt_phi_min,
                  m.belt_phi_max);
    return std::string(buf) + mb;
  }
  if (c.extended_si()) {
    const Region6Config& r = c.region6;
    char eb[400];
    std::snprintf(eb, sizeof eb,
                  "; extended single-ion: E_l from the flow speed inside 12 AU, kappa of the heliosheath PUIs %g; region 6: "
                  "kappas %g/%g/%g/%g (transmitted/reflected/heliosheath PUIs/SW), accelerated fraction %g with energy share %g, "
                  "kappa %g, from %s",
                  p.kappa_heliosheath_pui, r.kappa_transmitted, r.kappa_reflected, r.kappa_heliosheath_pui, r.kappa_sw,
                  r.accelerated_fraction, r.accelerated_energy_share, r.accelerated_kappa, r.accelerated_from.c_str());
    return std::string(buf) + eb;
  }
  return buf;
}

namespace {
void apply_override(TomlTable& root, const std::string& ov) {
  auto eq = ov.find('=');
  if (eq == std::string::npos || eq == 0) throw std::runtime_error("override '" + ov + "' is not table.key=value");
  std::string key = ov.substr(0, eq), val = ov.substr(eq + 1);
  std::vector<std::string> parts;
  std::size_t start = 0;
  for (std::size_t i = 0; i <= key.size(); ++i)
    if (i == key.size() || key[i] == '.') { parts.push_back(key.substr(start, i - start)); start = i + 1; }
  TomlTable* t = &root;
  for (std::size_t i = 0; i + 1 < parts.size(); ++i) t = &t->tables[parts[i]];
  TomlValue v;
  try {
    v = parse_toml("v = " + val, "override " + ov).values.at("v");
  } catch (const std::exception&) {  // bare word: a string
    v.kind = TomlValue::String;
    v.str = val;
  }
  t->values[parts.back()] = v;
}
}  // namespace

RunConfig load_config(const std::string& path, const std::vector<std::string>& overrides) {
  TomlTable root = parse_toml_file(path);
  for (auto& ov : overrides) apply_override(root, ov);
  RunConfig c = config_from_toml(root, path);
  // Relative input directories are taken relative to the configuration file.
  if (!c.inputs.dir.empty() && c.inputs.dir[0] != '/') {
    auto slash = path.find_last_of('/');
    std::string base = slash == std::string::npos ? "" : path.substr(0, slash);
    if (!base.empty()) c.inputs.dir = c.inputs.dir == "." ? base : base + "/" + c.inputs.dir;
  }
  PhysicsConfig& p = c.physics;
  if (!p.v2_trajectory_file.empty()) {
    std::string f = p.v2_trajectory_file;
    auto slash = path.find_last_of('/');
    if (f[0] != '/' && slash != std::string::npos) f = path.substr(0, slash) + "/" + f;
    for_each_line(f, [&](std::string_view text, long line) {
      auto first = text.find_first_not_of(" \t\r");
      if (first == std::string_view::npos || text[first] == '#') return;  // blank or comment
      double v[4];
      if (parse_numbers(text, v, 3) != 3)
        throw std::runtime_error(f + ":" + std::to_string(line) + ": expected r, theta, phi");
      if (!p.v2_r.empty() && !(v[0] > p.v2_r.back()))
        throw std::runtime_error(f + ":" + std::to_string(line) + ": r must increase");
      if (v[1] < 0 || v[1] > 180) throw std::runtime_error(f + ":" + std::to_string(line) + ": theta outside [0, 180]");
      p.v2_r.push_back(v[0]);
      p.v2_theta.push_back(v[1]);
      // unwrap phi so it can be interpolated across 0/360
      double ph = v[2];
      if (!p.v2_phi.empty()) {
        while (ph - p.v2_phi.back() > 180) ph -= 360;
        while (ph - p.v2_phi.back() < -180) ph += 360;
      }
      p.v2_phi.push_back(ph);
    });
    if (p.v2_r.size() < 2) throw std::runtime_error(f + ": a V2 trajectory needs at least two rows");
  }
  return c;
}

}  // namespace ena
