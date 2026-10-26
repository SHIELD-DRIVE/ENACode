#include "ena/h5io.hpp"

#include <hdf5.h>

#include <cmath>
#include <ctime>
#include <stdexcept>
#include <sys/stat.h>

namespace ena {

namespace {

void check(herr_t s, const std::string& what) {
  if (s < 0) throw std::runtime_error("HDF5 error: " + what);
}
hid_t check_id(hid_t id, const std::string& what) {
  if (id < 0) throw std::runtime_error("HDF5 error: " + what);
  return id;
}

void make_dirs(const std::string& dir) {
  std::string cur;
  for (std::size_t i = 0; i <= dir.size(); ++i) {
    if (i == dir.size() || dir[i] == '/') {
      if (!cur.empty()) mkdir(cur.c_str(), 0755);
    }
    if (i < dir.size()) cur += dir[i];
  }
}

// Fixed-length string array attribute (padded).
void write_string_attr(hid_t loc, const std::string& name, const std::vector<std::string>& vals, bool scalar) {
  std::size_t w = 1;
  for (auto& v : vals) w = std::max(w, v.size());
  std::vector<char> buf(vals.size() * w, ' ');
  for (std::size_t i = 0; i < vals.size(); ++i) std::copy(vals[i].begin(), vals[i].end(), buf.begin() + i * w);
  hid_t type = H5Tcopy(H5T_C_S1);
  H5Tset_size(type, w);
  H5Tset_strpad(type, H5T_STR_SPACEPAD);
  hsize_t n = vals.size();
  hid_t space = scalar ? H5Screate(H5S_SCALAR) : H5Screate_simple(1, &n, nullptr);
  hid_t a = check_id(H5Acreate2(loc, name.c_str(), type, space, H5P_DEFAULT, H5P_DEFAULT), "create attribute " + name);
  check(H5Awrite(a, type, buf.data()), "write attribute " + name);
  H5Aclose(a);
  H5Sclose(space);
  H5Tclose(type);
}

void write_num_attr(hid_t loc, const std::string& name, double v) {
  hid_t space = H5Screate(H5S_SCALAR);
  hid_t a = check_id(H5Acreate2(loc, name.c_str(), H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT), "attr " + name);
  H5Awrite(a, H5T_NATIVE_DOUBLE, &v);
  H5Aclose(a);
  H5Sclose(space);
}

double read_num_attr(hid_t loc, const std::string& name) {
  hid_t a = check_id(H5Aopen(loc, name.c_str(), H5P_DEFAULT), "missing attribute " + name);
  double v = 0;
  H5Aread(a, H5T_NATIVE_DOUBLE, &v);
  H5Aclose(a);
  return v;
}

void write_dataset(hid_t loc, const std::string& name, const void* data, hid_t type,
                   const std::vector<hsize_t>& dims, const std::vector<std::string>& dim_names, int deflate = 4) {
  hid_t space = dims.empty() ? H5Screate(H5S_SCALAR) : H5Screate_simple(int(dims.size()), dims.data(), nullptr);
  hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
  hsize_t total = 1;
  for (auto d : dims) total *= d;
  if (!dims.empty() && total >= 4096) {
    // Chunks of one leading-index slice keep chunks well under 4 GB.
    std::vector<hsize_t> chunk = dims;
    for (std::size_t i = 0; i + 1 < chunk.size() && chunk.size() > 1; ++i) {
      hsize_t rest = 1;
      for (std::size_t j = i + 1; j < chunk.size(); ++j) rest *= chunk[j];
      if (rest * 8 > (hsize_t(1) << 26)) chunk[i] = 1;  // > 64 MB per chunk: slice this dimension
      else break;
    }
    H5Pset_chunk(dcpl, int(chunk.size()), chunk.data());
    if (deflate > 0) {  // byte shuffle + deflate (gzip) level 1-9; 0: stored uncompressed
      H5Pset_shuffle(dcpl);
      H5Pset_deflate(dcpl, unsigned(deflate));
    }
  }
  hid_t ds = check_id(H5Dcreate2(loc, name.c_str(), type, space, H5P_DEFAULT, dcpl, H5P_DEFAULT), "create dataset " + name);
  check(H5Dwrite(ds, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, data), "write dataset " + name);
  if (!dim_names.empty()) write_string_attr(ds, "dim_names", dim_names, false);
  H5Dclose(ds);
  H5Pclose(dcpl);
  H5Sclose(space);
}

template <class T>
std::vector<T> read_dataset(hid_t loc, const std::string& name, hid_t memtype, std::vector<hsize_t>* dims_out = nullptr) {
  hid_t ds = check_id(H5Dopen2(loc, name.c_str(), H5P_DEFAULT), "missing dataset " + name);
  hid_t space = H5Dget_space(ds);
  int nd = H5Sget_simple_extent_ndims(space);
  std::vector<hsize_t> dims(std::max(nd, 0));
  if (nd > 0) H5Sget_simple_extent_dims(space, dims.data(), nullptr);
  hsize_t total = 1;
  for (auto d : dims) total *= d;
  std::vector<T> out(total);
  check(H5Dread(ds, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, out.data()), "read dataset " + name);
  H5Sclose(space);
  H5Dclose(ds);
  if (dims_out) *dims_out = dims;
  return out;
}

std::vector<std::string> grid_dims() { return {"theta", "phi", "r"}; }

bool has_link(hid_t loc, const std::string& path) {  // dataset or group at path (intermediate groups must exist)
  std::size_t pos = 0;
  while (true) {
    pos = path.find('/', pos + 1);
    const std::string part = path.substr(0, pos);
    if (H5Lexists(loc, part.c_str(), H5P_DEFAULT) <= 0) return false;
    if (pos == std::string::npos) return true;
  }
}

// model code stored in snapshot files: which fluids and neutral populations they hold
int model_code(const RunConfig& cfg) { return cfg.multi_ion() ? 1 : cfg.extended_si() ? 2 : 0; }
const char* model_name(int code) {
  return code == 1 ? "multi_ion_inca" : code == 2 ? "extended_si" : "four_population_ts";
}

}  // namespace

// ---------------------------------------------------------------------------
// CheckpointWriter
// ---------------------------------------------------------------------------
CheckpointWriter::CheckpointWriter(const std::string& dir, const std::string& name) {
  make_dirs(dir);
  path_ = (dir.empty() ? std::string(".") : dir) + "/" + name + ".h5";
  file_ = check_id(H5Fcreate(path_.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT), "create " + path_);
  hid_t root = H5Gopen2(file_, "/", H5P_DEFAULT);
  std::time_t now = std::time(nullptr);
  char stamp[64];
  std::strftime(stamp, sizeof stamp, "%a %b %d %H:%M:%S %Y", std::localtime(&now));
  write_string_attr(root, "checkpoint", {name}, true);
  write_string_attr(root, "created", {stamp}, true);
  write_string_attr(root, "writer", {"ena C++"}, true);
  write_string_attr(root, "layout_note",
                    {"dim_names lists dimensions in HDF5 (C, slowest-varying first) order."}, true);
  H5Gclose(root);
}

CheckpointWriter::~CheckpointWriter() {
  if (file_ >= 0) H5Fclose(file_);
}

void CheckpointWriter::write(const std::string& name, const void* data, long type_id,
                             const std::vector<unsigned long long>& dims, const std::vector<std::string>& dim_names) {
  std::vector<hsize_t> d(dims.begin(), dims.end());
  write_dataset(file_, name, data, hid_t(type_id), d, dim_names);
}

void CheckpointWriter::scalar(const std::string& n, double v) { write(n, &v, H5T_NATIVE_DOUBLE, {}, {}); }
void CheckpointWriter::scalar_float(const std::string& n, float v) { write(n, &v, H5T_NATIVE_FLOAT, {}, {}); }
void CheckpointWriter::scalar_int(const std::string& n, long v) { write(n, &v, H5T_NATIVE_LONG, {}, {}); }
void CheckpointWriter::string(const std::string& n, const std::string& v) {
  hid_t type = H5Tcopy(H5T_C_S1);
  H5Tset_size(type, std::max<std::size_t>(1, v.size()));
  write(n, v.c_str(), type, {}, {});
  H5Tclose(type);
}
void CheckpointWriter::field(const std::string& n, const Field3D<double>& f) {
  write(n, f.raw().data(), H5T_NATIVE_DOUBLE, {(unsigned long long)f.ntheta(), (unsigned long long)f.nphi(), (unsigned long long)f.nr()}, grid_dims());
}
void CheckpointWriter::field(const std::string& n, const Field3D<float>& f) {
  write(n, f.raw().data(), H5T_NATIVE_FLOAT, {(unsigned long long)f.ntheta(), (unsigned long long)f.nphi(), (unsigned long long)f.nr()}, grid_dims());
}
void CheckpointWriter::field(const std::string& n, const Field3D<std::int16_t>& f) {
  write(n, f.raw().data(), H5T_NATIVE_SHORT, {(unsigned long long)f.ntheta(), (unsigned long long)f.nphi(), (unsigned long long)f.nr()}, grid_dims());
}
void CheckpointWriter::field(const std::string& n, const Field2D<float>& f) {
  write(n, f.raw().data(), H5T_NATIVE_FLOAT, {(unsigned long long)f.ntheta(), (unsigned long long)f.nphi()}, {"theta", "phi"});
}
void CheckpointWriter::field(const std::string& n, const Field2D<double>& f) {
  write(n, f.raw().data(), H5T_NATIVE_DOUBLE, {(unsigned long long)f.ntheta(), (unsigned long long)f.nphi()}, {"theta", "phi"});
}
void CheckpointWriter::field_dims(const std::string& n, const Field3D<double>& f, const std::vector<std::string>& dims) {
  write(n, f.raw().data(), H5T_NATIVE_DOUBLE, {(unsigned long long)f.ntheta(), (unsigned long long)f.nphi(), (unsigned long long)f.nr()}, dims);
}
void CheckpointWriter::field4(const std::string& n, const std::vector<double>& data, const Grid& g, int nslot) {
  if (data.size() != g.size() * std::size_t(nslot)) throw std::runtime_error("field4: size mismatch for " + n);
  write(n, data.data(), H5T_NATIVE_DOUBLE,
        {(unsigned long long)g.ntheta, (unsigned long long)g.nphi, (unsigned long long)g.nr, (unsigned long long)nslot},
        {"theta", "phi", "r", "n" + std::to_string(nslot)});
}
void CheckpointWriter::vector(const std::string& n, const std::vector<double>& v, const std::string& dim) {
  write(n, v.data(), H5T_NATIVE_DOUBLE, {(unsigned long long)v.size()}, {dim});
}
void CheckpointWriter::vector(const std::string& n, const std::vector<long>& v, const std::string& dim) {
  write(n, v.data(), H5T_NATIVE_LONG, {(unsigned long long)v.size()}, {dim});
}

// ---------------------------------------------------------------------------
// Snapshot files
// ---------------------------------------------------------------------------
void write_snapshot_h5(const std::string& path, const RunConfig& cfg, const Snapshot& s, const StreamlineSet* sl,
                       int deflate) {
  // every dataset of the snapshot file is written with this compression level
  auto write_dataset = [deflate](hid_t loc, const std::string& name, const void* data, hid_t type,
                                 const std::vector<hsize_t>& dims, const std::vector<std::string>& dim_names) {
    ena::write_dataset(loc, name, data, type, dims, dim_names, deflate);
  };
  const Grid& g = cfg.grid;
  hid_t f = check_id(H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT), "create " + path);
  hid_t root = H5Gopen2(f, "/", H5P_DEFAULT);
  write_string_attr(root, "format", {"ena-snapshot"}, true);
  // format_version 2 (version 2.8): all ion fluids and all neutral populations with their bulk
  // velocities, so multi-ion and extended single-ion snapshots can be stored; version-1 files
  // (single-ion, n and T of the total H only) are still read.
  write_num_attr(root, "format_version", 2);
  write_num_attr(root, "populations_code", model_code(cfg));
  write_string_attr(root, "populations", {cfg.populations}, true);
  write_num_attr(root, "time", s.time);
  write_num_attr(root, "dtheta", g.dtheta);
  write_num_attr(root, "dphi", g.dphi);
  write_num_attr(root, "dr", g.dr);
  write_num_attr(root, "r_inner", g.r_inner);
  write_num_attr(root, "nr", g.nr);
  write_string_attr(root, "model_source", {cfg.model_source}, true);
  write_string_attr(root, "units", {"n: cm^-3, T: K (already halved for standard-model inputs), velocities: cm/s, B: as input"}, true);
  H5Gclose(root);

  auto stack = [&](auto getter, std::size_t count) {
    std::vector<double> buf(count * g.size());
    for (std::size_t m = 0; m < count; ++m) {
      const auto& src = getter(m).raw();
      std::copy(src.begin(), src.end(), buf.begin() + m * g.size());
    }
    return buf;
  };
  std::vector<hsize_t> d4p = {s.plasma.size(), hsize_t(g.ntheta), hsize_t(g.nphi), hsize_t(g.nr)};
  std::vector<std::string> n4p = {"fluid", "theta", "phi", "r"};
  hid_t gp = H5Gcreate2(f, "plasma", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  write_dataset(gp, "n", stack([&](std::size_t m) -> const Field3D<double>& { return s.plasma[m].n; }, s.plasma.size()).data(), H5T_NATIVE_DOUBLE, d4p, n4p);
  write_dataset(gp, "T", stack([&](std::size_t m) -> const Field3D<double>& { return s.plasma[m].T; }, s.plasma.size()).data(), H5T_NATIVE_DOUBLE, d4p, n4p);
  write_dataset(gp, "vr", stack([&](std::size_t m) -> const Field3D<double>& { return s.plasma[m].vr; }, s.plasma.size()).data(), H5T_NATIVE_DOUBLE, d4p, n4p);
  write_dataset(gp, "vtheta", stack([&](std::size_t m) -> const Field3D<double>& { return s.plasma[m].vtheta; }, s.plasma.size()).data(), H5T_NATIVE_DOUBLE, d4p, n4p);
  write_dataset(gp, "vphi", stack([&](std::size_t m) -> const Field3D<double>& { return s.plasma[m].vphi; }, s.plasma.size()).data(), H5T_NATIVE_DOUBLE, d4p, n4p);
  write_dataset(gp, "vt", stack([&](std::size_t m) -> const Field3D<double>& { return s.plasma[m].vt; }, s.plasma.size()).data(), H5T_NATIVE_DOUBLE, d4p, n4p);
  H5Gclose(gp);

  std::vector<hsize_t> d4n = {s.neutrals.size(), hsize_t(g.ntheta), hsize_t(g.nphi), hsize_t(g.nr)};
  hid_t gn = H5Gcreate2(f, "neutral", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  const std::vector<std::string> n4n = {"population", "theta", "phi", "r"};
  write_dataset(gn, "n", stack([&](std::size_t m) -> const Field3D<double>& { return s.neutrals[m].n; }, s.neutrals.size()).data(), H5T_NATIVE_DOUBLE, d4n, n4n);
  write_dataset(gn, "T", stack([&](std::size_t m) -> const Field3D<double>& { return s.neutrals[m].T; }, s.neutrals.size()).data(), H5T_NATIVE_DOUBLE, d4n, n4n);
  // bulk velocities of the neutral populations (multi-ion and extended single-ion files; as read, km/s)
  if (!s.neutrals.empty() && !s.neutrals[0].vr.empty()) {
    for (auto [name, mem] : {std::pair<const char*, Field3D<double> NeutralPopulation::*>{"vr", &NeutralPopulation::vr},
                             {"vt", &NeutralPopulation::vt}, {"vx", &NeutralPopulation::vx},
                             {"vy", &NeutralPopulation::vy}, {"vz", &NeutralPopulation::vz}})
      write_dataset(gn, name, stack([&](std::size_t m) -> const Field3D<double>& { return s.neutrals[m].*mem; }, s.neutrals.size()).data(), H5T_NATIVE_DOUBLE, d4n, n4n);
  }
  H5Gclose(gn);

  std::vector<hsize_t> d3 = {hsize_t(g.ntheta), hsize_t(g.nphi), hsize_t(g.nr)};
  write_dataset(f, "region", s.region.raw().data(), H5T_NATIVE_SHORT, d3, grid_dims());
  write_dataset(f, "Bx", s.Bx.raw().data(), H5T_NATIVE_DOUBLE, d3, grid_dims());
  write_dataset(f, "By", s.By.raw().data(), H5T_NATIVE_DOUBLE, d3, grid_dims());
  write_dataset(f, "Bz", s.Bz.raw().data(), H5T_NATIVE_DOUBLE, d3, grid_dims());
  if (!s.pbeta.empty()) write_dataset(f, "pbeta", s.pbeta.raw().data(), H5T_NATIVE_DOUBLE, d3, grid_dims());

  if (sl && sl->points()) {
    hid_t gs = H5Gcreate2(f, "streamlines", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    write_dataset(gs, "values", sl->values.data(), H5T_NATIVE_DOUBLE, {sl->points(), hsize_t(sl->ncols)}, {"point", "column"});
    std::vector<long long> off(sl->offsets.begin(), sl->offsets.end());
    write_dataset(gs, "offsets", off.data(), H5T_NATIVE_LLONG, {off.size()}, {"streamline_plus_1"});
    std::vector<int> gaps(sl->gap_lines.begin(), sl->gap_lines.end());
    write_dataset(gs, "gap_lines", gaps.data(), H5T_NATIVE_INT, {gaps.size()}, {"streamline"});
    write_string_attr(gs, "columns", cfg.inputs.streamlines.columns, false);
    if (!sl->variables.empty()) write_string_attr(gs, "tecplot_variables", sl->variables, false);
    H5Gclose(gs);
  }
  H5Fclose(f);
}

Snapshot read_snapshot_h5(const RunConfig& cfg, const std::string& path, const ReadOptions&) {
  const Grid& g = cfg.grid;
  hid_t f = check_id(H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT), "open " + path);
  auto same = [](double a, double b) { return std::fabs(a - b) <= 1e-9 * std::max(1.0, std::fabs(b)); };
  if (!same(read_num_attr(f, "dtheta"), g.dtheta) || !same(read_num_attr(f, "dphi"), g.dphi) ||
      !same(read_num_attr(f, "dr"), g.dr) || !same(read_num_attr(f, "r_inner"), g.r_inner) ||
      !same(read_num_attr(f, "nr"), g.nr)) {
    H5Fclose(f);
    throw std::runtime_error(path + ": grid stored in the file differs from [grid] in the configuration");
  }
  // the file must hold the populations of the configured model (format 1: single-ion only)
  const int code = H5Aexists(f, "populations_code") > 0 ? int(read_num_attr(f, "populations_code")) : 0;
  if (code != model_code(cfg)) {
    H5Fclose(f);
    throw std::runtime_error(path + ": snapshot file made for model.populations = \"" + model_name(code) +
                             "\", but the configuration uses \"" + cfg.populations + "\"");
  }
  Snapshot s;
  s.time = read_num_attr(f, "time");
  auto unstack = [&](const std::vector<double>& buf, std::size_t m, Field3D<double>& out) {
    out.resize(g);
    std::copy(buf.begin() + m * g.size(), buf.begin() + (m + 1) * g.size(), out.raw().begin());
  };
  std::vector<hsize_t> dims;
  auto pn = read_dataset<double>(f, "plasma/n", H5T_NATIVE_DOUBLE, &dims);
  std::size_t nf = dims.at(0);
  auto pT = read_dataset<double>(f, "plasma/T", H5T_NATIVE_DOUBLE);
  auto pvr = read_dataset<double>(f, "plasma/vr", H5T_NATIVE_DOUBLE);
  auto pvth = read_dataset<double>(f, "plasma/vtheta", H5T_NATIVE_DOUBLE);
  auto pvph = read_dataset<double>(f, "plasma/vphi", H5T_NATIVE_DOUBLE);
  auto pvt = read_dataset<double>(f, "plasma/vt", H5T_NATIVE_DOUBLE);
  s.plasma.resize(nf);
  for (std::size_t m = 0; m < nf; ++m) {
    unstack(pn, m, s.plasma[m].n); unstack(pT, m, s.plasma[m].T); unstack(pvr, m, s.plasma[m].vr);
    unstack(pvth, m, s.plasma[m].vtheta); unstack(pvph, m, s.plasma[m].vphi); unstack(pvt, m, s.plasma[m].vt);
  }
  auto nn = read_dataset<double>(f, "neutral/n", H5T_NATIVE_DOUBLE, &dims);
  auto nT = read_dataset<double>(f, "neutral/T", H5T_NATIVE_DOUBLE);
  s.neutrals.resize(dims.at(0));
  for (std::size_t m = 0; m < s.neutrals.size(); ++m) { unstack(nn, m, s.neutrals[m].n); unstack(nT, m, s.neutrals[m].T); }
  if (has_link(f, "neutral/vr")) {
    for (auto [name, mem] : {std::pair<const char*, Field3D<double> NeutralPopulation::*>{"neutral/vr", &NeutralPopulation::vr},
                             {"neutral/vt", &NeutralPopulation::vt}, {"neutral/vx", &NeutralPopulation::vx},
                             {"neutral/vy", &NeutralPopulation::vy}, {"neutral/vz", &NeutralPopulation::vz}}) {
      auto buf = read_dataset<double>(f, name, H5T_NATIVE_DOUBLE);
      for (std::size_t m = 0; m < s.neutrals.size(); ++m) unstack(buf, m, s.neutrals[m].*mem);
    }
  }
  if (has_link(f, "pbeta")) {
    s.pbeta.resize(g);
    s.pbeta.raw() = read_dataset<double>(f, "pbeta", H5T_NATIVE_DOUBLE);
  }
  auto reg = read_dataset<std::int16_t>(f, "region", H5T_NATIVE_SHORT);
  s.region.resize(g);
  s.region.raw() = reg;
  for (auto [name, fld] : {std::pair<const char*, Field3D<double>*>{"Bx", &s.Bx}, {"By", &s.By}, {"Bz", &s.Bz}}) {
    fld->resize(g);
    fld->raw() = read_dataset<double>(f, name, H5T_NATIVE_DOUBLE);
  }
  H5Fclose(f);
  return s;
}

StreamlineSet read_streamlines_h5(const std::string& path) {
  hid_t f = check_id(H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT), "open " + path);
  StreamlineSet sl;
  std::vector<hsize_t> dims;
  sl.values = read_dataset<double>(f, "streamlines/values", H5T_NATIVE_DOUBLE, &dims);
  sl.ncols = int(dims.at(1));
  auto off = read_dataset<long long>(f, "streamlines/offsets", H5T_NATIVE_LLONG);
  sl.offsets.assign(off.begin(), off.end());
  auto gaps = read_dataset<int>(f, "streamlines/gap_lines", H5T_NATIVE_INT);
  sl.gap_lines.assign(gaps.begin(), gaps.end());
  H5Fclose(f);
  return sl;
}

Field3D<double> read_grid_field_h5(const std::string& path, const std::string& name, const Grid& g) {
  hid_t f = check_id(H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT), "cannot open " + path);
  std::vector<hsize_t> dims;
  std::vector<double> v;
  try {
    v = read_dataset<double>(f, name, H5T_NATIVE_DOUBLE, &dims);
  } catch (...) {
    H5Fclose(f);
    throw;
  }
  H5Fclose(f);
  Field3D<double> out(g, 0.0);
  const hsize_t nt = g.ntheta, np = g.nphi, nr = g.nr;
  if (dims == std::vector<hsize_t>{nt, np, nr}) {
    out.raw() = v;
  } else if (dims == std::vector<hsize_t>{nr, np, nt}) {  // IDL dump: [r][phi][theta] in C order
    for (hsize_t k = 0; k < nr; ++k)
      for (hsize_t j = 0; j < np; ++j)
        for (hsize_t i = 0; i < nt; ++i) out(int(i), int(j), int(k)) = v[(k * np + j) * nt + i];
  } else {
    throw std::runtime_error(path + ": dataset " + name + " does not have the grid's shape");
  }
  return out;
}

}  // namespace ena
