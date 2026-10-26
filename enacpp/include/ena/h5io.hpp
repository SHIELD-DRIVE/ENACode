#pragma once
// HDF5 input/output: checkpoint dumps (same layout as dump_reference.pro, for
// validate.py) and converted snapshot files.

#include <cstdint>
#include <string>
#include <vector>

#include "ena/config.hpp"
#include "ena/grid.hpp"
#include "ena/snapshot.hpp"
#include "ena/streamlines.hpp"

namespace ena {

// One checkpoint file: <dir>/<name>.h5, one dataset per array. Every array
// gets a 'dim_names' attribute in HDF5 (C) order, e.g. ['theta','phi','r'].
class CheckpointWriter {
 public:
  CheckpointWriter(const std::string& dir, const std::string& name);
  ~CheckpointWriter();
  CheckpointWriter(const CheckpointWriter&) = delete;
  CheckpointWriter& operator=(const CheckpointWriter&) = delete;

  void scalar(const std::string& name, double v);
  void scalar_float(const std::string& name, float v);
  void scalar_int(const std::string& name, long v);
  void string(const std::string& name, const std::string& v);
  void field(const std::string& name, const Field3D<double>& f);
  void field(const std::string& name, const Field3D<float>& f);
  void field(const std::string& name, const Field3D<std::int16_t>& f);
  void field(const std::string& name, const Field2D<float>& f);
  void field(const std::string& name, const Field2D<double>& f);
  void field_dims(const std::string& name, const Field3D<double>& f, const std::vector<std::string>& dims);
  // [theta][phi][r][slot] array (slot fastest), e.g. per-species fields.
  void field4(const std::string& name, const std::vector<double>& data, const Grid& g, int nslot);
  void vector(const std::string& name, const std::vector<double>& v, const std::string& dim_name);
  void vector(const std::string& name, const std::vector<long>& v, const std::string& dim_name);

  const std::string& path() const { return path_; }

 private:
  void write(const std::string& name, const void* data, long type_id, const std::vector<unsigned long long>& dims,
             const std::vector<std::string>& dim_names);
  long file_ = -1;
  std::string path_;
};

// Snapshot files written by ena_convert: fields in internal units, plus the
// streamlines. Readers check the grid stored in the file against the config.
// `deflate`: gzip compression level 0-9 (0: uncompressed; smaller files, slower reads with higher levels).
void write_snapshot_h5(const std::string& path, const RunConfig& cfg, const Snapshot& s, const StreamlineSet* streamlines,
                       int deflate = 0);
Snapshot read_snapshot_h5(const RunConfig& cfg, const std::string& path, const ReadOptions& opt);
StreamlineSet read_streamlines_h5(const std::string& path);

// One grid-shaped dataset, in either C++ order [theta][phi][r] or the IDL
// dump order [r][phi][theta].
Field3D<double> read_grid_field_h5(const std::string& path, const std::string& name, const Grid& g);

}  // namespace ena
