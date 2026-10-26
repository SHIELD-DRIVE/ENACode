#pragma once
// Streamlines traced through the MHD solution (Tecplot ASCII, one ZONE per
// streamline, DATAPACKING=POINT).

#include <string>
#include <vector>

#include "ena/config.hpp"

namespace ena {

struct StreamlineSet {
  std::vector<std::string> variables;   // names from the Tecplot VARIABLES list
  int ncols = 0;
  std::vector<double> values;           // all points, row-major [point][column]
  std::vector<std::size_t> offsets;     // start point of each streamline; size = count + 1
  std::vector<int> gap_lines;           // non-data lines before each streamline (header for the first)

  std::size_t count() const { return offsets.empty() ? 0 : offsets.size() - 1; }
  std::size_t points() const { return ncols ? values.size() / ncols : 0; }
  double at(std::size_t point, int column) const { return values[point * ncols + column]; }
};

// Reads a Tecplot point-format file with any number of zones. If the file has
// no ZONE records, blank or non-numeric lines separate the streamlines.
StreamlineSet read_streamlines_tecplot(const std::string& path, int expected_columns);

}  // namespace ena
