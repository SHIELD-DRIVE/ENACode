#include "ena/streamlines.hpp"

#include <stdexcept>
#include <string_view>

#include "ena/ascii.hpp"

namespace ena {

namespace {
bool starts_with_number(std::string_view s) {
  for (char c : s) {
    if (c == ' ' || c == '\t' || c == '\r') continue;
    return (c >= '0' && c <= '9') || c == '-' || c == '+' || c == '.';
  }
  return false;
}
void collect_quoted(std::string_view s, std::vector<std::string>& out) {
  std::size_t p = 0;
  while ((p = s.find('"', p)) != std::string_view::npos) {
    std::size_t q = s.find('"', p + 1);
    if (q == std::string_view::npos) break;
    out.emplace_back(s.substr(p + 1, q - p - 1));
    p = q + 1;
  }
}
}  // namespace

StreamlineSet read_streamlines_tecplot(const std::string& path, int expected_columns) {
  StreamlineSet set;
  set.ncols = expected_columns;
  bool in_variables = false;
  bool in_data = false;   // currently inside a run of data rows
  int gap = 0;            // non-data lines since the last data row
  // numeric lines are parsed in parallel; this function sees them in file order
  for_each_parsed_line(path, expected_columns, starts_with_number,
                       [&](std::string_view s, const double* row, int n, long line) {
    if (n != kNotNumeric) {
      if (n != expected_columns)
        throw std::runtime_error(path + ":" + std::to_string(line) + ": expected " + std::to_string(expected_columns) +
                                 " numbers in a streamline row, found " + (n < 0 ? std::string("a non-numeric value") : std::to_string(n)));
      if (!in_data) {
        set.offsets.push_back(set.points());
        set.gap_lines.push_back(gap);
        in_data = true;
      }
      set.values.insert(set.values.end(), row, row + expected_columns);
      gap = 0;
      in_variables = false;
      return;
    }
    // Header / zone record / blank line.
    in_data = false;
    ++gap;
    std::string_view t = s;
    while (!t.empty() && (t.front() == ' ' || t.front() == '\t')) t.remove_prefix(1);
    if (t.substr(0, 9) == "VARIABLES") { in_variables = true; collect_quoted(t, set.variables); return; }
    if (in_variables && !t.empty() && t.front() == '"') { collect_quoted(t, set.variables); return; }
    in_variables = false;
  });
  if (set.points() == 0) throw std::runtime_error(path + ": no streamline data rows found");
  set.offsets.push_back(set.points());
  if (!set.variables.empty() && int(set.variables.size()) != expected_columns)
    throw std::runtime_error(path + ": VARIABLES lists " + std::to_string(set.variables.size()) +
                             " names but the column map has " + std::to_string(expected_columns));
  return set;
}

}  // namespace ena
