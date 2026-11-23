#pragma once
// Shared by ena_fieldtopo and ena_fieldtopo_mhd: writing a copy of a boundaries file with some region
// codes changed (every other character kept), the copy's name, and the --snapshots argument.

#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "ena/ascii.hpp"
#include "ena/config.hpp"

namespace ena_topo {

using namespace ena;

// "file.txt" + "_topo" -> "file_topo.txt" (the suffix goes before the extension of the file name)
inline std::string with_suffix(const std::string& path, const std::string& suffix) {
  const std::size_t slash = path.find_last_of('/');
  const std::size_t dot = path.find_last_of('.');
  if (dot == std::string::npos || (slash != std::string::npos && dot < slash)) return path + suffix;
  return path.substr(0, dot) + suffix + path.substr(dot);
}

// the region token rewritten from code `from` to code `to`, in the token's own number format where
// possible ("3.000000E+00" -> "2.000000E+00"), otherwise as a plain integer
inline std::string region_token(std::string_view tok, int from, int to) {
  std::string t(tok);
  const std::size_t p = (!t.empty() && (t[0] == '+' || t[0] == '-')) ? 1 : 0;
  if (p < t.size() && t[p] == char('0' + from) && from >= 0 && from <= 9 && to >= 0 && to <= 9) {
    t[p] = char('0' + to);
    char* end = nullptr;
    const double v = std::strtod(t.c_str(), &end);
    if (end && *end == '\0' && int(std::trunc(v)) == to) return t;
  }
  return std::to_string(to);
}

// Copies the boundaries file `in` to `out`, replacing the region code of every row whose cell has
// new_code >= 0; every other character is kept. Returns the number of rows changed.
inline long rewrite_boundaries(const RunConfig& cfg, const std::string& in, const std::string& out,
                        const std::vector<int>& new_code) {
  const ColumnMap& m = cfg.inputs.boundaries;
  const Grid& g = cfg.grid;
  const int ncols = int(m.columns.size());
  const int creg = m.column("region"), ci = m.column("i"), cj = m.column("j"), ck = m.column("k");
  if (creg < 0 || ci < 0 || cj < 0 || ck < 0)
    throw std::runtime_error("the boundaries column map needs region, i, j and k columns");
  const std::string tmp = out + ".part";
  FILE* fo = std::fopen(tmp.c_str(), "wb");
  if (!fo) throw std::runtime_error("cannot write " + tmp);
  long changed = 0;
  try {
    for_each_parsed_line(in, ncols, [](std::string_view) { return true; },
                         [&](std::string_view text, const double* v, int n, long line) {
      std::string row(text);
      if (n != 0) {  // blank lines are copied as they are
        if (n != ncols)
          throw std::runtime_error(in + ":" + std::to_string(line) + ": expected " + std::to_string(ncols) +
                                   " numbers, found " + (n < 0 ? std::string("a non-numeric value") : std::to_string(n)));
        const int i = int(std::lround(v[ci])), j = int(std::lround(v[cj])), k = int(std::lround(v[ck]));
        if (i < 0 || i >= g.ntheta || j < 0 || j >= g.nphi || k < 0 || k >= g.nr)
          throw std::runtime_error(in + ":" + std::to_string(line) + ": cell index outside the grid");
        const int to = new_code[g.index(i, j, k)];
        const int from = int(std::trunc(v[creg]));
        if (to >= 0 && to != from) {
          // find the region column's token: whitespace-separated, as the reader takes it
          std::size_t pos = 0, a = 0, b = 0;
          for (int c = 0; c <= creg; ++c) {
            while (pos < row.size() && std::isspace((unsigned char)row[pos])) ++pos;
            a = pos;
            while (pos < row.size() && !std::isspace((unsigned char)row[pos])) ++pos;
            b = pos;
          }
          row = row.substr(0, a) + region_token(std::string_view(row).substr(a, b - a), from, to) + row.substr(b);
          ++changed;
        }
      }
      row += '\n';
      if (std::fwrite(row.data(), 1, row.size(), fo) != row.size())
        throw std::runtime_error("error writing " + tmp + " (disk full?)");
    });
  } catch (...) {
    std::fclose(fo);
    std::remove(tmp.c_str());
    throw;
  }
  const bool ok = std::ferror(fo) == 0;
  if (std::fclose(fo) != 0 || !ok) throw std::runtime_error("error writing " + tmp + " (disk full?)");
  if (std::rename(tmp.c_str(), out.c_str()) != 0) throw std::runtime_error("cannot rename " + tmp);
  return changed;
}

// --snapshots: "all", "A-B" or "A,B,C" (indices into the configuration's snapshot list, from 0)
inline std::vector<std::size_t> parse_snapshots(const std::string& v, std::size_t count) {
  std::vector<std::size_t> out;
  if (v == "all") {
    for (std::size_t s = 0; s < count; ++s) out.push_back(s);
    return out;
  }
  std::stringstream ss(v);
  std::string tok;
  while (std::getline(ss, tok, ',')) {
    if (tok.empty()) continue;
    const std::size_t dash = tok.find('-');
    const std::size_t a = std::stoul(tok.substr(0, dash));
    const std::size_t b = dash == std::string::npos ? a : std::stoul(tok.substr(dash + 1));
    for (std::size_t s = a; s <= b; ++s) out.push_back(s);
  }
  for (std::size_t s : out)
    if (s >= count) throw std::runtime_error("--snapshots: index " + std::to_string(s) + " out of range (0-" +
                                             std::to_string(count - 1) + ")");
  return out;
}

}  // namespace ena_topo
