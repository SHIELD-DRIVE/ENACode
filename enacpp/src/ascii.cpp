#include "ena/ascii.hpp"

#include <charconv>
#include <cstdio>
#include <memory>
#include <stdexcept>

namespace ena {

int parse_numbers(std::string_view line, double* out, int max) {
  const char* p = line.data();
  const char* end = p + line.size();
  int n = 0;
  for (;;) {
    while (p < end && (*p == ' ' || *p == '\t' || *p == '\r' || *p == ',')) ++p;
    if (p >= end) break;
    if (n >= max) return max + 1;  // too many fields
    if (*p == '+') ++p;
    double v;
    auto res = std::from_chars(p, end, v);
    if (res.ec != std::errc()) {
      // Accept Fortran-style 'D' exponents and NaN/Inf spellings.
      const char* q = p;
      std::string tok;
      while (q < end && *q != ' ' && *q != '\t' && *q != '\r' && *q != ',') tok += (*q == 'D' || *q == 'd') ? 'E' : *q, ++q;
      char* e2 = nullptr;
      v = std::strtod(tok.c_str(), &e2);
      if (tok.empty() || *e2 != '\0') return -1;
      p = q;
    } else {
      p = res.ptr;
      if (p < end && *p != ' ' && *p != '\t' && *p != '\r' && *p != ',') return -1;
    }
    out[n++] = v;
  }
  return n;
}

void for_each_line(const std::string& path, const std::function<void(std::string_view, long)>& fn) {
  std::unique_ptr<FILE, int (*)(FILE*)> f(std::fopen(path.c_str(), "rb"), &std::fclose);
  if (!f) throw std::runtime_error("cannot open input file " + path);
  const std::size_t block = std::size_t(1) << 24;
  std::string buf;
  std::vector<char> chunk(block);
  long line = 0;
  for (;;) {
    std::size_t got = std::fread(chunk.data(), 1, block, f.get());
    if (got == 0) break;
    buf.append(chunk.data(), got);
    std::size_t start = 0;
    for (;;) {
      std::size_t nl = buf.find('\n', start);
      if (nl == std::string::npos) break;
      fn(std::string_view(buf.data() + start, nl - start), ++line);
      start = nl + 1;
    }
    buf.erase(0, start);
  }
  if (!buf.empty()) fn(std::string_view(buf), ++line);
}

void for_each_parsed_line(const std::string& path, int max_cols,
                          const std::function<bool(std::string_view)>& numeric,
                          const std::function<void(std::string_view, const double*, int, long)>& fn) {
  // The file is read in blocks of whole lines. The numeric lines of a block are
  // parsed in parallel (the costly part), then every line is handed to fn in
  // file order, so the result and the first error reported are those of a
  // serial read, for any number of threads.
  std::unique_ptr<FILE, int (*)(FILE*)> f(std::fopen(path.c_str(), "rb"), &std::fclose);
  if (!f) throw std::runtime_error("cannot open input file " + path);
  const std::size_t block = std::size_t(1) << 25;
  std::string buf;
  std::vector<char> chunk(block);
  std::vector<std::size_t> starts, ends;
  std::vector<double> vals;
  std::vector<int> counts;
  const std::size_t stride = std::size_t(max_cols) + 1;
  long line0 = 0;  // lines handed out before this block
  bool eof = false;
  while (!eof) {
    const std::size_t got = std::fread(chunk.data(), 1, block, f.get());
    if (got == 0) eof = true;
    else buf.append(chunk.data(), got);
    // complete lines in buf (at end of file, also the last line without '\n')
    starts.clear(); ends.clear();
    std::size_t s = 0;
    for (;;) {
      const std::size_t nl = buf.find('\n', s);
      if (nl == std::string::npos) break;
      starts.push_back(s); ends.push_back(nl); s = nl + 1;
    }
    if (eof && s < buf.size()) { starts.push_back(s); ends.push_back(buf.size()); s = buf.size(); }
    const long nlines = long(starts.size());
    vals.resize(std::size_t(nlines) * stride);
    counts.resize(std::size_t(nlines));
#pragma omp parallel for schedule(static)
    for (long i = 0; i < nlines; ++i) {
      const std::string_view text(buf.data() + starts[i], ends[i] - starts[i]);
      counts[i] = numeric(text) ? parse_numbers(text, &vals[std::size_t(i) * stride], max_cols) : kNotNumeric;
    }
    for (long i = 0; i < nlines; ++i)
      fn(std::string_view(buf.data() + starts[i], ends[i] - starts[i]), &vals[std::size_t(i) * stride], counts[i],
         line0 + i + 1);
    line0 += nlines;
    buf.erase(0, s);
  }
}

void for_each_row(const std::string& path, int ncols, const std::function<void(const double*, long)>& fn) {
  for_each_parsed_line(path, ncols, [](std::string_view) { return true; },
                       [&](std::string_view, const double* v, int n, long line) {
                         if (n == 0) return;  // blank line
                         if (n != ncols)
                           throw std::runtime_error(path + ":" + std::to_string(line) + ": expected " +
                                                    std::to_string(ncols) + " numbers, found " +
                                                    (n < 0 ? std::string("a non-numeric value") : std::to_string(n)));
                         fn(v, line);
                       });
}

}  // namespace ena
