#pragma once
// Fast reading of whitespace-separated numeric text files.

#include <cstddef>
#include <functional>
#include <string>
#include <string_view>
#include <vector>

namespace ena {

// Reads one number starting exactly at p (no blanks skipped, no leading '+'), as
// std::from_chars does: returns the position after it, or nullptr if p holds no number.
// Uses std::from_chars where the C++ library supports it for floating point (GCC 11 or
// newer), otherwise std::strtod on a copy of the token (GCC 9 and 10). Both round
// correctly, so the values are identical either way.
const char* read_double(const char* p, const char* end, double& v);

// Parses up to `max` numbers from a line into `out`. Returns the count parsed,
// or -1 if a token is not a number.
int parse_numbers(std::string_view line, double* out, int max);

// Calls `fn(line_text, line_number)` for every line of a text file, reading it
// in large blocks. Line numbers start at 1.
void for_each_line(const std::string& path,
                   const std::function<void(std::string_view, long)>& fn);

// Count passed to for_each_parsed_line's fn for a line that `numeric` rejected.
constexpr int kNotNumeric = -2;

// Calls `fn(line_text, values, count, line_number)` for every line of a text
// file, in order. Lines for which `numeric(line_text)` is true are parsed
// first (in parallel, by blocks of lines): count is parse_numbers' result
// (-1: a non-numeric token; max_cols + 1: too many numbers). Other lines get
// count kNotNumeric. `numeric` must be safe to call from several threads.
void for_each_parsed_line(const std::string& path, int max_cols,
                          const std::function<bool(std::string_view)>& numeric,
                          const std::function<void(std::string_view, const double*, int, long)>& fn);

// Calls `fn(values, line_number)` for every non-blank line, which must hold
// exactly `ncols` numbers; any other line is an error naming the line.
void for_each_row(const std::string& path, int ncols,
                  const std::function<void(const double*, long)>& fn);

}  // namespace ena
