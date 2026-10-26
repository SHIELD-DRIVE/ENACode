#pragma once
// Communication between MPI ranks (milestone 6). The time-dependent method uses
// only these operations. In a build without MPI (or with one rank) they are
// trivial: one rank holds everything and the exchange is a copy, so the same
// code path runs in every layout.
//
// Build with MPI: cmake -DENA_WITH_MPI=ON (needs an MPI library with mpicxx);
// run with e.g. mpiexec -n 8 ena_run run.toml --threads 16.

#include <cstddef>
#include <string>
#include <vector>

namespace ena {

class Comm {
 public:
  // The world communicator (MPI_COMM_WORLD), or a single rank without MPI.
  static Comm& world();
  int rank() const { return rank_; }
  int size() const { return size_; }
  bool root() const { return rank_ == 0; }
  static bool mpi_enabled();

  // Collective error agreement: every rank calls it with its own error message
  // ("" = no error). If any rank failed, every rank throws std::runtime_error
  // with the message of the lowest failing rank (so all ranks stop together
  // instead of waiting for each other forever).
  void agree(const std::string& error) const;

  // All-to-all exchange of variable-length double buffers: send[d] goes to rank
  // d; returns recv[s], the buffer received from rank s. `send` is consumed
  // (its memory is released), so the data exists at most twice per rank.
  std::vector<std::vector<double>> alltoallv(std::vector<std::vector<double>>&& send) const;

  // Gathers every rank's buffer on rank 0 (in rank order); other ranks get {}.
  std::vector<std::vector<double>> gather(const std::vector<double>& send) const;
  // Same for one text line per rank (e.g. log lines).
  std::vector<std::string> gather_text(const std::string& text) const;

  double max_double(double v) const;  // maximum over ranks, on every rank
  double sum_double(double v) const;  // sum over ranks, on every rank
  void barrier() const;

 private:
  Comm();
  int rank_ = 0, size_ = 1;
};

// Initialization and finalization (ena_run main); no-ops without MPI.
void comm_init(int* argc, char*** argv);
void comm_finalize();
// Aborts all ranks (after an error that not every rank sees).
[[noreturn]] void comm_abort(int code);

// Contiguous block partition of n items over `parts` parts: items
// [block_begin(n, parts, i), block_begin(n, parts, i + 1)) belong to part i;
// sizes differ by at most one.
inline long block_begin(long n, int parts, int i) {
  const long q = n / parts, r = n % parts;
  return long(i) * q + (i < r ? i : r);
}

}  // namespace ena
