#pragma once
// Shared-memory parallelism (OpenMP), milestone 5.
//
// Every parallel loop in the code distributes independent work items (lines of
// sight, grid cells, text lines) over threads, and each item writes only its
// own outputs; sums along a line of sight stay sequential inside one item. The
// results are therefore bit-identical for any number of threads, including a
// build without OpenMP.
//
// The number of threads is set with ena_run --threads N, or the OMP_NUM_THREADS
// environment variable (default: all cores).

#include <cstddef>
#include <exception>
#include <mutex>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace ena {

inline int max_threads() {
#ifdef _OPENMP
  return omp_get_max_threads();
#else
  return 1;
#endif
}

inline void set_threads(int n) {
#ifdef _OPENMP
  if (n > 0) omp_set_num_threads(n);
#else
  (void)n;
#endif
}

inline bool openmp_enabled() {
#ifdef _OPENMP
  return true;
#else
  return false;
#endif
}

// Exceptions must not leave an OpenMP parallel region. Each work item runs its
// body inside capture(); after the loop, rethrow() raises the exception of the
// lowest-numbered failing item, which is the one a serial run would have
// raised first.
class ParallelErrors {
 public:
  template <class F>
  void capture(std::size_t item, F&& body) {
    try {
      body();
    } catch (...) {
      std::lock_guard<std::mutex> lock(m_);
      if (!err_ || item < item_) { err_ = std::current_exception(); item_ = item; }
    }
  }
  void rethrow() const {
    if (err_) std::rethrow_exception(err_);
  }

 private:
  std::mutex m_;
  std::exception_ptr err_;
  std::size_t item_ = 0;
};

}  // namespace ena
