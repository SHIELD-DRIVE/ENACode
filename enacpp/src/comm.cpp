#include "ena/comm.hpp"

#include <algorithm>
#include <climits>
#include <cstdlib>
#include <stdexcept>

#ifdef ENA_HAVE_MPI
#include <mpi.h>
#endif

namespace ena {

Comm::Comm() {
#ifdef ENA_HAVE_MPI
  int init = 0;
  MPI_Initialized(&init);
  if (init) {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    MPI_Comm_size(MPI_COMM_WORLD, &size_);
  }
#endif
}

Comm& Comm::world() {
  static Comm c;
  return c;
}

bool Comm::mpi_enabled() {
#ifdef ENA_HAVE_MPI
  return true;
#else
  return false;
#endif
}

void comm_init(int* argc, char*** argv) {
#ifdef ENA_HAVE_MPI
  // Only the main thread calls MPI; OpenMP threads work between the calls.
  int provided = 0;
  MPI_Init_thread(argc, argv, MPI_THREAD_FUNNELED, &provided);
#else
  (void)argc;
  (void)argv;
#endif
  Comm::world();  // records rank and size
}

void comm_finalize() {
#ifdef ENA_HAVE_MPI
  int fin = 0;
  MPI_Finalized(&fin);
  if (!fin) MPI_Finalize();
#endif
}

void comm_abort(int code) {
#ifdef ENA_HAVE_MPI
  MPI_Abort(MPI_COMM_WORLD, code);
#endif
  std::exit(code);
}

void Comm::agree(const std::string& error) const {
#ifdef ENA_HAVE_MPI
  if (size_ > 1) {
    // lowest failing rank (size_ if none)
    int mine = error.empty() ? size_ : rank_, first = size_;
    MPI_Allreduce(&mine, &first, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    if (first == size_) return;
    int len = rank_ == first ? int(error.size()) : 0;
    MPI_Bcast(&len, 1, MPI_INT, first, MPI_COMM_WORLD);
    std::string msg(std::size_t(len), ' ');
    if (rank_ == first) msg = error;
    MPI_Bcast(msg.data(), len, MPI_CHAR, first, MPI_COMM_WORLD);
    throw std::runtime_error(first == 0 ? msg : msg + " (on MPI rank " + std::to_string(first) + ")");
  }
#endif
  if (!error.empty()) throw std::runtime_error(error);
}

std::vector<std::vector<double>> Comm::alltoallv(std::vector<std::vector<double>>&& send) const {
  if (int(send.size()) != size_) throw std::logic_error("Comm::alltoallv: one buffer per rank needed");
#ifdef ENA_HAVE_MPI
  if (size_ > 1) {
    // Point-to-point messages straight from and into the per-rank buffers (no packing
    // copies): first the lengths, then the data.
    std::vector<long long> scount(size_), rcount(size_);
    for (int d = 0; d < size_; ++d) scount[d] = (long long)send[d].size();
    MPI_Alltoall(scount.data(), 1, MPI_LONG_LONG, rcount.data(), 1, MPI_LONG_LONG, MPI_COMM_WORLD);
    std::vector<std::vector<double>> recv(size_);
    std::vector<MPI_Request> req;
    for (int src = 0; src < size_; ++src) {
      if (src == rank_ || rcount[src] == 0) continue;
      if (rcount[src] > INT_MAX) throw std::runtime_error("MPI message too large (use fewer energies per pass)");
      recv[src].resize(std::size_t(rcount[src]));
      req.emplace_back();
      MPI_Irecv(recv[src].data(), int(rcount[src]), MPI_DOUBLE, src, 7, MPI_COMM_WORLD, &req.back());
    }
    for (int d = 0; d < size_; ++d) {
      if (d == rank_ || scount[d] == 0) continue;
      if (scount[d] > INT_MAX) throw std::runtime_error("MPI message too large (use fewer energies per pass)");
      req.emplace_back();
      MPI_Isend(send[d].data(), int(scount[d]), MPI_DOUBLE, d, 7, MPI_COMM_WORLD, &req.back());
    }
    recv[rank_] = std::move(send[rank_]);
    MPI_Waitall(int(req.size()), req.data(), MPI_STATUSES_IGNORE);
    send.clear();
    return recv;
  }
#endif
  return std::move(send);
}

std::vector<std::vector<double>> Comm::gather(const std::vector<double>& send) const {
#ifdef ENA_HAVE_MPI
  if (size_ > 1) {
    if (send.size() > std::size_t(INT_MAX)) throw std::runtime_error("MPI gather buffer too large");
    int n = int(send.size());
    std::vector<int> counts(size_), disp(size_);
    MPI_Gather(&n, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
    long tot = 0;
    if (rank_ == 0)
      for (int s = 0; s < size_; ++s) {
        disp[s] = int(tot);
        tot += counts[s];
        if (tot > INT_MAX) throw std::runtime_error("MPI gather buffer too large");
      }
    std::vector<double> all(rank_ == 0 ? std::size_t(tot) : 0);
    MPI_Gatherv(send.data(), n, MPI_DOUBLE, all.data(), counts.data(), disp.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    std::vector<std::vector<double>> out;
    if (rank_ == 0)
      for (int s = 0; s < size_; ++s) out.emplace_back(all.begin() + disp[s], all.begin() + disp[s] + counts[s]);
    return out;
  }
#endif
  return {send};
}

std::vector<std::string> Comm::gather_text(const std::string& text) const {
  std::vector<double> buf(text.begin(), text.end());  // characters as doubles (short log lines)
  std::vector<std::string> out;
  for (const auto& b : gather(buf)) {
    std::string s;
    for (double c : b) s += char(c);
    out.push_back(s);
  }
  return out;
}

double Comm::max_double(double v) const {
#ifdef ENA_HAVE_MPI
  if (size_ > 1) {
    double out = v;
    MPI_Allreduce(&v, &out, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    return out;
  }
#endif
  return v;
}

double Comm::sum_double(double v) const {
#ifdef ENA_HAVE_MPI
  if (size_ > 1) {
    double out = v;
    MPI_Allreduce(&v, &out, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return out;
  }
#endif
  return v;
}

void Comm::barrier() const {
#ifdef ENA_HAVE_MPI
  if (size_ > 1) MPI_Barrier(MPI_COMM_WORLD);
#endif
}

}  // namespace ena
