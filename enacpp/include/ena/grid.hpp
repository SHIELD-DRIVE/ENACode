#pragma once
// Spherical grid and 3-D field storage.
//
// Grid points (not cell centres): theta_i = i*dtheta (i = 0 is the north pole),
// phi_j = j*dphi, r_k = r_inner + k*dr. Fields are stored contiguously in the
// order [theta][phi][r], so each line of sight (fixed theta, phi) is a
// contiguous run of nr values.

#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

namespace ena {

struct Grid {
  double dtheta = 6.0;   // deg
  double dphi = 6.0;     // deg
  double dr = 2.0;       // AU
  double r_inner = 2.0;  // AU, radius of index k = 0
  int ntheta = 30;
  int nphi = 60;
  int nr = 736;

  std::size_t size() const { return std::size_t(ntheta) * nphi * nr; }
  std::size_t index(int i, int j, int k) const {
    return (std::size_t(i) * nphi + j) * nr + k;
  }
  double theta(int i) const { return i * dtheta; }
  double phi(int j) const { return j * dphi; }
  double r(int k) const { return r_inner + k * dr; }
  // Radial index of radius r (nearest grid point).
  int k_of(double radius) const;

  void validate() const;
};

template <class T>
class Field3D {
 public:
  using value_type = T;
  Field3D() = default;
  explicit Field3D(const Grid& g, T fill = T()) { resize(g, fill); }
  void resize(const Grid& g, T fill = T()) {
    nt_ = g.ntheta; np_ = g.nphi; nr_ = g.nr;
    data_.assign(g.size(), fill);
  }
  T& operator()(int i, int j, int k) { return data_[(std::size_t(i) * np_ + j) * nr_ + k]; }
  const T& operator()(int i, int j, int k) const { return data_[(std::size_t(i) * np_ + j) * nr_ + k]; }
  T* los(int i, int j) { return &data_[(std::size_t(i) * np_ + j) * nr_]; }
  const T* los(int i, int j) const { return &data_[(std::size_t(i) * np_ + j) * nr_]; }
  std::vector<T>& raw() { return data_; }
  const std::vector<T>& raw() const { return data_; }
  bool empty() const { return data_.empty(); }
  int ntheta() const { return nt_; }
  int nphi() const { return np_; }
  int nr() const { return nr_; }

 private:
  int nt_ = 0, np_ = 0, nr_ = 0;
  std::vector<T> data_;
};

template <class T>
class Field2D {  // [theta][phi]
 public:
  Field2D() = default;
  explicit Field2D(const Grid& g, T fill = T()) { resize(g, fill); }
  void resize(const Grid& g, T fill = T()) {
    nt_ = g.ntheta; np_ = g.nphi;
    data_.assign(std::size_t(nt_) * np_, fill);
  }
  T& operator()(int i, int j) { return data_[std::size_t(i) * np_ + j]; }
  const T& operator()(int i, int j) const { return data_[std::size_t(i) * np_ + j]; }
  std::vector<T>& raw() { return data_; }
  const std::vector<T>& raw() const { return data_; }
  int ntheta() const { return nt_; }
  int nphi() const { return np_; }

 private:
  int nt_ = 0, np_ = 0;
  std::vector<T> data_;
};

}  // namespace ena
