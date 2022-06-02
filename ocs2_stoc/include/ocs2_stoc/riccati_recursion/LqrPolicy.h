#pragma once

#include <ocs2_core/Types.h>

namespace ocs2 {
namespace stoc {

struct LqrPolicy {
  matrix_t K;
  vector_t k;
  vector_t T;
  vector_t W;

  void resize(const size_t nx, const size_t nu) {
    K.resize(nu, nx);
    k.resize(nu);
    T.resize(nu);
    W.resize(nu);
  }

  void setZero() {
    K.setZero();
    k.setZero();
    T.setZero();
    W.setZero();
  }

  void swap(LqrPolicy& other) {
    K.swap(other.K);
    k.swap(other.k);
    T.swap(other.T);
    W.swap(other.W);
  }
};

} // namespace stoc
} // namespace ocs2