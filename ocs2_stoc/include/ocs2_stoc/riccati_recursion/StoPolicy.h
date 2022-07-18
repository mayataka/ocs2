#pragma once

#include <utility>
#include <ocs2_core/Types.h>

namespace ocs2 {
namespace stoc {

struct StoPolicy {
  vector_t dtsdx;
  scalar_t dtsdts;
  scalar_t dts0;

  void resize(size_t nx) {
    dtsdx.resize(nx);
  }

  void setZero() {
    dtsdx.setZero();
    dtsdts = 0.0;
    dts0 = 0.0;
  }

  void swap(StoPolicy& other) {
    dtsdx.swap(other.dtsdx);
    std::swap(dtsdts, other.dtsdts);
    std::swap(dts0, other.dts0);
  }
};

} // namespace stoc
} // namespace ocs2