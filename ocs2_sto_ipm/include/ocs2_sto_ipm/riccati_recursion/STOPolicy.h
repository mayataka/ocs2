#pragma once

#include <utility>
#include <ocs2_core/Types.h>

namespace ocs2 {
namespace sto_ipm {

struct STOPolicy {
  vector_t dtsdx;
  scalar_t dtsdts;
  scalar_t dts0;

  void resize(const size_t nx) {
    dtsdx.resize(nx);
  }

  void setZero() {
    dtsdx.setZero();
    dtsdts = 0.0;
    dts0 = 0.0;
  }

  void swap(STOPolicy& other) {
    dtsdx.swap(other.dtsdx);
    std::swap(dtsdts, other.dtsdts);
    std::swap(dts0, other.dts0);
  }
};

} // namespace sto_ipm
} // namespace ocs2