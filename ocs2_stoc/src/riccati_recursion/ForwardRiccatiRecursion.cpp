#include "ocs2_stoc/riccati_recursion/ForwardRiccatiRecursion.h"

namespace ocs2 {
namespace stoc {
namespace ForwardRiccatiRecursion {

void computeInput(const LqrPolicy& lqrPolicy, const vector_t& dx, vector_t& du, 
                  scalar_t dts, scalar_t dtsNext, bool sto, bool stoNext) {
  if (lqrPolicy.k.size() > 0) {
    du.noalias()  = lqrPolicy.K * dx;
    du.noalias() += lqrPolicy.k;
    if (sto) {
      du.noalias() += lqrPolicy.T * (dtsNext-dts);
      if (stoNext) {
        du.noalias() -= lqrPolicy.W * dtsNext;
      }
    }
  }
  else {
    du.resize(0);
  }
}


void computeState(const ipm::ModelData& modelData, 
                  const vector_t& dx, const vector_t& du, vector_t& dxNext, scalar_t dts, scalar_t dtsNext, bool sto) {
  const auto& dynamics = modelData.dynamics;
  dxNext = dynamics.f;
  dxNext.noalias() += dynamics.dfdx * dx;
  if (du.size() > 0) {
    dxNext.noalias() += dynamics.dfdu * du;
  }
  if (sto) {
    dxNext.noalias() += modelData.hamiltonian.dfdt * (dtsNext-dts);
  }
}


void computeCostate(const RiccatiRecursionData& riccati, const vector_t& dx, vector_t& dlmd, 
                    scalar_t dts, scalar_t dtsNext, bool sto, bool stoNext) {
  dlmd.noalias() = riccati.P * dx - riccati.s;
  if (sto) {
    dlmd.noalias() += riccati.Psi * (dtsNext-dts);
    if (stoNext) {
      dlmd.noalias() -= riccati.Phi * dtsNext;
    }
  }
}


scalar_t computeSwitchingTime(const StoPolicy& stoPolicy, const vector_t& dx, scalar_t dtsPrev, bool stoPrev) {
  scalar_t dts = stoPolicy.dtsdx.dot(dx) + stoPolicy.dts0;
  if (stoPrev) {
    dts += stoPolicy.dtsdts * dtsPrev;
  }
  return dts;
}

} // namespace ForwardRiccatiRecursion 
} // namespace stoc
} // namespace ocs2