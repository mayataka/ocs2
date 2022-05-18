#include <ocs2_sto_ipm/riccati_recursion/ForwardRiccatiRecursion.h>

#include <cmath>
#include <exception>

namespace ocs2 {
namespace sto_ipm {
namespace ForwardRiccatiRecursion {

void computeInput(const LQRPolicy& lqr_policy, const vector_t& dx, vector_t& du, 
                  const scalar_t dts, const scalar_t dts_next, 
                  const bool sto, const bool has_next_sto_phase) {
  du.noalias()  = lqr_policy.K * dx;
  du.noalias() += lqr_policy.k;
  if (sto) {
    du.noalias() += lqr_policy.T * (dts_next-dts);
    if (has_next_sto_phase) {
      du.noalias() -= lqr_policy.W * dts;
    }
  }
}


void computeState(const VectorFunctionLinearApproximationWrapper& dynamics, 
                  const vector_t& dx, const vector_t& du, vector_t& dx_next, 
                  const scalar_t dts, const scalar_t dts_next, const bool sto) {
  dx_next = dynamics.base.f;
  dx_next.noalias() += dynamics.base.dfdx * dx;
  dx_next.noalias() += dynamics.base.dfdu * du;
  if (sto) {
    dx_next.noalias() += dynamics.dfdt * (dts_next-dts);
  }
}


void computeCostate(const RiccatiRecursionData& riccati, const vector_t& dx,
                    vector_t& dlmd, const scalar_t dts, const scalar_t dts_next,  
                    const bool sto, const bool has_next_sto_phase) {
  dlmd.noalias() = riccati.P * dx - riccati.s;
  if (sto) {
    dlmd.noalias() += riccati.Psi * (dts_next-dts);
    if (has_next_sto_phase) {
      dlmd.noalias() -= riccati.Phi * dts_next;
    }
  }
}


double computeSwitchingTime(const STOPolicy& sto_policy, const vector_t& dx, 
                            const double dts_prev, const bool has_prev_sto_phase) {
  double dts = sto_policy.dtsdx.dot(dx) + sto_policy.dts0;
  if (has_prev_sto_phase) {
    dts += sto_policy.dtsdts * dts_prev;
  }
  return dts;
}

} // namespace ForwardRiccatiRecursion 
} // namespace sto_ipm
} // namespace ocs2