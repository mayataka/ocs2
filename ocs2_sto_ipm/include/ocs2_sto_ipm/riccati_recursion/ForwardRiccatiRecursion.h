#pragma once

#include <ocs2_sto_ipm/Types.h>
#include <ocs2_sto_ipm/riccati_recursion/RiccatiRecursionData.h>
#include <ocs2_sto_ipm/riccati_recursion/LQRPolicy.h>
#include <ocs2_sto_ipm/riccati_recursion/STOPolicy.h>

namespace ocs2 {
namespace sto_ipm {
namespace ForwardRiccatiRecursion {

void computeInput(const LQRPolicy& lqr_policy, const vector_t& dx, vector_t& du, 
                  const scalar_t dts=0.0, const scalar_t dts_next=0.0, 
                  const bool sto=false, const bool has_next_sto_phase=false);

void computeState(const VectorFunctionLinearApproximationWrapper& dynamics,
                  const vector_t& dx, const vector_t& du, vector_t& dx_next, 
                  const scalar_t dt=0.0, const scalar_t dts_next=0.0, 
                  const bool sto=false);

void computeCostate(const RiccatiRecursionData& riccati, const vector_t& dx, 
                    vector_t& dlmd, const scalar_t dts=0.0, 
                    const scalar_t dts_next=0.0, const bool sto=false, 
                    const bool has_next_sto_phase=false);

double computeSwitchingTime(const STOPolicy& sto_policy, const vector_t& dx, 
                            const double dts_prev, const bool has_prev_sto_phase);

} // namespace ForwardRiccatiRecursion 
} // namespace sto_ipm
} // namespace ocs2