#pragma once

#include <ocs2_core/Types.h>
#include <ocs2_ipm/model_data/Hamiltonian.h>
#include <ocs2_stoc/riccati_recursion/RiccatiRecursionData.h>
#include <ocs2_stoc/riccati_recursion/LqrPolicy.h>
#include <ocs2_stoc/riccati_recursion/StoPolicy.h>

namespace ocs2 {
namespace stoc {
namespace ForwardRiccatiRecursion {

void computeInput(const LqrPolicy& lqr_policy, const vector_t& dx, vector_t& du, 
                  const scalar_t dts=0.0, const scalar_t dts_next=0.0, 
                  const bool sto=false, const bool has_next_sto_phase=false);

void computeState(const VectorFunctionLinearApproximation& dynamics,
                  const ipm::Hamiltonian& hamiltonian,
                  const vector_t& dx, const vector_t& du, vector_t& dx_next, 
                  const scalar_t dt=0.0, const scalar_t dts_next=0.0, 
                  const bool sto=false);

void computeCostate(const RiccatiRecursionData& riccati, const vector_t& dx, 
                    vector_t& dlmd, const scalar_t dts=0.0, 
                    const scalar_t dts_next=0.0, const bool sto=false, 
                    const bool has_next_sto_phase=false);

double computeSwitchingTime(const StoPolicy& sto_policy, const vector_t& dx, 
                            const double dts_prev, const bool has_prev_sto_phase);

} // namespace ForwardRiccatiRecursion 
} // namespace stoc
} // namespace ocs2