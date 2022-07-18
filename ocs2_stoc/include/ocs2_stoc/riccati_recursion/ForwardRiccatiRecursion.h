#pragma once

#include <ocs2_core/Types.h>
#include <ocs2_ipm/model_data/Hamiltonian.h>
#include <ocs2_stoc/riccati_recursion/RiccatiRecursionData.h>
#include <ocs2_stoc/riccati_recursion/LqrPolicy.h>
#include <ocs2_stoc/riccati_recursion/StoPolicy.h>

namespace ocs2 {
namespace stoc {
namespace ForwardRiccatiRecursion {

void computeInput(const LqrPolicy& lqrPolicy, const vector_t& dx, vector_t& du, 
                  scalar_t dts=0.0, scalar_t dtsNext=0.0, bool sto=false, bool stoNext=false);

void computeState(const VectorFunctionLinearApproximation& dynamics, const ipm::Hamiltonian& hamiltonian,
                  const vector_t& dx, const vector_t& du, vector_t& dxNext, 
                  scalar_t dt=0.0, scalar_t dtsNext=0.0, bool sto=false);

void computeCostate(const RiccatiRecursionData& riccati, const vector_t& dx, vector_t& dlmd, 
                    scalar_t dts=0.0, scalar_t dtsNext=0.0, bool sto=false, bool stoNext=false);

scalar_t computeSwitchingTime(const StoPolicy& stoPolicy, const vector_t& dx, scalar_t dtsPrev, bool stoPrev);

} // namespace ForwardRiccatiRecursion 
} // namespace stoc
} // namespace ocs2