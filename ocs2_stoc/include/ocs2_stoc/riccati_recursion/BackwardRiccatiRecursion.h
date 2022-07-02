#pragma once

#include <ocs2_core/Types.h>
#include <ocs2_ipm/model_data/Hamiltonian.h>
#include <ocs2_stoc/riccati_recursion/RiccatiRecursionData.h>
#include <ocs2_stoc/riccati_recursion/LqrPolicy.h>
#include <ocs2_stoc/riccati_recursion/StoPolicy.h>

namespace ocs2 {
namespace stoc {

class BackwardRiccatiRecursion {
public:
  BackwardRiccatiRecursion(const size_t nx, const size_t nu, 
                           const scalar_t dts0_max=0.1);

  BackwardRiccatiRecursion();

  void setRegularization(const scalar_t dts0_max);

  void compute(const ScalarFunctionQuadraticApproximation& cost,
               RiccatiRecursionData& riccati);

  void compute(const RiccatiRecursionData& riccati_next,
               const VectorFunctionLinearApproximation& dynamics,
               ScalarFunctionQuadraticApproximation& cost,
               RiccatiRecursionData& riccati, LqrPolicy& lqr_policy);

  void compute(const RiccatiRecursionData& riccati_next,
               const VectorFunctionLinearApproximation& dynamics,
               ScalarFunctionQuadraticApproximation& cost, ipm::Hamiltonian& hamiltonian,
               RiccatiRecursionData& riccati, LqrPolicy& lqr_policy,
               const bool sto, const bool has_next_sto_phase);

  void phaseTransition(const RiccatiRecursionData& riccati, 
                       RiccatiRecursionData& riccati_pre_event, 
                       StoPolicy& sto_policy, 
                       const bool has_next_sto_phase) const;

private:
  size_t nx_, nu_;
  scalar_t dts0_max_, sgm_eps_;
  matrix_t AtP_, BtP_, GK_, Ginv_3_, Ginv_2_;
  vector_t Pf_;
  scalar_t Ginv_1_;
  Eigen::LLT<matrix_t> llt_;
};

} // namespace stoc
} // namespace ocs2