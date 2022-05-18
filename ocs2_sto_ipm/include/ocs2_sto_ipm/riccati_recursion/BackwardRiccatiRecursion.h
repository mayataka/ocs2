#pragma once

#include <ocs2_sto_ipm/Types.h>
#include <ocs2_sto_ipm/riccati_recursion/RiccatiRecursionData.h>
#include <ocs2_sto_ipm/riccati_recursion/LQRPolicy.h>
#include <ocs2_sto_ipm/riccati_recursion/STOPolicy.h>

namespace ocs2 {
namespace sto_ipm {

class BackwardRiccatiRecursion {
public:
  BackwardRiccatiRecursion(const size_t nx, const size_t nu, 
                           const scalar_t dts0_max=0.1);

  BackwardRiccatiRecursion();

  void setRegularization(const scalar_t dts0_max);

  void compute(const ScalarFunctionQuadraticApproximationWrapper& cost,
               RiccatiRecursionData& riccati);

  void compute(const RiccatiRecursionData& riccati_next,
               const VectorFunctionLinearApproximationWrapper& dynamics,
               ScalarFunctionQuadraticApproximationWrapper& cost,
               RiccatiRecursionData& riccati, LQRPolicy& lqr_policy);

  void compute(const RiccatiRecursionData& riccati_next,
               const VectorFunctionLinearApproximationWrapper& dynamics,
               ScalarFunctionQuadraticApproximationWrapper& cost,
               RiccatiRecursionData& riccati, LQRPolicy& lqr_policy,
               const bool sto, const bool has_next_sto_phase);

  void phaseTransition(const RiccatiRecursionData& riccati, 
                       RiccatiRecursionData& riccati_pre_event, 
                       STOPolicy& sto_policy, 
                       const bool has_next_sto_phase) const;

private:
  size_t nx_, nu_;
  scalar_t dts0_max_, sgm_eps_;
  matrix_t AtP_, BtP_, GK_, Ginv_3_, Ginv_2_;
  vector_t Pf_;
  scalar_t Ginv_1_;
  Eigen::LLT<matrix_t> llt_;
};

} // namespace sto_ipm
} // namespace ocs2