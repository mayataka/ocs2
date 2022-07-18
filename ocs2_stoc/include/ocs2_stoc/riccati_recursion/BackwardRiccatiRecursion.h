#pragma once

#include <ocs2_core/Types.h>
#include <ocs2_ipm/model_data/Hamiltonian.h>
#include <ocs2_stoc/riccati_recursion/RiccatiRecursionData.h>
#include <ocs2_stoc/riccati_recursion/LqrPolicy.h>
#include <ocs2_stoc/riccati_recursion/StoPolicy.h>

#include <Eigen/LU>

namespace ocs2 {
namespace stoc {

class BackwardRiccatiRecursion {
public:
  BackwardRiccatiRecursion(scalar_t dts0_max=0.1);

  void resize(size_t nx, size_t nu);

  void setRegularization(scalar_t dts0_max);

  static void compute(const ScalarFunctionQuadraticApproximation& cost, RiccatiRecursionData& riccati);

  void compute(const RiccatiRecursionData& riccatiNext, const VectorFunctionLinearApproximation& dynamics,
               ScalarFunctionQuadraticApproximation& cost, RiccatiRecursionData& riccati, LqrPolicy& lqrPolicy);

  void compute(const RiccatiRecursionData& riccatiNext, const VectorFunctionLinearApproximation& dynamics,
               ScalarFunctionQuadraticApproximation& cost, ipm::Hamiltonian& hamiltonian, 
               RiccatiRecursionData& riccati, LqrPolicy& lqrPolicy, const bool sto, const bool stoNext);

  void compute(RiccatiRecursionData& riccati, StoPolicy& stoPolicy, const bool stoNext) const;

private:
  scalar_t dts0_max_, sgm_eps_;
  matrix_t AtP_, BtP_, GK_, Ginv_3_, Ginv_2_;
  vector_t Pf_;
  scalar_t Ginv_1_;
  Eigen::LLT<matrix_t> llt_;
  Eigen::LDLT<matrix_t> ldlt_;
};

} // namespace stoc
} // namespace ocs2