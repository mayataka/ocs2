#pragma once

#include <ocs2_core/Types.h>
#include <ocs2_ipm/model_data/ModelData.h>
#include <ocs2_stoc/riccati_recursion/RiccatiRecursionData.h>
#include <ocs2_stoc/riccati_recursion/LqrPolicy.h>
#include <ocs2_stoc/riccati_recursion/StoPolicy.h>
#include <ocs2_stoc/riccati_recursion/RiccatiSolverMode.h>

#include <Eigen/LU>

namespace ocs2 {
namespace stoc {

class BackwardRiccatiRecursion {
public:
  BackwardRiccatiRecursion(RiccatiSolverMode riccatiSolverMode=RiccatiSolverMode::Robust, scalar_t switchingTimeTrustRegionRadius=0.1,
                           bool enableSwitchingTimeTrustRegion=true);

  void resize(size_t nx, size_t nu);

  void setRegularization(scalar_t switchingTimeTrustRegion, bool enableSwitchingTimeTrustRegion=true);

  static void computeFinal(const ipm::ModelData& modelData, RiccatiRecursionData& riccati);

  void computeIntermediate(const RiccatiRecursionData& riccatiNext, ipm::ModelData& modelData, 
                           RiccatiRecursionData& riccati, LqrPolicy& lqrPolicy, const bool sto, const bool stoNext);

  void computePreJump(const RiccatiRecursionData& riccatiNext, ipm::ModelData& modelData, 
                      RiccatiRecursionData& riccati, LqrPolicy& lqrPolicy, const bool sto, const bool stoNext);

  void computePostJump(const RiccatiRecursionData& riccatiNext, ipm::ModelData& modelData, 
                       RiccatiRecursionData& riccati, LqrPolicy& lqrPolicy, StoPolicy& stoPolicy, 
                       const bool sto, const bool stoNext);

  void modifyPostJump(RiccatiRecursionData& riccati, StoPolicy& stoPolicy, bool computeStoPolicy) const;

private:
  RiccatiSolverMode riccatiSolverMode_;
  scalar_t switchingTimeTrustRegionRadius_, minQuadraticCoeff_;
  bool enableSwitchingTimeTrustRegion_;
  matrix_t AtP_, BtP_, GK_, Ginv_4_, Ginv_3_, Ginv_2_;
  vector_t Pf_;
  scalar_t Ginv_1_;
  Eigen::LLT<matrix_t> llt_;
  Eigen::LDLT<matrix_t> ldlt_;

  void computeIntermediate(const RiccatiRecursionData& riccatiNext, ipm::ModelData& modelData, 
                           RiccatiRecursionData& riccati, LqrPolicy& lqrPolicy);
};

} // namespace stoc
} // namespace ocs2