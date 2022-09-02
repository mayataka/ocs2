#pragma once

#include <unordered_map>

#include <ocs2_core/Types.h>
#include <ocs2_ipm/model_data/ModelData.h>
#include <ocs2_stoc/riccati_recursion/RiccatiRecursionData.h>
#include <ocs2_stoc/riccati_recursion/LqrPolicy.h>
#include <ocs2_stoc/riccati_recursion/StoPolicy.h>
#include <ocs2_stoc/riccati_recursion/BackwardRiccatiRecursion.h>
#include <ocs2_stoc/riccati_recursion/ForwardRiccatiRecursion.h>
#include <ocs2_stoc/TimeDiscretization.h>

namespace ocs2 {
namespace stoc {

class RiccatiRecursion {
public:
  RiccatiRecursion(RiccatiSolverMode riccatiSolverMode=RiccatiSolverMode::Robust, scalar_t switchingTimeTrustRegion=0.1,
                   bool enableSwitchingTimeTrustRegion=true);

  void resize(size_t N);

  void backwardRecursion(const std::vector<Grid>& timeDiscretizationGrid, std::vector<ipm::ModelData>& modelData);

  void forwardRecursion(const std::vector<Grid>& timeDiscretizationGrid, const std::vector<ipm::ModelData>& modelData,
                        vector_array_t& stateTrajectory, vector_array_t& inputTrajectory, 
                        vector_array_t& costateTrajectory, scalar_array_t& switchingTimes);

  const std::vector<RiccatiRecursionData>& getRiccatiRecursionData() const { return riccatiData_; }

  const std::vector<LqrPolicy>& getLQRPolicies() const { return lqrPolicy_; }

private:
  std::vector<RiccatiRecursionData> riccatiData_;
  std::vector<LqrPolicy> lqrPolicy_;
  std::vector<StoPolicy> stoPolicy_;
  BackwardRiccatiRecursion backwardRecursion_;
};

} // namespace stoc
} // namespace ocs2