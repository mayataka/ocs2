#pragma once

#include <unordered_map>

#include <ocs2_core/Types.h>
#include <ocs2_ipm/model_data/ModelData.h>
#include <ocs2_stoc/riccati_recursion/RiccatiRecursionData.h>
#include <ocs2_stoc/riccati_recursion/LqrPolicy.h>
#include <ocs2_stoc/riccati_recursion/StoPolicy.h>
#include <ocs2_stoc/riccati_recursion/BackwardRiccatiRecursion.h>
#include <ocs2_stoc/riccati_recursion/ForwardRiccatiRecursion.h>
#include <ocs2_stoc/reference/DiscreteTimeModeSchedule.h>

namespace ocs2 {
namespace stoc {

class RiccatiRecursion {
public:
  RiccatiRecursion(const size_t nx, const size_t nu, const size_t N,
                   const size_t maxNumEvents, const scalar_t dts0_max=0.1);

  RiccatiRecursion();

  void backwardRecursion(const DiscreteTimeModeSchedule& modeSchedule, std::vector<ipm::ModelData>& modelData);

  void forwardRecursion(const DiscreteTimeModeSchedule& modeSchedule, const std::vector<ipm::ModelData>& modelData,
                        vector_array_t& stateTrajectory, vector_array_t& inputTrajectory, 
                        vector_array_t& costateTrajectory, scalar_array_t& switchingTimes);

  const std::vector<RiccatiRecursionData>& getRiccatiRecursionData() const { return riccati_; }

  const std::vector<LqrPolicy>& getLQRPolicies() const { return lqrPolicy_; }

private:
  size_t N_;
  std::vector<RiccatiRecursionData> riccati_, riccatiPreEvent_;
  std::vector<LqrPolicy> lqrPolicy_;
  std::vector<StoPolicy> stoPolicy_;
  BackwardRiccatiRecursion backwardRecursion_;
};

} // namespace stoc
} // namespace ocs2