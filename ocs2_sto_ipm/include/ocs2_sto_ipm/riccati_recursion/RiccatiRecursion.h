#pragma once

#include <unordered_map>

#include <ocs2_core/Types.h>
#include <ocs2_sto_ipm/riccati_recursion/RiccatiRecursionData.h>
#include <ocs2_sto_ipm/riccati_recursion/LQRPolicy.h>
#include <ocs2_sto_ipm/riccati_recursion/STOPolicy.h>
#include <ocs2_sto_ipm/riccati_recursion/BackwardRiccatiRecursion.h>
#include <ocs2_sto_ipm/riccati_recursion/ForwardRiccatiRecursion.h>
#include <ocs2_sto_ipm/reference/DiscreteTimeModeSchedule.h>

namespace ocs2 {
namespace sto_ipm {

class RiccatiRecursion {
public:
  RiccatiRecursion(const size_t nx, const size_t nu, const size_t N,
                   const size_t maxNumEvents, const scalar_t dts0_max=0.1);

  RiccatiRecursion();

  void backwardRecursion(const DiscreteTimeModeSchedule& modeSchedule,
                         const std::vector<VectorFunctionLinearApproximationWrapper>& dynamics,
                         std::vector<ScalarFunctionQuadraticApproximationWrapper>& cost);

  void forwardRecursion(const DiscreteTimeModeSchedule& modeSchedule, 
                        const std::vector<VectorFunctionLinearApproximationWrapper>& dynamics,
                        vector_array_t& stateTrajectory, vector_array_t& inputTrajectory, 
                        vector_array_t& costateTrajectory, scalar_array_t& switchingTimes);

  const std::vector<LQRPolicy>& getLQRPolicies() const;

private:
  size_t N_;
  std::vector<RiccatiRecursionData> riccati_, riccatiPreEvent_;
  std::vector<LQRPolicy> lqrPolicy_;
  std::vector<STOPolicy> stoPolicy_;
  BackwardRiccatiRecursion backwardRecursion_;
};

} // namespace sto_ipm
} // namespace ocs2