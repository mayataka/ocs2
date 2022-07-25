#include <iostream>

#include <ocs2_stoc/sto/LinearQuadraticApproximator.h>

namespace ocs2 {
namespace stoc {

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
void approximateStoProblem(const ipm::OptimalControlProblem& problem, scalar_t initTime, scalar_t finalTime, 
                           const ModeSchedule& stoModeSchedule, const ModeSchedule& referenceModeSchedule, StoModelData& stoModelData) {
  const auto& preComputation = *problem.preComputationPtr;

  stoModelData.dim = getNumValidSwitchingTimes(initTime, finalTime, referenceModeSchedule);

  // Cost
  stoModelData.stoCost = problem.stoCostPtr->getQuadraticApproximation(initTime, finalTime, stoModeSchedule, referenceModeSchedule, 
                                                                       preComputation);

  // Inequality constraints
  stoModelData.stoConstraint = problem.stoConstraintPtr->getLinearApproximation(initTime, finalTime, stoModeSchedule, referenceModeSchedule, 
                                                                                preComputation);
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
StoModelData approximateStoProblem(const ipm::OptimalControlProblem& problem, scalar_t initTime, scalar_t finalTime, 
                                   const ModeSchedule& stoModeSchedule, const ModeSchedule& referenceModeSchedule) {
  StoModelData stoModelData;
  approximateStoProblem(problem, initTime, finalTime, stoModeSchedule, referenceModeSchedule, stoModelData);
  return stoModelData;   
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
scalar_t computeCost(const ipm::OptimalControlProblem& problem, scalar_t initTime, scalar_t finalTime, 
                     const ModeSchedule& stoModeSchedule, const ModeSchedule& referenceModeSchedule) {
  const auto& preComputation = *problem.preComputationPtr;

  return problem.stoCostPtr->getValue(initTime, finalTime, stoModeSchedule, referenceModeSchedule, preComputation);
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
ScalarFunctionQuadraticApproximation approximateCost(const ipm::OptimalControlProblem& problem, scalar_t initTime, scalar_t finalTime, 
                                                     const ModeSchedule& stoModeSchedule, const ModeSchedule& referenceModeSchedule) {
  const auto& preComputation = *problem.preComputationPtr;

  return problem.stoCostPtr->getQuadraticApproximation(initTime, finalTime, stoModeSchedule, referenceModeSchedule, preComputation);
}

}  // namespace stoc
}  // namespace ocs2
