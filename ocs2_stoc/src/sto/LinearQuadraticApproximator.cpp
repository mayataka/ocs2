#include "ocs2_stoc/sto/LinearQuadraticApproximator.h"

#include <iostream>

namespace ocs2 {

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
void approximateStoProblem(const ipm::OptimalControlProblem& problem, scalar_t initTime, scalar_t finalTime, 
                           const ModeSchedule& referenceModeSchedule, const ModeSchedule& stoModeSchedule , StoModelData& stoModelData) {
  const auto& preComputation = *problem.preComputationPtr;

  stoModelData.dim = getNumValidSwitchingTimes(initTime, finalTime, referenceModeSchedule);

  // Cost
  stoModelData.stoCost = problem.stoCostPtr->getQuadraticApproximation(initTime, finalTime, referenceModeSchedule, stoModeSchedule, 
                                                                        preComputation);

  // Inequality constraints
  stoModelData.stoConstraint = problem.stoConstraintPtr->getLinearApproximation(initTime, finalTime, referenceModeSchedule, stoModeSchedule, 
                                                                                preComputation);
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
StoModelData approximateStoProblem(const ipm::OptimalControlProblem& problem, scalar_t initTime, scalar_t finalTime, 
                                   const ModeSchedule& referenceModeSchedule, const ModeSchedule& stoModeSchedule) {
  StoModelData stoModelData;
  approximateStoProblem(problem, initTime, finalTime, referenceModeSchedule, stoModeSchedule, stoModelData);
  return stoModelData;   
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
scalar_t computeCost(const ipm::OptimalControlProblem& problem, scalar_t initTime, scalar_t finalTime, 
                     const ModeSchedule& referenceModeSchedule, const ModeSchedule& stoModeSchedule) {
  const auto& preComputation = *problem.preComputationPtr;
  return problem.stoCostPtr->getValue(initTime, finalTime, referenceModeSchedule, stoModeSchedule, preComputation);
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
ScalarFunctionQuadraticApproximation approximateCost(const ipm::OptimalControlProblem& problem, scalar_t initTime, scalar_t finalTime, 
                                                     const ModeSchedule& stoModeSchedule, const ModeSchedule& referenceModeSchedule) {
  const auto& preComputation = *problem.preComputationPtr;
  return problem.stoCostPtr->getQuadraticApproximation(initTime, finalTime, referenceModeSchedule, stoModeSchedule, preComputation);
}

}  // namespace ocs2
