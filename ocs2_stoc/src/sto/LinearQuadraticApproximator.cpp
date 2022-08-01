#include <iostream>

#include <ocs2_stoc/sto/LinearQuadraticApproximator.h>

namespace ocs2 {

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
void approximateStoProblem(const ipm::OptimalControlProblem& problem, scalar_t initTime, scalar_t finalTime, 
                           const ModeSchedule& referenceModeSchedule, const ModeSchedule& stoModeSchedule , StoModelData& stoModelData) {
  const auto& preComputation = *problem.preComputationPtr;

  stoModelData.dim = getNumValidSwitchingTimes(initTime, finalTime, referenceModeSchedule);

  // Cost
  if (problem.stoCostPtr) 
    stoModelData.stoCost = problem.stoCostPtr->getQuadraticApproximation(initTime, finalTime, referenceModeSchedule, stoModeSchedule, 
                                                                         preComputation);
  else 
    stoModelData.stoCost = ScalarFunctionQuadraticApproximation(stoModelData.dim, 0);

  // Inequality constraints
  if (problem.stoConstraintPtr) 
    stoModelData.stoConstraint = problem.stoConstraintPtr->getLinearApproximation(initTime, finalTime, referenceModeSchedule, stoModeSchedule, 
                                                                                  preComputation);
  else 
    stoModelData.stoConstraint = VectorFunctionLinearApproximation(0, stoModelData.dim, 0);
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
  if (problem.stoCostPtr)
    return problem.stoCostPtr->getValue(initTime, finalTime, referenceModeSchedule, stoModeSchedule, preComputation);
  else 
    return 0.0;
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
ScalarFunctionQuadraticApproximation approximateCost(const ipm::OptimalControlProblem& problem, scalar_t initTime, scalar_t finalTime, 
                                                     const ModeSchedule& stoModeSchedule, const ModeSchedule& referenceModeSchedule) {
  const auto& preComputation = *problem.preComputationPtr;
  if (problem.stoCostPtr)
    return problem.stoCostPtr->getQuadraticApproximation(initTime, finalTime, referenceModeSchedule, stoModeSchedule, preComputation);
  else 
    ScalarFunctionQuadraticApproximation(getNumValidSwitchingTimes(initTime, finalTime, referenceModeSchedule), 0);
}

}  // namespace ocs2
