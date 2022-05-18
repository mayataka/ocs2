#include <ocs2_sto_ipm/stoc_problem/SwitchingTimeOptimalControlProblem.h>

namespace ocs2 {

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
SwitchingTimeOptimalControlProblem::SwitchingTimeOptimalControlProblem()
    : /* Cost */
      costPtr(new StateInputCostCollection),
      stateCostPtr(new StateCostCollection),
      preJumpCostPtr(new StateCostCollection),
      finalCostPtr(new StateCostCollection),
      /* Soft constraints */
      softConstraintPtr(new StateInputCostCollection),
      stateSoftConstraintPtr(new StateCostCollection),
      preJumpSoftConstraintPtr(new StateCostCollection),
      finalSoftConstraintPtr(new StateCostCollection),
      /* Inequality constraints */
      inequalityConstraintPtr(new StateInputConstraintCollection),
      stateInequalityConstraintPtr(new StateConstraintCollection),
      preJumpInequalityConstraintPtr(new StateConstraintCollection),
      finalInequalityConstraintPtr(new StateConstraintCollection),
      /* Equality constraints */
      equalityConstraintPtr(new StateInputConstraintCollection),
      stateEqualityConstraintPtr(new StateConstraintCollection),
      preJumpEqualityConstraintPtr(new StateConstraintCollection),
      finalEqualityConstraintPtr(new StateConstraintCollection),
      /* Lagrangians */
      equalityLagrangianPtr(new StateInputCostCollection),
      stateEqualityLagrangianPtr(new StateCostCollection),
      inequalityLagrangianPtr(new StateInputCostCollection),
      stateInequalityLagrangianPtr(new StateCostCollection),
      preJumpEqualityLagrangianPtr(new StateCostCollection),
      preJumpInequalityLagrangianPtr(new StateCostCollection),
      finalEqualityLagrangianPtr(new StateCostCollection),
      finalInequalityLagrangianPtr(new StateCostCollection),
      /* Misc. */
      preComputationPtr(new PreComputation),
      targetTrajectoriesPtr(nullptr) {}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
SwitchingTimeOptimalControlProblem::SwitchingTimeOptimalControlProblem(const SwitchingTimeOptimalControlProblem& other)
    : /* Cost */
      costPtr(other.costPtr->clone()),
      stateCostPtr(other.stateCostPtr->clone()),
      preJumpCostPtr(other.preJumpCostPtr->clone()),
      finalCostPtr(other.finalCostPtr->clone()),
      /* Soft constraints */
      softConstraintPtr(other.softConstraintPtr->clone()),
      stateSoftConstraintPtr(other.stateSoftConstraintPtr->clone()),
      preJumpSoftConstraintPtr(other.preJumpSoftConstraintPtr->clone()),
      finalSoftConstraintPtr(other.finalSoftConstraintPtr->clone()),
      /* Inequality constraints */
      inequalityConstraintPtr(other.inequalityConstraintPtr->clone()),
      stateInequalityConstraintPtr(other.stateInequalityConstraintPtr->clone()),
      preJumpInequalityConstraintPtr(other.preJumpInequalityConstraintPtr->clone()),
      finalInequalityConstraintPtr(other.finalInequalityConstraintPtr->clone()),
      /* Equality constraints */
      equalityConstraintPtr(other.equalityConstraintPtr->clone()),
      stateEqualityConstraintPtr(other.stateEqualityConstraintPtr->clone()),
      preJumpEqualityConstraintPtr(other.preJumpEqualityConstraintPtr->clone()),
      finalEqualityConstraintPtr(other.finalEqualityConstraintPtr->clone()),
      /* Lagrangians */
      equalityLagrangianPtr(other.equalityLagrangianPtr->clone()),
      stateEqualityLagrangianPtr(other.stateEqualityLagrangianPtr->clone()),
      inequalityLagrangianPtr(other.inequalityLagrangianPtr->clone()),
      stateInequalityLagrangianPtr(other.stateInequalityLagrangianPtr->clone()),
      preJumpEqualityLagrangianPtr(other.preJumpEqualityLagrangianPtr->clone()),
      preJumpInequalityLagrangianPtr(other.preJumpInequalityLagrangianPtr->clone()),
      finalEqualityLagrangianPtr(other.finalEqualityLagrangianPtr->clone()),
      finalInequalityLagrangianPtr(other.finalInequalityLagrangianPtr->clone()),
      /* Misc. */
      preComputationPtr(other.preComputationPtr->clone()),
      targetTrajectoriesPtr(other.targetTrajectoriesPtr) {
  if (other.dynamicsPtr != nullptr) {
    dynamicsPtr.reset(other.dynamicsPtr->clone());
  }
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
SwitchingTimeOptimalControlProblem& SwitchingTimeOptimalControlProblem::operator=(const SwitchingTimeOptimalControlProblem& rhs) {
  SwitchingTimeOptimalControlProblem tmp(rhs);
  swap(tmp);
  return *this;
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
void SwitchingTimeOptimalControlProblem::swap(SwitchingTimeOptimalControlProblem& other) noexcept {
  /* Cost */
  costPtr.swap(other.costPtr);
  stateCostPtr.swap(other.stateCostPtr);
  preJumpCostPtr.swap(other.preJumpCostPtr);
  finalCostPtr.swap(other.finalCostPtr);

  /* Soft constraints */
  softConstraintPtr.swap(other.softConstraintPtr);
  stateSoftConstraintPtr.swap(other.stateSoftConstraintPtr);
  preJumpSoftConstraintPtr.swap(other.preJumpSoftConstraintPtr);
  finalSoftConstraintPtr.swap(other.finalSoftConstraintPtr);

  /* Inequality constraints */
  inequalityConstraintPtr.swap(other.inequalityConstraintPtr);
  stateInequalityConstraintPtr.swap(other.stateInequalityConstraintPtr);
  preJumpInequalityConstraintPtr.swap(other.preJumpInequalityConstraintPtr);
  finalInequalityConstraintPtr.swap(other.finalInequalityConstraintPtr);

  /* Equality constraints */
  equalityConstraintPtr.swap(other.equalityConstraintPtr);
  stateEqualityConstraintPtr.swap(other.stateEqualityConstraintPtr);
  preJumpEqualityConstraintPtr.swap(other.preJumpEqualityConstraintPtr);
  finalEqualityConstraintPtr.swap(other.finalEqualityConstraintPtr);

  /* Lagrangians */
  equalityLagrangianPtr.swap(other.equalityLagrangianPtr);
  stateEqualityLagrangianPtr.swap(other.stateEqualityLagrangianPtr);
  inequalityLagrangianPtr.swap(other.inequalityLagrangianPtr);
  stateInequalityLagrangianPtr.swap(other.stateInequalityLagrangianPtr);
  preJumpEqualityLagrangianPtr.swap(other.preJumpEqualityLagrangianPtr);
  preJumpInequalityLagrangianPtr.swap(other.preJumpInequalityLagrangianPtr);
  finalEqualityLagrangianPtr.swap(other.finalEqualityLagrangianPtr);
  finalInequalityLagrangianPtr.swap(other.finalInequalityLagrangianPtr);

  /* Dynamics */
  dynamicsPtr.swap(other.dynamicsPtr);

  /* Misc. */
  preComputationPtr.swap(other.preComputationPtr);
  std::swap(targetTrajectoriesPtr, other.targetTrajectoriesPtr);
}

}  // namespace ocs2
