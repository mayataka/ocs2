#include "ocs2_ipm_oc/oc_problem/OptimalControlProblem.h"

namespace ocs2 {
namespace ipm {

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
OptimalControlProblem::OptimalControlProblem()
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
      equalityLagrangianPtr(new StateInputAugmentedLagrangianCollection),
      stateEqualityLagrangianPtr(new StateAugmentedLagrangianCollection),
      inequalityLagrangianPtr(new StateInputAugmentedLagrangianCollection),
      stateInequalityLagrangianPtr(new StateAugmentedLagrangianCollection),
      preJumpEqualityLagrangianPtr(new StateAugmentedLagrangianCollection),
      preJumpInequalityLagrangianPtr(new StateAugmentedLagrangianCollection),
      finalEqualityLagrangianPtr(new StateAugmentedLagrangianCollection),
      finalInequalityLagrangianPtr(new StateAugmentedLagrangianCollection),
      /* STO */
      stoCostPtr(new StoCostCollection),
      stoConstraintPtr(new StoConstraintCollection),
      /* Misc. */
      preComputationPtr(new PreComputation),
      targetTrajectoriesPtr(nullptr) {}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
OptimalControlProblem::OptimalControlProblem(const ::ocs2::OptimalControlProblem& other)
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
      inequalityConstraintPtr(new StateInputConstraintCollection),
      stateInequalityConstraintPtr(new StateConstraintCollection),
      preJumpInequalityConstraintPtr(new StateConstraintCollection),
      finalInequalityConstraintPtr(new StateConstraintCollection),
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
      /* STO */
      stoCostPtr(new StoCostCollection),
      stoConstraintPtr(new StoConstraintCollection),
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
OptimalControlProblem::OptimalControlProblem(const OptimalControlProblem& other)
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
      /* STO */
      stoCostPtr(other.stoCostPtr->clone()),
      stoConstraintPtr(other.stoConstraintPtr->clone()),
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
OptimalControlProblem& OptimalControlProblem::operator=(const OptimalControlProblem& rhs) {
  OptimalControlProblem tmp(rhs);
  swap(tmp);
  return *this;
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
void OptimalControlProblem::swap(OptimalControlProblem& other) noexcept {
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

  /* STO */
  stoCostPtr.swap(other.stoCostPtr);
  stoConstraintPtr.swap(other.stoConstraintPtr);

  /* Misc. */
  preComputationPtr.swap(other.preComputationPtr);
  std::swap(targetTrajectoriesPtr, other.targetTrajectoriesPtr);
}

::ocs2::OptimalControlProblem OptimalControlProblem::toOcs2OptimalControlProblem() const {
  ::ocs2::OptimalControlProblem ocs2OptimalControlProblem;
  /* Cost */
  ocs2OptimalControlProblem.costPtr.reset(costPtr->clone());
  ocs2OptimalControlProblem.stateCostPtr.reset(stateCostPtr->clone());
  ocs2OptimalControlProblem.preJumpCostPtr.reset(preJumpCostPtr->clone());
  ocs2OptimalControlProblem.finalCostPtr.reset(finalCostPtr->clone());
  /* Soft constraints */
  ocs2OptimalControlProblem.softConstraintPtr.reset(softConstraintPtr->clone());
  ocs2OptimalControlProblem.stateSoftConstraintPtr.reset(stateSoftConstraintPtr->clone());
  ocs2OptimalControlProblem.preJumpSoftConstraintPtr.reset(preJumpSoftConstraintPtr->clone());
  ocs2OptimalControlProblem.finalSoftConstraintPtr.reset(finalSoftConstraintPtr->clone());
  // TODO: convert inequality constraints to Augmented Lagrangians
  /* Inequality constraints */
  // inequalityConstraintPtr->clone();
  // stateInequalityConstraintPtr->clone();
  // preJumpInequalityConstraintPtr->clone();
  // finalInequalityConstraintPtr->clone();
  /* Equality constraints */
  ocs2OptimalControlProblem.equalityConstraintPtr.reset(equalityConstraintPtr->clone());
  ocs2OptimalControlProblem.stateEqualityConstraintPtr.reset(stateEqualityConstraintPtr->clone());
  ocs2OptimalControlProblem.preJumpEqualityConstraintPtr.reset(preJumpEqualityConstraintPtr->clone());
  ocs2OptimalControlProblem.finalEqualityConstraintPtr.reset(finalEqualityConstraintPtr->clone());
  /* Lagrangians */
  ocs2OptimalControlProblem.equalityLagrangianPtr.reset(equalityLagrangianPtr->clone());
  ocs2OptimalControlProblem.stateEqualityLagrangianPtr.reset(stateEqualityLagrangianPtr->clone());
  ocs2OptimalControlProblem.inequalityLagrangianPtr.reset(inequalityLagrangianPtr->clone());
  ocs2OptimalControlProblem.stateInequalityLagrangianPtr.reset(stateInequalityLagrangianPtr->clone());
  ocs2OptimalControlProblem.preJumpEqualityLagrangianPtr.reset(preJumpEqualityLagrangianPtr->clone());
  ocs2OptimalControlProblem.preJumpInequalityLagrangianPtr.reset(preJumpInequalityLagrangianPtr->clone());
  ocs2OptimalControlProblem.finalEqualityLagrangianPtr.reset(finalEqualityLagrangianPtr->clone());
  ocs2OptimalControlProblem.finalInequalityLagrangianPtr.reset(finalInequalityLagrangianPtr->clone());
  /** System dynamics pointer */
  if (dynamicsPtr != nullptr) {
    ocs2OptimalControlProblem.dynamicsPtr.reset(dynamicsPtr->clone());
  }
  /* Misc. */
  ocs2OptimalControlProblem.preComputationPtr.reset(preComputationPtr->clone());
  ocs2OptimalControlProblem.targetTrajectoriesPtr = targetTrajectoriesPtr;
  return ocs2OptimalControlProblem;
}

}  // namespace ipm
}  // namespace ocs2
