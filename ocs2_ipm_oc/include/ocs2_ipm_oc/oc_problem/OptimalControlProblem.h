#pragma once

#include <memory>

#include <ocs2_core/PreComputation.h>
#include <ocs2_core/Types.h>
#include <ocs2_core/augmented_lagrangian/StateAugmentedLagrangianCollection.h>
#include <ocs2_core/augmented_lagrangian/StateInputAugmentedLagrangianCollection.h>
#include <ocs2_core/constraint/StateConstraintCollection.h>
#include <ocs2_core/constraint/StateInputConstraintCollection.h>
#include <ocs2_core/cost/StateCostCollection.h>
#include <ocs2_core/cost/StateInputCostCollection.h>
#include <ocs2_core/dynamics/SystemDynamicsBase.h>
#include <ocs2_core/reference/TargetTrajectories.h>
#include <ocs2_oc/oc_problem/OptimalControlProblem.h>
#include <ocs2_sto/cost/StoCostCollection.h>
#include <ocs2_sto/constraint/StoConstraintCollection.h>

namespace ocs2 {
namespace ipm {

/** Optimal Control Problem including Switching Time Optimization definition */
struct OptimalControlProblem {
  /* Cost */
  /** Intermediate cost */
  std::unique_ptr<StateInputCostCollection> costPtr;
  /** Intermediate state-only cost */
  std::unique_ptr<StateCostCollection> stateCostPtr;
  /** Pre-jump cost */
  std::unique_ptr<StateCostCollection> preJumpCostPtr;
  /** Final cost */
  std::unique_ptr<StateCostCollection> finalCostPtr;

  /* Soft constraints */
  /** Intermediate soft constraint penalty */
  std::unique_ptr<StateInputCostCollection> softConstraintPtr;
  /** Intermediate state-only soft constraint penalty */
  std::unique_ptr<StateCostCollection> stateSoftConstraintPtr;
  /** Pre-jump soft constraint penalty */
  std::unique_ptr<StateCostCollection> preJumpSoftConstraintPtr;
  /** Final soft constraint penalty */
  std::unique_ptr<StateCostCollection> finalSoftConstraintPtr;

  /* Inequality constraints */
  /** Intermediate inequality constraints */
  std::unique_ptr<StateInputConstraintCollection> inequalityConstraintPtr;
  /** Intermediate state-only inequality constraints */
  std::unique_ptr<StateConstraintCollection> stateInequalityConstraintPtr;
  /** Pre-jump inequality constraints */
  std::unique_ptr<StateConstraintCollection> preJumpInequalityConstraintPtr;
  /** Final inequality constraints */
  std::unique_ptr<StateConstraintCollection> finalInequalityConstraintPtr;

  /* Equality constraints */
  /** Intermediate equality constraints, full row rank w.r.t. inputs */
  std::unique_ptr<StateInputConstraintCollection> equalityConstraintPtr;
  /** Intermediate state-only equality constraints */
  std::unique_ptr<StateConstraintCollection> stateEqualityConstraintPtr;
  /** Pre-jump equality constraints */
  std::unique_ptr<StateConstraintCollection> preJumpEqualityConstraintPtr;
  /** Final equality constraints */
  std::unique_ptr<StateConstraintCollection> finalEqualityConstraintPtr;

  /* Lagrangians */
  /** Lagrangian for intermediate equality constraints */
  std::unique_ptr<StateInputAugmentedLagrangianCollection> equalityLagrangianPtr;
  /** Lagrangian for intermediate state-only equality constraints */
  std::unique_ptr<StateAugmentedLagrangianCollection> stateEqualityLagrangianPtr;
  /** Lagrangian for intermediate inequality constraints */
  std::unique_ptr<StateInputAugmentedLagrangianCollection> inequalityLagrangianPtr;
  /** Lagrangian for intermediate state-only inequality constraints */
  std::unique_ptr<StateAugmentedLagrangianCollection> stateInequalityLagrangianPtr;
  /** Lagrangian for pre-jump equality constraints */
  std::unique_ptr<StateAugmentedLagrangianCollection> preJumpEqualityLagrangianPtr;
  /** Lagrangian for pre-jump inequality constraints */
  std::unique_ptr<StateAugmentedLagrangianCollection> preJumpInequalityLagrangianPtr;
  /** Lagrangian for final equality constraints */
  std::unique_ptr<StateAugmentedLagrangianCollection> finalEqualityLagrangianPtr;
  /** Lagrangian for final inequality constraints */
  std::unique_ptr<StateAugmentedLagrangianCollection> finalInequalityLagrangianPtr;

  /* Dynamics */
  /** System dynamics pointer */
  std::unique_ptr<SystemDynamicsBase> dynamicsPtr;

  /* STO cost */
  std::unique_ptr<StoCostCollection> stoCostPtr;

  /* STO constraints */
  std::unique_ptr<StoConstraintCollection> stoConstraintPtr;

  /* Misc. */
  /** The pre-computation module */
  std::unique_ptr<PreComputation> preComputationPtr;

  /** The cost desired trajectories (will be substitute by ReferenceManager) */
  const TargetTrajectories* targetTrajectoriesPtr;

  /** Default constructor */
  OptimalControlProblem();

  /** Constructor from ::ocs2::OptimalControlProblem */
  OptimalControlProblem(const ::ocs2::OptimalControlProblem& optimalControlProblem);

  /** Default destructor */
  ~OptimalControlProblem() = default;

  /** Copy constructor */
  OptimalControlProblem(const OptimalControlProblem& other);

  /** Copy assignment */
  OptimalControlProblem& operator=(const OptimalControlProblem& rhs);

  /** Move constructor */
  OptimalControlProblem(OptimalControlProblem&& other) noexcept = default;

  /** Move assignment */
  OptimalControlProblem& operator=(OptimalControlProblem&& rhs) noexcept = default;

  /** Swap */
  void swap(OptimalControlProblem& other) noexcept;

  ::ocs2::OptimalControlProblem toOcs2OptimalControlProblem() const;
};

}  // namespace ipm
}  // namespace ocs2
