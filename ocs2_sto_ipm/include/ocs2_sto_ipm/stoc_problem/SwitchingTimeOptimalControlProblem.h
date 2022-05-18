#pragma once

#include <memory>

#include <ocs2_core/PreComputation.h>
#include <ocs2_core/Types.h>
//  #include <ocs2_core/augmented_lagrangian/StateAugmentedLagrangianCollection.h>
//  #include <ocs2_core/augmented_lagrangian/StateInputAugmentedLagrangianCollection.h>
#include <ocs2_core/constraint/StateConstraintCollection.h>
#include <ocs2_core/constraint/StateInputConstraintCollection.h>
#include <ocs2_core/cost/StateCostCollection.h>
#include <ocs2_core/cost/StateInputCostCollection.h>
#include <ocs2_core/dynamics/SystemDynamicsBase.h>
#include <ocs2_core/reference/TargetTrajectories.h>

namespace ocs2 {

/** Optimal Control Problem including Switching Time Optimization definition */
struct SwitchingTimeOptimalControlProblem {
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
  /** Intermediate equality constraints, full row rank w.r.t. inputs */
  std::unique_ptr<StateInputConstraintCollection> inequalityConstraintPtr;
  /** Intermediate state-only equality constraints */
  std::unique_ptr<StateConstraintCollection> stateInequalityConstraintPtr;
  /** Pre-jump equality constraints */
  std::unique_ptr<StateConstraintCollection> preJumpInequalityConstraintPtr;
  /** Final equality constraints */
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
  std::unique_ptr<StateInputCostCollection> equalityLagrangianPtr;
  /** Lagrangian for intermediate state-only equality constraints */
  std::unique_ptr<StateCostCollection> stateEqualityLagrangianPtr;
  /** Lagrangian for intermediate inequality constraints */
  std::unique_ptr<StateInputCostCollection> inequalityLagrangianPtr;
  /** Lagrangian for intermediate state-only inequality constraints */
  std::unique_ptr<StateCostCollection> stateInequalityLagrangianPtr;
  /** Lagrangian for pre-jump equality constraints */
  std::unique_ptr<StateCostCollection> preJumpEqualityLagrangianPtr;
  /** Lagrangian for pre-jump inequality constraints */
  std::unique_ptr<StateCostCollection> preJumpInequalityLagrangianPtr;
  /** Lagrangian for final equality constraints */
  std::unique_ptr<StateCostCollection> finalEqualityLagrangianPtr;
  /** Lagrangian for final inequality constraints */
  std::unique_ptr<StateCostCollection> finalInequalityLagrangianPtr;

  /* Dynamics */
  /** System dynamics pointer */
  std::unique_ptr<SystemDynamicsBase> dynamicsPtr;

  /* Misc. */
  /** The pre-computation module */
  std::unique_ptr<PreComputation> preComputationPtr;

  /** The cost desired trajectories (will be substitute by ReferenceManager) */
  const TargetTrajectories* targetTrajectoriesPtr;

  /** Default constructor */
  SwitchingTimeOptimalControlProblem();

  /** Default destructor */
  ~SwitchingTimeOptimalControlProblem() = default;

  /** Copy constructor */
  SwitchingTimeOptimalControlProblem(const SwitchingTimeOptimalControlProblem& other);

  /** Copy assignment */
  SwitchingTimeOptimalControlProblem& operator=(const SwitchingTimeOptimalControlProblem& rhs);

  /** Move constructor */
  SwitchingTimeOptimalControlProblem(SwitchingTimeOptimalControlProblem&& other) noexcept = default;

  /** Move assignment */
  SwitchingTimeOptimalControlProblem& operator=(SwitchingTimeOptimalControlProblem&& rhs) noexcept = default;

  /** Swap */
  void swap(SwitchingTimeOptimalControlProblem& other) noexcept;
};

}  // namespace ocs2
