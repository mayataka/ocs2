#include "ocs2_ipm_oc/approximate_model/LinearQuadraticApproximator.h"

#include <ocs2_core/misc/LinearAlgebra.h>

#include <iostream>

namespace ocs2 {
namespace ipm {

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
void approximateIntermediateLQ(OptimalControlProblem& problem, const scalar_t time, const vector_t& state, const vector_t& input,
                               ModelData& modelData, bool enableStateOnlyIneqConstraint) {
  const auto& targetTrajectories = *problem.targetTrajectoriesPtr;
  auto& preComputation = *problem.preComputationPtr;
  constexpr auto request = Request::Cost + Request::SoftConstraint + Request::Constraint + Request::Dynamics + Request::Approximation;
  preComputation.request(request, time, state, input);

  modelData.time = time;
  modelData.stateDim = state.rows();
  modelData.inputDim = input.rows();
  modelData.dynamicsBias.setZero(state.rows());

  // Dynamics
  modelData.dynamicsCovariance = problem.dynamicsPtr->dynamicsCovariance(time, state, input);
  modelData.dynamics = problem.dynamicsPtr->linearApproximation(time, state, input, preComputation);

  // Cost
  modelData.cost = approximateCost(problem, time, state, input);

  // Inequality constraints. State-only inequality constraints should be disabled at the initial time of the horizon.
  if (enableStateOnlyIneqConstraint) {
    modelData.stateIneqConstraint = problem.stateInequalityConstraintPtr->getLinearApproximation(time, state, preComputation);
  }
  modelData.stateInputIneqConstraint = problem.inequalityConstraintPtr->getLinearApproximation(time, state, input, preComputation);

  // Equality constraints
  modelData.stateEqConstraint = problem.stateEqualityConstraintPtr->getLinearApproximation(time, state, preComputation);
  modelData.stateInputEqConstraint = problem.equalityConstraintPtr->getLinearApproximation(time, state, input, preComputation);

  // TODO: add Augmented Lagrangian terms
  // // Lagrangians
  // if (!problem.stateEqualityLagrangianPtr->empty()) {
  //   auto approx = problem.stateEqualityLagrangianPtr->getQuadraticApproximation(time, state, multipliers.stateEq, preComputation);
  //   modelData.cost.f += approx.f;
  //   modelData.cost.dfdx += approx.dfdx;
  //   modelData.cost.dfdxx += approx.dfdxx;
  // }
  // if (!problem.stateInequalityLagrangianPtr->empty()) {
  //   auto approx = problem.stateInequalityLagrangianPtr->getQuadraticApproximation(time, state, multipliers.stateIneq, preComputation);
  //   modelData.cost.f += approx.f;
  //   modelData.cost.dfdx += approx.dfdx;
  //   modelData.cost.dfdxx += approx.dfdxx;
  // }
  // if (!problem.equalityLagrangianPtr->empty()) {
  //   modelData.cost +=
  //       problem.equalityLagrangianPtr->getQuadraticApproximation(time, state, input, multipliers.stateInputEq, preComputation);
  // }
  // if (!problem.inequalityLagrangianPtr->empty()) {
  //   modelData.cost +=
  //       problem.inequalityLagrangianPtr->getQuadraticApproximation(time, state, input, multipliers.stateInputIneq, preComputation);
  // }
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
void approximatePreJumpLQ(OptimalControlProblem& problem, const scalar_t& time, const vector_t& state, ModelData& modelData) {
  const auto& targetTrajectories = *problem.targetTrajectoriesPtr;
  auto& preComputation = *problem.preComputationPtr;
  constexpr auto request = Request::Cost + Request::SoftConstraint + Request::Constraint + Request::Dynamics + Request::Approximation;
  preComputation.requestPreJump(request, time, state);

  modelData.time = time;
  modelData.stateDim = state.rows();
  modelData.inputDim = 0;
  modelData.dynamicsBias.setZero(state.rows());

  // Jump map
  modelData.dynamics = problem.dynamicsPtr->jumpMapLinearApproximation(time, state, preComputation);

  // Pre-jump cost
  modelData.cost = approximateEventCost(problem, time, state);

  // state inequality constraints
  modelData.stateIneqConstraint = problem.preJumpInequalityConstraintPtr->getLinearApproximation(time, state, preComputation);

  // state equality constraint
  modelData.stateEqConstraint = problem.preJumpEqualityConstraintPtr->getLinearApproximation(time, state, preComputation);

  // TODO: add Augmented Lagrangian terms
  // // Lagrangians
  // if (!problem.preJumpEqualityLagrangianPtr->empty()) {
  //   auto approx = problem.preJumpEqualityLagrangianPtr->getQuadraticApproximation(time, state, multipliers.stateEq, preComputation);
  //   modelData.cost.f += approx.f;
  //   modelData.cost.dfdx += approx.dfdx;
  //   modelData.cost.dfdxx += approx.dfdxx;
  // }
  // if (!problem.preJumpInequalityLagrangianPtr->empty()) {
  //   auto approx = problem.preJumpInequalityLagrangianPtr->getQuadraticApproximation(time, state, multipliers.stateIneq, preComputation);
  //   modelData.cost.f += approx.f;
  //   modelData.cost.dfdx += approx.dfdx;
  //   modelData.cost.dfdxx += approx.dfdxx;
  // }
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
void approximateFinalLQ(OptimalControlProblem& problem, const scalar_t& time, const vector_t& state, ModelData& modelData) {
  const auto& targetTrajectories = *problem.targetTrajectoriesPtr;
  auto& preComputation = *problem.preComputationPtr;
  constexpr auto request = Request::Cost + Request::SoftConstraint + Request::Constraint + Request::Approximation;
  preComputation.requestFinal(request, time, state);

  modelData.time = time;
  modelData.stateDim = state.rows();
  modelData.inputDim = 0;
  modelData.dynamicsBias = vector_t();

  // Dynamics
  modelData.dynamics = VectorFunctionLinearApproximation();

  // state inequality constraints
  modelData.stateIneqConstraint = problem.finalInequalityConstraintPtr->getLinearApproximation(time, state, preComputation);

  // state equality constraint
  modelData.stateEqConstraint = problem.finalEqualityConstraintPtr->getLinearApproximation(time, state, preComputation);

  // Final cost
  modelData.cost = approximateFinalCost(problem, time, state);

  // TODO: add Augmented Lagrangian terms
  // // Lagrangians
  // if (!problem.finalEqualityLagrangianPtr->empty()) {
  //   auto approx = problem.finalEqualityLagrangianPtr->getQuadraticApproximation(time, state, multipliers.stateEq, preComputation);
  //   modelData.cost.f += approx.f;
  //   modelData.cost.dfdx += approx.dfdx;
  //   modelData.cost.dfdxx += approx.dfdxx;
  // }
  // if (!problem.finalInequalityLagrangianPtr->empty()) {
  //   auto approx = problem.finalInequalityLagrangianPtr->getQuadraticApproximation(time, state, multipliers.stateIneq, preComputation);
  //   modelData.cost.f += approx.f;
  //   modelData.cost.dfdx += approx.dfdx;
  //   modelData.cost.dfdxx += approx.dfdxx;
  // }
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
scalar_t computeCost(const OptimalControlProblem& problem, const scalar_t& time, const vector_t& state, const vector_t& input) {
  const auto& targetTrajectories = *problem.targetTrajectoriesPtr;
  const auto& preComputation = *problem.preComputationPtr;

  // Compute and sum all costs
  auto cost = problem.costPtr->getValue(time, state, input, targetTrajectories, preComputation);
  cost += problem.softConstraintPtr->getValue(time, state, input, targetTrajectories, preComputation);
  cost += problem.stateCostPtr->getValue(time, state, targetTrajectories, preComputation);
  cost += problem.stateSoftConstraintPtr->getValue(time, state, targetTrajectories, preComputation);

  return cost;
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
ScalarFunctionQuadraticApproximation approximateCost(const OptimalControlProblem& problem, const scalar_t& time, const vector_t& state,
                                                     const vector_t& input) {
  const auto& targetTrajectories = *problem.targetTrajectoriesPtr;
  const auto& preComputation = *problem.preComputationPtr;

  // get the state-input cost approximations
  auto cost = problem.costPtr->getQuadraticApproximation(time, state, input, targetTrajectories, preComputation);

  if (!problem.softConstraintPtr->empty()) {
    cost += problem.softConstraintPtr->getQuadraticApproximation(time, state, input, targetTrajectories, preComputation);
  }

  // get the state only cost approximations
  if (!problem.stateCostPtr->empty()) {
    auto stateCost = problem.stateCostPtr->getQuadraticApproximation(time, state, targetTrajectories, preComputation);
    cost.f += stateCost.f;
    cost.dfdx.noalias() += stateCost.dfdx;
    cost.dfdxx.noalias() += stateCost.dfdxx;
  }

  if (!problem.stateSoftConstraintPtr->empty()) {
    auto stateCost = problem.stateSoftConstraintPtr->getQuadraticApproximation(time, state, targetTrajectories, preComputation);
    cost.f += stateCost.f;
    cost.dfdx.noalias() += stateCost.dfdx;
    cost.dfdxx.noalias() += stateCost.dfdxx;
  }

  return cost;
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
scalar_t computeEventCost(const OptimalControlProblem& problem, const scalar_t& time, const vector_t& state) {
  const auto& targetTrajectories = *problem.targetTrajectoriesPtr;
  const auto& preComputation = *problem.preComputationPtr;

  auto cost = problem.preJumpCostPtr->getValue(time, state, targetTrajectories, preComputation);
  cost += problem.preJumpSoftConstraintPtr->getValue(time, state, targetTrajectories, preComputation);

  return cost;
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
ScalarFunctionQuadraticApproximation approximateEventCost(const OptimalControlProblem& problem, const scalar_t& time,
                                                          const vector_t& state) {
  const auto& targetTrajectories = *problem.targetTrajectoriesPtr;
  const auto& preComputation = *problem.preComputationPtr;

  auto cost = problem.preJumpCostPtr->getQuadraticApproximation(time, state, targetTrajectories, preComputation);
  if (!problem.preJumpSoftConstraintPtr->empty()) {
    cost += problem.preJumpSoftConstraintPtr->getQuadraticApproximation(time, state, targetTrajectories, preComputation);
  }

  return cost;
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
scalar_t computeFinalCost(const OptimalControlProblem& problem, const scalar_t& time, const vector_t& state) {
  const auto& targetTrajectories = *problem.targetTrajectoriesPtr;
  const auto& preComputation = *problem.preComputationPtr;

  auto cost = problem.finalCostPtr->getValue(time, state, targetTrajectories, preComputation);
  cost += problem.finalSoftConstraintPtr->getValue(time, state, targetTrajectories, preComputation);

  return cost;
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
ScalarFunctionQuadraticApproximation approximateFinalCost(const OptimalControlProblem& problem, const scalar_t& time,
                                                          const vector_t& state) {
  const auto& targetTrajectories = *problem.targetTrajectoriesPtr;
  const auto& preComputation = *problem.preComputationPtr;

  auto cost = problem.finalCostPtr->getQuadraticApproximation(time, state, targetTrajectories, preComputation);
  if (!problem.finalSoftConstraintPtr->empty()) {
    cost += problem.finalSoftConstraintPtr->getQuadraticApproximation(time, state, targetTrajectories, preComputation);
  }

  return cost;
}

}  // namespace ipm
}  // namespace ocs2
