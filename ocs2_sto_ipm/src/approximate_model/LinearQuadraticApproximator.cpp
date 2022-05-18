#include <iostream>

#include <ocs2_core/misc/LinearAlgebra.h>
#include <ocs2_oc/approximate_model/LinearQuadraticApproximator.h>
#include <ocs2_sto_ipm/approximate_model/LinearQuadraticApproximator.h>

namespace ocs2 {
namespace sto_ipm {

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
void approximateIntermediateLQ(SwitchingTimeOptimalControlProblem& problem, const scalar_t time, const vector_t& state, const vector_t& input,
                               StoModelData& modelData) {
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
  modelData.dynamics.base = problem.dynamicsPtr->linearApproximation(time, state, input, preComputation);

  // Cost
  modelData.cost.base = ocs2::approximateCost(problem, time, state, input);

  // Inequality constraints
  modelData.stateIneqConstraint.base = problem.stateInequalityConstraintPtr->getLinearApproximation(time, state, preComputation);
  modelData.stateInputIneqConstraint.base = problem.inequalityConstraintPtr->getLinearApproximation(time, state, input, preComputation);

  // Equality constraints
  modelData.stateEqConstraint.base = problem.stateEqualityConstraintPtr->getLinearApproximation(time, state, preComputation);
  modelData.stateInputEqConstraint.base = problem.equalityConstraintPtr->getLinearApproximation(time, state, input, preComputation);

  // Lagrangians
  if (!problem.stateEqualityLagrangianPtr->empty()) {
    auto approx = problem.stateEqualityLagrangianPtr->getQuadraticApproximation(time, state, targetTrajectories, preComputation);
    modelData.cost.base.f.noalias() += approx.f;
    modelData.cost.base.dfdx.noalias() += approx.dfdx;
    modelData.cost.base.dfdxx.noalias() += approx.dfdxx;
  }
  if (!problem.stateInequalityLagrangianPtr->empty()) {
    auto approx = problem.stateInequalityLagrangianPtr->getQuadraticApproximation(time, state, targetTrajectories, preComputation);
    modelData.cost.base.f.noalias() += approx.f;
    modelData.cost.base.dfdx.noalias() += approx.dfdx;
    modelData.cost.base.dfdxx.noalias() += approx.dfdxx;
  }
  if (!problem.equalityLagrangianPtr->empty()) {
    modelData.cost.base += problem.equalityLagrangianPtr->getQuadraticApproximation(time, state, input, targetTrajectories, preComputation);
  }
  if (!problem.inequalityLagrangianPtr->empty()) {
    modelData.cost.base += problem.inequalityLagrangianPtr->getQuadraticApproximation(time, state, input, targetTrajectories, preComputation);
  }
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
void approximatePreJumpLQ(SwitchingTimeOptimalControlProblem& problem, const scalar_t& time, const vector_t& state, StoModeData& modelData) {
  const auto& targetTrajectories = *problem.targetTrajectoriesPtr;
  auto& preComputation = *problem.preComputationPtr;
  constexpr auto request = Request::Cost + Request::SoftConstraint + Request::Constraint + Request::Dynamics + Request::Approximation;
  preComputation.requestPreJump(request, time, state);

  modelData.time = time;
  modelData.stateDim = state.rows();
  modelData.inputDim = 0;
  modelData.dynamicsBias.setZero(state.rows());

  // Jump map
  modelData.dynamics.base = problem.dynamicsPtr->jumpMapLinearApproximation(time, state, preComputation);

  // Pre-jump cost
  modelData.cost.base = approximateEventCost(problem, time, state);

  // state inequality constraint
  modelData.stateIneqConstraint.base = problem.preJumpInequalityConstraintPtr->getLinearApproximation(time, state, preComputation);

  // state equality constraint
  modelData.stateEqConstraint.base = problem.preJumpEqualityConstraintPtr->getLinearApproximation(time, state, preComputation);

  // Lagrangians
  if (!problem.preJumpEqualityLagrangianPtr->empty()) {
    auto approx = problem.preJumpEqualityLagrangianPtr->getQuadraticApproximation(time, state, targetTrajectories, preComputation);
    modelData.cost.f.noalias() += approx.f;
    modelData.cost.dfdx.noalias() += approx.dfdx;
    modelData.cost.dfdxx.noalias() += approx.dfdxx;
  }
  if (!problem.preJumpInequalityLagrangianPtr->empty()) {
    auto approx = problem.preJumpInequalityLagrangianPtr->getQuadraticApproximation(time, state, targetTrajectories, preComputation);
    modelData.cost.f.noalias() += approx.f;
    modelData.cost.dfdx.noalias() += approx.dfdx;
    modelData.cost.dfdxx.noalias() += approx.dfdxx;
  }
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
void approximateFinalLQ(SwitchingTimeOptimalControlProblem& problem, const scalar_t& time, const vector_t& state, StoModeData& modelData) {
  const auto& targetTrajectories = *problem.targetTrajectoriesPtr;
  auto& preComputation = *problem.preComputationPtr;
  constexpr auto request = Request::Cost + Request::SoftConstraint + Request::Constraint + Request::Approximation;
  preComputation.requestFinal(request, time, state);

  modelData.time = time;
  modelData.stateDim = state.rows();
  modelData.inputDim = 0;
  modelData.dynamicsBias = vector_t();

  // Dynamics
  modelData.dynamics.base = VectorFunctionLinearApproximation();

  // state inequality constraint
  modelData.stateIneqConstraint.base = problem.finalInequalityConstraintPtr->getLinearApproximation(time, state, preComputation);

  // state equality constraint
  modelData.stateEqConstraint.base = problem.finalEqualityConstraintPtr->getLinearApproximation(time, state, preComputation);

  // Final cost
  modelData.cost.base = approximateFinalCost(problem, time, state);

  // Lagrangians
  if (!problem.finalEqualityLagrangianPtr->empty()) {
    auto approx = problem.finalEqualityLagrangianPtr->getQuadraticApproximation(time, state, targetTrajectories, preComputation);
    modelData.cost.base.f.noalias() += approx.f;
    modelData.cost.base.dfdx.noalias() += approx.dfdx;
    modelData.cost.base.dfdxx.noalias() += approx.dfdxx;
  }
  if (!problem.finalInequalityLagrangianPtr->empty()) {
    auto approx = problem.finalInequalityLagrangianPtr->getQuadraticApproximation(time, state, targetTrajectories, preComputation);
    modelData.cost.base.f.noalias() += approx.f;
    modelData.cost.base.dfdx.noalias() += approx.dfdx;
    modelData.cost.base.dfdxx.noalias() += approx.dfdxx;
  }
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
scalar_t computeCost(const SwitchingTimeOptimalControlProblem& problem, const scalar_t& time, const vector_t& state, const vector_t& input) {
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
ScalarFunctionQuadraticApproximation approximateCost(const SwitchingTimeOptimalControlProblem& problem, const scalar_t& time, 
                                                     const vector_t& state, const vector_t& input) {
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
    cost.dfdx += stateCost.dfdx;
    cost.dfdxx += stateCost.dfdxx;
  }

  if (!problem.stateSoftConstraintPtr->empty()) {
    auto stateCost = problem.stateSoftConstraintPtr->getQuadraticApproximation(time, state, targetTrajectories, preComputation);
    cost.f += stateCost.f;
    cost.dfdx += stateCost.dfdx;
    cost.dfdxx += stateCost.dfdxx;
  }

  return cost;
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
scalar_t computeEventCost(const SwitchingTimeOptimalControlProblem& problem, const scalar_t& time, const vector_t& state) {
  const auto& targetTrajectories = *problem.targetTrajectoriesPtr;
  const auto& preComputation = *problem.preComputationPtr;

  auto cost = problem.preJumpCostPtr->getValue(time, state, targetTrajectories, preComputation);
  cost += problem.preJumpSoftConstraintPtr->getValue(time, state, targetTrajectories, preComputation);

  return cost;
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
ScalarFunctionQuadraticApproximation approximateEventCost(const SwitchingTimeOptimalControlProblem& problem, const scalar_t& time,
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
scalar_t computeFinalCost(const SwitchingTimeOptimalControlProblem& problem, const scalar_t& time, const vector_t& state) {
  const auto& targetTrajectories = *problem.targetTrajectoriesPtr;
  const auto& preComputation = *problem.preComputationPtr;

  auto cost = problem.finalCostPtr->getValue(time, state, targetTrajectories, preComputation);
  cost += problem.finalSoftConstraintPtr->getValue(time, state, targetTrajectories, preComputation);

  return cost;
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
ScalarFunctionQuadraticApproximation approximateFinalCost(const SwitchingTimeOptimalControlProblem& problem, 
                                                          const scalar_t& time, const vector_t& state) {
  const auto& targetTrajectories = *problem.targetTrajectoriesPtr;
  const auto& preComputation = *problem.preComputationPtr;

  auto cost = problem.finalCostPtr->getQuadraticApproximation(time, state, targetTrajectories, preComputation);
  if (!problem.finalSoftConstraintPtr->empty()) {
    cost += problem.finalSoftConstraintPtr->getQuadraticApproximation(time, state, targetTrajectories, preComputation);
  }

  return cost;
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
Metrics computeIntermediateMetrics(SwitchingTimeOptimalControlProblem& problem, const scalar_t time, const vector_t& state, const vector_t& input) {
  auto& preComputation = *problem.preComputationPtr;
  const auto& targetTrajectories = *problem.targetTrajectoriesPtr;

  Metrics metrics;

  // Cost
  metrics.cost = computeCost(problem, time, state, input);

  // Equality constraints
  metrics.stateEqConstraint = problem.stateEqualityConstraintPtr->getValue(time, state, preComputation);
  metrics.stateInputEqConstraint = problem.equalityConstraintPtr->getValue(time, state, input, preComputation);

  // Lagrangians
  metrics.stateEqLagrangian = problem.stateEqualityLagrangianPtr->getValue(time, state, targetTrajectories, preComputation);
  metrics.stateIneqLagrangian = problem.stateInequalityLagrangianPtr->getValue(time, state, targetTrajectories, preComputation);
  metrics.stateInputEqLagrangian = problem.equalityLagrangianPtr->getValue(time, state, input, targetTrajectories, preComputation);
  metrics.stateInputIneqLagrangian = problem.inequalityLagrangianPtr->getValue(time, state, input, targetTrajectories, preComputation);

  return metrics;
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
Metrics computePreJumpMetrics(SwitchingTimeOptimalControlProblem& problem, const scalar_t time, const vector_t& state) {
  auto& preComputation = *problem.preComputationPtr;
  const auto& targetTrajectories = *problem.targetTrajectoriesPtr;

  Metrics metrics;

  // Cost
  metrics.cost = computeEventCost(problem, time, state);

  // Equality constraint
  metrics.stateEqConstraint = problem.preJumpEqualityConstraintPtr->getValue(time, state, preComputation);

  // Lagrangians
  metrics.stateEqLagrangian = problem.preJumpEqualityLagrangianPtr->getValue(time, state, targetTrajectories, preComputation);
  metrics.stateIneqLagrangian = problem.preJumpInequalityLagrangianPtr->getValue(time, state, targetTrajectories, preComputation);

  return metrics;
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
Metrics computeFinalMetrics(SwitchingTimeOptimalControlProblem& problem, const scalar_t time, const vector_t& state) {
  auto& preComputation = *problem.preComputationPtr;
  const auto& targetTrajectories = *problem.targetTrajectoriesPtr;

  Metrics metrics;

  // Cost
  metrics.cost = computeFinalCost(problem, time, state);

  // Equality constraint
  metrics.stateEqConstraint = problem.finalEqualityConstraintPtr->getValue(time, state, preComputation);

  // Lagrangians
  metrics.stateEqLagrangian = problem.finalEqualityLagrangianPtr->getValue(time, state, targetTrajectories, preComputation);
  metrics.stateIneqLagrangian = problem.finalInequalityLagrangianPtr->getValue(time, state, targetTrajectories, preComputation);

  return metrics;
}

}  // namespace sto_ipm
}  // namespace ocs2
