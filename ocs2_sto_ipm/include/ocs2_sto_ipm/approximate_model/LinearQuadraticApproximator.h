#pragma once

#include <memory>
#include <string>

#include <ocs2_core/PreComputation.h>
#include <ocs2_core/Types.h>
#include <ocs2_oc/oc_data/Metrics.h>
#include <ocs2_sto_ipm/stoc_problem/SwitchingTimeOptimalControlProblem.h>
#include <ocs2_sto_ipm/model_data/StoModelData.h>

namespace ocs2 {
namespace sto_ipm {

/**
 * Calculates an LQ approximate of the constrained optimal control problem at a given time, state, and input.
 *
 * @param [in] problem: The optimal control problem
 * @param [in] time: The current time.
 * @param [in] state: The current state.
 * @param [in] input: The current input.
 * @param [out] modelData: The output data model.
 */
void approximateIntermediateLQ(SwitchingTimeOptimalControlProblem& problem, const scalar_t time, const vector_t& state, const vector_t& input,
                               ModelData& modelData);

/**
 * Calculates an LQ approximate of the constrained optimal control problem at a given time, state, and input.
 *
 * @param [in] problem: The optimal control problem
 * @param [in] time: The current time.
 * @param [in] state: The current state.
 * @param [in] input: The current input.
 * @return The output data model.
 */
inline ModelData approximateIntermediateLQ(SwitchingTimeOptimalControlProblem& problem, const scalar_t time, const vector_t& state,
                                           const vector_t& input) {
  ModelData md;
  approximateIntermediateLQ(problem, time, state, input, md);
  return md;
}

/**
 * Calculates an LQ approximate of the constrained optimal control problem at a jump event time.
 *
 * @param [in] problem: The optimal control problem
 * @param [in] time: The current time.
 * @param [in] state: The current state.
 * @param [out] modelData: The output data model.
 */
void approximatePreJumpLQ(SwitchingTimeOptimalControlProblem& problem, const scalar_t& time, const vector_t& state, ModelData& modelData);

/**
 * Calculates an LQ approximate of the constrained optimal control problem at a jump event time.
 *
 * @param [in] problem: The optimal control problem
 * @param [in] time: The current time.
 * @param [in] state: The current state.
 * @return The output data model.
 */
inline ModelData approximatePreJumpLQ(SwitchingTimeOptimalControlProblem& problem, const scalar_t& time, const vector_t& state) {
  ModelData md;
  approximatePreJumpLQ(problem, time, state, md);
  return md;
}

/**
 * Calculates an LQ approximate of the constrained optimal control problem at final time.
 *
 * @param [in] problem: The optimal control problem
 * @param [in] time: The current time.
 * @param [in] state: The current state.
 * @param [out] modelData: The output data model.
 */
void approximateFinalLQ(SwitchingTimeOptimalControlProblem& problem, const scalar_t& time, const vector_t& state, ModelData& modelData);

/**
 * Calculates an LQ approximate of the constrained optimal control problem at final time.
 *
 * @param [in] problem: The optimal control problem
 * @param [in] time: The current time.
 * @param [in] state: The current state.
 * @return The output data model.
 */
inline ModelData approximateFinalLQ(SwitchingTimeOptimalControlProblem& problem, const scalar_t& time, const vector_t& state) {
  ModelData md;
  approximateFinalLQ(problem, time, state, md);
  return md;
}

/**
 * Compute the total intermediate cost (i.e. cost + softConstraints). It is assumed that the precomputation request is already made.
 */
scalar_t computeCost(const SwitchingTimeOptimalControlProblem& problem, const scalar_t& time, const vector_t& state, const vector_t& input);

/**
 * Compute the quadratic approximation of the total intermediate cost (i.e. cost + softConstraints). It is assumed that the precomputation
 * request is already made.
 */
ScalarFunctionQuadraticApproximation approximateCost(const SwitchingTimeOptimalControlProblem& problem, const scalar_t& time, const vector_t& state,
                                                     const vector_t& input);

/**
 * Compute the total preJump cost (i.e. cost + softConstraints). It is assumed that the precomputation request is already made.
 */
scalar_t computeEventCost(const SwitchingTimeOptimalControlProblem& problem, const scalar_t& time, const vector_t& state);

/**
 * Compute the quadratic approximation of the total preJump cost (i.e. cost + softConstraints). It is assumed that the precomputation
 * request is already made.
 */
ScalarFunctionQuadraticApproximation approximateEventCost(const SwitchingTimeOptimalControlProblem& problem, const scalar_t& time,
                                                          const vector_t& state);

/**
 * Compute the total final cost (i.e. cost + softConstraints). It is assumed that the precomputation request is already made.
 */
scalar_t computeFinalCost(const SwitchingTimeOptimalControlProblem& problem, const scalar_t& time, const vector_t& state);

/**
 * Compute the quadratic approximation of the total final cost (i.e. cost + softConstraints). It is assumed that the precomputation
 * request is already made.
 */
ScalarFunctionQuadraticApproximation approximateFinalCost(const SwitchingTimeOptimalControlProblem& problem, const scalar_t& time,
                                                          const vector_t& state);

/**
 * Compute the intermediate-time metrics (i.e. cost, softConstraints, and constraints).
 *
 * @note It is assumed that the precomputation request is already made.
 * problem.preComputationPtr->request(Request::Cost + Request::Constraint + Request::SoftConstraint, t, x, u)
 */
Metrics computeIntermediateMetrics(SwitchingTimeOptimalControlProblem& problem, const scalar_t time, const vector_t& state, const vector_t& input);

/**
 * Compute the event-time metrics based on pre-jump state value (i.e. cost, softConstraints, and constraints).
 *
 * @note It is assumed that the precomputation request is already made.
 * problem.preComputationPtr->requestPreJump(Request::Cost + Request::Constraint + Request::SoftConstraint, t, x)
 */
Metrics computePreJumpMetrics(SwitchingTimeOptimalControlProblem& problem, const scalar_t time, const vector_t& state);

/**
 * Compute the final-time metrics (i.e. cost, softConstraints, and constraints).
 *
 * @note It is assumed that the precomputation request is already made.
 * problem.preComputationPtr->requestFinal(Request::Cost + Request::Constraint + Request::SoftConstraint, t, x)
 */
Metrics computeFinalMetrics(SwitchingTimeOptimalControlProblem& problem, const scalar_t time, const vector_t& state);

}  // namespace sto_ipm
}  // namespace ocs2
