#pragma once

#include <memory>
#include <string>

#include <ocs2_core/PreComputation.h>
#include <ocs2_core/Types.h>
#include <ocs2_ipm/model_data/ModelData.h>
#include <ocs2_ipm/oc_data/Metrics.h>
#include <ocs2_ipm/oc_data/DualVariable.h>
#include <ocs2_ipm/oc_problem/OptimalControlProblem.h>

namespace ocs2 {
namespace ipm {

/**
 * Calculates an LQ approximate of the constrained optimal control problem at a given time, state, and input.
 *
 * @param [in] problem: The optimal control problem
 * @param [in] time: The current time.
 * @param [in] state: The current state.
 * @param [in] input: The current input.
 * @param [out] modelData: The output data model.
 */
void discretizeIntermediateLQ(const scalar_t dt, const vector_t& state, const vector_t& state_next,
                              const DualVariable& dual, const DualVariable& dual_next, ModelData& modelData);

/**
 * Calculates an LQ approximate of the constrained optimal control problem at a jump event time.
 *
 * @param [in] problem: The optimal control problem
 * @param [in] time: The current time.
 * @param [in] state: The current state.
 * @param [out] modelData: The output data model.
 */
void discretizePreJumpLQ(const vector_t& state, const vector_t& state_next,
                         const DualVariable& dual, const DualVariable& dual_next, ModelData& modelData);

/**
 * Calculates an LQ approximate of the constrained optimal control problem at final time.
 *
 * @param [in] problem: The optimal control problem
 * @param [in] time: The current time.
 * @param [in] state: The current state.
 * @param [out] modelData: The output data model.
 */
void discretizeFinalLQ(const DualVariable& dual, ModelData& modelData);

/**
 * Compute the intermediate-time metrics (i.e. cost, softConstraints, and constraints).
 *
 * @note It is assumed that the precomputation request is already made.
 * problem.preComputationPtr->request(Request::Cost + Request::Constraint + Request::SoftConstraint, t, x, u)
 */
void computeIntermediateDiscretizedMetrics(const scalar_t dt, const vector_t& state, const vector_t& state_next, 
                                           const DualVariable& dual, const ModelData& modelData, Metrics& metrics);

/**
 * Compute the event-time metrics based on pre-jump state value (i.e. cost, softConstraints, and constraints).
 *
 * @note It is assumed that the precomputation request is already made.
 * problem.preComputationPtr->requestPreJump(Request::Cost + Request::Constraint + Request::SoftConstraint, t, x)
 */
void computePreJumpDiscretizedMetrics(const vector_t& state, const vector_t& state_next, 
                                      const DualVariable& dual, const ModelData& modelData, Metrics& metrics);

/**
 * Compute the event-time metrics based on pre-jump state value (i.e. cost, softConstraints, and constraints).
 *
 * @note It is assumed that the precomputation request is already made.
 * problem.preComputationPtr->requestPreJump(Request::Cost + Request::Constraint + Request::SoftConstraint, t, x)
 */
void computeFinalDiscretizedMetrics(const vector_t& state, const DualVariable& dual, const ModelData& modelData, Metrics& metrics);

}  // namespace ipm
}  // namespace ocs2
