#pragma once

#include <ocs2_core/PreComputation.h>
#include <ocs2_core/Types.h>
#include <ocs2_ipm/model_data/ModelData.h>

namespace ocs2 {
namespace ipm {

/**
 * Discretizes a continuous-time LQ approximation of the constrained optimal control problem into a discrete-time counterpart. It is assumed
 * that the continuous-time LQ approximateion is already computed.
 *
 * @param [in] dt: The time step of the discretization.
 * @param [in] state: The current state.
 * @param [in] stateNext: The next state.
 * @param [in] costate: The current costate.
 * @param [in] costateNext: The next costate.
 * @param [in, out] modelData: The input continous-time data model. Also an output as discrete-time data model.
 */
void discretizeIntermediateLQ(const scalar_t dt, const vector_t& state, const vector_t& stateNext, const vector_t& costate, 
                              const vector_t& costateNext, ModelData& modelData);

/**
 * Discretizes a continuous-time LQ approximation of the constrained optimal control problem at a jump event time into a discrete-time 
 * counterpart. It is assumed that the continuous-time LQ approximateion is already computed.
 *
 * @param [in] state: The current state.
 * @param [in] stateNext: The next state.
 * @param [in] costate: The current costate.
 * @param [in] costateNext: The next costate.
 * @param [in, out] modelData: The input continous-time data model. Also an output as discrete-time data model.
 */
void discretizePreJumpLQ(const vector_t& state, const vector_t& stateNext, const vector_t& costate, 
                         const vector_t& costateNext, ModelData& modelData);

/**
 * Discretizes a continuous-time LQ approximation of the constrained optimal control problem at a final time into a discrete-time 
 * counterpart. It is assumed that the continuous-time LQ approximateion is already computed.
 *
 * @param [in] costate: The current costate.
 * @param [in, out] modelData: The input continous-time data model. Also an output as discrete-time data model.
 */
void discretizeFinalLQ(const vector_t& costate, ModelData& modelData);

}  // namespace ipm
}  // namespace ocs2
