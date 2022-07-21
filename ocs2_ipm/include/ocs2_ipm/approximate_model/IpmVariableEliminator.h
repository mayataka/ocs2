#pragma once

#include <ocs2_core/Types.h>
#include <ocs2_ipm/model_data/ModelData.h>
#include <ocs2_ipm/oc_data/IpmVariables.h>
#include <ocs2_ipm/oc_data/IpmData.h>

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
void eliminateIpmVariablesIntermediateLQ(const IpmVariables& ipmVariables, ModelData& modelData, IpmData& ipmData, scalar_t barrierParam);

IpmData eliminateIpmVariablesIntermediateLQ(const IpmVariables& ipmVariables, ModelData& modelData, scalar_t barrierParam);

/**
 * Calculates an LQ approximate of the constrained optimal control problem at a jump event time.
 *
 * @param [in] problem: The optimal control problem
 * @param [in] time: The current time.
 * @param [in] state: The current state.
 * @param [out] modelData: The output data model.
 */
void eliminateIpmVariablesPreJumpLQ(const IpmVariables& ipmVariables, ModelData& modelData, IpmData& ipmData, scalar_t barrierParam);

IpmData eliminateIpmVariablesPreJumpLQ(const IpmVariables& ipmVariables, ModelData& modelData, scalar_t barrierParam);

/**
 * Calculates an LQ approximate of the constrained optimal control problem at final time.
 *
 * @param [in] problem: The optimal control problem
 * @param [in] time: The current time.
 * @param [in] state: The current state.
 * @param [out] modelData: The output data model.
 */
void eliminateIpmVariablesFinalLQ(const IpmVariables& ipmVariables, ModelData& modelData, IpmData& ipmData, scalar_t barrierParam);

IpmData eliminateIpmVariablesFinalLQ(const IpmVariables& ipmVariables, ModelData& modelData, scalar_t barrierParam);

}  // namespace ipm
}  // namespace ocs2
