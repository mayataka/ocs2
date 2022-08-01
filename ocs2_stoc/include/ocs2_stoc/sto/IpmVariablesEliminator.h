#pragma once

#include <ocs2_core/Types.h>
#include <ocs2_sto/model_data/StoModelData.h>
#include <ocs2_ipm/oc_data/IpmVariables.h>
#include <ocs2_ipm/oc_data/IpmData.h>

namespace ocs2 {

/**
 * Calculates an LQ approximate of the constrained optimal control problem at a given time, state, and input.
 *
 * @param [in] problem: The optimal control problem
 * @param [in] time: The current time.
 * @param [in] state: The current state.
 * @param [in] input: The current input.
 * @param [out] modelData: The output data model.
 */
void eliminateIpmVariablesSTO(const ipm::IpmVariables& ipmVariables, StoModelData& stoModelData, ipm::IpmData& ipmData, 
                              scalar_t barrierParam);

ipm::IpmData eliminateIpmVariablesSTO(const ipm::IpmVariables& ipmVariables, StoModelData& stoModelData, scalar_t barrierParam);

}  // namespace ocs2
