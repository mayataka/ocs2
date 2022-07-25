#pragma once

#include <ocs2_core/Types.h>
#include <ocs2_sto/model_data/StoModelData.h>
#include <ocs2_ipm/core/SlackDual.h>
#include <ocs2_ipm/core/InteriorPointMethodData.h>
#include <ocs2_ipm/core/InteriorPointMethod.h>

namespace ocs2 {
namespace stoc {

/**
 * Calculates an LQ approximate of the constrained optimal control problem at a given time, state, and input.
 *
 * @param [in] problem: The optimal control problem
 * @param [in] time: The current time.
 * @param [in] state: The current state.
 * @param [in] input: The current input.
 * @param [out] modelData: The output data model.
 */
void eliminateIpmVariablesSTO(const ipm::SlackDual& ipmVariables, StoModelData& stoModelData, ipm::InteriorPointMethodData& ipmData, 
                              scalar_t barrierParam);

ipm::InteriorPointMethodData eliminateIpmVariablesSTO(const ipm::SlackDual& ipmVariables, StoModelData& stoModelData, scalar_t barrierParam);

}  // namespace stoc
}  // namespace ocs2
