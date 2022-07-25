#pragma once

#include <ocs2_core/Types.h>
#include <ocs2_sto/model_data/StoModelData.h>
#include <ocs2_ipm/core/SlackDual.h>
#include <ocs2_ipm/core/InteriorPointMethodData.h>
#include <ocs2_ipm/core/InteriorPointMethod.h>

namespace ocs2 {
namespace stoc {

void retriveStoIpmVariablesDirection(const StoModelData& stoModelData, const ipm::InteriorPointMethodData& ipmData, 
                                     const ipm::SlackDual& ipmVariables, const vector_t& dts, ipm::SlackDualDirection& ipmVariablesDirection);

ipm::SlackDualDirection retriveStoIpmVariablesDirection(const StoModelData& stoModelData, const ipm::InteriorPointMethodData& ipmData, 
                                                        const ipm::SlackDual& ipmVariables, const vector_t& dts);

scalar_t stoPrimalStepSize(const ipm::InteriorPointMethodData& ipmData, const ipm::SlackDual& ipmVariables, 
                           const ipm::SlackDualDirection& ipmVariablesDirection, scalar_t marginRate=0.995);

scalar_t stoDualStepSize(const ipm::InteriorPointMethodData& ipmData, const ipm::SlackDual& ipmVariables, 
                         const ipm::SlackDualDirection& ipmVariablesDirection, scalar_t marginRate=0.995);

}  // namespace stoc
}  // namespace ocs2
