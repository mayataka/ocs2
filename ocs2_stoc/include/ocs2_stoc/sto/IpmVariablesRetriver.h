#pragma once

#include <ocs2_core/Types.h>
#include <ocs2_sto/model_data/StoModelData.h>
#include <ocs2_ipm/oc_data/IpmVariables.h>
#include <ocs2_ipm/oc_data/IpmData.h>

namespace ocs2 {

void retriveStoIpmVariablesDirection(const StoModelData& stoModelData, const ipm::IpmData& ipmData, const ipm::IpmVariables& ipmVariables, 
                                     const vector_t& dts, ipm::IpmVariablesDirection& ipmVariablesDirection);

ipm::IpmVariablesDirection retriveStoIpmVariablesDirection(const StoModelData& stoModelData, const ipm::IpmData& ipmData, 
                                                           const ipm::IpmVariables& ipmVariables, const vector_t& dts);

scalar_t stoPrimalStepSize(const ipm::IpmData& ipmData, const ipm::IpmVariables& ipmVariables, 
                           const ipm::IpmVariablesDirection& ipmVariablesDirection, scalar_t marginRate=0.995);

scalar_t stoDualStepSize(const ipm::IpmData& ipmData, const ipm::IpmVariables& ipmVariables, 
                         const ipm::IpmVariablesDirection& ipmVariablesDirection, scalar_t marginRate=0.995);

}  // namespace ocs2
