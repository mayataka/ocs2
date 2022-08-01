#pragma once

#include <ocs2_core/Types.h>
#include <ocs2_ipm/oc_data/IpmVariables.h>
#include <ocs2_ipm/oc_data/IpmData.h>
#include <ocs2_sto/model_data/StoModelData.h>

namespace ocs2 {

void initIpmVariables(const StoModelData& stoModelData, ipm::IpmVariables& ipmVariables, scalar_t barrier);

ipm::IpmVariables initIpmVariables(const StoModelData& stoModelData, scalar_t barrier);

void initIpmVariablesDirection(const StoModelData& stoModelData, ipm::IpmVariablesDirection& ipmVariablesDirection);

ipm::IpmVariablesDirection initIpmVariablesDirection(const StoModelData& stoModelData);

}  // namespace ocs2
