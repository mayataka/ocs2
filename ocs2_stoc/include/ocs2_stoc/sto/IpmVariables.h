#pragma once

#include <ocs2_core/Types.h>
#include <ocs2_ipm/core/SlackDual.h>
#include <ocs2_ipm/core/SlackDualDirection.h>
#include <ocs2_sto/model_data/StoModelData.h>

namespace ocs2 {
namespace stoc {

void initIpmVariables(const StoModelData& stoModelData, ipm::SlackDual& ipmVariables, scalar_t barrier);

ipm::SlackDual initIpmVariables(const StoModelData& stoModelData, scalar_t barrier);

void initIpmVariablesDirection(const StoModelData& stoModelData, ipm::SlackDualDirection& ipmVariablesDirection);

void initIpmVariablesDirection(const ipm::SlackDual& ipmVariables, ipm::SlackDualDirection& ipmVariablesDirection);

ipm::SlackDualDirection initIpmVariablesDirection(const StoModelData& stoModelData);

ipm::SlackDualDirection initIpmVariablesDirection(const ipm::SlackDual& ipmVariables);

void updateIpmVariables(ipm::SlackDual& ipmVariables, const ipm::SlackDualDirection& ipmVariablesDirection,
                        scalar_t primalStepSize, scalar_t dualStepSize);

}  // namespace stoc
}  // namespace ocs2
