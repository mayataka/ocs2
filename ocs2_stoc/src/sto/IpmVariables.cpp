#include <ocs2_stoc/sto/IpmVariables.h>


namespace ocs2 {

void initIpmVariables(const StoModelData& stoModelData, ipm::IpmVariables& ipmVariables, scalar_t barrier) {
  ipm::initSlackDual(stoModelData.stoConstraint, ipmVariables.slackDualStateIneqConstraint, barrier);
}

ipm::IpmVariables initIpmVariables(const StoModelData& stoModelData, scalar_t barrier) {
  ipm::IpmVariables ipmVariables;
  initIpmVariables(stoModelData, ipmVariables, barrier);
  return ipmVariables;
}

void initIpmVariablesDirection(const StoModelData& stoModelData, ipm::IpmVariablesDirection& ipmVariablesDirection) {
  ipmVariablesDirection.slackDualDirectionStateIneqConstraint.resize(stoModelData.stoConstraint.f.size());
}

ipm::IpmVariablesDirection initIpmVariablesDirection(const StoModelData& stoModelData) { 
  ipm::IpmVariablesDirection ipmVariablesDirection;
  initIpmVariablesDirection(stoModelData, ipmVariablesDirection);
  return ipmVariablesDirection;
}

}  // namespace ocs2
