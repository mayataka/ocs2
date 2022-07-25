#include <ocs2_stoc/sto/IpmVariables.h>
#include <ocs2_ipm/core/InteriorPointMethod.h>

#include <cassert>

namespace ocs2 {
namespace stoc {

void initIpmVariables(const StoModelData& stoModelData, ipm::SlackDual& ipmVariables, scalar_t barrier) {
  ipm::initSlackDual(stoModelData.stoConstraint, ipmVariables, barrier);
}

ipm::SlackDual initIpmVariables(const StoModelData& stoModelData, scalar_t barrier) {
  ipm::SlackDual ipmVariables;
  ipm::initSlackDual(stoModelData.stoConstraint, ipmVariables, barrier);
  return ipmVariables;
}

void initIpmVariablesDirection(const StoModelData& stoModelData, ipm::SlackDualDirection& ipmVariablesDirection) {
  ipmVariablesDirection.resize(stoModelData.stoConstraint.f.size());
}

void initIpmVariablesDirection(const ipm::SlackDual& ipmVariables, ipm::SlackDualDirection& ipmVariablesDirection) {
  ipmVariablesDirection.resize(ipmVariables.slack.size());
}

ipm::SlackDualDirection initIpmVariablesDirection(const StoModelData& stoModelData) { 
  ipm::SlackDualDirection ipmVariablesDirection;
  initIpmVariablesDirection(stoModelData, ipmVariablesDirection);
  return ipmVariablesDirection;
}

ipm::SlackDualDirection initIpmVariablesDirection(const ipm::SlackDual& ipmVariables) { 
  ipm::SlackDualDirection ipmVariablesDirection;
  initIpmVariablesDirection(ipmVariables, ipmVariablesDirection);
  return ipmVariablesDirection;
}

void updateIpmVariables(ipm::SlackDual& ipmVariables, const ipm::SlackDualDirection& ipmVariablesDirection,
                        scalar_t primalStepSize, scalar_t dualStepSize) {
  if (ipmVariables.slack.size() > 0) {
    assert(ipmVariables.slack.size() == ipmVariablesDirection.slackDirection.size());
    ipmVariables.slack.noalias() += primalStepSize * ipmVariablesDirection.slackDirection; 
    ipmVariables.dual.noalias() += dualStepSize * ipmVariablesDirection.dualDirection; 
  }
}

}  // namespace stoc
}  // namespace ocs2
