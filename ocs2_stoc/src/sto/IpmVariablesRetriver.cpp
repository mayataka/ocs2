#include <ocs2_stoc/sto/IpmVariablesRetriver.h>
#include <ocs2_ipm/core/InteriorPointMethod.h>

#include <cassert>

namespace ocs2 {
namespace stoc {

void retriveStoIpmVariablesDirection(const StoModelData& stoModelData, const ipm::InteriorPointMethodData& ipmData, 
                                     const ipm::SlackDual& ipmVariables, const vector_t& dts, ipm::SlackDualDirection& ipmVariablesDirection) {
  ipmVariablesDirection.resize(stoModelData.stoConstraint.f.size());
  ipm::expandSlackDual(stoModelData.stoConstraint, ipmVariables, ipmData, dts, ipmVariablesDirection);
}

ipm::SlackDualDirection retriveStoIpmVariablesDirection(const StoModelData& stoModelData, const ipm::InteriorPointMethodData& ipmData, 
                                                        const ipm::SlackDual& ipmVariables, const vector_t& dts) {
  ipm::SlackDualDirection ipmVariablesDirection;
  retriveStoIpmVariablesDirection(stoModelData, ipmData, ipmVariables, dts, ipmVariablesDirection);
  return ipmVariablesDirection;
}

scalar_t stoPrimalStepSize(const ipm::InteriorPointMethodData& ipmData, const ipm::SlackDual& ipmVariables, 
                           const ipm::SlackDualDirection& ipmVariablesDirection, scalar_t marginRate) {
  assert(marginRate > 0.0);
  assert(marginRate < 1.0);
  return ipm::fractionToBoundaryPrimalStepSize(ipmVariables, ipmVariablesDirection, ipmData, marginRate);
}

scalar_t stoDualStepSize(const ipm::InteriorPointMethodData& ipmData, const ipm::SlackDual& ipmVariables, 
                         const ipm::SlackDualDirection& ipmVariablesDirection, scalar_t marginRate) {
  assert(marginRate > 0.0);
  assert(marginRate < 1.0);
  return ipm::fractionToBoundaryDualStepSize(ipmVariables, ipmVariablesDirection, ipmData, marginRate);
}

}  // namespace stoc
}  // namespace ocs2
