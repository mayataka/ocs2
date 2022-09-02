#include "ocs2_stoc/sto/IpmVariablesRetriver.h"

#include <ocs2_ipm/core/InteriorPointMethod.h>

#include <cassert>
#include <iostream>

namespace ocs2 {

void retriveStoIpmVariablesDirection(const StoModelData& stoModelData, const ipm::IpmData& ipmData, const ipm::IpmVariables& ipmVariables, 
                                     const vector_t& dts, ipm::IpmVariablesDirection& ipmVariablesDirection) {
  ipm::initIpmVariablesDirection(ipmVariables, ipmVariablesDirection);
  ipm::expandSlackDual(stoModelData.stoConstraint, ipmVariables.slackDualStateIneqConstraint, ipmData.dataStateIneqConstraint, dts, 
                       ipmVariablesDirection.slackDualDirectionStateIneqConstraint);
}

ipm::IpmVariablesDirection retriveStoIpmVariablesDirection(const StoModelData& stoModelData, const ipm::IpmData& ipmData, 
                                                           const ipm::IpmVariables& ipmVariables, const vector_t& dts) {
  ipm::IpmVariablesDirection ipmVariablesDirection;
  retriveStoIpmVariablesDirection(stoModelData, ipmData, ipmVariables, dts, ipmVariablesDirection);
  return ipmVariablesDirection;
}

scalar_t stoPrimalStepSize(const ipm::IpmData& ipmData, const ipm::IpmVariables& ipmVariables, 
                           const ipm::IpmVariablesDirection& ipmVariablesDirection, scalar_t marginRate) {
  assert(marginRate > 0.0);
  assert(marginRate < 1.0);
  return ipm::fractionToBoundaryPrimalStepSize(ipmVariables.slackDualStateIneqConstraint, 
                                               ipmVariablesDirection.slackDualDirectionStateIneqConstraint,  
                                               ipmData.dataStateIneqConstraint, marginRate);
}

scalar_t stoDualStepSize(const ipm::IpmData& ipmData, const ipm::IpmVariables& ipmVariables, 
                         const ipm::IpmVariablesDirection& ipmVariablesDirection, scalar_t marginRate) {
  assert(marginRate > 0.0);
  assert(marginRate < 1.0);
  return ipm::fractionToBoundaryDualStepSize(ipmVariables.slackDualStateIneqConstraint, 
                                             ipmVariablesDirection.slackDualDirectionStateIneqConstraint,  
                                             ipmData.dataStateIneqConstraint, marginRate);
}

}  // namespace ocs2
