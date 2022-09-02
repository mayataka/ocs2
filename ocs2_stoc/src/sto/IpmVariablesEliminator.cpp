#include "ocs2_stoc/sto/IpmVariablesEliminator.h"

#include <ocs2_ipm/core/InteriorPointMethod.h>

namespace ocs2 {

void eliminateIpmVariablesSTO(const ipm::IpmVariables& ipmVariables, StoModelData& stoModelData, ipm::IpmData& ipmData, 
                              scalar_t barrierParam) {
  ipm::initInteriorPointMethodData(stoModelData.stoConstraint, ipmData.dataStateIneqConstraint);
  ipm::evalPerturbedResidual(stoModelData.stoConstraint, ipmVariables.slackDualStateIneqConstraint, ipmData.dataStateIneqConstraint, 
                             stoModelData.stoCost, barrierParam, true, false);
  ipm::condenseSlackDual(stoModelData.stoConstraint, ipmVariables.slackDualStateIneqConstraint, ipmData.dataStateIneqConstraint, 
                         stoModelData.stoCost, true, false);
}

ipm::IpmData eliminateIpmVariablesSTO(const ipm::IpmVariables& ipmVariables, StoModelData& stoModelData, scalar_t barrierParam) {
  ipm::IpmData ipmData;
  eliminateIpmVariablesSTO(ipmVariables, stoModelData, ipmData, barrierParam);
  return ipmData;
}

}  // namespace ocs2
