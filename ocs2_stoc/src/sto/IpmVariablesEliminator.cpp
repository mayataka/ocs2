#include <ocs2_stoc/sto/IpmVariablesEliminator.h>
#include <ocs2_ipm/core/InteriorPointMethod.h>


namespace ocs2 {
namespace ipm {

void eliminateIpmVariablesSTO(const ipm::SlackDual& ipmVariables, StoModelData& stoModelData, ipm::InteriorPointMethodData& ipmData, 
                              scalar_t barrierParam) {
  ipm::initInteriorPointMethodData(stoModelData.stoConstraint, ipmData);
  evalPerturbedResidual(stoModelData.stoConstraint, ipmVariables, ipmData, stoModelData.stoCost, barrierParam, true, false);
  condenseSlackDual(stoModelData.stoConstraint, ipmVariables, ipmData, stoModelData.stoCost, true, false);
}

ipm::InteriorPointMethodData eliminateIpmVariablesSTO(const ipm::SlackDual& ipmVariables, StoModelData& stoModelData, 
                                                      scalar_t barrierParam) {
  ipm::InteriorPointMethodData ipmData;
  eliminateIpmVariablesSTO(ipmVariables, stoModelData, ipmData, barrierParam);
  return ipmData;
}

}  // namespace ipm
}  // namespace ocs2
