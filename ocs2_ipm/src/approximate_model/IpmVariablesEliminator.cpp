#include <ocs2_ipm/approximate_model/IpmVariablesEliminator.h>
#include <ocs2_ipm/core/InteriorPointMethod.h>


namespace ocs2 {
namespace ipm {

void eliminateIpmVariablesIntermediateLQ(const IpmVariables& ipmVariables, ModelData& modelData, IpmData& ipmData, scalar_t barrierParam) {
  initIpmData(modelData, ipmData);
  // state inequality constraints
  evalPerturbedResidual(modelData.stateIneqConstraint, ipmVariables.slackDualStateIneqConstraint, 
                        ipmData.dataStateIneqConstraint, modelData.cost, barrierParam, true, false);
  condenseSlackDual(modelData.stateIneqConstraint, ipmVariables.slackDualStateIneqConstraint, 
                    ipmData.dataStateIneqConstraint, modelData.cost, true, false);
  // state-input inequality constraints
  evalPerturbedResidual(modelData.stateInputIneqConstraint, ipmVariables.slackDualStateInputIneqConstraint, 
                        ipmData.dataStateInputIneqConstraint, modelData.cost, barrierParam, true, true);
  condenseSlackDual(modelData.stateInputIneqConstraint, ipmVariables.slackDualStateInputIneqConstraint, 
                    ipmData.dataStateInputIneqConstraint, modelData.cost, true, true);
}

IpmData eliminateIpmVariablesIntermediateLQ(const IpmVariables& ipmVariables, ModelData& modelData, scalar_t barrierParam) {
  IpmData ipmData;
  eliminateIpmVariablesIntermediateLQ(ipmVariables, modelData, ipmData, barrierParam);
  return ipmData;
}

void eliminateIpmVariablesPreJumpLQ(const IpmVariables& ipmVariables, ModelData& modelData, IpmData& ipmData, scalar_t barrierParam) {
  initIpmData(modelData, ipmData);
  // state inequality constraints
  evalPerturbedResidual(modelData.stateIneqConstraint, ipmVariables.slackDualStateIneqConstraint, 
                        ipmData.dataStateIneqConstraint, modelData.cost, barrierParam, true, false);
  condenseSlackDual(modelData.stateIneqConstraint, ipmVariables.slackDualStateIneqConstraint, 
                    ipmData.dataStateIneqConstraint, modelData.cost, true, false);
}

IpmData eliminateIpmVariablesPreJumpLQ(const IpmVariables& ipmVariables, ModelData& modelData, scalar_t barrierParam) {
  IpmData ipmData;
  eliminateIpmVariablesPreJumpLQ(ipmVariables, modelData, ipmData, barrierParam);
  return ipmData;
}

void eliminateIpmVariablesFinalLQ(const IpmVariables& ipmVariables, ModelData& modelData, IpmData& ipmData, scalar_t barrierParam) {
  eliminateIpmVariablesPreJumpLQ(ipmVariables, modelData, ipmData, barrierParam);
}

IpmData eliminateIpmVariablesFinalLQ(const IpmVariables& ipmVariables, ModelData& modelData, scalar_t barrierParam) {
  IpmData ipmData;
  eliminateIpmVariablesFinalLQ(ipmVariables, modelData, ipmData, barrierParam);
  return ipmData;
}

}  // namespace ipm
}  // namespace ocs2
