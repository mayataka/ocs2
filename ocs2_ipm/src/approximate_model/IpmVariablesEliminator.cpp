#include <ocs2_ipm/approximate_model/IpmVariableEliminator.h>
#include <ocs2_ipm/core/InteriorPointMethod.h>


namespace ocs2 {
namespace ipm {

void eliminateIpmVariablesIntermediateLQ(const IpmVariables& ipmVariables, ModelData& modelData, IpmData& ipmData) {
  initIpmData(modelData, ipmData);
  // state inequality constraints
  evalPerturbedResidual(modelData.stateIneqConstraint, ipmVariables.slackDualStateIneqConstraint, 
                        ipmData.dataStateIneqConstraint, modelData.cost, true, false);
  condenseSlackDual(modelData.stateIneqConstraint, ipmVariables.slackDualStateIneqConstraint, 
                    ipmData.dataStateIneqConstraint, modelData.cost, true, false);
  // state-input inequality constraints
  evalPerturbedResidual(modelData.stateInputIneqConstraint, ipmVariables.slackDualStateInputIneqConstraint, 
                        ipmData.dataStateInputIneqConstraint, modelData.cost, true, true);
  condenseSlackDual(modelData.stateInputIneqConstraint, ipmVariables.slackDualStateInputIneqConstraint, 
                    ipmData.dataStateInputIneqConstraint, modelData.cost, true, true);
}

IpmData eliminateIpmVariablesIntermediateLQ(const IpmVariables& ipmVariables, ModelData& modelData) {
  IpmData ipmData;
  eliminateIpmVariablesIntermediateLQ(ipmVariables, modelData, ipmData);
  return ipmData;
}

void eliminateIpmVariablesPreJumpLQ(const IpmVariables& ipmVariables, ModelData& modelData, IpmData& ipmData) {
  initIpmData(modelData, ipmData);
  // state inequality constraints
  evalPerturbedResidual(modelData.stateIneqConstraint, ipmVariables.slackDualStateIneqConstraint, 
                        ipmData.dataStateIneqConstraint, modelData.cost, true, false);
  condenseSlackDual(modelData.stateIneqConstraint, ipmVariables.slackDualStateIneqConstraint, 
                    ipmData.dataStateIneqConstraint, modelData.cost, true, false);
}

IpmData eliminateIpmVariablesPreJumpLQ(const IpmVariables& ipmVariables, ModelData& modelData) {
  IpmData ipmData;
  eliminateIpmVariablesPreJumpLQ(ipmVariables, modelData, ipmData);
  return ipmData;
}

void eliminateIpmVariablesFinalLQ(const IpmVariables& ipmVariables, ModelData& modelData, IpmData& ipmData) {
  eliminateIpmVariablesPreJumpLQ(ipmVariables, modelData, ipmData);
}

IpmData eliminateIpmVariablesFinalLQ(const IpmVariables& ipmVariables, ModelData& modelData) {
  IpmData ipmData;
  eliminateIpmVariablesFinalLQ(ipmVariables, modelData, ipmData);
  return ipmData;
}

}  // namespace ipm
}  // namespace ocs2
