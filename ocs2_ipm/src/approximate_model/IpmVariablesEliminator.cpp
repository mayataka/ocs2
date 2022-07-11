#include <ocs2_ipm/approximate_model/IpmVariableEliminator.h>
#include <ocs2_ipm/core/InteriorPointMethod.h>

namespace ocs2 {
namespace ipm {

void eliminateIpmVariablesIntermediateLQ(const IpmVariables& ipmVariables, ModelData& modelData, IpmData& ipmData) {
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

void eliminateIpmVariablesPreJumpLQ(const IpmVariables& ipmVariables, ModelData& modelData, IpmData& ipmData) {
  // state inequality constraints
  evalPerturbedResidual(modelData.stateIneqConstraint, ipmVariables.slackDualStateIneqConstraint, 
                        ipmData.dataStateIneqConstraint, modelData.cost, true, false);
  condenseSlackDual(modelData.stateIneqConstraint, ipmVariables.slackDualStateIneqConstraint, 
                    ipmData.dataStateIneqConstraint, modelData.cost, true, false);
}

void eliminateIpmVariablesFinalLQ(const IpmVariables& ipmVariables, ModelData& modelData, IpmData& ipmData) {
  eliminateIpmVariablesPreJumpLQ(ipmVariables, modelData, ipmData);
}

}  // namespace ipm
}  // namespace ocs2
