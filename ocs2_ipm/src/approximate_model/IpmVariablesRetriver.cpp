#include <ocs2_ipm/approximate_model/IpmVariablesRetriver.h>
#include <ocs2_ipm/core/InteriorPointMethod.h>

namespace ocs2 {
namespace ipm {

void retriveIntermediateIpmVariablesDirection(const ModelData& modelData, const IpmData& ipmData, const IpmVariables& ipmVariables,
                                              const vector_t& dx, const vector_t& du, IpmVariablesDirection& ipmVariablesDirection) {
  // state inequality constraints
  expandSlackDual(modelData.stateIneqConstraint, ipmVariables.slackDualStateIneqConstraint, 
                  ipmData.dataStateIneqConstraint, dx, 
                  ipmVariablesDirection.slackDualDirectionStateIneqConstraint);
  // state-input inequality constraints
  expandSlackDual(modelData.stateInputIneqConstraint, ipmVariables.slackDualStateInputIneqConstraint, 
                  ipmData.dataStateInputIneqConstraint, dx, du, 
                  ipmVariablesDirection.slackDualDirectionStateInputIneqConstraint);
}

void retrivePreJumpIpmVariablesDirection(const ModelData& modelData, const IpmData& ipmData, const IpmVariables& ipmVariables,
                                         const vector_t& dx, IpmVariablesDirection& ipmVariablesDirection) {
  expandSlackDual(modelData.stateIneqConstraint, ipmVariables.slackDualStateIneqConstraint, 
                  ipmData.dataStateIneqConstraint, dx, 
                  ipmVariablesDirection.slackDualDirectionStateIneqConstraint);
}

void retriveFinalIpmVariablesDirection(const ModelData& modelData, const IpmData& ipmData, const vector_t& dx, 
                                       const IpmVariables& ipmVariables, IpmVariablesDirection& ipmVariablesDirection) {
  retrivePreJumpIpmVariablesDirection(modelData, ipmData, ipmVariables, dx, ipmVariablesDirection);
}

scalar_t intermediatePrimalStepSize(const ModelData& modelData, const IpmData& ipmData, const IpmVariables& ipmVariables, 
                                    const IpmVariablesDirection& ipmVariablesDirection, scalar_t marginRate) {
  return std::min(
      fractionToBoundaryPrimalStepSize(ipmVariables.slackDualStateIneqConstraint, 
                                       ipmVariablesDirection.slackDualDirectionStateIneqConstraint,  
                                       ipmData.dataStateIneqConstraint, marginRate),
      fractionToBoundaryPrimalStepSize(ipmVariables.slackDualStateInputIneqConstraint, 
                                       ipmVariablesDirection.slackDualDirectionStateInputIneqConstraint,  
                                       ipmData.dataStateInputIneqConstraint, marginRate));
}

scalar_t intermediateDualStepSize(const ModelData& modelData, const IpmData& ipmData, const IpmVariables& ipmVariables, 
                                  const IpmVariablesDirection& ipmVariablesDirection, scalar_t marginRate) {
  return std::min(
      fractionToBoundaryDualStepSize(ipmVariables.slackDualStateIneqConstraint, 
                                     ipmVariablesDirection.slackDualDirectionStateIneqConstraint,  
                                     ipmData.dataStateIneqConstraint, marginRate),
      fractionToBoundaryDualStepSize(ipmVariables.slackDualStateInputIneqConstraint, 
                                     ipmVariablesDirection.slackDualDirectionStateInputIneqConstraint,  
                                     ipmData.dataStateInputIneqConstraint, marginRate));
}


scalar_t preJumpPrimalStepSize(const ModelData& modelData, const IpmData& ipmData, const IpmVariables& ipmVariables, 
                               const IpmVariablesDirection& ipmVariablesDirection, scalar_t marginRate) {
  return fractionToBoundaryPrimalStepSize(ipmVariables.slackDualStateIneqConstraint, 
                                          ipmVariablesDirection.slackDualDirectionStateIneqConstraint,  
                                          ipmData.dataStateIneqConstraint, marginRate);
}

scalar_t preJumpDualStepSize(const ModelData& modelData, const IpmData& ipmData, const IpmVariables& ipmVariables, 
                             const IpmVariablesDirection& ipmVariablesDirection, scalar_t marginRate) {
  return fractionToBoundaryDualStepSize(ipmVariables.slackDualStateIneqConstraint, 
                                        ipmVariablesDirection.slackDualDirectionStateIneqConstraint,  
                                        ipmData.dataStateIneqConstraint, marginRate);
}

scalar_t finalPrimalStepSize(const ModelData& modelData, const IpmData& ipmData, const IpmVariables& ipmVariables, 
                             const IpmVariablesDirection& ipmVariablesDirection, scalar_t marginRate) {
  return preJumpPrimalStepSize(modelData, ipmData, ipmVariables, ipmVariablesDirection, marginRate);
}

scalar_t finalDualStepSize(const ModelData& modelData, const IpmData& ipmData, const IpmVariables& ipmVariables, 
                           const IpmVariablesDirection& ipmVariablesDirection, scalar_t marginRate) {
  return preJumpDualStepSize(modelData, ipmData, ipmVariables, ipmVariablesDirection, marginRate);
}

}  // namespace ipm
}  // namespace ocs2
