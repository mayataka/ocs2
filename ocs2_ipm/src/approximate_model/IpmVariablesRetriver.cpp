#include <ocs2_ipm/approximate_model/IpmVariablesRetriver.h>
#include <ocs2_ipm/core/InteriorPointMethod.h>

#include <cassert>

namespace ocs2 {
namespace ipm {

void retriveIntermediateIpmVariablesDirection(const ModelData& modelData, const IpmData& ipmData, const IpmVariables& ipmVariables,
                                              const vector_t& dx, const vector_t& du, IpmVariablesDirection& ipmVariablesDirection) {
  initIpmVariablesDirection(ipmVariables, ipmVariablesDirection);
  // state inequality constraints
  expandSlackDual(modelData.stateIneqConstraint, ipmVariables.slackDualStateIneqConstraint, 
                  ipmData.dataStateIneqConstraint, dx, 
                  ipmVariablesDirection.slackDualDirectionStateIneqConstraint);
  // state-input inequality constraints
  expandSlackDual(modelData.stateInputIneqConstraint, ipmVariables.slackDualStateInputIneqConstraint, 
                  ipmData.dataStateInputIneqConstraint, dx, du, 
                  ipmVariablesDirection.slackDualDirectionStateInputIneqConstraint);
}

IpmVariablesDirection retriveIntermediateIpmVariablesDirection(const ModelData& modelData, const IpmData& ipmData, const IpmVariables& ipmVariables,
                                                               const vector_t& dx, const vector_t& du) {
  IpmVariablesDirection ipmVariablesDirection;
  retriveIntermediateIpmVariablesDirection(modelData, ipmData, ipmVariables, dx, du, ipmVariablesDirection);
  return ipmVariablesDirection;
}

void retrivePreJumpIpmVariablesDirection(const ModelData& modelData, const IpmData& ipmData, const IpmVariables& ipmVariables,
                                         const vector_t& dx, IpmVariablesDirection& ipmVariablesDirection) {
  initIpmVariablesDirection(ipmVariables, ipmVariablesDirection);
  expandSlackDual(modelData.stateIneqConstraint, ipmVariables.slackDualStateIneqConstraint, 
                  ipmData.dataStateIneqConstraint, dx, 
                  ipmVariablesDirection.slackDualDirectionStateIneqConstraint);
}

IpmVariablesDirection retrivePreJumpIpmVariablesDirection(const ModelData& modelData, const IpmData& ipmData, const IpmVariables& ipmVariables,
                                                          const vector_t& dx) {
  IpmVariablesDirection ipmVariablesDirection;
  retrivePreJumpIpmVariablesDirection(modelData, ipmData, ipmVariables, dx, ipmVariablesDirection);
  return ipmVariablesDirection;
}

void retriveFinalIpmVariablesDirection(const ModelData& modelData, const IpmData& ipmData, const IpmVariables& ipmVariables, 
                                       const vector_t& dx, IpmVariablesDirection& ipmVariablesDirection) {
  initIpmVariablesDirection(ipmVariables, ipmVariablesDirection);
  retrivePreJumpIpmVariablesDirection(modelData, ipmData, ipmVariables, dx, ipmVariablesDirection);
}

IpmVariablesDirection retriveFinalIpmVariablesDirection(const ModelData& modelData, const IpmData& ipmData, const IpmVariables& ipmVariables,
                                                        const vector_t& dx) {
  IpmVariablesDirection ipmVariablesDirection;
  retriveFinalIpmVariablesDirection(modelData, ipmData, ipmVariables, dx, ipmVariablesDirection);
  return ipmVariablesDirection;
}

scalar_t intermediatePrimalStepSize(const IpmData& ipmData, const IpmVariables& ipmVariables, const IpmVariablesDirection& ipmVariablesDirection, 
                                    scalar_t marginRate) {
  return std::min(
      fractionToBoundaryPrimalStepSize(ipmVariables.slackDualStateIneqConstraint, 
                                       ipmVariablesDirection.slackDualDirectionStateIneqConstraint,  
                                       ipmData.dataStateIneqConstraint, marginRate),
      fractionToBoundaryPrimalStepSize(ipmVariables.slackDualStateInputIneqConstraint, 
                                       ipmVariablesDirection.slackDualDirectionStateInputIneqConstraint,  
                                       ipmData.dataStateInputIneqConstraint, marginRate));
}

scalar_t intermediateDualStepSize(const IpmData& ipmData, const IpmVariables& ipmVariables, const IpmVariablesDirection& ipmVariablesDirection, 
                                  scalar_t marginRate) {
  assert(marginRate > 0.0);
  assert(marginRate < 1.0);
  return std::min(
      fractionToBoundaryDualStepSize(ipmVariables.slackDualStateIneqConstraint, 
                                     ipmVariablesDirection.slackDualDirectionStateIneqConstraint,  
                                     ipmData.dataStateIneqConstraint, marginRate),
      fractionToBoundaryDualStepSize(ipmVariables.slackDualStateInputIneqConstraint, 
                                     ipmVariablesDirection.slackDualDirectionStateInputIneqConstraint,  
                                     ipmData.dataStateInputIneqConstraint, marginRate));
}

scalar_t preJumpPrimalStepSize(const IpmData& ipmData, const IpmVariables& ipmVariables, const IpmVariablesDirection& ipmVariablesDirection, 
                               scalar_t marginRate) {
  assert(marginRate > 0.0);
  assert(marginRate < 1.0);
  return fractionToBoundaryPrimalStepSize(ipmVariables.slackDualStateIneqConstraint, 
                                          ipmVariablesDirection.slackDualDirectionStateIneqConstraint,  
                                          ipmData.dataStateIneqConstraint, marginRate);
}

scalar_t preJumpDualStepSize(const IpmData& ipmData, const IpmVariables& ipmVariables, const IpmVariablesDirection& ipmVariablesDirection, 
                             scalar_t marginRate) {
  assert(marginRate > 0.0);
  assert(marginRate < 1.0);
  return fractionToBoundaryDualStepSize(ipmVariables.slackDualStateIneqConstraint, 
                                        ipmVariablesDirection.slackDualDirectionStateIneqConstraint,  
                                        ipmData.dataStateIneqConstraint, marginRate);
}

scalar_t finalPrimalStepSize(const IpmData& ipmData, const IpmVariables& ipmVariables, const IpmVariablesDirection& ipmVariablesDirection, 
                             scalar_t marginRate) {
  assert(marginRate > 0.0);
  assert(marginRate < 1.0);
  return preJumpPrimalStepSize(ipmData, ipmVariables, ipmVariablesDirection, marginRate);
}

scalar_t finalDualStepSize(const IpmData& ipmData, const IpmVariables& ipmVariables, const IpmVariablesDirection& ipmVariablesDirection, 
                           scalar_t marginRate) {
  assert(marginRate > 0.0);
  assert(marginRate < 1.0);
  return preJumpDualStepSize(ipmData, ipmVariables, ipmVariablesDirection, marginRate);
}

}  // namespace ipm
}  // namespace ocs2
