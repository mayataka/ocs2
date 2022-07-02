#include <ocs2_ipm/approximate_model/DualVariableReconstructor.h>
#include <ocs2_ipm/ipm/InteriorPointMethod.h>

namespace ocs2 {
namespace ipm {

void expandIntermediateDualDirection(const ModelData& modelData, const DualVariable& dualVariable,
                                     const vector_t& dx, const vector_t& du, DualVariableDirection& dualDirection) {
  expandSlackDual(modelData.stateIneqConstraint, dualVariable.slackDualStateIneqConstraint, 
                  modelData.ipmDataStateIneqConstraint, dx, 
                  dualDirection.slackDualDirectionStateIneqConstraint);
  expandSlackDual(modelData.stateInputIneqConstraint, dualVariable.slackDualStateInputIneqConstraint, 
                  modelData.ipmDataStateInputIneqConstraint, dx, du, 
                  dualDirection.slackDualDirectionStateInputIneqConstraint);
}

void expandPreJumpDualDirection(const ModelData& modelData, const DualVariable& dualVariable,
                                const vector_t& dx, DualVariableDirection& dualDirection) {
  expandSlackDual(modelData.stateIneqConstraint, dualVariable.slackDualStateIneqConstraint, 
                  modelData.ipmDataStateIneqConstraint, dx, 
                  dualDirection.slackDualDirectionStateIneqConstraint);
}

void expandFinalDualDirection(const ModelData& modelData, const vector_t& dx, 
                              const DualVariable& dualVariable, DualVariableDirection& dualDirection) {
  expandPreJumpDualDirection(modelData, dualVariable, dx, dualDirection);
}

scalar_t intermediatePrimalStepSize(const ModelData& modelData, const DualVariable& dualVariable, 
                                    const DualVariableDirection& dualDirection) {
  return std::min(
      fractionToBoundaryPrimalStepSize(dualVariable.slackDualStateIneqConstraint, 
                                       dualDirection.slackDualDirectionStateIneqConstraint,  
                                       modelData.ipmDataStateIneqConstraint),
      fractionToBoundaryPrimalStepSize(dualVariable.slackDualStateInputIneqConstraint, 
                                       dualDirection.slackDualDirectionStateInputIneqConstraint,  
                                       modelData.ipmDataStateInputIneqConstraint));
}

scalar_t intermediateDualStepSize(const ModelData& modelData, const DualVariable& dualVariable, 
                                  const DualVariableDirection& dualDirection) {
  return std::min(
      fractionToBoundaryDualStepSize(dualVariable.slackDualStateIneqConstraint, 
                                     dualDirection.slackDualDirectionStateIneqConstraint,  
                                     modelData.ipmDataStateIneqConstraint),
      fractionToBoundaryDualStepSize(dualVariable.slackDualStateInputIneqConstraint, 
                                     dualDirection.slackDualDirectionStateInputIneqConstraint,  
                                     modelData.ipmDataStateInputIneqConstraint));
}


scalar_t preJumpPrimalStepSize(const ModelData& modelData, const DualVariable& dualVariable, 
                               const DualVariableDirection& dualDirection) {
  return fractionToBoundaryPrimalStepSize(dualVariable.slackDualStateIneqConstraint, 
                                          dualDirection.slackDualDirectionStateIneqConstraint,  
                                          modelData.ipmDataStateIneqConstraint);
}

scalar_t preJumpDualStepSize(const ModelData& modelData, const DualVariable& dualVariable, 
                             const DualVariableDirection& dualDirection) {
  return fractionToBoundaryDualStepSize(dualVariable.slackDualStateIneqConstraint, 
                                        dualDirection.slackDualDirectionStateIneqConstraint,  
                                        modelData.ipmDataStateIneqConstraint);
}

scalar_t finalPrimalStepSize(const ModelData& modelData, const DualVariable& dualVariable, 
                             const DualVariableDirection& dualDirection) {
  return preJumpPrimalStepSize(modelData, dualVariable, dualDirection);
}

scalar_t finalDualStepSize(const ModelData& modelData, const DualVariable& dualVariable, 
                           const DualVariableDirection& dualDirection) {
  return preJumpDualStepSize(modelData, dualVariable, dualDirection);
}

}  // namespace ipm
}  // namespace ocs2
