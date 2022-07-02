#include <ocs2_ipm/approximate_model/DualVariableReconstructor.h>
#include <ocs2_ipm/ipm/InteriorPointMethod.h>

namespace ocs2 {
namespace ipm {

void expandIntermediateDualVariables(const ModelData& modelData, const DualVariable& dualVariable,
                                     const vector_t& dx, const vector_t& du, DualVariableDirection& dualDirection) {
  expandSlackDual(modelData.stateIneqConstraint, dualVariable.slackDualStateIneqConstraint, 
                  modelData.ipmDataStateIneqConstraint, dx, 
                  dualDirection.slackDualDirectionStateIneqConstraint);
  expandSlackDual(modelData.stateInputIneqConstraint, dualVariable.slackDualStateInputIneqConstraint, 
                  modelData.ipmDataStateInputIneqConstraint, dx, du, 
                  dualDirection.slackDualDirectionStateInputIneqConstraint);
}

void expandPreJumpDualVariables(const ModelData& modelData, const DualVariable& dualVariable,
                                const vector_t& dx, DualVariableDirection& dualDirection) {
  expandSlackDual(modelData.stateIneqConstraint, dualVariable.slackDualStateIneqConstraint, 
                  modelData.ipmDataStateIneqConstraint, dx, 
                  dualDirection.slackDualDirectionStateIneqConstraint);
}

void expandFinalDualVariables(const ModelData& modelData, const vector_t& dx, 
                              const DualVariable& dualVariable, DualVariableDirection& dualDirection) {
  expandPreJumpDualVariables(modelData, dualVariable, dx, dualDirection);
}

scalar_t primalStepSizeIntermediateVariables(const ModelData& modelData, const DualVariable& dualVariable, 
                                             const DualVariableDirection& dualDirection) {
  return std::min(
      fractionToBoundaryPrimalStepSize(dualVariable.slackDualStateIneqConstraint, 
                                       dualDirection.slackDualDirectionStateIneqConstraint,  
                                       modelData.ipmDataStateIneqConstraint),
      fractionToBoundaryPrimalStepSize(dualVariable.slackDualStateInputIneqConstraint, 
                                       dualDirection.slackDualDirectionStateInputIneqConstraint,  
                                       modelData.ipmDataStateInputIneqConstraint));
}

scalar_t dualStepSizeIntermediateVariables(const ModelData& modelData, const DualVariable& dualVariable, 
                                           const DualVariableDirection& dualDirection) {
  return std::min(
      fractionToBoundaryDualStepSize(dualVariable.slackDualStateIneqConstraint, 
                                     dualDirection.slackDualDirectionStateIneqConstraint,  
                                     modelData.ipmDataStateIneqConstraint),
      fractionToBoundaryDualStepSize(dualVariable.slackDualStateInputIneqConstraint, 
                                     dualDirection.slackDualDirectionStateInputIneqConstraint,  
                                     modelData.ipmDataStateInputIneqConstraint));
}


scalar_t primalStepSizePreJumpVariables(const ModelData& modelData, const DualVariable& dualVariable, 
                                        const DualVariableDirection& dualDirection) {
  return fractionToBoundaryPrimalStepSize(dualVariable.slackDualStateIneqConstraint, 
                                          dualDirection.slackDualDirectionStateIneqConstraint,  
                                          modelData.ipmDataStateIneqConstraint);
}

scalar_t dualStepSizePreJumpVariables(const ModelData& modelData, const DualVariable& dualVariable, 
                                      const DualVariableDirection& dualDirection) {
  return fractionToBoundaryDualStepSize(dualVariable.slackDualStateIneqConstraint, 
                                        dualDirection.slackDualDirectionStateIneqConstraint,  
                                        modelData.ipmDataStateIneqConstraint);
}

scalar_t primalStepSizeFinalVariables(const ModelData& modelData, const DualVariable& dualVariable, 
                                      const DualVariableDirection& dualDirection) {
  return primalStepSizePreJumpVariables(modelData, dualVariable, dualDirection);
}

scalar_t dualStepSizeFinalVariables(const ModelData& modelData, const DualVariable& dualVariable, 
                                    const DualVariableDirection& dualDirection) {
  return dualStepSizePreJumpVariables(modelData, dualVariable, dualDirection);
}

}  // namespace ipm
}  // namespace ocs2
