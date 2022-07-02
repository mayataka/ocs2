#pragma once

#include <ocs2_core/Types.h>
#include <ocs2_ipm/ipm/SlackDual.h>
#include <ocs2_ipm/model_data/ModelData.h>
#include <ocs2_ipm/oc_data/DualVariable.h>

namespace ocs2 {
namespace ipm {

void expandIntermediateDualVariables(const ModelData& modelData, const DualVariable& dualVariable,
                                     const vector_t& dx, const vector_t& du, DualVariableDirection& dualDirection);

void expandPreJumpDualVariables(const ModelData& modelData, const DualVariable& dualVariable,
                                const vector_t& dx, DualVariableDirection& dualDirection);

void expandFinalDualVariables(const ModelData& modelData, const DualVariable& dualVariable,
                              const vector_t& dx, DualVariableDirection& dualDirection);

scalar_t primalStepSizeIntermediateVariables(const ModelData& modelData, const DualVariable& dualVariable, 
                                             const DualVariableDirection& dualDirection);

scalar_t dualStepSizeIntermediateVariables(const ModelData& modelData, const DualVariable& dualVariable, 
                                           const DualVariableDirection& dualDirection);

scalar_t primalStepSizePreJumpVariables(const ModelData& modelData, const DualVariable& dualVariable, 
                                        const DualVariableDirection& dualDirection);

scalar_t dualStepSizePreJumpVariables(const ModelData& modelData, const DualVariable& dualVariable, 
                                      const DualVariableDirection& dualDirection);

scalar_t primalStepSizeFinalVariables(const ModelData& modelData, const DualVariable& dualVariable, 
                                      const DualVariableDirection& dualDirection);

scalar_t dualStepSizeFinalVariables(const ModelData& modelData, const DualVariable& dualVariable, 
                                    const DualVariableDirection& dualDirection);

}  // namespace ipm
}  // namespace ocs2
