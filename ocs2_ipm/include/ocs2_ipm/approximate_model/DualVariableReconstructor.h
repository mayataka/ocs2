#pragma once

#include <ocs2_core/Types.h>
#include <ocs2_ipm/ipm/SlackDual.h>
#include <ocs2_ipm/model_data/ModelData.h>
#include <ocs2_ipm/oc_data/DualVariable.h>

namespace ocs2 {
namespace ipm {

void expandIntermediateDualDirection(const ModelData& modelData, const DualVariable& dualVariable,
                                     const vector_t& dx, const vector_t& du, DualVariableDirection& dualDirection);

void expandPreJumpDualDirection(const ModelData& modelData, const DualVariable& dualVariable,
                                const vector_t& dx, DualVariableDirection& dualDirection);

void expandFinalDualDirection(const ModelData& modelData, const DualVariable& dualVariable,
                              const vector_t& dx, DualVariableDirection& dualDirection);

scalar_t intermediatePrimalStepSize(const ModelData& modelData, const DualVariable& dualVariable, 
                                    const DualVariableDirection& dualDirection);

scalar_t intermediateDualStepSize(const ModelData& modelData, const DualVariable& dualVariable, 
                                  const DualVariableDirection& dualDirection);

scalar_t preJumpPrimalStepSize(const ModelData& modelData, const DualVariable& dualVariable, 
                               const DualVariableDirection& dualDirection);

scalar_t preJumpDualStepSize(const ModelData& modelData, const DualVariable& dualVariable, 
                             const DualVariableDirection& dualDirection);

scalar_t finalPrimalStepSize(const ModelData& modelData, const DualVariable& dualVariable, 
                             const DualVariableDirection& dualDirection);

scalar_t finalDualStepSize(const ModelData& modelData, const DualVariable& dualVariable, 
                           const DualVariableDirection& dualDirection);

}  // namespace ipm
}  // namespace ocs2
