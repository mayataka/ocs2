#pragma once

#include <ocs2_core/Types.h>
#include <ocs2_ipm/core/SlackDual.h>
#include <ocs2_ipm/model_data/ModelData.h>
#include <ocs2_ipm/oc_data/IpmVariables.h>
#include <ocs2_ipm/oc_data/IpmData.h>

namespace ocs2 {
namespace ipm {

void retriveIntermediateIpmVariablesDirection(const ModelData& modelData, const IpmData& ipmData, const IpmVariables& ipmVariables,
                                              const vector_t& dx, const vector_t& du, IpmVariablesDirection& ipmVariablesDirection);

void retrivePreJumpIpmVariablesDirection(const ModelData& modelData, const IpmData& ipmData, const IpmVariables& ipmVariables,
                                         const vector_t& dx, IpmVariablesDirection& ipmVariablesDirection);

void retriveFinalIpmVariablesDirection(const ModelData& modelData, const IpmData& ipmData, const IpmVariables& ipmVariables,
                                       const vector_t& dx, IpmVariablesDirection& ipmVariablesDirection);

scalar_t intermediatePrimalStepSize(const ModelData& modelData, const IpmData& ipmData, const IpmVariables& ipmVariables, 
                                    const IpmVariablesDirection& ipmVariablesDirection, scalar_t marginRate=0.995);

scalar_t intermediateDualStepSize(const ModelData& modelData, const IpmData& ipmData, const IpmVariables& ipmVariables, 
                                  const IpmVariablesDirection& ipmVariablesDirection, scalar_t marginRate=0.995);

scalar_t preJumpPrimalStepSize(const ModelData& modelData, const IpmData& ipmData, const IpmVariables& ipmVariables, 
                               const IpmVariablesDirection& ipmVariablesDirection, scalar_t marginRate=0.995);

scalar_t preJumpDualStepSize(const ModelData& modelData, const IpmData& ipmData, const IpmVariables& ipmVariables, 
                             const IpmVariablesDirection& ipmVariablesDirection, scalar_t marginRate=0.995);

scalar_t finalPrimalStepSize(const ModelData& modelData, const IpmData& ipmData, const IpmVariables& ipmVariables, 
                             const IpmVariablesDirection& ipmVariablesDirection, scalar_t marginRate=0.995);

scalar_t finalDualStepSize(const ModelData& modelData, const IpmData& ipmData, const IpmVariables& ipmVariables, 
                           const IpmVariablesDirection& ipmVariablesDirection, scalar_t marginRate=0.995);

}  // namespace ipm
}  // namespace ocs2
