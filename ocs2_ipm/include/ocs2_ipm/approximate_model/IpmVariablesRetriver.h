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

IpmVariablesDirection retriveIntermediateIpmVariablesDirection(const ModelData& modelData, const IpmData& ipmData, const IpmVariables& ipmVariables,
                                                               const vector_t& dx, const vector_t& du);

void retrivePreJumpIpmVariablesDirection(const ModelData& modelData, const IpmData& ipmData, const IpmVariables& ipmVariables,
                                         const vector_t& dx, IpmVariablesDirection& ipmVariablesDirection);

IpmVariablesDirection retrivePreJumpIpmVariablesDirection(const ModelData& modelData, const IpmData& ipmData, const IpmVariables& ipmVariables,
                                                          const vector_t& dx);

void retriveFinalIpmVariablesDirection(const ModelData& modelData, const IpmData& ipmData, const IpmVariables& ipmVariables,
                                       const vector_t& dx, IpmVariablesDirection& ipmVariablesDirection);

IpmVariablesDirection retriveFinalIpmVariablesDirection(const ModelData& modelData, const IpmData& ipmData, const IpmVariables& ipmVariables,
                                                        const vector_t& dx);

scalar_t intermediatePrimalStepSize(const IpmData& ipmData, const IpmVariables& ipmVariables, const IpmVariablesDirection& ipmVariablesDirection, 
                                    scalar_t marginRate=0.995);

scalar_t intermediateDualStepSize(const IpmData& ipmData, const IpmVariables& ipmVariables, const IpmVariablesDirection& ipmVariablesDirection, 
                                  scalar_t marginRate=0.995);

scalar_t preJumpPrimalStepSize(const IpmData& ipmData, const IpmVariables& ipmVariables, const IpmVariablesDirection& ipmVariablesDirection, 
                               scalar_t marginRate=0.995);

scalar_t preJumpDualStepSize(const IpmData& ipmData, const IpmVariables& ipmVariables, const IpmVariablesDirection& ipmVariablesDirection, 
                             scalar_t marginRate=0.995);

scalar_t finalPrimalStepSize(const IpmData& ipmData, const IpmVariables& ipmVariables, const IpmVariablesDirection& ipmVariablesDirection, 
                             scalar_t marginRate=0.995);

scalar_t finalDualStepSize(const IpmData& ipmData, const IpmVariables& ipmVariables, const IpmVariablesDirection& ipmVariablesDirection, 
                           scalar_t marginRate=0.995);

}  // namespace ipm
}  // namespace ocs2
