#pragma once

#include <ocs2_core/Types.h>
#include <ocs2_ipm/core/SlackDual.h>
#include <ocs2_ipm/model_data/ModelData.h>

#include "ocs2_ipm_oc/oc_data/IpmVariables.h"
#include "ocs2_ipm_oc/oc_data/IpmData.h"

namespace ocs2 {
namespace ipm {

/**
 * Retrives the Newton direction of IPM-related variables from the directions of primal variables.
 *
 * @param [in] modelData: The discrete-time LQ approximation model data.
 * @param [in] ipmData: The IPM-related data.
 * @param [in] ipmVariables: The IPM-related variables (a collection of slack and dual variables).
 * @param [in] dx: The state direction.
 * @param [in] du: The input direction.
 * @param [out] ipmVariablesDirection: The direction of the IPM-related variables.
 */
void retriveIntermediateIpmVariablesDirection(const ModelData& modelData, const IpmData& ipmData, const IpmVariables& ipmVariables,
                                              const vector_t& dx, const vector_t& du, IpmVariablesDirection& ipmVariablesDirection);

/**
 * Retrives the Newton direction of IPM-related variables from the directions of primal variables.
 *
 * @param [in] modelData: The discrete-time LQ approximation model data.
 * @param [in] ipmData: The IPM-related data.
 * @param [in] ipmVariables: The IPM-related variables (a collection of slack and dual variables).
 * @param [in] dx: The state direction.
 * @param [in] du: The input direction.
 * @return The direction of the IPM-related variables.
 */
IpmVariablesDirection retriveIntermediateIpmVariablesDirection(const ModelData& modelData, const IpmData& ipmData, const IpmVariables& ipmVariables,
                                                               const vector_t& dx, const vector_t& du);

/**
 * Retrives the Newton direction of IPM-related variables from the directions of primal variables at a jump event time.
 *
 * @param [in] modelData: The discrete-time LQ approximation model data.
 * @param [in] ipmData: The IPM-related data.
 * @param [in] ipmVariables: The IPM-related variables (a collection of slack and dual variables).
 * @param [in] dx: The state direction.
 * @param [out] ipmVariablesDirection: The direction of the IPM-related variables.
 */
void retrivePreJumpIpmVariablesDirection(const ModelData& modelData, const IpmData& ipmData, const IpmVariables& ipmVariables,
                                         const vector_t& dx, IpmVariablesDirection& ipmVariablesDirection);

/**
 * Retrives the Newton direction of IPM-related variables from the directions of primal variables at a jump event time.
 *
 * @param [in] modelData: The discrete-time LQ approximation model data.
 * @param [in] ipmData: The IPM-related data.
 * @param [in] ipmVariables: The IPM-related variables (a collection of slack and dual variables).
 * @param [in] dx: The state direction.
 * @return The direction of the IPM-related variables.
 */
IpmVariablesDirection retrivePreJumpIpmVariablesDirection(const ModelData& modelData, const IpmData& ipmData, const IpmVariables& ipmVariables,
                                                          const vector_t& dx);

/**
 * Retrives the Newton direction of IPM-related variables from the directions of primal variables at final time.
 *
 * @param [in] modelData: The discrete-time LQ approximation model data.
 * @param [in] ipmData: The IPM-related data.
 * @param [in] ipmVariables: The IPM-related variables (a collection of slack and dual variables).
 * @param [in] dx: The state direction.
 * @param [out] ipmVariablesDirection: The direction of the IPM-related variables.
 */
void retriveFinalIpmVariablesDirection(const ModelData& modelData, const IpmData& ipmData, const IpmVariables& ipmVariables,
                                       const vector_t& dx, IpmVariablesDirection& ipmVariablesDirection);

/**
 * Retrives the Newton direction of IPM-related variables from the directions of primal variables at final time.
 *
 * @param [in] modelData: The discrete-time LQ approximation model data.
 * @param [in] ipmData: The IPM-related data.
 * @param [in] ipmVariables: The IPM-related variables (a collection of slack and dual variables).
 * @param [in] dx: The state direction.
 * @return The direction of the IPM-related variables.
 */
IpmVariablesDirection retriveFinalIpmVariablesDirection(const ModelData& modelData, const IpmData& ipmData, const IpmVariables& ipmVariables,
                                                        const vector_t& dx);

/**
 * Computes the primal step size via the fraction-to-boundary rule.
 *
 * @param [in] ipmData: The IPM-related data.
 * @param [in] ipmVariables: The IPM-related variables (a collection of slack and dual variables).
 * @param [in] ipmVariablesDirection: The direction of the IPM-related variables.
 * @param [in] marginRate: Margin rate of the fraction-to-boundary rule. Must be between 0 to 1.0. Default is 0.995.
 * @return The primal step size.
 */
scalar_t intermediatePrimalStepSize(const IpmData& ipmData, const IpmVariables& ipmVariables, const IpmVariablesDirection& ipmVariablesDirection, 
                                    scalar_t marginRate=0.995);

/**
 * Computes the dual step size via the fraction-to-boundary rule.
 *
 * @param [in] ipmData: The IPM-related data.
 * @param [in] ipmVariables: The IPM-related variables (a collection of slack and dual variables).
 * @param [in] ipmVariablesDirection: The direction of the IPM-related variables.
 * @param [in] marginRate: Margin rate of the fraction-to-boundary rule. Must be between 0 to 1.0. Default is 0.995.
 * @return The dual step size.
 */
scalar_t intermediateDualStepSize(const IpmData& ipmData, const IpmVariables& ipmVariables, const IpmVariablesDirection& ipmVariablesDirection, 
                                  scalar_t marginRate=0.995);

/**
 * Computes the primal step size via the fraction-to-boundary rule.
 *
 * @param [in] ipmData: The IPM-related data.
 * @param [in] ipmVariables: The IPM-related variables (a collection of slack and dual variables).
 * @param [in] ipmVariablesDirection: The direction of the IPM-related variables.
 * @param [in] marginRate: Margin rate of the fraction-to-boundary rule. Must be between 0 to 1.0. Default is 0.995.
 * @return The primal step size.
 */
scalar_t preJumpPrimalStepSize(const IpmData& ipmData, const IpmVariables& ipmVariables, const IpmVariablesDirection& ipmVariablesDirection, 
                               scalar_t marginRate=0.995);

/**
 * Computes the dual step size via the fraction-to-boundary rule.
 *
 * @param [in] ipmData: The IPM-related data.
 * @param [in] ipmVariables: The IPM-related variables (a collection of slack and dual variables).
 * @param [in] ipmVariablesDirection: The direction of the IPM-related variables.
 * @param [in] marginRate: Margin rate of the fraction-to-boundary rule. Must be between 0 to 1.0. Default is 0.995.
 * @return The dual step size.
 */
scalar_t preJumpDualStepSize(const IpmData& ipmData, const IpmVariables& ipmVariables, const IpmVariablesDirection& ipmVariablesDirection, 
                             scalar_t marginRate=0.995);

/**
 * Computes the primal step size via the fraction-to-boundary rule.
 *
 * @param [in] ipmData: The IPM-related data.
 * @param [in] ipmVariables: The IPM-related variables (a collection of slack and dual variables).
 * @param [in] ipmVariablesDirection: The direction of the IPM-related variables.
 * @param [in] marginRate: Margin rate of the fraction-to-boundary rule. Must be between 0 to 1.0. Default is 0.995.
 * @return The primal step size.
 */
scalar_t finalPrimalStepSize(const IpmData& ipmData, const IpmVariables& ipmVariables, const IpmVariablesDirection& ipmVariablesDirection, 
                             scalar_t marginRate=0.995);

/**
 * Computes the dual step size via the fraction-to-boundary rule.
 *
 * @param [in] ipmData: The IPM-related data.
 * @param [in] ipmVariables: The IPM-related variables (a collection of slack and dual variables).
 * @param [in] ipmVariablesDirection: The direction of the IPM-related variables.
 * @param [in] marginRate: Margin rate of the fraction-to-boundary rule. Must be between 0 to 1.0. Default is 0.995.
 * @return The dual step size.
 */
scalar_t finalDualStepSize(const IpmData& ipmData, const IpmVariables& ipmVariables, const IpmVariablesDirection& ipmVariablesDirection, 
                           scalar_t marginRate=0.995);

}  // namespace ipm
}  // namespace ocs2
