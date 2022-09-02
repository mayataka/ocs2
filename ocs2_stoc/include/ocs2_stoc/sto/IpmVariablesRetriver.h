#pragma once

#include <ocs2_core/Types.h>
#include <ocs2_sto/model_data/StoModelData.h>
#include <ocs2_ipm_oc/oc_data/IpmVariables.h>
#include <ocs2_ipm_oc/oc_data/IpmData.h>

namespace ocs2 {

/**
 * Retrives the Newton direction of IPM-related variables from the directions of primal variables.
 *
 * @param [in] stoModelData: The quadratic approximation model data.
 * @param [in] ipmData: The IPM-related data.
 * @param [in] ipmVariables: The IPM-related variables (a collection of slack and dual variables).
 * @param [in] dts: The switching time directions.
 * @param [out] ipmVariablesDirection: The direction of the IPM-related variables.
 */
void retriveStoIpmVariablesDirection(const StoModelData& stoModelData, const ipm::IpmData& ipmData, const ipm::IpmVariables& ipmVariables, 
                                     const vector_t& dts, ipm::IpmVariablesDirection& ipmVariablesDirection);

/**
 * Retrives the Newton direction of IPM-related variables from the directions of primal variables.
 *
 * @param [in] stoModelData: The quadratic approximation model data.
 * @param [in] ipmData: The IPM-related data.
 * @param [in] ipmVariables: The IPM-related variables (a collection of slack and dual variables).
 * @param [in] dts: The switching time directions.
 * @return The direction of the IPM-related variables.
 */
ipm::IpmVariablesDirection retriveStoIpmVariablesDirection(const StoModelData& stoModelData, const ipm::IpmData& ipmData, 
                                                           const ipm::IpmVariables& ipmVariables, const vector_t& dts);

/**
 * Computes the primal step size via the fraction-to-boundary rule.
 *
 * @param [in] ipmData: The IPM-related data.
 * @param [in] ipmVariables: The IPM-related variables (a collection of slack and dual variables).
 * @param [in] ipmVariablesDirection: The direction of the IPM-related variables.
 * @param [in] marginRate: Margin rate of the fraction-to-boundary rule. Must be between 0 to 1.0. Default is 0.995.
 * @return The primal step size.
 */
scalar_t stoPrimalStepSize(const ipm::IpmData& ipmData, const ipm::IpmVariables& ipmVariables, 
                           const ipm::IpmVariablesDirection& ipmVariablesDirection, scalar_t marginRate=0.995);

/**
 * Computes the dual step size via the fraction-to-boundary rule.
 *
 * @param [in] ipmData: The IPM-related data.
 * @param [in] ipmVariables: The IPM-related variables (a collection of slack and dual variables).
 * @param [in] ipmVariablesDirection: The direction of the IPM-related variables.
 * @param [in] marginRate: Margin rate of the fraction-to-boundary rule. Must be between 0 to 1.0. Default is 0.995.
 * @return The dual step size.
 */
scalar_t stoDualStepSize(const ipm::IpmData& ipmData, const ipm::IpmVariables& ipmVariables, 
                         const ipm::IpmVariablesDirection& ipmVariablesDirection, scalar_t marginRate=0.995);

}  // namespace ocs2
