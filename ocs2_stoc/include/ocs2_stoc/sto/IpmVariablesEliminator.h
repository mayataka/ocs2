#pragma once

#include <ocs2_core/Types.h>
#include <ocs2_sto/model_data/StoModelData.h>
#include <ocs2_ipm_oc/oc_data/IpmVariables.h>
#include <ocs2_ipm_oc/oc_data/IpmData.h>

namespace ocs2 {

/**
 * Eliminates the IPM-related variables and inequality constraints from a quadratic approximation of the STO problem.
 *
 * @param [in] ipmVariables: The IPM-related variables (a collection of slack and dual variables).
 * @param [in, out] stoModelData: The quadratic approximation model data.
 * @param [out] ipmData: The IPM-related data.
 * @param [out] barrierParam: The barrier parameter of the IPM.
 */
void eliminateIpmVariablesSTO(const ipm::IpmVariables& ipmVariables, StoModelData& stoModelData, ipm::IpmData& ipmData, 
                              scalar_t barrierParam);

/**
 * Eliminates the IPM-related variables and inequality constraints from a quadratic approximation of the STO problem.
 *
 * @param [in] ipmVariables: The IPM-related variables (a collection of slack and dual variables).
 * @param [in, out] stoModelData: The quadratic approximation model data.
 * @param [out] barrierParam: The barrier parameter of the IPM.
 * @return The IPM-related data.
 */
ipm::IpmData eliminateIpmVariablesSTO(const ipm::IpmVariables& ipmVariables, StoModelData& stoModelData, scalar_t barrierParam);

}  // namespace ocs2
