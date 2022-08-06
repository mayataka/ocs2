#pragma once

#include <ocs2_core/Types.h>
#include <ocs2_ipm/oc_data/IpmVariables.h>
#include <ocs2_ipm/oc_data/IpmData.h>
#include <ocs2_sto/model_data/StoModelData.h>

namespace ocs2 {

/**
 * Initializes the IPM variables (slack and dual variables).
 *
 * @param [in] stoModelData: The STO model data.
 * @param [out] ipmVariables: The ipm variables.
 * @param [in] barrier: The barrier parameter of the IPM.
 */
void initIpmVariables(const StoModelData& stoModelData, ipm::IpmVariables& ipmVariables, scalar_t barrier);

/**
 * Initializes the IPM variables (slack and dual variables).
 *
 * @param [in] stoModelData: The STO model data.
 * @param [in] barrier: The barrier parameter of the IPM.
 * @return The ipm variables.
 */
ipm::IpmVariables initIpmVariables(const StoModelData& stoModelData, scalar_t barrier);

/**
 * Initializes the IPM variables direction (directions of the slack and dual variables).
 *
 * @param [in] stoModelData: The STO model data.
 * @param [out] ipmVariablesDirection: The ipm variables direction.
 */
void initIpmVariablesDirection(const StoModelData& stoModelData, ipm::IpmVariablesDirection& ipmVariablesDirection);

/**
 * Initializes the IPM variables direction (directions of the slack and dual variables).
 *
 * @param [in] stoModelData: The STO model data.
 * @return The ipm variables direction.
 */
ipm::IpmVariablesDirection initIpmVariablesDirection(const StoModelData& stoModelData);

}  // namespace ocs2
